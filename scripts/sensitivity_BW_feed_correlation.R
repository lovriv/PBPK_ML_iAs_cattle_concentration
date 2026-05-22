# =============================================================================
# sensitivity_BW_feed_correlation.R
#
# Question (Prof. Wu): does imposing a correlation between CATTLE body weight
# and CATTLE feed intake change the transfer factors (TF) or the derived MRL?
#
# This script reuses the exact PBPK model definitions from
# PBPK_MRL_cattle_Adult.R (lines 1:457, i.e. everything BEFORE the 10000-run
# population Monte Carlo) and runs a small, fast experiment:
#
#   For a fixed cattle body weight, vary the daily feed intake (FR_feed) and
#   make the ingested iAs dose proportional to it (dose = FR_feed x C_feed).
#   Recompute TF = C_tissue / (C_feed x FR_feed) at steady state.
#
# If TF is invariant to FR_feed (for a given BW), then correlating FR_feed with
# BW cannot change the TF distribution -- and since cattle BW does not appear in
# the MRL exposure equation (Eq. 7), the MRL is unaffected as well.
# =============================================================================

# Robustly locate the main script relative to this file
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
here <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
main <- file.path(here, "PBPK_MRL_cattle_Adult.R")
stopifnot("Cannot find main pipeline script" = file.exists(main))

# --- Load model definitions only (lines 1:457, before the population MC) ------
src <- readLines(main, encoding = "UTF-8", warn = FALSE)
cut <- grep("^## A7\\. Population PBPK setup", src)[1] - 1L
cat(sprintf("Sourcing model definitions: lines 1:%d of %s\n", cut, basename(main)))
eval(parse(text = paste(src[seq_len(cut)], collapse = "\n")), envir = globalenv())

MW_As <- unname(parms["MW"])
C_FEED <- unname(parms["C_feed_iAs_ug_kg"])      # 103.3 ug/kg
# iAs species split in the Hung 2021 control feed (ug per kg DM)
conc <- c(AsIII = 14.3, AsV = 89.0, MMA = 14.4, DMA = 1.8)

# --- One PBPK solve at a given (BW, feed intake) ------------------------------
solve_tf <- function(BW, FR) {
  p <- parms
  p["BW"]           <- BW
  p["PdoseC_AsIII"] <- FR * conc["AsIII"]   # dose proportional to feed intake
  p["PdoseC_AsV"]   <- FR * conc["AsV"]
  p["PdoseC_MMA"]   <- FR * conc["MMA"]
  p["PdoseC_DMA"]   <- FR * conc["DMA"]
  allp <- c(p, compute_derived(p))
  out  <- as.data.frame(deSolve::ode(y = y0, times = times,
                                     func = pbpk_cattle_ode,
                                     parms = allp, method = "lsoda"))
  ss <- out[nrow(out), ]
  Cmus_uM <- (ss$AMT_AsIII_muscle + ss$AMT_AsV_muscle) / allp["V_muscle"]
  Cliv_uM <- (ss$AMT_AsIII_liv    + ss$AMT_AsV_liv)    / allp["V_liv"]
  Ckid_uM <- (ss$AMT_AsIII_kid    + ss$AMT_AsV_kid)    / allp["V_kid"]
  # TF = (C_tissue ug/kg) / (C_feed ug/kg * FR kg/day)   -> d/kg
  data.frame(
    BW = BW, FR_feed = FR,
    Cmus_uM = unname(Cmus_uM), Cliv_uM = unname(Cliv_uM), Ckid_uM = unname(Ckid_uM),
    TF_muscle = unname(Cmus_uM * MW_As) / (C_FEED * FR),
    TF_liver  = unname(Cliv_uM * MW_As) / (C_FEED * FR),
    TF_kidney = unname(Ckid_uM * MW_As) / (C_FEED * FR)
  )
}

# --- Experiment: 3 body weights x 5 feed-intake levels ------------------------
BWs <- c(540, 621, 700)
FRs <- c(6, 8, 9.42, 11, 13)            # kg DM/day  (baseline = 9.42)
grid <- expand.grid(BW = BWs, FR = FRs)
res <- do.call(rbind, Map(solve_tf, grid$BW, grid$FR))
res <- res[order(res$BW, res$FR_feed), ]

cat("\n================ TF vs feed intake, by body weight ================\n")
print(round(res[, c("BW","FR_feed","TF_muscle","TF_liver","TF_kidney")], 7), row.names = FALSE)

# --- Quantify: spread of TF across feed-intake levels, within each BW ---------
cat("\n================ Coefficient of variation of TF across the 5 feed-intake levels (per BW) ================\n")
agg <- do.call(rbind, lapply(split(res, res$BW), function(d) {
  data.frame(BW = d$BW[1],
             TFmus_cv_pct = 100 * sd(d$TF_muscle) / mean(d$TF_muscle),
             TFliv_cv_pct = 100 * sd(d$TF_liver)  / mean(d$TF_liver),
             TFkid_cv_pct = 100 * sd(d$TF_kidney) / mean(d$TF_kidney))
}))
print(format(agg, digits = 3), row.names = FALSE)

# --- Linearity check: steady-state tissue conc vs Michaelis-Menten Km ---------
Km <- unname(parms["Km_AsIII_to_MMA_liv"])
cat(sprintf("\nLinearity check: max steady-state tissue conc = %.4f umol/L vs Km = %.0f umol/L (ratio ~ %.0e)\n",
            max(res$Cmus_uM, res$Cliv_uM, res$Ckid_uM), Km,
            max(res$Cmus_uM, res$Cliv_uM, res$Ckid_uM) / Km))

cat("\nINTERPRETATION:\n")
cat("  If TF CV across feed-intake levels is ~0, TF is independent of feed intake.\n")
cat("  Then any BW<->feed-intake correlation leaves TF (and hence MRL) unchanged.\n")
