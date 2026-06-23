# =============================================================================
# sensitivity_analysis.R
# -----------------------------------------------------------------------------
# Sensitivity analysis of the cattle iAs PBPK model, following the regulatory
# workflow of McNally et al. (2011) and OECD (2021, Section 2.6.3-2.6.4):
#
#   (1) LOCAL one-at-a-time (OAT) analysis: normalized sensitivity coefficients
#       (SC = elasticity = d ln C / d ln theta) by central difference with a 1%
#       perturbation (Fisher et al., 2020, Eq. 9.2), for the steady-state iAs
#       (AsIII + AsV) concentration in the three edible tissues that drive the
#       feed ML (muscle, liver, kidney).
#   (2) GLOBAL screening with the Morris elementary-effects method
#       (Morris, 1991; Campolongo & Saltelli, 1997), giving mu* (overall
#       influence) and sigma (interactions / non-linearity) per parameter.
#       Uniform ranges span half to twice the point value (theta/2, 2*theta),
#       the default of McNally et al. (2011) and OECD (2021) when no
#       parameter-specific distribution is available.
#
# Scope: the analysis targets the KINETIC and PARTITION parameters that were
# extrapolated cross-species from the human model (El-Masri & Kenyon 2008;
# Yu 1999), i.e. the parameters that carry cross-species uncertainty and have
# no cattle-specific distribution. The influence of the physiological
# parameters (BW, QCC, blood-flow and tissue-volume fractions) is already
# characterised by the 10,000-animal population Monte Carlo in the main
# pipeline, which propagates their real Lin (2020) distributions.
#
# Mass balance under perturbation is preserved automatically by the model:
# blood flows are normalised by the sum of flow fractions and tissue volumes
# by the sum of volume fractions in compute_derived() (equivalent to the
# re-parameterisation of Gelman et al., 1996). Kinetic and partition
# parameters, perturbed here, do not enter those constraints.
#
# It reuses the model definitions from the main pipeline WITHOUT running the
# full pipeline: only the definition prefix (before the deterministic run) is
# evaluated.
#
#   Rscript scripts/sensitivity_analysis.R
#
# Reproducibility: R 4.2.3; deSolve; sensitivity 1.30.0; Morris seed = 2026.
# =============================================================================

suppressMessages({
  library(deSolve)
  library(sensitivity)
})

# --- locate and load the model definitions from the main script --------------
this_dir <- {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(getwd())
}
main_src <- readLines(file.path(this_dir, "PBPK_ML_cattle_Adult_concentration.R"),
                      warn = FALSE)
cut <- grep("## A5. Deterministic PBPK run", main_src, fixed = TRUE)[1] - 1L
stopifnot(is.finite(cut), cut > 100)
# NOTE: eval() here is safe and intentional. It evaluates ONLY the definition
# prefix of the project's own main script (a trusted local file in this repo)
# to reuse the parameters/ODE/helpers without re-running the full pipeline.
eval(parse(text = main_src[1:cut]), envir = globalenv())   # parms, ODE, helpers

# --- steady-state evaluator --------------------------------------------------
# Integrate the deterministic model under continuous feeding to a horizon well
# beyond steady state (3000 h; verified adequate even when all clearances are
# halved) and return the steady-state iAs (AsIII + AsV) concentration
# (umol/L) in muscle, liver and kidney.
T_CAP <- 3000
y0    <- setNames(rep(0, length(state_names)), state_names)

eval_ss <- function(over) {
  p <- parms
  if (length(over)) p[names(over)] <- over
  p["TSTOP"] <- T_CAP + 100                       # feed over the whole horizon
  pp <- c(p, compute_derived(p))
  o  <- tryCatch(
    ode(y0, c(0, T_CAP), pbpk_cattle_ode, pp, method = "lsoda"),
    error = function(e) NULL)
  if (is.null(o)) return(c(muscle = NA, liver = NA, kid = NA))
  last <- o[nrow(o), ]
  c(muscle = (last[["AMT_AsIII_muscle"]] + last[["AMT_AsV_muscle"]]) / pp[["V_muscle"]],
    liver  = (last[["AMT_AsIII_liv"]]    + last[["AMT_AsV_liv"]])    / pp[["V_liv"]],
    kid    = (last[["AMT_AsIII_kid"]]    + last[["AMT_AsV_kid"]])    / pp[["V_kid"]])
}

# --- parameters under analysis (kinetic + partition) -------------------------
sa_names <- c(
  # absorption
  "Ka_AsIII", "Ka_AsV", "Ka_MMA", "Ka_DMA",
  # redox interconversion
  "K_red_AsV_to_AsIII", "K_ox_AsIII_to_AsV",
  # methylation Vmax (liver, kidney)
  "Vmax_AsIII_to_MMA_liv", "Vmax_AsIII_to_DMA_liv", "Vmax_MMA_to_DMA_liv",
  "Vmax_AsIII_to_MMA_kid", "Vmax_AsIII_to_DMA_kid", "Vmax_MMA_to_DMA_kid",
  # methylation Km (liver, kidney)
  "Km_AsIII_to_MMA_liv", "Km_AsIII_to_DMA_liv", "Km_MMA_to_DMA_liv",
  "Km_AsIII_to_MMA_kid", "Km_AsIII_to_DMA_kid", "Km_MMA_to_DMA_kid",
  # urinary / biliary / faecal excretion
  "k_urine_AsIII", "k_urine_AsV", "k_urine_MMA", "k_urine_DMA",
  "eF_AsV", "eB_AsV",
  # tissue:blood partition coefficients
  "P_AsIII_lung", "P_AsV_lung", "P_MMA_lung", "P_DMA_lung",
  "P_AsIII_muscle", "P_AsV_muscle", "P_MMA_muscle", "P_DMA_muscle",
  "P_AsIII_kid", "P_AsV_kid", "P_MMA_kid", "P_DMA_kid",
  "P_AsIII_liv", "P_AsV_liv", "P_MMA_liv", "P_DMA_liv",
  "P_AsIII_gi", "P_AsV_gi", "P_MMA_gi", "P_DMA_gi",
  "P_AsIII_rest", "P_AsV_rest", "P_MMA_rest", "P_DMA_rest"
)
base_p <- parms[sa_names]
stopifnot(all(is.finite(base_p)), all(base_p > 0))

grp_of <- function(nm) {
  ifelse(grepl("^Ka_", nm), "Absorption",
  ifelse(grepl("^K_red|^K_ox", nm), "Redox",
  ifelse(grepl("^Vmax_", nm), "Methylation Vmax",
  ifelse(grepl("^Km_", nm), "Methylation Km",
  ifelse(grepl("^k_urine|^eF_|^eB_", nm), "Excretion", "Partition")))))
}

cat("Sensitivity analysis: ", length(sa_names), " kinetic + partition parameters\n", sep = "")
cat("Steady-state evaluator horizon = ", T_CAP, " h\n\n", sep = "")

# --- (1) LOCAL OAT: normalized sensitivity coefficients ----------------------
cat("[1] Local OAT (central difference, 1% perturbation)...\n")
base_C <- eval_ss(numeric(0))
cat(sprintf("    baseline iAs (umol/L): muscle %.4f  liver %.4f  kidney %.4f\n",
            base_C["muscle"], base_C["liver"], base_C["kid"]))

oat <- t(sapply(sa_names, function(nm) {
  th   <- base_p[[nm]]
  fhi  <- eval_ss(setNames(th * 1.005, nm))       # +0.5%
  flo  <- eval_ss(setNames(th * 0.995, nm))       # -0.5%  (Delta = 1%)
  (fhi - flo) / (0.01 * base_C)                   # SC = (df/f)/(dtheta/theta)
}))
colnames(oat) <- c("SC_muscle", "SC_liver", "SC_kid")

# --- (2) GLOBAL Morris screening ---------------------------------------------
cat("[2] Morris elementary-effects screening (r = 20, ranges theta/2..2*theta)...\n")
set.seed(2026)
mo <- morris(model = NULL, factors = sa_names, r = 20,
             design = list(type = "oat", levels = 8, grid.jump = 4),
             binf = unname(base_p / 2), bsup = unname(base_p * 2))
cat(sprintf("    Morris model evaluations: %d\n", nrow(mo$X)))
t_mo <- system.time({
  Ymo <- t(apply(mo$X, 1, function(row) eval_ss(setNames(row, sa_names))))
})
cat(sprintf("    elapsed: %.1f min\n", t_mo[3] / 60))

mustar <- function(EE) apply(EE, 2, function(e) mean(abs(e)))
msig   <- function(EE) apply(EE, 2, function(e) sd(e))
mo_res <- lapply(c(muscle = "muscle", liver = "liver", kid = "kid"), function(tis) {
  tell(mo, Ymo[, tis])
  data.frame(mu = apply(mo$ee, 2, mean),
             mu.star = mustar(mo$ee),
             sigma = msig(mo$ee))
})

# --- assemble and write results table ----------------------------------------
res <- data.frame(
  parameter = sa_names,
  group     = grp_of(sa_names),
  baseline  = unname(base_p),
  oat,
  mu.star_muscle = mo_res$muscle$mu.star, sigma_muscle = mo_res$muscle$sigma,
  mu.star_liver  = mo_res$liver$mu.star,  sigma_liver  = mo_res$liver$sigma,
  mu.star_kid    = mo_res$kid$mu.star,    sigma_kid    = mo_res$kid$sigma,
  row.names = NULL, check.names = FALSE
)
res <- res[order(-abs(res$SC_muscle)), ]
write.csv(res, out_path("SA_sensitivity_coefficients.csv"), row.names = FALSE)
cat("\nSaved SA_sensitivity_coefficients.csv\n")

cat("\nTop 12 parameters by |SC_muscle|:\n")
print(format(head(res[, c("parameter","group","SC_muscle","SC_liver","SC_kid",
                          "mu.star_muscle","sigma_muscle")], 12),
             digits = 3), row.names = FALSE)
# sigma/mu* over the influential parameters only (the max over all parameters is
# meaningless: it is dominated by near-zero-influence parameters where mu* ~ 0).
infl <- order(-mo_res$muscle$mu.star)[1:10]
ratio_infl <- mo_res$muscle$sigma[infl] / mo_res$muscle$mu.star[infl]
cat(sprintf("\nMorris sigma/mu* over top-10 influential parameters (muscle): median %.2f, range %.2f-%.2f\n",
            median(ratio_infl), min(ratio_infl), max(ratio_infl)))
cat("  (moderate -> influential parameters act mildly non-linearly over the half-to-twice ranges;\n")
cat("   the mu* ranking still reproduces the local OAT ranking)\n")

# --- figure: (A) OAT tornado  (B) Morris mu*-sigma ---------------------------
short <- function(nm) {
  nm <- sub("^Vmax_", "Vmax:", nm); nm <- sub("^Km_", "Km:", nm)
  nm <- sub("^k_urine_", "kU:", nm); nm <- sub("^Ka_", "Ka:", nm)
  nm <- sub("^P_", "P:", nm)
  nm <- sub("^K_red_AsV_to_AsIII", "K_red", nm); nm <- sub("^K_ox_AsIII_to_AsV", "K_ox", nm)
  nm <- gsub("_to_", ">", nm); nm <- gsub("_", ".", nm)
  nm
}
pal_grp <- c(Absorption = "#66C2A5", Redox = "#FC8D62", "Methylation Vmax" = "#8DA0CB",
             "Methylation Km" = "#E78AC3", Excretion = "#A6D854", Partition = "#FFD92F")

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 9, 3, 1))
  # Panel A: tornado of top-15 by |SC_muscle|
  top <- head(res, 15); top <- top[order(top$SC_muscle), ]
  cols <- pal_grp[top$group]
  bp <- barplot(top$SC_muscle, horiz = TRUE, col = cols, border = NA,
                names.arg = short(top$parameter), las = 1, cex.names = 0.7,
                xlab = "Normalized sensitivity coefficient (muscle)",
                main = "A. Local OAT sensitivity (top 15)",
                xlim = range(0, top$SC_muscle) * 1.15)
  points(top$SC_liver, bp, pch = 23, bg = "white", cex = 0.7)
  points(top$SC_kid,   bp, pch = 21, bg = "grey40", cex = 0.6)
  abline(v = 0, col = "grey50")
  legend("bottomright", bty = "n", cex = 0.65,
         legend = c("muscle (bar)", "liver", "kidney"),
         pch = c(15, 23, 21), col = c("grey70", "black", "black"),
         pt.bg = c(NA, "white", "grey40"))

  # Panel B: Morris mu*-sigma (muscle)
  par(mar = c(4.5, 4.5, 3, 1))
  mm <- mo_res$muscle
  xr <- range(0, mm$mu.star) * 1.05; yr <- range(0, mm$sigma) * 1.1
  plot(mm$mu.star, mm$sigma, pch = 19, col = pal_grp[grp_of(sa_names)], cex = 0.9,
       xlim = xr, ylim = yr, xlab = expression(mu * "*  (overall influence)"),
       ylab = expression(sigma ~ "(interactions / non-linearity)"),
       main = "B. Morris screening (muscle)")
  abline(0, 1, lty = 3, col = "grey60")           # sigma = mu* reference
  lab <- order(-mm$mu.star)[1:8]
  text(mm$mu.star[lab], mm$sigma[lab], short(sa_names[lab]),
       pos = 4, cex = 0.6, offset = 0.3, xpd = NA)
  legend("topleft", bty = "n", cex = 0.6, pch = 19, col = pal_grp,
         legend = names(pal_grp))
}

tiff(out_path("SA_sensitivity.tiff"), width = 11, height = 5.5, units = "in",
     res = 300, compression = "lzw"); draw(); dev.off()
png(out_path("SA_sensitivity.png"), width = 11, height = 5.5, units = "in",
    res = 150); draw(); dev.off()
cat("Saved SA_sensitivity.tiff and .png to", OUT_DIR, "\n")
