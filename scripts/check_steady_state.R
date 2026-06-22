# =============================================================================
# check_steady_state.R
# Time-to-steady-state diagnostics for the population PBPK output.
#
# Reads the deterministic time course written by the main pipeline
# (default: ../output_concentration/A_pbpk_deterministic.csv) and reports,
# per edible tissue (iAs = AsIII + AsV):
#   (1) time to reach 90/95/99/99.9% of the plateau (value at the final time)
#   (2) rate-based time: when the relative rate |dC/dt| / C_plateau drops
#       below 1.0, 0.5, 0.1 %/day, plus the residual rate at the final time.
#
# Usage:
#   Rscript check_steady_state.R [path/to/A_pbpk_deterministic.csv]
# (run the main pipeline first to generate the CSV.)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else
  file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
            "..", "output_concentration", "A_pbpk_deterministic.csv")
stopifnot("deterministic CSV not found" = file.exists(f))

df    <- read.csv(f)
hours <- if ("time" %in% names(df)) df$time else df$Day * 24
tis   <- c(muscle = "C_iAs_muscle", liver = "C_iAs_liv", kidney = "C_iAs_kid")

cat(sprintf("Source: %s\nFinal time = %.0f h (%.1f d), %d points\n\n",
            f, max(hours), max(hours) / 24, nrow(df)))

## (1) Fraction-of-plateau criterion -----------------------------------------
pct <- c(0.90, 0.95, 0.99, 0.999)
cat("(1) Time to reach a fraction of plateau (hours):\n")
cat(sprintf("%-7s %12s", "tissue", "C_plateau"))
for (p in pct) cat(sprintf("  t%.1f%%", p * 100)); cat("\n")
for (nm in names(tis)) {
  C <- df[[tis[[nm]]]]; Css <- C[length(C)]
  cat(sprintf("%-7s %12.5f", nm, Css))
  for (p in pct) cat(sprintf("  %6.0f", hours[which(C >= p * Css)[1]]))
  cat("\n")
}

## (2) Rate-based criterion ---------------------------------------------------
thr <- c(1.0, 0.5, 0.1)   # %/day relative to plateau
cat("\n(2) Rate-based: time when |dC/dt| / C_plateau < threshold (hours):\n")
cat(sprintf("%-7s", "tissue"))
for (p in thr) cat(sprintf("  t<%.1f%%/d", p)); cat("   rate_end(%/d)\n")
for (nm in names(tis)) {
  C    <- df[[tis[[nm]]]]; Css <- C[length(C)]
  dCdt <- diff(C) / diff(hours)
  tmid <- (hours[-1] + hours[-length(hours)]) / 2
  rel  <- abs(dCdt) * 24 / Css * 100
  cat(sprintf("%-7s", nm))
  for (p in thr) {
    idx <- which(rel < p)[1]
    cat(sprintf("  %8.0f", if (is.na(idx)) NA_real_ else tmid[idx]))
  }
  cat(sprintf("   %10.3f\n", rel[length(rel)]))
}
cat("\nPractical steady state ~ t99 (criterion 1) or t(<1%/day) (criterion 2).\n")
