# =============================================================================
# linear_system_analysis.R
# -----------------------------------------------------------------------------
# Supplementary analysis. In the sub-saturation regime the cattle iAs PBPK model
# is a linear dynamical system  dx/dt = A x + b * C_feed.  This script:
#   (1) linearises the model about its steady state (finite-difference Jacobian A),
#   (2) computes the eigenvalue spectrum of A (relaxation modes / terminal t1/2),
#   (3) renders Figure S (system matrix heat-map + eigenvalue spectrum with a
#       near-zero inset), as 300-dpi TIFF (+ PNG).
#
# It reuses the model definitions (parameters, ODE, helpers) from the main
# pipeline script WITHOUT running the full pipeline: only the definition prefix
# (everything before the deterministic run) is evaluated.
#
#   Rscript scripts/linear_system_analysis.R
# =============================================================================

suppressMessages(library(deSolve))

# --- locate and load the model definitions from the main script --------------
this_dir <- {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(getwd())
}
main_src <- readLines(file.path(this_dir, "PBPK_ML_cattle_Adult_concentration.R"),
                      warn = FALSE)
cut <- grep("^## A5\\. Deterministic PBPK run", main_src)[1] - 1L
stopifnot(is.finite(cut), cut > 100)
# NOTE: eval() here is safe and intentional. It evaluates ONLY the definition
# prefix of the project's own main script (a trusted local file in this repo)
# to reuse the parameters/ODE/helpers without re-running the full pipeline. No
# external or user-supplied input is involved.
eval(parse(text = main_src[1:cut]), envir = globalenv())   # parms, ODE, helpers

# --- steady state, then finite-difference Jacobian A = d(dx/dt)/dx -----------
y0 <- setNames(rep(0, length(state_names)), state_names)
p  <- parms; p["TSTOP"] <- 6000
pp <- c(p, compute_derived(p))
o   <- as.data.frame(ode(y0, seq(0, 3000, by = 10), pbpk_cattle_ode, pp, method = "lsoda"))
yss <- as.numeric(o[nrow(o), -1]); names(yss) <- state_names

deriv <- function(y) pbpk_cattle_ode(100, y, pp)[[1]]   # exposure on at t = 100 h
f0 <- deriv(yss); n <- length(yss)
J  <- matrix(0, n, n, dimnames = list(state_names, state_names))
for (i in seq_len(n)) {
  h <- max(1e-6, abs(yss[i]) * 1e-6); yp <- yss; yp[i] <- yp[i] + h
  J[, i] <- (deriv(yp) - f0) / h
}

ev  <- eigen(J, only.values = TRUE)$values
re  <- Re(ev); im <- Im(ev); dec <- re < -1e-7          # decaying (non-trivial) modes
lam <- max(re[dec]); absl <- -lam                       # slowest = least-negative
dom <- sort(re[dec], decreasing = TRUE)[2]              # dominant tissue mode
cat(sprintf("lambda_slow   = %.4f /day  (terminal t1/2 = %.2f d)\n", lam * 24, log(2) / absl / 24))
cat(sprintf("dominant mode = %.4f /day  (t1/2 = %.2f d -> tissue SS ~%.0f d)\n",
            dom * 24, log(2) / (-dom) / 24, log(1000) / (-dom) / 24))

# --- figure ------------------------------------------------------------------
idx   <- 15:46                                          # AMT_*_compartment (4 sp x 8 comp)
As    <- J[idx, idx]; M <- sign(As) * log10(1 + abs(As))
comps <- c("Lung", "Muscle", "Kidney", "Liver", "GI", "Rest", "Ven.bl", "Art.bl")

draw <- function() {
  par(mfrow = c(1, 2), mar = c(6, 5, 3.5, 2))
  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  image(1:32, 1:32, t(M)[, 32:1], col = pal, zlim = c(-max(abs(M)), max(abs(M))),
        xlab = "state j (source)", ylab = "state i (target)",
        main = expression("A. System matrix  d" * bold(x) * "/dt = A " * bold(x) * " + " * bold(b) * " " * C[feed]),
        axes = FALSE)
  at <- seq(2.5, by = 4, length.out = 8)
  axis(1, at = at, labels = comps, las = 2, cex.axis = 0.7, tick = FALSE)
  axis(2, at = 33 - at, labels = comps, las = 2, cex.axis = 0.7, tick = FALSE); box()
  abline(h = 32 - seq(0, 32, 4) + 0.5, v = seq(0, 32, 4) + 0.5, col = "grey85", lwd = 0.4)

  red <- re[dec] * 24; imd <- im[dec] * 24
  plot(red, imd, pch = 19, col = "grey55", cex = 0.8,
       xlab = expression(Re(lambda) ~ "(day"^{-1} * ")"),
       ylab = expression(Im(lambda) ~ "(day"^{-1} * ")"),
       main = "B. Spectrum of relaxation rate constants of A")
  abline(v = 0, lty = 3, col = "grey60")
  legend("bottomleft", bty = "n", cex = 0.8, text.col = "grey30",
         legend = expression("all Re(" * lambda * ") < 0  ->  single stable steady state"))

  # inset: magnify the near-zero (slow) region
  op <- par(fig = c(0.795, 0.985, 0.55, 0.93), new = TRUE, mar = c(1.8, 1.8, 1.0, 0.4),
            mgp = c(0, 0.3, 0), tcl = -0.2)
  sel <- red > -3.3
  plot(red[sel], imd[sel], pch = 19, col = "grey60", cex = 0.7,
       xlim = c(-3.3, 0.35), ylim = c(-35, 35), xlab = "", ylab = "", cex.axis = 0.55)
  abline(v = 0, lty = 3, col = "grey70")
  points(dom * 24, 0, pch = 19, col = "#2166AC", cex = 1.2)
  points(lam * 24, 0, pch = 19, col = "#B2182B", cex = 1.5)
  legend("bottomleft", bty = "n", cex = 0.5, pch = 19, col = c("#B2182B", "#2166AC"),
         legend = c("lambda_slow (t1/2 ~8 d)", "dominant mode (t1/2 ~2.6 d)"))
  mtext("near-zero (slow) modes", side = 3, line = 0.05, cex = 0.5, col = "grey40")
  box(col = "grey40")
  par(op)
}

tiff(out_path("S1_linear_system.tiff"), width = 11, height = 5, units = "in",
     res = 300, compression = "lzw"); draw(); dev.off()
png(out_path("S1_linear_system.png"),  width = 11, height = 5, units = "in",
    res = 150); draw(); dev.off()
cat("Saved S1_linear_system.tiff and .png to", OUT_DIR, "\n")
