# Cattle BW <-> FR_feed correlation panel for Figure 3 (Panel C)
# Reproduces the per-animal samples used by the main pipeline (same seeds:
# 4729 for BW, 4730 for the correlated FR_feed) and plots the scatter with a
# linear fit, matching the style of the consumer-side panels (A, B).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Truncated-normal helper (same as PBPK_MRL_cattle_Adult.R)
rtnorm <- function(n, mean, sd, lo, hi) {
  p_lo <- pnorm(lo, mean, sd); p_hi <- pnorm(hi, mean, sd)
  qnorm(runif(n, p_lo, p_hi), mean, sd)
}

n <- 10000L

# BW: same call as the first rtnorm in parm_samples
set.seed(4729)
BW <- rtnorm(n, mean = 621, sd = 6.21, lo = 602.37, hi = 639.63)

# FR_feed: Gaussian z-score combination with target r = 0.5 (Waegeneers 2011)
fr_mean <- 9.42; fr_sd <- 0.94
fr_lo <- fr_mean - 3 * fr_sd; fr_hi <- fr_mean + 3 * fr_sd
set.seed(4730)
z_bw <- (BW - 621) / 6.21
eps  <- rnorm(n)
z_fr <- 0.5 * z_bw + sqrt(1 - 0.5^2) * eps
fr   <- pmin(pmax(fr_mean + z_fr * fr_sd, fr_lo), fr_hi)

r_real <- cor(BW, fr)
cat(sprintf("Realized cattle r(BW, FR_feed) = %.3f (target 0.50)\n", r_real))

df <- data.frame(BW = BW, FR = fr)

p <- ggplot(df, aes(BW, FR)) +
  geom_point(alpha = 0.10, size = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.3f (target = 0.50)", r_real),
           hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold") +
  labs(title = "Cattle  -  BW vs FR_feed (population PBPK)",
       subtitle = sprintf("n = %s virtual cattle; red line = linear fit",
                          format(n, big.mark = ",")),
       x = "Cattle body weight (kg)",
       y = expression(paste("FR"["feed"], " (kg DM/day)"))) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
out_dir <- normalizePath(file.path(script_dir, "..", "output"), mustWork = FALSE)
if (!dir.exists(out_dir)) out_dir <- "D:/repos/PBPK_MRL_iAs_cattle/output"

ggsave(file.path(out_dir, "Figure3C_Cattle_BW_FR_correlation.png"),
       p, width = 5, height = 4, dpi = 300)
ggsave(file.path(out_dir, "Figure3C_Cattle_BW_FR_correlation.tiff"),
       p, width = 5, height = 4, dpi = 300, compression = "lzw")

cat("Saved:\n  ", file.path(out_dir, "Figure3C_Cattle_BW_FR_correlation.png"), "\n  ",
    file.path(out_dir, "Figure3C_Cattle_BW_FR_correlation.tiff"), "\n")
