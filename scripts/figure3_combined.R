# Figure 3: combined 3-panel image
#   (A) consumer BW vs dietary meat intake
#   (B) consumer BW vs dietary other-offal intake
#   (C) cattle BW vs FR_feed (population PBPK)
# Each panel reproduces the per-animal sampling used by the pipeline
# (target Pearson r = 0.50, Gaussian z-score combination).

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(readxl)
  library(truncnorm)
  library(dplyr)
})

# ---- Paths -----------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
xlsx       <- file.path(repo_root, "data", "food_intake_iAs.xlsx")
out_dir    <- file.path(repo_root, "output")

# ---- Read consumer-side parameters ----------------------------------------
fi_data <- read_excel(xlsx, sheet = "food intake")
bw_data <- read_excel(xlsx, sheet = "bodyweight")

bw_adult <- bw_data |> filter(agegroup == "Adult")
bw_m  <- as.numeric(bw_adult$mean[1]); bw_s <- as.numeric(bw_adult$sd[1])
bw_lo <- as.numeric(bw_adult$min[1]);  bw_hi <- as.numeric(bw_adult$max[1])

get_fi <- function(tis, pop = "general public") {
  r <- fi_data |> filter(agegroup == "Adult", population == pop, tissue == tis)
  if (nrow(r) == 0) {
    # try alternate naming
    r <- fi_data |> filter(agegroup == "Adult", tissue == tis) |> slice(1)
  }
  list(mean = as.numeric(r$mean[1]), sd = as.numeric(r$sd[1]),
       min  = as.numeric(r$min[1]),  max = as.numeric(r$max[1]))
}

fi_meat <- get_fi("meat")
fi_oth  <- get_fi("others")

# ---- Helper: correlated truncated-normal pair (matching pipeline) ---------
sample_corr <- function(n, bw_m, bw_s, bw_lo, bw_hi,
                        fi_m, fi_s, fi_lo, fi_hi, r = 0.5) {
  BW <- rtruncnorm(n, a = bw_lo, b = bw_hi, mean = bw_m, sd = bw_s)
  z_bw <- (BW - bw_m) / bw_s
  eps  <- rnorm(n)
  z_fi <- r * z_bw + sqrt(1 - r^2) * eps
  fi_raw <- fi_m + z_fi * fi_s
  FI <- pmin(pmax(fi_raw, fi_lo), fi_hi)
  data.frame(BW = BW, FI = FI)
}

# ---- Panels A and B: consumer BW vs dietary intake (general public) -------
n <- 10000L
set.seed(123)
df_A <- sample_corr(n, bw_m, bw_s, bw_lo, bw_hi,
                    fi_meat$mean, fi_meat$sd, fi_meat$min, fi_meat$max, 0.5)
df_B <- sample_corr(n, bw_m, bw_s, bw_lo, bw_hi,
                    fi_oth$mean, fi_oth$sd, fi_oth$min, fi_oth$max, 0.5)
r_A <- cor(df_A$BW, df_A$FI); r_B <- cor(df_B$BW, df_B$FI)

# ---- Panel C: cattle BW vs FR_feed (same seeds as main pipeline) ----------
rtnorm <- function(n, mean, sd, lo, hi) {
  p_lo <- pnorm(lo, mean, sd); p_hi <- pnorm(hi, mean, sd)
  qnorm(runif(n, p_lo, p_hi), mean, sd)
}
set.seed(4729)
BW_c <- rtnorm(n, mean = 621, sd = 6.21, lo = 602.37, hi = 639.63)
set.seed(4730)
z_bw <- (BW_c - 621) / 6.21
eps  <- rnorm(n)
z_fr <- 0.5 * z_bw + sqrt(0.75) * eps
fr_lo <- 9.42 - 3 * 0.94; fr_hi <- 9.42 + 3 * 0.94
FR_c <- pmin(pmax(9.42 + z_fr * 0.94, fr_lo), fr_hi)
df_C <- data.frame(BW = BW_c, FI = FR_c)
r_C  <- cor(df_C$BW, df_C$FI)

cat(sprintf("Realized r: A(meat)=%.3f  B(others)=%.3f  C(cattle)=%.3f\n",
            r_A, r_B, r_C))

# ---- Plot helper (consistent style across panels) ------------------------
mk_panel <- function(d, r_real, tag, title, xlab, ylab) {
  ggplot(d, aes(BW, FI)) +
    geom_point(alpha = 0.10, size = 0.4, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.3f (target = 0.50)", r_real),
             hjust = -0.08, vjust = 1.5, size = 3.6, fontface = "bold") +
    labs(title = title, tag = tag, x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.tag      = element_text(face = "bold", size = 14),
          plot.tag.position = c(0.02, 0.98))
}

pA <- mk_panel(df_A, r_A, "A", "Consumer - meat",
               "Consumer body weight (kg)", "Dietary intake (g/day)")
pB <- mk_panel(df_B, r_B, "B", "Consumer - other offal",
               "Consumer body weight (kg)", "Dietary intake (g/day)")
pC <- mk_panel(df_C, r_C, "C", "Cattle - FR_feed",
               "Cattle body weight (kg)",
               expression(paste("FR"["feed"], " (kg DM/day)")))

combined <- arrangeGrob(pA, pB, pC, ncol = 3)

ggsave(file.path(out_dir, "Figure3_combined.png"),
       combined, width = 12, height = 4, dpi = 300)
ggsave(file.path(out_dir, "Figure3_combined.tiff"),
       combined, width = 12, height = 4, dpi = 300, compression = "lzw")

cat("Saved:\n  ", file.path(out_dir, "Figure3_combined.png"), "\n  ",
    file.path(out_dir, "Figure3_combined.tiff"), "\n")
