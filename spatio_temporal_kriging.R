# CM52068 Spatio-Temporal Analytics
# Coursework: Spatio-Temporal Kriging
# Author: Mo Altrabsheh
# Email:  ma2968@bath.ac.uk
# ============================================================

# Create output directory for figures
if (!dir.exists("figures")) dir.create("figures")

# ============================================================
# TASK 1: Loading and Inspecting the Data
# ============================================================

# (a) Load packages and data
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(patchwork)

noaa_df <- load("data/NOAA_df_1990.rda")
noaa_df <- get(noaa_df)

# Filter for Tmax, July 1993
july1993 <- noaa_df |>
  filter(year == 1993, month == 7, proc == "Tmax")

# (b) Basic summary statistics
cat("Number of spatial locations (stations):", length(unique(july1993$id)), "\n")
cat("Number of time steps (days):           ", length(unique(july1993$date)), "\n")
cat("Number of observations:                ", nrow(july1993), "\n")
cat("Number of missing Tmax observations:   ", sum(is.na(july1993$z)), "\n")

lonmin <- min(july1993$lon)
lonmax <- max(july1993$lon)
latmin <- min(july1993$lat)
latmax <- max(july1993$lat)
cat("Longitude range:", lonmin, "to", lonmax, "\n")
cat("Latitude range: ", latmin, "to", latmax, "\n")

# (c) Visualisations

# Histogram of Tmax values
ggplot() +
  geom_histogram(
    data  = july1993,
    aes(x = z),
    bins  = 25,
    fill  = "#4A90C4",
    color = "white",
    alpha = 0.85
  ) +
  labs(x = "Tmax (°F)", y = "Frequency") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 11)
  )

ggsave("figures/tmax_histogram.png", width = 7, height = 4.5, dpi = 300)

# Time series of Tmax for each station (grey) with median overlay (blue)
ggplot(july1993, aes(x = day, y = z)) +
  geom_line(aes(group = id), alpha = 0.15, colour = "grey70", linewidth = 0.3) +
  stat_summary(fun = median, geom = "line", colour = "#4A90C4", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 31, by = 5)) +
  labs(x = "Day (July 1993)", y = "Tmax (°F)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 11)
  )

ggsave("figures/tmax_timeseries.png", width = 7, height = 4.5, dpi = 300)

# Spatial map of station locations coloured by mean Tmax
world   <- ne_countries(scale = "medium", returnclass = "sf")
pad_lon <- 4
pad_lat <- 4

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "white") +
  geom_point(
    data  = july1993,
    aes(x = lon, y = lat, color = z),
    size  = 3,
    alpha = 0.8
  ) +
  scale_color_viridis_c(option = "plasma", name = "Tmax (°F)") +
  coord_sf(
    xlim = c(lonmin - pad_lon, lonmax + pad_lon),
    ylim = c(latmin - pad_lat, latmax + pad_lat)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 10)
  )

ggsave("figures/tmax_station_map.png", width = 9, height = 8, dpi = 300)

# (d) Remove missing observations and reindex as D
july1993 <- july1993[!is.na(july1993$z), ]
cat("Observations after removing NAs (m):", nrow(july1993), "\n")

# ============================================================
# TASK 2: OLS Regression
# ============================================================

ols_model <- lm(z ~ lon + lat + day, data = july1993)
summary(ols_model)

# ============================================================
# TASK 3: Estimating the Covariance Function via MLE
# ============================================================

# OLS residuals
residuals_ols <- resid(ols_model)

# Subsample for computational tractability
set.seed(26)
subset_idx  <- sample(seq_len(nrow(july1993)), 2000)
july1993_sub <- july1993[subset_idx, ]
resid_sub    <- residuals_ols[subset_idx]

# Gaussian spatio-temporal covariance function
# theta = (sigma2, a, b, c, tau2)
cov_func <- function(theta, lon, lat, day) {
  sigma2  <- theta[1]
  a       <- theta[2]
  b       <- theta[3]
  cc      <- theta[4]
  tau2    <- theta[5]

  lon_mat <- outer(lon, lon, "-")
  lat_mat <- outer(lat, lat, "-")
  day_mat <- outer(day, day, "-")

  sigma2 * exp(
    -(lon_mat^2 / a^2) -
      (lat_mat^2 / b^2) -
      (day_mat^2 / cc^2)
  ) + tau2 * diag(length(lon))
}

# Negative log-likelihood via Cholesky decomposition
neg_log_lik <- function(theta, lon, lat, day, residuals) {
  C <- cov_func(theta, lon, lat, day)
  L <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(L)) return(1e10)

  log_det  <- 2 * sum(log(diag(L)))
  Linv_r   <- backsolve(L, residuals, transpose = TRUE)
  0.5 * (length(residuals) * log(2 * pi) + log_det + sum(Linv_r^2))
}

# MLE optimisation
result <- optim(
  par    = c(16, 5, 3, 1, 4),
  fn     = neg_log_lik,
  lon    = july1993_sub$lon,
  lat    = july1993_sub$lat,
  day    = july1993_sub$day,
  residuals = resid_sub,
  method = "L-BFGS-B",
  lower  = rep(0.01, 5)
)

cat("Estimated covariance parameters (sigma2, a, b, c, tau2):\n")
print(round(result$par, 4))
cat("Convergence code (0 = success):", result$convergence, "\n")

theta <- result$par

# ============================================================
# TASK 4: Building the Prediction Grid
# ============================================================

lon_seq <- seq(-100, -80, length.out = 20)
lat_seq <- seq(32,   46,  length.out = 20)
days    <- c(1, 6, 11, 16, 21, 26)

G <- expand.grid(lon = lon_seq, lat = lat_seq, day = days)

cat("Prediction grid summary:\n")
cat("  Spatial locations:", 20 * 20, "\n")
cat("  Time steps:       ", length(days), "\n")
cat("  Total grid points:", nrow(G), "\n")

# ============================================================
# TASK 5: Universal Kriging
# ============================================================

# Covariance matrix of the observed subsample
cov_matrix <- cov_func(theta,
                       july1993_sub$lon,
                       july1993_sub$lat,
                       july1993_sub$day)

# Cholesky-based inverse
L    <- chol(cov_matrix)
Cinv <- chol2inv(L)

X <- cbind(1, july1993_sub$lon, july1993_sub$lat, july1993_sub$day)
Z <- july1993_sub$z

# GLS coefficient estimates
beta_gls  <- solve(t(X) %*% Cinv %*% X) %*% t(X) %*% Cinv %*% Z
resid_gls <- Z - X %*% beta_gls
Cinv_resid <- Cinv %*% resid_gls

# Cross-covariance matrix between observed points and grid (vectorised)
lon_diff2 <- outer(july1993_sub$lon, G$lon, "-")^2 / theta[2]^2
lat_diff2 <- outer(july1993_sub$lat, G$lat, "-")^2 / theta[3]^2
day_diff2 <- outer(july1993_sub$day, G$day, "-")^2 / theta[4]^2
C0 <- theta[1] * exp(-lon_diff2 - lat_diff2 - day_diff2)

X0       <- cbind(1, G$lon, G$lat, G$day)
Cinv_C0  <- Cinv %*% C0
A        <- solve(t(X) %*% Cinv %*% X)
M        <- t(X0) - t(X) %*% Cinv_C0

# Kriging predictions and variances
G$pred <- as.vector(X0 %*% beta_gls) + as.vector(t(C0) %*% Cinv_resid)
G$var  <- as.numeric((theta[1] + theta[5]) - colSums(C0 * Cinv_C0) +
                       colSums(M * (A %*% M)))
G$se   <- sqrt(G$var)

# Plot: kriging predictions
pred_plot <- ggplot(G, aes(x = lon, y = lat, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Tmax (°F)") +
  facet_wrap(~day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

# Plot: kriging standard errors
se_plot <- ggplot(G, aes(x = lon, y = lat, fill = se)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "SE (°F)") +
  facet_wrap(~day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

ggsave("figures/task5_predictions.png",   pred_plot, width = 14, height = 6,  dpi = 300)
ggsave("figures/task5_se.png",            se_plot,   width = 14, height = 6,  dpi = 300)
ggsave("figures/task5_combined.png",      pred_plot / se_plot,
       width = 14, height = 12, dpi = 300)

# ============================================================
# TASK 6: Comparison with OLS (Trend-Only) Predictor
# ============================================================

# OLS predictions at grid points
G$pred_ols <- as.vector(X0 %*% coef(ols_model))

# Side-by-side comparison plot
G_long <- G |>
  pivot_longer(
    cols      = c(pred, pred_ols),
    names_to  = "method",
    values_to = "temperature"
  ) |>
  mutate(method = ifelse(method == "pred", "Kriging", "OLS"))

ggplot(G_long, aes(x = lon, y = lat, fill = temperature)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Tmax (°F)") +
  facet_grid(method ~ day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    strip.text   = element_text(size = 11)
  )

ggsave("figures/task6_comparison.png", width = 14, height = 8, dpi = 300)

# Held-out RMSE evaluation (20% test split)
set.seed(26)
test_idx  <- sample(seq_len(nrow(july1993)), size = floor(0.2 * nrow(july1993)))
test_set  <- july1993[test_idx, ]
X0_test   <- cbind(1, test_set$lon, test_set$lat, test_set$day)

# OLS RMSE
ols_test_preds <- as.vector(X0_test %*% coef(ols_model))
ols_rmse       <- sqrt(mean((test_set$z - ols_test_preds)^2))

# Kriging RMSE
lon_diff2_t <- outer(july1993_sub$lon, test_set$lon, "-")^2 / theta[2]^2
lat_diff2_t <- outer(july1993_sub$lat, test_set$lat, "-")^2 / theta[3]^2
day_diff2_t <- outer(july1993_sub$day, test_set$day, "-")^2 / theta[4]^2
C0_test     <- theta[1] * exp(-lon_diff2_t - lat_diff2_t - day_diff2_t)

krig_test_preds <- as.vector(X0_test %*% beta_gls) +
  as.vector(t(C0_test) %*% Cinv_resid)
krig_rmse <- sqrt(mean((test_set$z - krig_test_preds)^2))

cat("OLS test RMSE:    ", round(ols_rmse,  3), "F\n")
cat("Kriging test RMSE:", round(krig_rmse, 3), "F\n")
cat("Improvement:      ", round((1 - krig_rmse / ols_rmse) * 100, 1), "%\n")