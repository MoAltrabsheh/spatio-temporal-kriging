# CM52068 Spatio-temporal analytics

# AUTHOR: Mo Altrabsheh
# Email: ma2968@bath.ac.uk

####################################################

# TASK 1: loading the data

# (a) loading packages
library(tidyverse)
library(ggplot2)

noaa_df <- load("data/NOAA_df_1990.rda")
noaa_df <- get(noaa_df)

# filtering the data for July 1993 and Tmax
july1993 <- noaa_df %>%
  filter(year == 1993, month == 7, proc == "Tmax")


# reporting basic statistics
cat("The number of spatial locations: ", length(unique(july1993$id)), "\n")
cat("The number of time steps: ", length(unique(july1993$date)), "\n")
cat("The number of observations: ", nrow(july1993), "\n")
cat("The number of missing observations Tmax: ",
    sum(is.na(july1993$z)), "\n")

lonmin <- min(july1993$lon)
lonmax <- max(july1993$lon)
latmin <- min(july1993$lat)
latmax <- max(july1993$lat)

cat("lon range:", lonmin, "to", lonmax, "\n")
cat("lat range:", latmin, "to", latmax, "\n")

# (b) visualising the data

# histogram of Tmax values
ggplot() +
  geom_histogram(data = july1993, aes(x = z),
                 bins = 25,
                 fill = "#4A90C4",
                 color = "white",
                 alpha = 0.85) +
  labs(x = "Tmax (°F)",
       y = "Frequency") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )
ggsave("tmax_histogram.png", width = 7, height = 4.5, dpi = 300)


# time series of Tmax for each station
ggplot(july1993, aes(x = day, y = z)) +
  geom_line(aes(group = id), alpha = 0.15, colour = "grey70", linewidth = 0.3) +
  stat_summary(fun = median, geom = "line", colour = "#4A90C4", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 31, by = 5)) +
  labs(x = "Day (July 1993)",
       y = "Tmax (°F)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )
ggsave("tmax_timeseries.png", width = 7, height = 4.5, dpi = 300)


# plotting the spatial distribution of Tmax on a map
library(sf)
library(rnaturalearth)

# get map
world <- ne_countries(scale = "medium", returnclass = "sf")
pad_lon <- 4
pad_lat <- 4

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "white") +
  geom_point(data = july1993,
             aes(x = lon, y = lat, color = z),
             size = 3, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Tmax (°F)") +
  coord_sf(
    xlim = c(lonmin - pad_lon, lonmax + pad_lon),
    ylim = c(latmin - pad_lat, latmax + pad_lat)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
ggsave("tmax_station_map.png", width = 9, height = 8, dpi = 300)

# (d) missing data

# There are no missing values in the Tmax variable
sum(!is.na(july1993$z)) # 4122
sum(is.na(july1993$z)) # 0



# TASK 2: fitting a regression model by ordinary least square (OLS)

ols_model <- lm(z ~ lon + lat + day, data = july1993)
summary(ols_model)


# TASK 3: estimating covariance function from known observations

# defining the residuals of the OLS model
residuals <- resid(ols_model)


# using a subset of the data for computational efficiency

set.seed(26)
subset_idx <- sample(1:nrow(july1993), 2000)

july1993_sub <- july1993[subset_idx, ]
resid_sub <- residuals[subset_idx]

# defining covariance function
cov_func <- function(theta, lon, lat, day) {
  sigma2 <- theta[1]
  a <- theta[2]
  b <- theta[3]
  cc <- theta[4]
  tau2 <- theta[5]

  lon_mat <- outer(lon, lon, "-")
  lat_mat <- outer(lat, lat, "-")
  day_mat <- outer(day, day, "-")

  cov <- sigma2 * exp(
    -(lon_mat^2 / a^2) -
      (lat_mat^2 / b^2) -
      (day_mat^2 / cc^2)
  ) + tau2 * diag(length(lon))

  return(cov)
}

# defining negative log-likelihood function using Cholesky decomposition
neg_log_lik <- function(theta, lon, lat, day, residuals) {
  C <- cov_func(theta, lon, lat, day)

  L <- tryCatch(chol(C), error = function(e) return(NULL))
  if (is.null(L)) return(1e10)

  log_det <- 2 * sum(log(diag(L)))
  Linv_r <- backsolve(L, residuals, transpose = TRUE)

  nll <- 0.5 * (length(residuals) * log(2 * pi) + log_det + sum(Linv_r^2))
  return(nll)
}

# running the optimisation
result <- optim(
  par = c(16, 5, 3, 1, 4),
  fn = neg_log_lik,
  lon = july1993_sub$lon,
  lat = july1993_sub$lat,
  day = july1993_sub$day,
  residuals = resid_sub,
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, 0.01, 0.01, 0.01)
)

result$par # estimated parameters for theta
result$convergence # should be 0 for successful convergence

# TASK 4: building a prediction grid

lon_seq <- seq(-100, -80, length.out = 20)
lat_seq <- seq(32, 46, length.out = 20)

days <- c(1, 6, 11, 16, 21, 26)

G <- expand.grid(lon = lon_seq, lat = lat_seq, day = days)

nrow(G) # 20 x 20 x 6 = 2400


# TASK 5: performing covariance-based universal kriging

theta <- result$par

# defining the covariance matrix for the observed data
cov_matrix <- cov_func(theta,
                       july1993_sub$lon,
                       july1993_sub$lat,
                       july1993_sub$day)

# using Cholesky decomposition to compute the inverse of the covariance matrix
L <- chol(cov_matrix)
Cinv <- chol2inv(L)

X <- cbind(1, july1993_sub$lon, july1993_sub$lat, july1993_sub$day)
Z <- july1993_sub$z

# computing the GLS estimates of the regression coefficients
beta_gls <- solve(t(X) %*% Cinv %*% X) %*% t(X) %*% Cinv %*% Z
resid_gls <- Z - X %*% beta_gls
Cinv_resid <- Cinv %*% resid_gls

# computing the covariance between the grid points and the observed data
lon_diff2 <- outer(july1993_sub$lon, G$lon, "-")^2 / theta[2]^2
lat_diff2 <- outer(july1993_sub$lat, G$lat, "-")^2 / theta[3]^2
day_diff2 <- outer(july1993_sub$day, G$day, "-")^2 / theta[4]^2
C0 <- theta[1] * exp(-lon_diff2 - lat_diff2 - day_diff2)

X0 <- cbind(1, G$lon, G$lat, G$day)

# computing the kriging predictions and variances
Cinv_C0 <- Cinv %*% C0
A <- solve(t(X) %*% Cinv %*% X)
M <- t(X0) - t(X) %*% Cinv_C0

G$pred <- as.vector(X0 %*% beta_gls) + as.vector(t(C0) %*% Cinv_resid)
G$var  <- as.numeric((theta[1] + theta[5]) - colSums(C0 * Cinv_C0) + colSums(M * (A %*% M)))
G$se   <- sqrt(G$var)


# visualising the predictions and standard errors
library(patchwork)
ggplot(G, aes(x = lon, y = lat, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~day) +
  labs(
    title = "Kriging Predictions of Tmax (July 1993)",
    x = "Longitude",
    y = "Latitude",
    fill = "Tmax"
  ) +
  theme_minimal()
  
ggsave("predictions_map.png", width = 12, height = 8, dpi = 300)

# visualising the standard errors
ggplot(G, aes(x = lon, y = lat, fill = se)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(~day) +
  labs(
    title = "Kriging Standard Errors (July 1993)",
    x = "Longitude",
    y = "Latitude",
    fill = "SE"
  ) +
  theme_minimal()
ggsave("se_map.png", width = 12, height = 8, dpi = 300)

# visualising the predictions
pred_plot <- ggplot(G, aes(x = lon, y = lat, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Tmax (°F)") +
  facet_wrap(~day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

se_plot <- ggplot(G, aes(x = lon, y = lat, fill = se)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "SE (°F)") +
  facet_wrap(~day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )


combined <- pred_plot | se_plot  # stacked vertically
ggsave("task5_combined.png", width = 14, height = 12, dpi = 300)

# TASK 6: comparining with OLS predictions

# OLS predictions at grid points
G$pred_ols <- cbind(1, G$lon, G$lat, G$day) %*% coef(ols_model)

# side by side predictions
G_long <- G %>%
  pivot_longer(cols = c(pred, pred_ols), names_to = "method", values_to = "temperature") %>%
  mutate(method = ifelse(method == "pred", "Kriging", "OLS"))

ggplot(G_long, aes(x = lon, y = lat, fill = temperature)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Tmax (°F)") +
  facet_grid(method ~ day, labeller = labeller(day = function(x) paste("Day", x))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11)
  )

ggsave("task6_comparison.png", width = 14, height = 12, dpi = 300)

set.seed(26)
test_idx <- sample(1:nrow(july1993), size = floor(0.2 * nrow(july1993)))
test_set <- july1993[test_idx, ]

# OLS RMSE on test set
ols_test_preds <- cbind(1,
                        test_set$lon,
                        test_set$lat,
                        test_set$day) %*% coef(ols_model)
ols_rmse <- sqrt(mean((test_set$z - ols_test_preds)^2))

# Kriging RMSE on test set
lon_diff2_t <- outer(july1993_sub$lon, test_set$lon, "-")^2 / theta[2]^2
lat_diff2_t <- outer(july1993_sub$lat, test_set$lat, "-")^2 / theta[3]^2
day_diff2_t <- outer(july1993_sub$day, test_set$day, "-")^2 / theta[4]^2

C0_test <- theta[1] * exp(-lon_diff2_t - lat_diff2_t - day_diff2_t)

X0_test <- cbind(1, test_set$lon, test_set$lat, test_set$day)

krig_test_preds <- as.vector(X0_test %*% beta_gls) +
  as.vector(t(C0_test) %*% Cinv_resid)

krig_rmse <- sqrt(mean((test_set$z - krig_test_preds)^2))

cat("OLS test RMSE:    ", round(ols_rmse, 3), "°F\n")
cat("Kriging test RMSE:", round(krig_rmse, 3), "°F\n")
cat("Improvement:      ", round((1 - krig_rmse/ols_rmse) * 100, 1), "%\n")
