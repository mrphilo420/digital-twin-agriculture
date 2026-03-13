# ==============================================================================
# MEMORY‑OPTIMIZED DIGITAL TWIN OF AGRICULTURAL SOIL
# Multivariate Lattice Model with 3D Volumetrics, Anisotropic Transport,
# Michaelis‑Menten Kinetics, Economic Translation, and EnKF Data Assimilation
# Optimized for large grids (e.g., 1000×1000) with reduced RAM footprint.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PACKAGE SETUP
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "reshape2", "MASS", "xtable", "viridis", "imager", "spdep")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(reshape2)
library(MASS)
library(xtable)
library(viridis)
library(imager)
library(spdep)      # for efficient Moran's I (sparse matrices)

# Optional: use float for half memory (experimental)
# library(float)
# use_float <- FALSE   # set to TRUE to convert arrays to float

# ------------------------------------------------------------------------------
# 1. OPTIMIZED PHYSICS: FFT‑BASED ANISOTROPIC SMOOTHING (unchanged)
# ------------------------------------------------------------------------------
apply_anisotropic_smoothing_fast <- function(mat, sigma_x, sigma_y, global_theta = 0) {
  if (max(sigma_x, sigma_y) <= 0) return(mat)
  radius <- ceiling(3 * max(sigma_x, sigma_y))
  x_seq <- seq(-radius, radius); y_seq <- seq(-radius, radius)
  grid <- expand.grid(dx = x_seq, dy = y_seq)
  cos_t <- cos(global_theta); sin_t <- sin(global_theta)
  x_rot <- grid$dx * cos_t + grid$dy * sin_t
  y_rot <- -grid$dx * sin_t + grid$dy * cos_t
  kernel_vals <- exp(-(x_rot^2) / (2 * sigma_x^2) - (y_rot^2) / (2 * sigma_y^2))
  kernel_mat <- matrix(kernel_vals, nrow = length(x_seq), ncol = length(y_seq))
  
  img <- as.cimg(mat)
  k_img <- as.cimg(kernel_mat)
  num <- imager::convolve(img, k_img, dirichlet = TRUE)
  ones_img <- as.cimg(matrix(1, nrow = nrow(mat), ncol = ncol(mat)))
  den <- imager::convolve(ones_img, k_img, dirichlet = TRUE)
  out <- as.matrix(num / den)
  return(out)
}

# ------------------------------------------------------------------------------
# 2. MEMORY‑SAFE ENSEMBLE KALMAN FILTER (unchanged, but reduce ensemble size later)
# ------------------------------------------------------------------------------
ensemble_kalman_update_fast <- function(ensemble, obs, H_operator_func, R_diag, inflation = 1.05) {
  n_ens <- ncol(ensemble)
  ens_mean <- rowMeans(ensemble)
  A <- (ensemble - ens_mean) * sqrt(inflation)
  HA <- H_operator_func(A)
  R_inv_HA <- HA / R_diag
  D <- t(HA) %*% R_inv_HA + diag(n_ens - 1, n_ens)
  D_inv <- solve(D)
  obs_pert <- matrix(obs, nrow = length(obs), ncol = n_ens) +
    matrix(rnorm(length(obs) * n_ens, 0, sqrt(R_diag)), ncol = n_ens)
  Innovations <- obs_pert - H_operator_func(ensemble)
  update_step <- A %*% (D_inv %*% (t(R_inv_HA) %*% Innovations))
  posterior_ensemble <- ensemble + update_step
  return(rowMeans(posterior_ensemble))
}

# ------------------------------------------------------------------------------
# 3. MAIN SIMULATION DRIVER (3D DIGITAL TWIN) – MEMORY‑OPTIMIZED
#    Stores only current and next time slices (2 time steps).
# ------------------------------------------------------------------------------
run_digital_twin_3D_opt <- function(rotation_seq, sigma_x, sigma_y, theta,
                                    alpha_map_3D, forces_dict, km_vals, leach_rates,
                                    root_dist, dt = 1.0) {
  X <- dim(alpha_map_3D)[1]
  Y <- dim(alpha_map_3D)[2]
  Z <- dim(alpha_map_3D)[3]
  T_steps <- length(rotation_seq) + 1
  C <- 3
  
  # Allocate only two time slices: current (t) and next (t+1)
  # Use a list: current = array(1.0, dim = c(X, Y, Z, C))
  current <- array(1.0, dim = c(X, Y, Z, C))
  next_state <- array(0, dim = c(X, Y, Z, C))
  
  for (yr in 1:length(rotation_seq)) {
    crop <- rotation_seq[yr]
    force_base <- forces_dict[[crop]]
    rho <- root_dist[[crop]]
    
    # Copy current to next as starting point? Actually we build next from scratch.
    # We'll compute next_state and then swap.
    for (c_idx in 1:C) {
      for (z in 1:Z) {
        current_S <- current[, , z, c_idx]   # from previous year
        f_max <- force_base[c_idx] * rho[z]
        km <- km_vals[c_idx]
        
        # Reaction
        if (f_max >= 0) {
          new_layer <- current_S * (1 + f_max * dt * alpha_map_3D[, , z])
        } else {
          uptake_rate <- -f_max * alpha_map_3D[, , z]
          new_layer <- current_S * exp(-uptake_rate * dt / (km + current_S))
        }
        
        # Horizontal smoothing
        new_layer <- apply_anisotropic_smoothing_fast(new_layer, sigma_x, sigma_y, theta)
        
        # Store temporarily (will be used for leaching)
        # We'll accumulate into next_state[,,z,c_idx] later after leaching across layers
        # But leaching needs all layers of the same nutrient together.
        # So we need to keep a temporary 3D array for this nutrient across depths.
        # Let's create a temporary 3D array for this nutrient.
        if (c_idx == 1 && z == 1) {
          # Initialize temp array for this nutrient
          temp_nutrient <- array(0, dim = c(X, Y, Z))
        }
        temp_nutrient[, , z] <- new_layer
      } # end depth loop
      
      # Now temp_nutrient holds the post‑smoothing values for nutrient c_idx.
      # Apply leaching cascade
      for (z in 1:(Z-1)) {
        downward_flux <- temp_nutrient[, , z] * leach_rates[c_idx] * dt
        temp_nutrient[, , z] <- temp_nutrient[, , z] - downward_flux
        temp_nutrient[, , z+1] <- temp_nutrient[, , z+1] + downward_flux
      }
      # Groundwater loss from deepest layer
      temp_nutrient[, , Z] <- temp_nutrient[, , Z] * (1 - leach_rates[c_idx] * dt)
      
      # Store final result in next_state
      for (z in 1:Z) {
        next_state[, , z, c_idx] <- temp_nutrient[, , z]
      }
    } # end nutrient loop
    
    # Swap: current <- next_state, and reset next_state to zeros
    current <- next_state
    next_state <- array(0, dim = c(X, Y, Z, C))
  }
  
  # Return final state (after all years)
  return(list(Final_State = current))
}

# ------------------------------------------------------------------------------
# 4. PARAMETER DEFINITIONS (same as before, but can scale up)
# ------------------------------------------------------------------------------
set.seed(42)
X <- 200; Y <- 200; Z <- 3

# 4.1 Soil buffering map
alpha_2d <- matrix(runif(X*Y, 0.3, 1.0), X, Y)
alpha_2d[1:(X/2), 1:(Y/2)] <- 0.9
alpha_2d[(X/2+1):X, (Y/2+1):Y] <- 0.5
alpha_2d <- apply_anisotropic_smoothing_fast(alpha_2d, 2.0, 2.0, 0)

alpha_3D <- array(dim = c(X, Y, Z))
alpha_3D[,,1] <- alpha_2d
alpha_3D[,,2] <- pmax(0.1, alpha_2d - 0.2)
alpha_3D[,,3] <- pmax(0.1, alpha_2d - 0.4)

# 4.2 Crop forces
forces <- list(
  Corn    = c(-0.6, -0.2, -0.2),
  Soybean = c( 0.2, -0.1, -0.1),
  Wheat   = c(-0.2, -0.4, -0.1)
)

# 4.3 Half‑saturation constants
km_vals <- c(N = 0.3, P = 0.2, K = 0.25)

# 4.4 Leaching rates
leach_rates <- c(N = 0.15, P = 0.01, K = 0.05)

# 4.5 Root distributions
root_dist <- list(
  Corn    = c(0.6, 0.3, 0.1),
  Soybean = c(0.7, 0.2, 0.1),
  Wheat   = c(0.4, 0.4, 0.2)
)

# ------------------------------------------------------------------------------
# 5. SCENARIO EXECUTION (using optimized driver)
# ------------------------------------------------------------------------------
cat("\n=== Running Optimized Digital Twin Simulations ===\n")
cat("Grid size:", X, "x", Y, "x", Z, "=", X*Y*Z, "voxels\n")

cat("Scenario 0: Baseline rotation (Corn‑Soybean‑Wheat)...\n")
res_base <- run_digital_twin_3D_opt(
  rotation_seq = c("Corn", "Soybean", "Wheat"),
  sigma_x = 2.0, sigma_y = 0.5, theta = 0,
  alpha_map_3D = alpha_3D,
  forces_dict = forces,
  km_vals = km_vals,
  leach_rates = leach_rates,
  root_dist = root_dist
)

cat("Scenario S4: Continuous corn...\n")
res_cont <- run_digital_twin_3D_opt(
  rotation_seq = c("Corn", "Corn", "Corn"),
  sigma_x = 2.0, sigma_y = 0.5, theta = 0,
  alpha_map_3D = alpha_3D,
  forces_dict = forces,
  km_vals = km_vals,
  leach_rates = leach_rates,
  root_dist = root_dist
)

# ------------------------------------------------------------------------------
# 6. AGRONOMIC AND ECONOMIC TRANSLATION
# ------------------------------------------------------------------------------
mitscherlich_baule <- function(S_eff, crop_type, Y_max) {
  k_coefs <- list(
    Corn    = c(N = 2.5, P = 1.8, K = 1.2),
    Soybean = c(N = 1.5, P = 2.2, K = 1.0),
    Wheat   = c(N = 2.0, P = 2.5, K = 1.3)
  )
  k <- k_coefs[[crop_type]]
  yield_frac <- apply(S_eff, 1, function(x) prod(1 - exp(-k * x)))
  return(Y_max * yield_frac)
}

root_corn <- root_dist[["Corn"]]

S_base <- res_base$Final_State
S_cont <- res_cont$Final_State

depth_weighted_avg <- function(S_array, root_wt) {
  eff <- S_array[,,1] * root_wt[1] +
    S_array[,,2] * root_wt[2] +
    S_array[,,3] * root_wt[3]
  return(as.vector(eff))
}

eff_N_base <- depth_weighted_avg(S_base[,,,1], root_corn)
eff_P_base <- depth_weighted_avg(S_base[,,,2], root_corn)
eff_K_base <- depth_weighted_avg(S_base[,,,3], root_corn)

eff_N_cont <- depth_weighted_avg(S_cont[,,,1], root_corn)
eff_P_cont <- depth_weighted_avg(S_cont[,,,2], root_corn)
eff_K_cont <- depth_weighted_avg(S_cont[,,,3], root_corn)

Y_max <- 11.0; price_per_ton <- 180
yield_base <- mitscherlich_baule(cbind(eff_N_base, eff_P_base, eff_K_base), "Corn", Y_max)
yield_cont <- mitscherlich_baule(cbind(eff_N_cont, eff_P_cont, eff_K_cont), "Corn", Y_max)
loss_base <- (Y_max - yield_base) * price_per_ton
loss_cont <- (Y_max - yield_cont) * price_per_ton

cat("\nMean economic loss (baseline):", mean(loss_base, na.rm = TRUE), "USD/ha")
cat("\nMean economic loss (continuous corn):", mean(loss_cont, na.rm = TRUE), "USD/ha\n")

# ------------------------------------------------------------------------------
# 7. DATA ASSIMILATION DEMO (reduced ensemble size)
# ------------------------------------------------------------------------------
cat("\n=== Ensemble Kalman Filter Demonstration ===\n")
true_state <- as.vector(S_base[,,1,1])
n_state <- length(true_state)
n_ens <- 15   # reduced from 30 to save memory

ensemble <- matrix(rnorm(n_state * n_ens, mean = true_state, sd = 0.1),
                   nrow = n_state, ncol = n_ens)
obs_indices <- sample(1:n_state, size = floor(0.1 * n_state))
obs <- true_state[obs_indices] + rnorm(length(obs_indices), 0, 0.05)
R_diag <- 0.05^2

H_operator <- function(ens) ens[obs_indices, , drop = FALSE]

posterior_mean <- ensemble_kalman_update_fast(
  ensemble = ensemble,
  obs = obs,
  H_operator_func = H_operator,
  R_diag = R_diag,
  inflation = 1.05
)

prior_mean <- rowMeans(ensemble)
rmse_prior <- sqrt(mean((prior_mean - true_state)^2))
rmse_post  <- sqrt(mean((posterior_mean - true_state)^2))
cat(sprintf("Prior RMSE: %.4f\n", rmse_prior))
cat(sprintf("Posterior RMSE: %.4f\n", rmse_post))
cat(sprintf("RMSE reduction: %.1f%%\n", (1 - rmse_post/rmse_prior)*100))

# ------------------------------------------------------------------------------
# 8. SUMMARY STATISTICS (with sparse Moran's I)
# ------------------------------------------------------------------------------
cat("\n=== Generating Summary Statistics ===\n")

stress_base <- sqrt((S_base[,,1,1] - 1)^2 + (S_base[,,2,1] - 1)^2 + (S_base[,,3,1] - 1)^2)
stress_cont <- sqrt((S_cont[,,1,1] - 1)^2 + (S_cont[,,2,1] - 1)^2 + (S_cont[,,3,1] - 1)^2)

# Efficient Moran's I using spdep (sparse weights)
compute_moran_sparse <- function(mat) {
  vals <- as.vector(mat)
  coords <- expand.grid(x = 1:nrow(mat), y = 1:ncol(mat))
  nb <- knn2nb(knearneigh(coords, k = 8))   # 8 nearest neighbours
  listw <- nb2listw(nb)
  moran <- moran.test(vals, listw, randomisation = FALSE)
  return(c(moran$estimate[1], moran$p.value))
}

moran_base <- compute_moran_sparse(stress_base)
moran_cont <- compute_moran_sparse(stress_cont)

summary_df <- data.frame(
  Scenario = c("Baseline (Corn-Soy-Wheat)", "Continuous Corn"),
  Mean_Stress = c(mean(stress_base), mean(stress_cont)),
  Max_Stress = c(max(stress_base), max(stress_cont)),
  CV_Stress = c(sd(stress_base)/mean(stress_base), sd(stress_cont)/mean(stress_cont)),
  Mean_Econ_Loss_USD = c(mean(loss_base), mean(loss_cont)),
  Moran_I = c(moran_base[1], moran_cont[1])
)

print(summary_df, digits = 3)

# LaTeX table
latex_table <- xtable(summary_df,
                      caption = "Summary statistics for baseline and continuous corn scenarios.",
                      label = "tab:summary",
                      digits = 3)
print(latex_table, include.rownames = FALSE, booktabs = TRUE)


# ------------------------------------------------------------------------------
# 9. PUBLICATION FIGURES
# ------------------------------------------------------------------------------
cat("\n=== Generating Figures ===\n")

# Helper to convert matrix to data frame for ggplot
matrix_to_df <- function(mat, name) {
  df <- expand.grid(X = 1:nrow(mat), Y = 1:ncol(mat))
  df$Value <- as.vector(mat)
  df$Metric <- name
  return(df)
}

# --- FIGURE 1: 3D Volumetrics ---
df_top <- matrix_to_df(S_base[,,1,1], "Topsoil N (0-20 cm)")
df_sub <- matrix_to_df(S_base[,,3,1], "Subsoil N (40-60 cm)")
df_vol <- rbind(df_top, df_sub)

p1 <- ggplot(df_vol, aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  facet_wrap(~Metric) +
  scale_fill_viridis_c(option = "mako", name = "N Status") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(title = "3D Volumetric Mining: Topsoil vs. Subsoil Nitrogen")

ggsave("Figure1_3D_Volumetrics.png", plot = p1, width = 10, height = 5, dpi = 300)


# --- FIGURE 2: Economic Loss ---
df_loss <- matrix_to_df(matrix(loss_base, nrow = X, ncol = Y), "Economic Loss (USD/ha)")
p2 <- ggplot(df_loss, aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = "USD / ha") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank()) +
  labs(title = "Economic Translation",
       subtitle = "Yield loss penalty for subsequent Corn crop")

ggsave("Figure2_Economic_Loss.png", plot = p2, width = 6, height = 5, dpi = 300)


# --- FIGURE 3: EnKF RMSE Reduction ---
# Note: rmse_prior and rmse_post are calculated in Section 7 (Data Assimilation Demo)
rmse_data <- data.frame(
  State = factor(c("Prior (Model Only)", "Posterior (Data Assimilated)"), 
                 levels = c("Prior (Model Only)", "Posterior (Data Assimilated)")),
  RMSE = c(rmse_prior, rmse_post)
)

p3 <- ggplot(rmse_data, aes(x = State, y = RMSE, fill = State)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", size = 0.5) +
  scale_fill_manual(values = c("Prior (Model Only)" = "#E07A5F", 
                               "Posterior (Data Assimilated)" = "#3D405B")) +
  geom_text(aes(label = sprintf("%.4f", RMSE)), vjust = -0.5, size = 5, fontface = "bold") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Digital Twin Self-Correction via EnKF",
    subtitle = "Root Mean Square Error (RMSE) reduction after assimilating spatial observations",
    y = "Spatial RMSE (Normalized Concentration)"
  )

ggsave("Figure3_EnKF_RMSE.png", plot = p3, width = 7, height = 5, dpi = 300)

cat("Figures saved: Figure1_3D_Volumetrics.png, Figure2_Economic_Loss.png, Figure3_EnKF_RMSE.png\n")

cat("\n=== Digital Twin Pipeline Complete (Memory Optimized) ===\n")