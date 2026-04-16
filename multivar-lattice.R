# =============================================================================
# COMPLETE WORKING DIGITAL TWIN SIMULATION SCRIPT
# =============================================================================

rm(list = ls())
gc()
set.seed(42)

# =============================================================================
# PACKAGES
# =============================================================================

for (pkg in c("ggplot2", "viridis")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# PARAMETERS
# =============================================================================

N_x <- 100
N_y <- 100
N_cells <- N_x * N_y
Delta <- 0.01
T_max <- 3
n_channels <- 3
sigma_kernel <- 1.2

# Crop forces [N, P, K] - as numeric vectors
crop_forces <- list(
  corn = c(-0.60, -0.20, -0.20),
  soybean = c(+0.20, -0.10, -0.10),
  wheat = c(-0.20, -0.40, -0.10)
)

# =============================================================================
# GAUSSIAN SMOOTHING FUNCTION
# =============================================================================

gaussian_smooth_matrix <- function(mat, sigma) {
  n_x <- nrow(mat)
  n_y <- ncol(mat)
  
  kernel_size <- max(1, ceiling(3 * sigma))
  if (kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
  
  k_half <- (kernel_size - 1) / 2
  k_seq <- (-k_half):k_half
  kernel_1d <- dnorm(k_seq, mean = 0, sd = sigma)
  kernel_1d <- kernel_1d / sum(kernel_1d)
  kernel_2d <- outer(kernel_1d, kernel_1d)
  
  result <- matrix(NA, n_x, n_y)
  
  for (i in 1:n_x) {
    for (j in 1:n_y) {
      total <- 0
      w_sum <- 0
      
      for (ki in 1:kernel_size) {
        for (kj in 1:kernel_size) {
          # Map to matrix indices with reflection at boundaries
          ii <- i + (ki - k_half - 1)
          jj <- j + (kj - k_half - 1)
          
          # Reflect at boundaries
          if (ii < 1) ii <- 2 - ii
          if (ii > n_x) ii <- 2 * n_x - ii
          if (jj < 1) jj <- 2 - jj
          if (jj > n_y) jj <- 2 * n_y - jj
          
          # Clamp
          ii <- max(1, min(n_x, ii))
          jj <- max(1, min(n_y, jj))
          
          weight <- kernel_2d[ki, kj]
          total <- total + weight * mat[ii, jj]
          w_sum <- w_sum + weight
        }
      }
      
      result[i, j] <- total / w_sum
    }
  }
  
  return(result)
}

# =============================================================================
# CREATE BUFFERING CAPACITY MAP
# =============================================================================

create_buffering_map <- function(N_x, N_y, sigma = 1.2) {
  
  x_coords <- seq(0, 1, length.out = N_x)
  y_coords <- seq(0, 1, length.out = N_y)
  
  # Create quadrant-based alpha
  alpha_raw <- matrix(NA, nrow = N_x, ncol = N_y)
  
  for (i in 1:N_x) {
    for (j in 1:N_y) {
      x <- x_coords[i]
      y <- y_coords[j]
      
      if (x < 0.5 && y > 0.5) {
        alpha_raw[i, j] <- 1.0   # Sandy
      } else if (x > 0.5 && y < 0.5) {
        alpha_raw[i, j] <- 0.2   # Clay
      } else {
        alpha_raw[i, j] <- 0.5   # Loam
      }
    }
  }
  
  # Smooth
  alpha_smooth <- gaussian_smooth_matrix(alpha_raw, sigma)
  alpha_smooth <- pmax(pmin(alpha_smooth, 1.0), 0.2)
  
  return(list(
    raw = as.vector(alpha_raw),
    smoothed = as.vector(alpha_smooth),
    matrix_raw = alpha_raw,
    matrix_smoothed = alpha_smooth,
    x = x_coords,
    y = y_coords,
    n_x = N_x,
    n_y = N_y
  ))
}

# =============================================================================
# CALCULATE MORAN'S I
# =============================================================================

calc_morans_I <- function(z, n_x, n_y) {
  
  z <- as.vector(z)
  n <- length(z)
  z_bar <- mean(z)
  z_c <- z - z_bar
  
  numerator <- 0
  weight_sum <- 0
  
  for (i in 1:n_x) {
    for (j in 1:n_y) {
      idx <- (j - 1) * n_x + i
      
      # Rook neighbors
      neighbors <- c()
      if (j > 1) neighbors <- c(neighbors, idx - n_x)
      if (j < n_y) neighbors <- c(neighbors, idx + n_x)
      if (i > 1) neighbors <- c(neighbors, idx - 1)
      if (i < n_x) neighbors <- c(neighbors, idx + 1)
      
      for (n_idx in neighbors) {
        numerator <- numerator + z_c[idx] * z_c[n_idx]
        weight_sum <- weight_sum + 1
      }
    }
  }
  
  denominator <- sum(z_c^2)
  
  if (denominator > 0) {
    I <- (n / weight_sum) * (numerator / denominator)
  } else {
    I <- 0
  }
  
  return(list(I = I, p_value = NA))
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

run_simulation <- function(alpha_map, rotation, T_max = 3, verbose = TRUE) {
  
  if (verbose) cat("Starting simulation:", paste(rotation, collapse = " -> "), "\n")
  
  n_x <- alpha_map$n_x
  n_y <- alpha_map$n_y
  alpha <- alpha_map$smoothed
  
  # State array: S[x, y, time, channel]
  # time: 1 = initial, 2 = year 1, etc.
  S <- array(1.0, dim = c(n_x, n_y, T_max + 1, n_channels))
  D <- matrix(NA, nrow = n_x * n_y, ncol = T_max)
  
  for (t in 1:T_max) {
    
    if (verbose) cat("  Year", t, ":", rotation[t], "\n")
    
    force <- crop_forces[[rotation[t]]]
    prev_t_idx <- t
    curr_t_idx <- t + 1
    
    # Apply crop forces
    for (c in 1:n_channels) {
      for (i in 1:n_x) {
        for (j in 1:n_y) {
          idx <- (j - 1) * n_x + i
          S[i, j, curr_t_idx, c] <- S[i, j, prev_t_idx, c] * (1 + force[c] * alpha[idx])
        }
      }
    }
    
    # Smooth each channel
    for (c in 1:n_channels) {
      S[, , curr_t_idx, c] <- gaussian_smooth_matrix(S[, , curr_t_idx, c], sigma_kernel)
    }
    
    # Calculate stress
    for (i in 1:n_x) {
      for (j in 1:n_y) {
        idx <- (j - 1) * n_x + i
        current_state <- S[i, j, curr_t_idx, ]
        initial_state <- S[i, j, 1, ]
        D[idx, t] <- sqrt(sum((current_state - initial_state)^2))
      }
    }
    
    if (verbose) cat("    Mean stress:", round(mean(D[, t]), 4), "\n")
  }
  
  # Calculate final stress weights
  final_t_idx <- T_max + 1
  
  # Convert to vectors
  N_final <- as.vector(S[, , final_t_idx, 1])
  P_final <- as.vector(S[, , final_t_idx, 2])
  K_final <- as.vector(S[, , final_t_idx, 3])
  
  # Variance of squared deviations
  var_N <- var((N_final - 1)^2)
  var_P <- var((P_final - 1)^2)
  var_K <- var((K_final - 1)^2)
  
  total_var <- var_N + var_P + var_K
  
  if (total_var > 0) {
    w_N <- var_N / total_var
    w_P <- var_P / total_var
    w_K <- var_K / total_var
  } else {
    w_N <- w_P <- w_K <- 1/3
  }
  
  # Moran's I
  moran_N <- calc_morans_I(N_final, n_x, n_y)
  moran_alpha <- calc_morans_I(alpha, n_x, n_y)
  moran_D <- calc_morans_I(D[, T_max], n_x, n_y)
  
  # Statistics
  final_stress <- D[, T_max]
  mean_stress <- mean(final_stress)
  max_stress <- max(final_stress)
  stress_cv <- sd(final_stress) / mean_stress
  
  # Variance ratio
  var_D <- var(final_stress)
  var_alpha <- var(alpha)
  R_ratio <- var_D / var_alpha
  
  if (verbose) {
    cat("\n  === Final Results ===\n")
    cat("  Mean stress:", round(mean_stress, 3), "\n")
    cat("  Max stress:", round(max_stress, 3), "\n")
    cat("  Moran's I:", round(moran_N$I, 3), "\n")
    cat("  R ratio:", round(R_ratio, 3), "\n")
    cat("  Weights: N =", round(w_N, 3), 
        "P =", round(w_P, 3), 
        "K =", round(w_K, 3), "\n\n")
  }
  
  return(list(
    S = S,
    D = D,
    stress_weights = c(N = w_N, P = w_P, K = w_K),
    moran_I = moran_N$I,
    moran_I_alpha = moran_alpha$I,
    mean_stress = mean_stress,
    max_stress = max_stress,
    stress_cv = stress_cv,
    R_ratio = R_ratio,
    alpha_map = alpha_map,
    rotation = rotation
  ))
}

# =============================================================================
# RUN SIMULATIONS
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("RUNNING DIGITAL TWIN SIMULATIONS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Step 1: Create buffering map
cat("Step 1: Creating buffering capacity map...\n")
alpha_map <- create_buffering_map(N_x, N_y, sigma_kernel)
cat("  Done! Alpha range:", round(range(alpha_map$smoothed), 3), "\n\n")

# Step 2: Run baseline rotation
cat("Step 2: Running BASELINE ROTATION (Corn -> Soybean -> Wheat)...\n")
baseline_result <- run_simulation(
  alpha_map, 
  rotation = c("corn", "soybean", "wheat"), 
  T_max = T_max,
  verbose = TRUE
)

# Step 3: Run continuous corn
cat("Step 3: Running CONTINUOUS CORN...\n")
continuous_result <- run_simulation(
  alpha_map,
  rotation = c("corn", "corn", "corn"),
  T_max = T_max,
  verbose = TRUE
)

# =============================================================================
# PRINT COMPARISON TABLE
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("SIMULATION RESULTS COMPARISON\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat(sprintf("%-25s %12s %15s\n", "Metric", "Baseline", "Continuous Corn"))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat(sprintf("%-25s %12.3f %15.3f\n", "Mean Stress D", 
            baseline_result$mean_stress, continuous_result$mean_stress))
cat(sprintf("%-25s %12.3f %15.3f\n", "Max Stress", 
            baseline_result$max_stress, continuous_result$max_stress))
cat(sprintf("%-25s %12.3f %15.3f\n", "Stress CV", 
            baseline_result$stress_cv, continuous_result$stress_cv))
cat(sprintf("%-25s %12.3f %15.3f\n", "Moran's I", 
            baseline_result$moran_I, continuous_result$moran_I))
cat(sprintf("%-25s %12.3f %15.3f\n", "R ratio (Lemma 2)", 
            baseline_result$R_ratio, continuous_result$R_ratio))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat("Stress Weights:\n")
cat(sprintf("  %-23s %11.1f%% %14.1f%%\n", "Nitrogen (N)", 
            baseline_result$stress_weights["N"] * 100, 
            continuous_result$stress_weights["N"] * 100))
cat(sprintf("  %-23s %11.1f%% %14.1f%%\n", "Phosphorus (P)", 
            baseline_result$stress_weights["P"] * 100, 
            continuous_result$stress_weights["P"] * 100))
cat(sprintf("  %-23s %11.1f%% %14.1f%%\n", "Potassium (K)", 
            baseline_result$stress_weights["K"] * 100, 
            continuous_result$stress_weights["K"] * 100))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat(sprintf("Dominant stressor        %12s %15s\n",
            ifelse(baseline_result$stress_weights["N"] > 0.5, "N-dominated", "N ≈ P"),
            "N-dominated"))
cat(paste(rep("=", 60), collapse = ""), "\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

if (!dir.exists("results")) dir.create("results")

saveRDS(baseline_result, "results/baseline_rotation.rds")
saveRDS(continuous_result, "results/continuous_corn.rds")
saveRDS(alpha_map, "results/alpha_map.rds")

cat("\nResults saved to results/ directory\n")
