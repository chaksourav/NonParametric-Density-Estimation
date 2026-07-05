# ==============================================================================
# Non-Parametric Density Estimation
# File    : 05_glivenko_cantelli.R
# Authors : Sourav Chakraborty , Saswata Seal ,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Empirical Validation of the Glivenko-Cantelli Theorem for KDE.
#           Tracks V_n = sup_x |f_hat_n(x) - f(x)| as n grows to 50,000.
#           Also demonstrates failure cases:
#             - Uniform continuity violated (Exp, Uniform)
#             - Bandwidth h_n = n^{-1/2} (series diverges)
#             - Kernels of unbounded variation (Sine, Damped Hyper-Oscillating)
# ==============================================================================

source("00_setup.R")

library(ggplot2)

# ==============================================================================
# Glivenko-Cantelli Theorem Verification
#  Uniform Convergence of KDE to True Density f(x)
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# Verification for Single Densities
# ------------------------------------------------------------------------------
UC_GL_verify <- function(n_vector, 
                         kernel = "gaussian", 
                         bw_method = function(n) n^(-1/5), 
                         dist_name = "norm", 
                         x_range = c(-4, 4), 
                         ..., dname) {
  
  r_func <- match.fun(paste0("r", dist_name))
  d_func <- match.fun(paste0("d", dist_name))
  
  grid_points <- seq(x_range[1], x_range[2], length.out = 1000)
  true_vals <- d_func(grid_points, ...)
  
  sup_norms <- sapply(n_vector, function(n) {
    data_samples <- r_func(n, ...)
    h_n <- bw_method(n)
    # Ensure density is evaluated exactly on the grid_points
    kde_fit <- density(data_samples, bw = h_n, kernel = kernel, 
                       n = 1000, from = x_range[1], to = x_range[2])
    max(abs(kde_fit$y - true_vals))
  })
  
  plot_df <- data.frame(n = n_vector, sup = sup_norms)
  
  ggplot(plot_df, aes(x = n, y = sup)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1, alpha = 0.2) +
    geom_line(color = "#2E7D32", linewidth = 1) + 
    labs(title = paste("Glivenko Cantelli Theorem Verification:", kernel, "kernel"),
         subtitle = paste0("Distribution: ", dname, " | Bandwidth hn = n^(-1/5)"),
         x = "Sample Size (n)", 
         y = "sup |fn(x) - f(x)|" )+
    theme_bw() + 
    theme(
      panel.background = element_rect(fill = "#F1F8E9", color = NA),
      plot.background = element_rect(fill = "#F1F8E9", color = NA),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(color = "#1B5E20", face = "bold")
    )
}

# ------------------------------------------------------------------------------
# 2. Verification for Mixture Densities
# ------------------------------------------------------------------------------
UC_GL_two_mixture <- function(n_vector, 
                              p1,                
                              dist1_name = "norm", params1 = list(), 
                              dist2_name = "norm", params2 = list(), 
                              kernel = "gaussian", 
                              bw_method = function(n) n^(-1/5), 
                              x_range = c(-5, 10), dname) {
  
  r_func1 <- match.fun(paste0("r", dist1_name))
  d_func1 <- match.fun(paste0("d", dist1_name))
  r_func2 <- match.fun(paste0("r", dist2_name))
  d_func2 <- match.fun(paste0("d", dist2_name))
  
  true_mixture_density <- function(x) {
    p1 * do.call(d_func1, c(list(x = x), params1)) + 
      (1 - p1) * do.call(d_func2, c(list(x = x), params2))
  }
  
  grid_points <- seq(x_range[1], x_range[2], length.out = 1000)
  true_vals <- true_mixture_density(grid_points)
  
  sup_norms <- sapply(n_vector, function(n) {
    n1 <- rbinom(1, n, p1)
    samples <- c(do.call(r_func1, c(list(n = n1), params1)), 
                 do.call(r_func2, c(list(n = n - n1), params2)))
    
    h_n <- bw_method(n)
    kde_fit <- density(samples, bw = h_n, kernel = kernel, 
                       n = 1000, from = x_range[1], to = x_range[2])
    
    # Using approx to ensure KDE and True values align perfectly on the grid
    kde_interp <- approx(kde_fit$x, kde_fit$y, xout = grid_points)$y
    max(abs(kde_interp - true_vals), na.rm = TRUE)
  })
  
  plot_df <- data.frame(n = n_vector, sup = sup_norms)
  
  ggplot(plot_df, aes(x = n, y = sup)) +
    geom_hline(yintercept = 0, color = "#1B5E20", linewidth = 1.5, alpha = 0.15) +
    geom_line(color = "#1B5E20", linewidth = 1) + 
    labs(title = paste("Glivenko Cantelli Theorem Verification:", kernel, "kernel"),
         subtitle = paste0("Distribution: ", dname, " | Bandwidth hn = n^(-1/5)"),
         x = "Sample Size (n)", 
         y = "sup|fn(x)-f(x)|") +
    theme_bw() + 
    theme(
      panel.background = element_rect(fill = "#E8F5E9", color = NA),
      plot.background = element_rect(fill = "#E8F5E9", color = NA),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(color = "#1B5E20", size = 12, face = "bold")
    )
}

# ------------------------------------------------------------------------------
# 3. Execution & Verification
# ------------------------------------------------------------------------------

n_sequence <- seq(100, 50000, by = 200)

# Single Density: Normal N(2, 1)
UC_GL_verify(n_vector = n_sequence, 
             kernel = "epanechnikov", 
             bw_method = function(n) n^(-1/5), 
             dist_name = "norm", mean = 2, sd = 1, dname = "N(2,1)")

# Single Density: Exponential Exp(1)
UC_GL_verify(n_vector = n_sequence, 
             kernel = "triangular", 
             bw_method = function(n) n^(-1/5), 
             dist_name = "exp", rate = 1, 
             x_range = c(0, 5), dname = "Exp(1)")

# Mixture: Normal + Cauchy
UC_GL_two_mixture(
  n_vector = n_sequence,
  p1 = 0.5, 
  dist1_name = "norm", params1 = list(mean = -2, sd = 1),
  dist2_name = "cauchy", params2 = list(location = 2, scale = 1),
  x_range = c(-8, 8), dname = "0.5*N(-2,1) + 0.5*Cauchy(2,1)"
)

# Setup sample size sequence
n_seq <- seq(10, 50000, by =50 )

# --- Normal N(0,1) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "norm", 
             mean = 0, sd = 1, x_range = c(-4, 4), dname = "Normal N(0,1)")

# --- Exponential Exp(1) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "exp", 
             rate = 1, x_range = c(0, 5), dname = "Exponential Exp(1)")

# --- Double Exponential (Laplace) L(0,1) ---
# 1. Define Laplace Density (Double Exponential)
dlaplace <- function(x, m = 0, s = 1) {
  (1 / (2 * s)) * exp(-abs(x - m) / s)
}

# 2. Define Laplace Random Generator
rlaplace <- function(n, m = 0, s = 1) {
  u <- runif(n, -0.5, 0.5)
  m - s * sign(u) * log(1 - 2 * abs(u))
}
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "laplace", 
             m = 0, s = 1, x_range = c(-5, 5), dname = "Laplace L(0,1)")

# --- Uniform U(-5,5) ---
UC_GL_verify(n_vector = n_seq, kernel = "epanechnikov", dist_name = "unif", 
             min = -5, max = 5, x_range = c(-6, 6), dname = "Uniform U(-5,5)")

# --- Beta(2,3) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "beta", 
             shape1 = 2, shape2 = 3, x_range = c(0, 1), dname = "Beta(2,3)")

# --- Gamma(2,2) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "gamma", 
             shape = 2, rate = 2, x_range = c(0, 6), dname = "Gamma(2,2)")

# --- Weibull(2,1) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "weibull", 
             shape = 2, scale = 1, x_range = c(0, 4), dname = "Weibull(2,1)")

# --- Cauchy(0,1) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "cauchy", 
             location = 0, scale = 1, x_range = c(-10, 10), dname = "Cauchy(0,1)")

# --- Logistic(0,1) ---
UC_GL_verify(n_vector = n_seq, kernel = "gaussian", dist_name = "logis", 
             location = 0, scale = 1, x_range = c(-6, 6), dname = "Logistic(0,1)")

# --- Mixture 1: 0.5*N(-2,1) + 0.5*N(2,1) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.5, 
                  dist1_name = "norm", params1 = list(mean = -2, sd = 1),
                  dist2_name = "norm", params2 = list(mean = 2, sd = 1),
                  x_range = c(-6, 6), dname = "Mixture 1: 0.5*N(-2,1) + 0.5*N(2,1)")

# --- Mixture 2: 0.3*N(0,1) + 0.7*N(4,0.5) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.3, 
                  dist1_name = "norm", params1 = list(mean = 0, sd = 1),
                  dist2_name = "norm", params2 = list(mean = 4, sd = 0.5),
                  x_range = c(-3, 7), dname = "Mixture 2: 0.3*N(0,1) + 0.7*N(4,0.5)")

# --- Mixture 4: 0.6*N(0,1) + 0.4*Exp(1) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.6, 
                  dist1_name = "norm", params1 = list(mean = 0, sd = 1),
                  dist2_name = "exp", params2 = list(rate = 1),
                  x_range = c(-3, 6), dname = "Mixture 4: 0.6*N(0,1) + 0.4*Exp(1)")

# --- Mixture 5: 0.5*N(5,1) + 0.5*Gamma(2,1) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.5, 
                  dist1_name = "norm", params1 = list(mean = 5, sd = 1),
                  dist2_name = "gamma", params2 = list(shape = 2, rate = 1),
                  x_range = c(0, 10), dname = "Mixture 5: 0.5*N(5,1) + 0.5*Gamma(2,1)")

# --- Mixture 7: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.8, 
                  dist1_name = "logis", params1 = list(location = 0, scale = 1),
                  dist2_name = "cauchy", params2 = list(location = 0, scale = 1),
                  x_range = c(-10, 10), dname = "Mixture 7: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1)")

# --- Mixture 9: 0.5*Normal(2,1) + 0.5*Cauchy(0,1) ---
UC_GL_two_mixture(n_vector = n_seq, p1 = 0.5, 
                  dist1_name = "norm", params1 = list(mean = 2, sd = 1),
                  dist2_name = "cauchy", params2 = list(location = 0, scale = 1),
                  x_range = c(-10, 10), dname = "Mixture 9: 0.5*Normal(2,1) + 0.5*Cauchy(0,1)")

# =====================================================================
# 1. THE KERNEL FUNCTIONS (Not of Bounded Variation)
# =====================================================================

# The Oscillating Sine Kernel
sine_osc_kernel <- function(u) {
  ifelse(abs(u) > 0 & abs(u) <= 1, sin(1 / u), 0)
}

#The Infinite Square-Wave Kernel
square_wave_kernel <- function(u) {
  abs_u <- abs(u)
  k <- floor(1 / abs_u)
  ifelse(abs_u > 0 & abs_u <= 1, (-1)^k, 0)
}

#The Damped Hyper-Oscillating Kernel
damped_osc_kernel <- function(u) {
  ifelse(abs(u) > 0 & abs(u) <= 1, sqrt(abs(u)) * cos(1 / u^2), 0)
}


# =====================================================================
# 2. THE SIMULATION FUNCTION
# =====================================================================
UC_GL_verify_custom <- function(n_vector, 
                                dist_name = "norm", 
                                bw_method = function(n) n^(-1/5), 
                                x_range = c(-4, 4), 
                                custom_kernel = NULL,
                                custom_kernel_name = "Custom",
                                ..., dname, h_func) {
  
  r_func <- match.fun(paste0("r", dist_name))
  d_func <- match.fun(paste0("d", dist_name))
  
  grid_points <- seq(x_range[1], x_range[2], length.out = 1000)
  true_vals <- d_func(grid_points, ...)
  
  sup_norms <- sapply(n_vector, function(n) {
    data_samples <- r_func(n, ...)
    h_n <- bw_method(n)
    
    # If a custom kernel is provided, we calculate KDE manually
    if (!is.null(custom_kernel)) {
      kde_y <- sapply(grid_points, function(x) {
        # Standard KDE formula: (1/(n*h)) * sum( K( (x - Xi) / h ) )
        mean(custom_kernel((x - data_samples) / h_n)) / h_n
      })
    } else {
      # Fallback to standard Gaussian if no custom kernel
      kde_fit <- density(data_samples, bw = h_n, kernel = "gaussian", 
                         n = 1000, from = x_range[1], to = x_range[2])
      kde_y <- kde_fit$y
    }
    
    max(abs(kde_y - true_vals))
  })
  
  plot_df <- data.frame(n = n_vector, sup = sup_norms)
  
  # Determine the title based on whether a custom kernel is used
  title_kernel <- ifelse(is.null(custom_kernel), "Gaussian", custom_kernel_name)
  
  ggplot(plot_df, aes(x = n, y = sup)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1, alpha = 0.2) +
    geom_line(color = "#D32F2F", linewidth = 1) + # Changed to Red to indicate failure
    labs(title = paste("Glivenko-Cantelli Theorem Failure:", title_kernel, "kernel"),
         subtitle = paste0("Distribution: ", dname, " | Bandwidth h(n) = ", h_func),
         x = "Sample Size (n)", 
         y = "sup |fn(x) - f(x)|" )+
    theme_bw() + 
    theme(
      panel.background = element_rect(fill = "#FFEBEE", color = NA), # Light red background
      plot.background = element_rect(fill = "#FFEBEE", color = NA),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(color = "#B71C1C", face = "bold") # Dark red title
    )
}

#  kernels of unbounded variation 

n_sequence <- seq(100, 50000, by = 200)

# Oscillating Sine Kernel on Standard Normal N(0,1)
UC_GL_verify_custom(n_vector = n_sequence,
                    dist_name = "norm", mean = 0, sd = 1,
                    x_range = c(-4, 4),
                    custom_kernel = sine_osc_kernel,
                    custom_kernel_name = "Oscillating Sine",
                    dname = "Standard Normal N(0,1)",
                    h_func = "n^(-1/5)")

# Damped Hyper-Oscillating Kernel on Standard Normal N(0,1)
UC_GL_verify_custom(n_vector = n_sequence,
                    dist_name = "norm", mean = 0, sd = 1,
                    x_range = c(-4, 4),
                    custom_kernel = damped_osc_kernel,
                    custom_kernel_name = "Damped Hyper-Oscillating",
                    dname = "Standard Normal N(0,1)",
                    h_func = "n^(-1/5)")
