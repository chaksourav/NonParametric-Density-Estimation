# ==============================================================================
# Non-Parametric Density Estimation
# File    : 06_kernel_density_estimation.R
# Authors : Sourav Chakraborty , Saswata Seal ,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Kernel Density Estimation (KDE) — visual analysis covering:
#           - Standard univariate distributions (Normal, Exp, Beta, Uniform,
#             Gamma, Laplace, Logistic, Cauchy, Weibull)
#           - Mixture distributions (9 types)
#           - Effect of bandwidth h (h=0.1 vs h=0.7)
#           - Effect of kernel choice (Gaussian, Epanechnikov, Rectangular,
#             Triangular) — all produce nearly identical estimates
#           - Beta kernel for boundary-corrected estimation on (0,1)
# ==============================================================================

source("00_setup.R")

#Kernel Density Estimations
#
#
#===================================================================================
set.seed(123)

plot_kde <- function(n, h, kernel = "gaussian", x_values, dist_name = "norm", ..., Name = "name") {
  r_func <- match.fun(paste0("r", dist_name))
  d_func <- match.fun(paste0("d", dist_name))
  
  params <- list(...)
  param_label <- if(length(params) > 0) paste(names(params), params, sep = "=", collapse = ", ") else "Default"
  
  samples <- r_func(n, ...)
  
  true_df <- data.frame(x = x_values, y = d_func(x_values, ...), type = "True Density")
  
  if (kernel == "rectangular_01") {
    kde_y  <- rect_kde_01(samples, x_values, h)
    kde_df <- data.frame(x = x_values, y = kde_y, type = "KDE Estimation")
  } else {
    kde_res <- density(samples, kernel = kernel, bw = h, n = 1024)
    kde_df  <- data.frame(x = kde_res$x, y = kde_res$y, type = "KDE Estimation")
  }
  
  ggplot() +
    geom_area(data = kde_df, aes(x = x, y = y), fill = "#2E7D32", alpha = 0.07) +
    geom_rug(data = data.frame(x = samples), aes(x = x), color = "#1B5E20", alpha = 0.3, size = 0.4) +
    geom_line(data = kde_df, aes(x = x, y = y, color = "KDE Estimation"), size = 1.2) +
    geom_line(data = true_df, aes(x = x, y = y, color = "True Density"), size = 1, linetype = "dashed") +
    labs(
      title = paste("Kernel Density Estimation:", Name, "Distribution"),
      subtitle = paste0("Parameters: | ", param_label, " |\n",
                        "Kernel: ", kernel, " | Bandwidth (h): ", round(as.numeric(h), 4), " | n: ", n),
      x = "Random Variable (X)", y = "Density f(x)", color = "Legend"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("KDE Estimation" = "#2E7D32", "True Density" = "black")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.5))) +
    theme(
      plot.background = element_rect(fill = "#F9FBF9", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 14, color = "#1B5E20"),
      plot.subtitle = element_text(size = 10, color = "grey20", lineheight = 1.2),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey85", size = 0.5),
      axis.title = element_text(face = "bold")
    )
}

plot_mixture_kde <- function(n, h, kernel = "gaussian", x_values, p1, 
                             dist1_name = "norm", params1 = list(), 
                             dist2_name = "norm", params2 = list(),
                             d1_lab = "Dist 1", d2_lab = "Dist 2") {
  
  
  r_func1 <- match.fun(paste0("r", dist1_name)); d_func1 <- match.fun(paste0("d", dist1_name))
  r_func2 <- match.fun(paste0("r", dist2_name)); d_func2 <- match.fun(paste0("d", dist2_name))
  
  
  n1 <- rbinom(1, n, p1)
  samples <- c(do.call(r_func1, c(list(n = n1), params1)), 
               do.call(r_func2, c(list(n = n - n1), params2)))
  
  
  true_y <- p1 * do.call(d_func1, c(list(x = x_values), params1)) + 
    (1 - p1) * do.call(d_func2, c(list(x = x_values), params2))
  
  true_df <- data.frame(x = x_values, y = true_y, type = "True Density")
  
  kde_res <- density(samples, kernel = kernel, bw = h, n = 1024)
  kde_df <- data.frame(x = kde_res$x, y = kde_res$y, type = "KDE Estimation")
  
  format_p <- function(p) if(length(p)==0) "Default" else paste(names(p), p, sep = "=", collapse = ", ")
  label_subtitle <- paste0("Dist1: ", d1_lab, " | Dist2: ", d2_lab,  ")| p1=",p1,"\n",
                           "Kernel: ", kernel, " | h: ", round(as.numeric(h), 4), " | n: ", n)
  
  ggplot() +
    geom_area(data = kde_df, aes(x = x, y = y), fill = "#0D47A1", alpha = 0.07) +
    geom_rug(data = data.frame(x = samples), aes(x = x), color = "#0D47A1", alpha = 0.3) +
    
    geom_line(data = kde_df, aes(x = x, y = y, color = "KDE Estimation"), size = 1.2) +
    geom_line(data = true_df, aes(x = x, y = y, color = "True Density"), size = 1, linetype = "dashed") +
    
    labs(title = "Density Estimation: Mixture Distribution", 
         subtitle = label_subtitle,
         x = "Random Variable (X)", 
         y = "Density f(x)", 
         color = "Legend") +
    
    theme_minimal() +
    
    scale_color_manual(values = c("KDE Estimation" = "#0D47A1", 
                                  "True Density" = "black")) +
    
    guides(color = guide_legend(ncol = 1, override.aes = list(linetype = c("solid", "dashed"), size = 1))) +
    
    theme(
      plot.background = element_rect(fill = "#F8F9FA", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", color = "#0D47A1"),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey85", size = 0.5),
      axis.title = element_text(face = "bold")
    )
}
# Global Grid for Visualization
x_grid <- seq(-8, 8, length.out = 1000)
n_val  <- 500
h_val  <- 0.3# Fixed bandwidth for comparison



# --- Normal N(0,1) ---
plot_kde(n = n_val, h = h_val, kernel = "gaussian", x_values = x_grid, 
         dist_name = "norm", mean = 0, sd = 1, Name = "Normal")

plot_kde(n = n_val, h = h_val, kernel = "epanechnikov", x_values = x_grid, 
         dist_name = "norm", mean = 0, sd = 1, Name = "Normal")

plot_kde(n = n_val, h = 0.1, kernel = "gaussian", x_values = x_grid, 
         dist_name = "norm", mean = 0, sd = 1, Name = "Normal")

plot_kde(n = n_val, h = 0.7, kernel = "gaussian", x_values = x_grid, 
         dist_name = "norm", mean = 0, sd = 1, Name = "Normal")

# --- Exponential Exp(1) ---
plot_kde(n = n_val, h = 0.3, kernel = "epanechnikov", x_values = seq(0, 6, 0.01), 
         dist_name = "exp", rate = 1, Name = "Exponential")

# --- Double Exponential (Laplace) L(0,1) ---

plot_kde(n = n_val, h = 0.3, kernel = "epanechnikov", x_values = x_grid, 
         dist_name = "laplace", m = 0, s = 1, Name = "Double Exponential")

# --- Uniform U(-5,5) ---
plot_kde(n = n_val, h = 0.3, kernel = "gaussian", x_values = seq(-7, 7, 0.01), 
         dist_name = "unif", min = -5, max = 5, Name = "Uniform")



rect_kde_01 <- function(samples, x_values, h) {
  sapply(x_values, function(xi) {
    u <- (xi - samples) / h
    mean(u >= 0 & u <= 1) / h
  })
}

plot_kde(n = n_val, h = 0.3, kernel = "rectangular_01",
         x_values = seq(0, 1, 0.001),
         dist_name = "beta", shape1 = 2, shape2 = 3,
         Name = "Beta")


# --- Gamma(2,2) ---
plot_kde(n = n_val, h = 0.3, kernel = "gaussian", x_values = seq(0, 8, 0.01), 
         dist_name = "gamma", shape = 2, rate = 2, Name = "Gamma")

# --- Weibull(2,1) ---
plot_kde(n = n_val, h = 0.3, kernel = "epanechnikov", x_values = seq(0, 5, 0.01), 
         dist_name = "weibull", shape = 2, scale = 1, Name = "Weibull")

# --- Cauchy(0,1) ---
plot_kde(n = n_val, h = 0.3, kernel = "gaussian", x_values = seq(-10, 10, 0.01), 
         dist_name = "cauchy", location = 0, scale = 1, Name = "Cauchy")

# --- Logistic(0,1) ---
plot_kde(n = n_val, h = 0.3, kernel = "triangular", x_values = x_grid, 
         dist_name = "logis", location = 0, scale = 1, Name = "Logistic")


# --- Mixture 1: 0.3*N(-2,1) + 0.7*N(2,1) ---
plot_mixture_kde(n = n_val, h = 0.3,kernel = "triangular", x_values = x_grid, p1 = 0.3,
                 dist1_name = "norm", params1 = list(mean = -2, sd = 1),
                 dist2_name = "norm", params2 = list(mean = 2, sd = 1),
                 d1_lab = "N(-2,1)", d2_lab = "N(2,1)")

plot_mixture_kde(n = n_val, h = 0.3,kernel = "gaussian", x_values = x_grid, p1 = 0.3,
                 dist1_name = "norm", params1 = list(mean = -2, sd = 1),
                 dist2_name = "norm", params2 = list(mean = 2, sd = 1),
                 d1_lab = "N(-2,1)", d2_lab = "N(2,1)")

plot_mixture_kde(n = n_val, h = 0.3,kernel = "epanechnikov", x_values = x_grid, p1 = 0.3,
                 dist1_name = "norm", params1 = list(mean = -2, sd = 1),
                 dist2_name = "norm", params2 = list(mean = 2, sd = 1),
                 d1_lab = "N(-2,1)", d2_lab = "N(2,1)")

# --- Mixture 2: 0.3*N(0,1) + 0.7*N(4,0.5) ---
plot_mixture_kde(n = n_val, h = 0.3, kernel = "rectangular",x_values = x_grid, p1 = 0.3,
                 dist1_name = "norm", params1 = list(mean = 0, sd = 1),
                 dist2_name = "norm", params2 = list(mean = 4, sd = 0.5),
                 d1_lab = "N(0,1)", d2_lab = "N(4,0.5)")

# --- Mixture 4: 0.6*N(0,1) + 0.4*Exp(1) ---
plot_mixture_kde(n = 500, h = 0.3, kernel = "epanechnikov",x_values = seq(-3, 6, 0.01), p1 = 0.6,
                 dist1_name = "norm", params1 = list(mean = 0, sd = 1),
                 dist2_name = "exp", params2 = list(rate = 1),
                 d1_lab = "Normal", d2_lab = "Exponential")

# --- Mixture 5: 0.5*N(5,1) + 0.5*Gamma(2,1) ---
plot_mixture_kde(n = 500, h = 0.3,kernel = "gaussian", x_values = seq(0, 12, 0.01), p1 = 0.5,
                 dist1_name = "norm", params1 = list(mean = 5, sd = 1),
                 dist2_name = "gamma", params2 = list(shape = 2, rate = 1),
                 d1_lab = "Normal", d2_lab = "Gamma")

# --- Mixture 6: 0.7*Beta(2,5) + 0.3*U(0,1) ---
plot_mixture_kde(n = 500, h = 0.3, kernel = "triangular",x_values = seq(0, 1, 0.001), p1 = 0.7,
                 dist1_name = "beta", params1 = list(shape1 = 2, shape2 = 5),
                 dist2_name = "unif", params2 = list(min = 0, max = 1),
                 d1_lab = "Beta", d2_lab = "Uniform")

# --- Mixture 7: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1) ---
plot_mixture_kde(n = 500, h = 0.3, kernel = "epanechnikov",x_values = seq(-10, 10, 0.1), p1 = 0.8,
                 dist1_name = "logis", params1 = list(location = 0, scale = 1),
                 dist2_name = "cauchy", params2 = list(location = 0, scale = 1),
                 d1_lab = "Logistic", d2_lab = "Cauchy")

# --- Mixture 8: 0.5*Weibull(2,1) + 0.5*Gamma(3,2) ---
plot_mixture_kde(n = 500, h = 0.3, kernel = "triangular",x_values = seq(0, 8, 0.01), p1 = 0.5,
                 dist1_name = "weibull", params1 = list(shape = 2, scale = 1),
                 dist2_name = "gamma", params2 = list(shape = 3, rate = 2),
                 d1_lab = "Weibull", d2_lab = "Gamma")

# --- Mixture 9: 0.5*Normal(2,1) + 0.5*Cauchy(0,1) ---
plot_mixture_kde(n = 500, h = 0.3, kernel = "gaussian",x_values = seq(-10, 10, 0.1), p1 = 0.3,
                 dist1_name = "norm", params1 = list(mean = 2, sd = 1),
                 dist2_name = "cauchy", params2 = list(location = 0, scale = 1),
                 d1_lab = "Normal", d2_lab = "Cauchy")





## only for uniform with beta kernal

#  Beta kernel helper
beta_kde <- function(samples, x_values, h) {
  a <- min(samples) - h
  b <- max(samples) + h
  z <- (samples - a) / (b - a)
  xz <- (x_values - a) / (b - a)
  
  kde_y <- sapply(xz, function(xi) {
    xi <- max(1e-4, min(1 - 1e-4, xi))
    alpha <- xi / h + 1
    beta_par <- (1 - xi) / h + 1
    mean(dbeta(z, shape1 = alpha, shape2 = beta_par))
  })
  
  kde_y / (b - a)
}

#  plot_kde with beta kernel
plot_kde <- function(n, h, kernel = "gaussian", x_values, dist_name = "norm", ..., Name = "name") {
  r_func <- match.fun(paste0("r", dist_name))
  d_func <- match.fun(paste0("d", dist_name))
  
  params <- list(...)
  param_label <- if(length(params) > 0) paste(names(params), params, sep = "=", collapse = ", ") else "Default"
  
  samples <- r_func(n, ...)
  true_df <- data.frame(x = x_values, y = d_func(x_values, ...), type = "True Density")
  
  if (kernel == "beta") {
    kde_y  <- beta_kde(samples, x_values, h)
    kde_df <- data.frame(x = x_values, y = kde_y, type = "KDE Estimation")
  } else {
    kde_res <- density(samples, kernel = kernel, bw = h, n = 1024)
    kde_df  <- data.frame(x = kde_res$x, y = kde_res$y, type = "KDE Estimation")
  }
  
  ggplot() +
    geom_area(data = kde_df, aes(x = x, y = y), fill = "#2E7D32", alpha = 0.07) +
    geom_rug(data = data.frame(x = samples), aes(x = x), color = "#1B5E20", alpha = 0.3, size = 0.4) +
    geom_line(data = kde_df, aes(x = x, y = y, color = "KDE Estimation"), size = 1.2) +
    geom_line(data = true_df, aes(x = x, y = y, color = "True Density"), size = 1, linetype = "dashed") +
    labs(
      title = paste("Kernel Density Estimation:", Name, "Distribution"),
      subtitle = paste0("Parameters: | ", param_label, " |\n",
                        "Kernel: ", kernel, " | Bandwidth (h): ", round(as.numeric(h), 4), " | n: ", n),
      x = "Random Variable (X)", y = "Density f(x)", color = "Legend"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("KDE Estimation" = "#2E7D32", "True Density" = "black")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.5))) +
    theme(
      plot.background  = element_rect(fill = "#F9FBF9", color = NA),
      panel.background = element_rect(fill = "white",   color = NA),
      plot.title       = element_text(face = "bold", size = 14, color = "#1B5E20"),
      plot.subtitle    = element_text(size = 10, color = "grey20", lineheight = 1.2),
      legend.position  = "right",
      legend.background = element_rect(fill = "white", color = "grey85", size = 0.5),
      axis.title       = element_text(face = "bold")
    )
}

n_val <- 500

plot_kde(n = n_val, h = 0.06, kernel = "beta",
         x_values = seq(0, 1, 0.001),
         dist_name = "unif", min = 0, max = 1,
         Name = "Uniform")




plot_mixture_kde <- function(n, h, kernel = "gaussian", x_values, p1, 
                             dist1_name = "norm", params1 = list(), 
                             dist2_name = "norm", params2 = list(),
                             d1_lab = "Dist 1", d2_lab = "Dist 2",
                             xlim = NULL) {          # ── ADD xlim argument
  
  r_func1 <- match.fun(paste0("r", dist1_name)); d_func1 <- match.fun(paste0("d", dist1_name))
  r_func2 <- match.fun(paste0("r", dist2_name)); d_func2 <- match.fun(paste0("d", dist2_name))
  
  n1 <- rbinom(1, n, p1)
  samples <- c(do.call(r_func1, c(list(n = n1), params1)), 
               do.call(r_func2, c(list(n = n - n1), params2)))
  
  true_y <- p1 * do.call(d_func1, c(list(x = x_values), params1)) + 
    (1 - p1) * do.call(d_func2, c(list(x = x_values), params2))
  
  true_df <- data.frame(x = x_values, y = true_y, type = "True Density")
  kde_res  <- density(samples, kernel = kernel, bw = h, n = 1024)
  kde_df   <- data.frame(x = kde_res$x, y = kde_res$y, type = "KDE Estimation")
  
  format_p     <- function(p) if(length(p)==0) "Default" else paste(names(p), p, sep = "=", collapse = ", ")
  label_subtitle <- paste0("Dist1: ", d1_lab, " | Dist2: ", d2_lab, ")| p1=", p1, "\n",
                           "Kernel: ", kernel, " | h: ", round(as.numeric(h), 4), " | n: ", n)
  
  p <- ggplot() +
    geom_area(data = kde_df, aes(x = x, y = y), fill = "#0D47A1", alpha = 0.07) +
    geom_rug(data = data.frame(x = samples), aes(x = x), color = "#0D47A1", alpha = 0.3) +
    geom_line(data = kde_df,  aes(x = x, y = y, color = "KDE Estimation"), size = 1.2) +
    geom_line(data = true_df, aes(x = x, y = y, color = "True Density"),   size = 1, linetype = "dashed") +
    labs(title    = "Density Estimation: Mixture Distribution",
         subtitle = label_subtitle,
         x = "Random Variable (X)", y = "Density f(x)", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("KDE Estimation" = "#0D47A1", "True Density" = "black")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(linetype = c("solid", "dashed"), size = 1))) +
    theme(
      plot.background  = element_rect(fill = "#F8F9FA", color = NA),
      panel.background = element_rect(fill = "white",   color = NA),
      plot.title       = element_text(face = "bold", color = "#0D47A1"),
      legend.position  = "right",
      legend.background = element_rect(fill = "white", color = "grey85", size = 0.5),
      axis.title       = element_text(face = "bold")
    )
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)   # coord_cartesian clips without removing data
  }
  
  p
}


plot_mixture_kde(n = 500, h = 0.3, kernel = "gaussian",
                 x_values = seq(-10, 10, 0.01), p1 = 0.2,
                 dist1_name = "norm",  params1 = list(mean = 0, sd = 1),
                 dist2_name = "cauchy", params2 = list(location = 0, scale = 1),
                 d1_lab = "Normal", d2_lab = "Cauchy",
                 xlim = c(-10, 10))       
