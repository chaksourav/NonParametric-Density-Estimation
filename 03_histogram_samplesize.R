# ==============================================================================
# Non-Parametric Density Estimation
# File    : 03_histogram_samplesize.R
# Authors : Sourav Chakraborty , Saswata Seal ,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Histogram with variable sample sizes and fixed origin Xo and
#           bandwidth h.
#           Fixed: h = 0.3, Xo = 0.
#           Varied: n in {10, 20, 30, 70, 100, 250, 500, 1000}
# ==============================================================================

source("00_setup.R")

#Histogram with variable sample sizes and fixed origin Xo and bandwidth h
## take n vector of size 8
## we can send as input in function a list of arrays of variable size denoting the samples from different distributions
## here data is variable storing vectors of variables size representing the samples of different sizes from
##given distribution

plot_variable_n <- function(n_vector, h, Xo, true_pdf, rfun, dist_name) {
  
  plot_list <- list()
  
  for (i in 1:length(n_vector)) {
    n      <- n_vector[i]
    sample <- rfun(n)                    # generate sample inside function
    
    min_val <- min(sample) - h
    max_val <- max(sample) + h
    breaks  <- seq(Xo + h * floor((min_val - Xo) / h),
                   Xo + h * ceiling((max_val - Xo) / h) + h,
                   by = h)
    
    p <- ggplot(data.frame(x = sample), aes(x = x)) +
      geom_histogram(breaks = breaks,
                     fill = "#4e79a7",
                     color = "white",
                     aes(y = after_stat(density)),
                     alpha = 0.7) +
      stat_function(fun = true_pdf, color = "red", linewidth = 1) +
      labs(title = paste("n =", n),
           x = NULL, y = "Density") +
      theme_minimal(base_size = 10)
    
    plot_list[[i]] <- p
  }
  
  combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]) /
    (plot_list[[5]] | plot_list[[6]] | plot_list[[7]] | plot_list[[8]])
  
  combined + plot_annotation(
    title    = paste("Effect of sample size on Histogram for", dist_name),
    subtitle = paste("Bandwidth h =", h, "| Origin Xo =", Xo)
  )
}
n_vector <- c(10, 20, 30, 70, 100, 250, 500, 1000)
h_val <- 0.3
Xo <-0

# --- Normal N(0,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = dnorm,
                rfun      = function(n) rnorm(n, 0, 1),
                dist_name = "Normal N(0,1)")

# --- Exponential Exp(1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = dexp,
                rfun      = function(n) rexp(n, 1),
                dist_name = "Exponential Exp(1)")

# --- Double Exponential L(0,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = ddoublex,
                rfun      = function(n) rdoublex(n, 0, 1),
                dist_name = "Double Exponential L(0,1)")

# --- Uniform U(-5,5) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = function(x) dunif(x, -5, 5),
                rfun      = function(n) runif(n, -5, 5),
                dist_name = "Uniform U(-5,5)")

# --- Beta(2,3) ---
plot_variable_n(n_vector, h_val, Xo = 0,
                true_pdf  = function(x) dbeta(x, 2, 3),
                rfun      = function(n) rbeta(n, 2, 3),
                dist_name = "Beta(2,3)")

# --- Gamma(2,2) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = function(x) dgamma(x, 2, 2),
                rfun      = function(n) rgamma(n, 2, 2),
                dist_name = "Gamma(2,2)")

# --- Weibull(2,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = function(x) dweibull(x, 2, 1),
                rfun      = function(n) rweibull(n, 2, 1),
                dist_name = "Weibull(shape=2, scale=1)")

# --- Cauchy(0,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = function(x) dcauchy(x, 0, 1),
                rfun      = function(n) rcauchy(n, 0, 1),
                dist_name = "Cauchy(0,1)")

# --- Logistic(0,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = function(x) dlogis(x, 0, 1),
                rfun      = function(n) rlogis(n, 0, 1),
                dist_name = "Logistic(0,1)")

# MIXTURE DISTRIBUTIONS 

# --- Mixture 1: 0.5*N(-2,1) + 0.5*N(2,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix1$dmix,
                rfun      = mix1$rmix,
                dist_name = "Mixture: 0.5*N(-2,1) + 0.5*N(2,1)")

# --- Mixture 2: 0.3*N(0,1) + 0.7*N(4,0.5) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix2$dmix,
                rfun      = mix2$rmix,
                dist_name = "Mixture: 0.3*N(0,1) + 0.7*N(4,0.5)")

# --- Mixture 3: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix3$dmix,
                rfun      = mix3$rmix,
                dist_name = "Mixture: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5)")

# --- Mixture 4: 0.6*N(0,1) + 0.4*Exp(1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix4$dmix,
                rfun      = mix4$rmix,
                dist_name = "Mixture: 0.6*N(0,1) + 0.4*Exp(1)")

# --- Mixture 5: 0.5*N(5,1) + 0.5*Gamma(2,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix5$dmix,
                rfun      = mix5$rmix,
                dist_name = "Mixture: 0.5*N(5,1) + 0.5*Gamma(2,1)")

# --- Mixture 6: 0.7*Beta(2,5) + 0.3*U(0,1) ---
plot_variable_n(n_vector, h = h_val, Xo = 0,
                true_pdf  = mix6$dmix,
                rfun      = mix6$rmix,
                dist_name = "Mixture: 0.7*Beta(2,5) + 0.3*U(0,1)")




# --- Mixture 8: 0.5*Weibull(2,1) + 0.5*Gamma(3,2) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix8$dmix,
                rfun      = mix8$rmix,
                dist_name = "Mixture: 0.5*Weibull(2,1) + 0.5*Gamma(3,2)")



# --- Mixture 7: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1) ---

# only for cauchy
plot_variable_n <- function(n_vector, h, Xo, true_pdf, rfun, dist_name, x_limits = c(-10, 10)) {
  
  plot_list <- list()
  
  for (i in 1:length(n_vector)) {
    n      <- n_vector[i]
    sample <- rfun(n)                    # generate sample inside function
    
    min_val <- min(sample) - h
    max_val <- max(sample) + h
    breaks  <- seq(Xo + h * floor((min_val - Xo) / h),
                   Xo + h * ceiling((max_val - Xo) / h) + h,
                   by = h)
    
    p <- ggplot(data.frame(x = sample), aes(x = x)) +
      geom_histogram(breaks = breaks,
                     fill = "#4e79a7",
                     color = "white",
                     aes(y = after_stat(density)),
                     alpha = 0.7) +
      stat_function(fun = true_pdf, color = "red", linewidth = 1) +
      
      # ADDED: Zoom in on the x-axis without clipping data calculations
      coord_cartesian(xlim = x_limits) + 
      
      labs(title = paste("n =", n),
           x = NULL, y = "Density") +
      theme_minimal(base_size = 10)
    
    plot_list[[i]] <- p
  }
  
  combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]) /
    (plot_list[[5]] | plot_list[[6]] | plot_list[[7]] | plot_list[[8]])
  
  combined + plot_annotation(
    title    = paste("Effect of sample size on Histogram for", dist_name),
    subtitle = paste("Bandwidth h =", h, "| Origin Xo =", Xo)
  )
}

plot_variable_n(n_vector, h = h_val, Xo,
                true_pdf  = mix7$dmix,
                rfun      = mix7$rmix,
                dist_name = "Mixture: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1)")



# --- Mixture 9: 0.5*Normal(0,1) + 0.5*Cauchy(0,1) ---
plot_variable_n(n_vector, h_val, Xo,
                true_pdf  = mix9$dmix,
                rfun      = mix9$rmix,
                dist_name = "Mixture: 0.5*Normal(0,1) + 0.5*Cauchy(0,1)")

