# ==============================================================================
# Non-Parametric Density Estimation
# File    : 02_histogram_bandwidth.R
# Authors : Sourav Chakraborty , Saswata Seal ,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Histogram with variable bandwidth h and fixed origin Xo and
#           sample size n.
#           Fixed: origin Xo = 0, sample size n = 50.
#           Varied: h in seq(0.05, 1, length.out = 8)
# ==============================================================================

source("00_setup.R")

#Histogram with variable bandwidth h and fixed origin Xo and sample size n
# take h_vector of size 8

plot_variable_h <- function(data, h_vector, Xo, true_pdf, dist_name) {
  
  plot_list <- list()
  
  for (i in 1:length(h_vector)) {
    h <- h_vector[i]
    
    # Calculate breaks starting from Xo that cover the data range for different h
    min_val <- min(data)
    max_val <- max(data)
    breaks <- seq(Xo + h * floor((min_val - Xo) / h), 
                  Xo + h * ceiling((max_val - Xo) / h) + h, 
                  by = h)
    
    p <- ggplot(data.frame(x = data), aes(x = x)) +
      
      geom_histogram(breaks = breaks, 
                     fill = "steelblue", 
                     color = "white", 
                     aes(y = after_stat(density)), 
                     alpha = 0.7) +
      
      stat_function(fun = true_pdf, color = "red", linewidth = 1) +
      labs(title = paste("Bandwidth =", h), 
           x = NULL, y = "Density") +
      theme_minimal(base_size = 10)
    
    plot_list[[i]] <- p
  }
  
  
  combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]) / (plot_list[[5]] | plot_list[[6]] | plot_list[[7]] | plot_list[[8]] )
  
  combined + plot_annotation(
    title = paste("Effect of Variable Bandwidth on Histogram for", dist_name),
    subtitle = paste("Sample size n =", length(data), "| origin =", Xo)
  )
}

x0_vals = 0
h_val = c(0.05, 0.19, 0.34, 0.48, 0.62, 0.76, 0.91, 1.00)
n = 50

set.seed(123)

#Normal Distribution mean 0 and sd 1
plot_variable_h(rnorm(n, mean = 0, sd = 1),
                 h_val, x0_vals, dnorm, "Normal Distribution N(0,1)")

#Exponential Distribution with mean 1
plot_variable_h(rexp(n, rate = 1),
                 h_val, x0_vals, dexp, "Exponential Distribution Exp(1)")

#Double exponential with mu = 0 and lambda = 1
plot_variable_h(rdoublex(n,mu = 0,lambda = 1),
                 h_val, x0_vals, ddoublex, "Double Exponential Distribution L(0,1)")

# Uniform U(-5,5)
plot_variable_h(runif(n, -5, 5), h_val, x0_vals,
                 function(x) dunif(x, -5, 5),
                 "Uniform U(-5,5)")

# --- Beta(2,3) --- support (0,1)
plot_variable_h(rbeta(n, 2, 3), h_val, x0_vals,
                 function(x) dbeta(x, 2, 3),
                 "Beta(2,3)")

# --- Gamma(2,2) ---
plot_variable_h(rgamma(n, shape = 2, rate = 2), h_val, x0_vals,
                 function(x) dgamma(x, shape = 2, rate = 2),
                 "Gamma(2,2)")

# --- Weibull(2,1) --- shape=2, scale=1
plot_variable_h(rweibull(n, shape = 2, scale = 1), h_val, x0_vals,
                 function(x) dweibull(x, shape = 2, scale = 1),
                 "Weibull(shape=2, scale=1)")

# --- Cauchy(0,1) --- 
plot_variable_h(rcauchy(n, 0, 1), h_val, x0_vals,
                 function(x) dcauchy(x, 0, 1),
                 "Cauchy(0,1)")

# --- Logistic(0,1) ---
plot_variable_h(rlogis(n, 0, 1), h_val, x0_vals,
                 function(x) dlogis(x, 0, 1),
                 "Logistic(0,1)")

# Mixture of Distribution

# 1. Normal Mixture: 0.5*N(-2,1) + 0.5*N(2,1)
plot_variable_h(mix1$rmix(n), h_val, x0_vals,
                 mix1$dmix,
                 "Mixture: 0.5*N(-2,1) + 0.5*N(2,1)")

# 2. 0.3*N(0,1) + 0.7*N(4,0.5) 
plot_variable_h(mix2$rmix(n), h_val, x0_vals,
                 mix2$dmix,
                 "Mixture: 0.3*N(0,1) + 0.7*N(4,0.5)")

# 3. Three-component Normal: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5)
plot_variable_h(mix3$rmix(n), h_val, x0_vals,
                 mix3$dmix,
                 "Mixture: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5)")

# 4. Normal + Exponential: 0.6*N(0,1) + 0.4*Exp(1) 
plot_variable_h(mix4$rmix(n), h_val, x0_vals,
                 mix4$dmix,
                 "Mixture: 0.6*N(0,1) + 0.4*Exp(1)")

# 5. Normal + Gamma: 0.5*N(5,1) + 0.5*Gamma(2,1)
plot_variable_h(mix5$rmix(n), h_val, x0_vals,
                 mix5$dmix,
                 "Mixture: 0.5*N(5,1) + 0.5*Gamma(2,1)")

# 6. Beta + Uniform: 0.7*Beta(2,5) + 0.3*U(0,1) 
plot_variable_h(mix6$rmix(n),h_val,
                 x0_vals,
                 mix6$dmix,
                 "Mixture: 0.7*Beta(2,5) + 0.3*U(0,1)")

#  7. Logistic + Cauchy: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1) 
plot_variable_h(mix7$rmix(n), h_val,
                 x0_vals,
                 mix7$dmix,
                 "Mixture: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1)")

#  8. Weibull + Gamma: 0.5*Weibull(2,1) + 0.5*Gamma(3,2) 
plot_variable_h(mix8$rmix(n), h_val, x0_vals,
                 mix8$dmix,
                 "Mixture: 0.5*Weibull(2,1) + 0.5*Gamma(3,2)")


# 9. Normal + Cauchy: 0.5*Normal(0,1) + 0.5*Cauchy(0,1)
plot_variable_h(mix9$rmix(n), h_val, x0_vals,
                 mix9$dmix,
                 "Mixture: 0.5*Normal(0,1) + 0.5*Cauchy(0,1)")


