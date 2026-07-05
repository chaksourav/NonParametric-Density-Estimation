# ==============================================================================
# Non-Parametric Density Estimation
# File    : 01_histogram_origin.R
# Authors : Sourav Chakraborty , Saswata Seal ,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Histogram with variable origin Xo and fixed h and sample size n
#           for a set of distribution functions.
#           Fixed: bandwidth h = 0.5, sample size n = 50.
#           Varied: origin Xo in {0, 0.2, 0.4, 1, -0.15, -0.35, -0.45, -2}
# ==============================================================================

source("00_setup.R")

# Histogram with variable origin Xo and fixed h and sample size n for set of distribution functions
# Let us fix bandwidth =0.5 and number of samples to be  50 for variable origin Xo

#Given below is the code for plot of histograms with variable Xo which takes as input set of random samples
#,bandwidth h , a vector of choices of origin , true probability distribution function and Name of the Distribution Function

set.seed(123)

plot_variable_x0 <- function(data, h, x0_vector, true_pdf, dist_name) {
  
  
  plot_list <- list()
  
  for (i in 1:length(x0_vector)) {
    x0 <- x0_vector[i]
    
    # Calculate breaks starting from x0 that cover the data range
    min_val <- min(data)-h
    max_val <- max(data)+h
    breaks <- seq(x0 + h * floor((min_val - x0) / h), 
                  x0 + h * ceiling((max_val - x0) / h) + h, 
                  by = h)
    
    
    p <- ggplot(data.frame(x = data), aes(x = x)) +
      
      geom_histogram(breaks = breaks, 
                     fill = "green", 
                     color = "white", 
                     aes(y = after_stat(density)), 
                     alpha = 0.7) +
      
      stat_function(fun = true_pdf, color = "black", linewidth = 1) +
      labs(title = paste("Xo =", x0), 
           x = "Bins", y = "Density") +
      theme_minimal(base_size = 10)
    
    plot_list[[i]] <- p
  }
  
  # Combined the plots into a 2x3 grid using patchwork
  combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]) / (plot_list[[5]] | plot_list[[6]] | plot_list[[7]] | plot_list[[8]] )
  
  combined + plot_annotation(
    title = paste("Effect of variable origin on Histogram for", dist_name),
    subtitle = paste("Sample size n =", length(data), "| Bandwidth h =", h)
  )
}



n <- 50
h_val <- 0.5
x0_vals <- c(0,0.2,0.4,1,-0.15,-0.35, -0.45, -2)

set.seed(123)



#Normal Distribution mean 0 and sd 1
plot_variable_x0(rnorm(n, mean = 0, sd = 1),
                 h_val, x0_vals, dnorm, "Normal Distribution N(0,1)")

#Exponential Distribution with mean 1
plot_variable_x0(rexp(n, rate = 1),
                 h_val, x0_vals, dexp, "Exponential Distribution Exp(1)")

#Double exponential with mu = 0 and lambda = 1
plot_variable_x0(rdoublex(n,mu = 0,lambda = 1),
                 h_val, x0_vals, ddoublex, "Double Exponential Distribution L(0,1)")

# Uniform U(-5,5)
plot_variable_x0(runif(n, -5, 5), h_val, x0_vals,
                 function(x) dunif(x, -5, 5),
                 "Uniform U(-5,5)")

# --- Beta(2,3) --- support (0,1)
plot_variable_x0(rbeta(n, 2, 3), h_val, x0_vals,
                 function(x) dbeta(x, 2, 3),
                 "Beta(2,3)")

# --- Gamma(2,2) ---
plot_variable_x0(rgamma(n, shape = 2, rate = 2), h_val, x0_vals,
                 function(x) dgamma(x, shape = 2, rate = 2),
                 "Gamma(2,2)")

# --- Weibull(2,1) --- shape=2, scale=1
plot_variable_x0(rweibull(n, shape = 2, scale = 1), h_val, x0_vals,
                 function(x) dweibull(x, shape = 2, scale = 1),
                 "Weibull(shape=2, scale=1)")

# --- Cauchy(0,1) --- 
plot_variable_x0(rcauchy(n, 0, 1), h_val, x0_vals,
                 function(x) dcauchy(x, 0, 1),
                 "Cauchy(0,1)")

# --- Logistic(0,1) ---
plot_variable_x0(rlogis(n, 0, 1), h_val, x0_vals,
                 function(x) dlogis(x, 0, 1),
                 "Logistic(0,1)")

# Mixture of Distribution

plot_variable_x0(mix1$rmix(n), h_val, x0_vals,
                 mix1$dmix,
                 "Mixture: 0.5*N(-2,1) + 0.5*N(2,1)")

plot_variable_x0(mix2$rmix(n), h_val, x0_vals,
                 mix2$dmix,
                 "Mixture: 0.3*N(0,1) + 0.7*N(4,0.5)")

plot_variable_x0(mix3$rmix(n), h_val, x0_vals,
                 mix3$dmix,
                 "Mixture: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5)")

plot_variable_x0(mix4$rmix(n), h_val, x0_vals,
                 mix4$dmix,
                 "Mixture: 0.6*N(0,1) + 0.4*Exp(1)")

plot_variable_x0(mix5$rmix(n), h_val, x0_vals,
                 mix5$dmix,
                 "Mixture: 0.5*N(5,1) + 0.5*Gamma(2,1)")

plot_variable_x0(mix6$rmix(n),h_val,
                 x0_vals,
                 mix6$dmix,
                 "Mixture: 0.7*Beta(2,5) + 0.3*U(0,1)")

plot_variable_x0(mix7$rmix(n), h_val,
                 x0_vals,
                 mix7$dmix,
                 "Mixture: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1)")

plot_variable_x0(mix8$rmix(n), h_val, x0_vals,
                 mix8$dmix,
                 "Mixture: 0.5*Weibull(2,1) + 0.5*Gamma(3,2)")

plot_variable_x0(mix9$rmix(n), h_val, x0_vals,
                 mix9$dmix,
                 "Mixture: 0.5*Normal(0,1) + 0.5*Cauchy(0,1)")
