# ==============================================================================
# Non-Parametric Density Estimation
# File    : 04_histogram_variablebins.R
# Authors : Sourav Chakraborty , Saswata Seal,
#           Piwal Abhishek Satish 
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Histogram with variable bin widths (adaptive / random grids).
#           Three binning schemes compared for Normal, Cauchy and Exponential:
#             Scheme 1 (h_breaks) : structured adaptive grid
#             Scheme 2 (x)        : random breaks from U(-10,10)
#             Scheme 3 (x1)       : random breaks from U(0,5)
# ==============================================================================

source("00_setup.R")

## variable binwidth
plot_variable_bins <- function(n, breaks_vec, rfun, true_pdf, dist_name, color_fill = "steelblue") {
  
  # Generate the data using the provided random function
  data <- rfun(n)
  
  # Create the plot
  p <- ggplot(data.frame(x = data), aes(x = x)) +
    # Use 'breaks' for unequal bin widths and 'density' for the y-axis
    geom_histogram(breaks = breaks_vec, 
                   fill = color_fill, 
                   color = "black", 
                   aes(y = after_stat(density)), 
                   alpha = 0.8) +
    # Overlay the true PDF
    # We evaluate over 1000 points within the range of our breaks
    stat_function(fun = true_pdf, color = "red", linewidth = 1, 
                  xlim = c(min(breaks_vec), max(breaks_vec)), n = 1000) +
    # Lock the viewport to the specific range of the breaks
    coord_cartesian(xlim = c(min(breaks_vec), max(breaks_vec))) +
    labs(title = paste(dist_name, "with Variable Bin Widths"),
         subtitle = paste("Breaks:", paste(round(breaks_vec, 2), collapse = ", ")),
         x = "x", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  print(p)
}


# --- Example 1: Normal Distribution ---
# Defining bins that are narrow at the center (0) and wide at tails
set.seed(123)

x=runif(15,-10,10)
h_breaks <- c(-3,-1.5,-1,-0.5,-0.2,0,.2,0.5,1,1.5,3)
x1=runif(15,0,5)
plot_variable_bins(n = 50, 
                   breaks_vec = h_breaks, 
                   rfun = function(n) rnorm(n, 0, 1), 
                   true_pdf = function(x) dnorm(x, 0, 1), 
                   dist_name = "N(0,1)", 
                   color_fill = "springgreen3")
plot_variable_bins(n = 50, 
                   breaks_vec = x, 
                   rfun = function(n) rnorm(n, 0, 1), 
                   true_pdf = function(x) dnorm(x, 0, 1), 
                   dist_name = "N(0,1)", 
                   color_fill = "springgreen3")
plot_variable_bins(n = 50, 
                   breaks_vec = x1, 
                   rfun = function(n) rnorm(n, 0, 1), 
                   true_pdf = function(x) dnorm(x, 0, 1), 
                   dist_name = "N(0,1)", 
                   color_fill = "springgreen3")

# --- Example 2: Cauchy Distribution ---
# Cauchy needs wide breaks in the tails to handle extreme values


plot_variable_bins(n = 50, 
                   breaks_vec = h_breaks, 
                   rfun = function(n) rcauchy(n, 0, 1), 
                   true_pdf = function(x) dcauchy(x, 0, 1), 
                   dist_name = "Cauchy(0,1)", 
                   color_fill = "cyan3")

plot_variable_bins(n = 50, 
                   breaks_vec = x, 
                   rfun = function(n) rcauchy(n, 0, 1), 
                   true_pdf = function(x) dcauchy(x, 0, 1), 
                   dist_name = "Cauchy(0,1)", 
                   color_fill = "cyan3")

plot_variable_bins(n = 50, 
                   breaks_vec = x1, 
                   rfun = function(n) rcauchy(n, 0, 1), 
                   true_pdf = function(x) dcauchy(x, 0, 1), 
                   dist_name = "Cauchy(0,1)", 
                   color_fill = "cyan3")

# --- Example 3: Exponential Distribution ---
# Bins starting at 0 and getting wider

plot_variable_bins(n = 50, 
                   breaks_vec = h_breaks, 
                   rfun = function(n) rexp(n, 1), 
                   true_pdf = function(x) dexp(x, 1), 
                   dist_name = "Exp(1)", 
                   color_fill = "blue")

plot_variable_bins(n = 50, 
                   breaks_vec = x, 
                   rfun = function(n) rexp(n, 1), 
                   true_pdf = function(x) dexp(x, 1), 
                   dist_name = "Exp(1)", 
                   color_fill = "blue")

plot_variable_bins(n = 50, 
                   breaks_vec = x1, 
                   rfun = function(n) rexp(n, 1), 
                   true_pdf = function(x) dexp(x, 1), 
                   dist_name = "Exp(1)", 
                   color_fill = "blue")



