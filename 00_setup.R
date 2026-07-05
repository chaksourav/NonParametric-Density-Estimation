# ==============================================================================
# Non-Parametric Density Estimation
# File    : 00_setup.R
# Authors : Sourav Chakraborty, Saswata Seal,
#           Piwal Abhishek Satish
# Course  : Non Parametric Inference and Sequential Models, ISI
#
# Purpose : Installs / loads packages, defines shared utilities (make_mixture,
#           dlaplace / rlaplace) and builds mixture objects mix1–mix9 that are
#           used across all other scripts.  Source this file first.
# ==============================================================================

rm(list = ls())

# Package Installation & Loading 

required_packages <- c("ggplot2","patchwork","smoothmest")

install_and_load <- function(packages) {
  
  # Identify packages that are not yet installed
  not_installed <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # Install missing packages
  if (length(not_installed) > 0) {
    message("Installing missing packages: ", paste(not_installed, collapse = ", "))
    install.packages(not_installed, dependencies = TRUE)
  } else {
    message("All packages already installed.")
  }
  
  # Load all packages and track failures
  load_status <- sapply(packages, function(pkg) {
    success <- requireNamespace(pkg, quietly = TRUE) && 
      tryCatch({
        library(pkg, character.only = TRUE)
        TRUE
      }, error = function(e) FALSE)
    success
  })
  
  # Report loading results
  loaded     <- names(load_status[load_status  == TRUE])
  failed     <- names(load_status[load_status  == FALSE])
  
  if (length(loaded) > 0)
    message("Successfully loaded : ", paste(loaded, collapse = ", "))
  if (length(failed) > 0)
    warning("Failed to load     : ", paste(failed, collapse = ", "))
  
  invisible(load_status)
}

install_and_load(required_packages)

set.seed(123)

# Laplace (Double Exponential) helpers 
# Defined here so both 05_glivenko_cantelli.R and 06_kde.R can use them.

dlaplace <- function(x, m = 0, s = 1) {
  (1 / (2 * s)) * exp(-abs(x - m) / s)
}

rlaplace <- function(n, m = 0, s = 1) {
  u <- runif(n, -0.5, 0.5)
  m - s * sign(u) * log(1 - 2 * abs(u))
}


# Function to create mixture density and sampler
make_mixture <- function(components) {
  # components is a list of list(weight, rfun, dfun)
  # weights should sum to 1
  
  # Mixture density
  dmix <- function(x) {
    dens <- numeric(length(x))
    for (comp in components) {
      dens <- dens + comp$weight * comp$dfun(x)
    }
    dens
  }
  
  # Mixture sampler
  rmix <- function(n) {
    # sample component indices according to weights
    weights <- sapply(components, function(c) c$weight)
    idx <- sample(seq_along(components), n, replace = TRUE, prob = weights)
    x <- numeric(n)
    for (k in seq_along(components)) {
      nk <- sum(idx == k)
      if (nk > 0) x[idx == k] <- components[[k]]$rfun(nk)
    }
    x
  }
  
  list(dmix = dmix, rmix = rmix)
}

# Mixture distributions 

# 1. Normal Mixture: 0.5*N(-2,1) + 0.5*N(2,1) 
mix1 <- make_mixture(list(
  list(weight = 0.5,
       rfun = function(n) rnorm(n, -2, 1),
       dfun = function(x) dnorm(x, -2, 1)),
  list(weight = 0.5,
       rfun = function(n) rnorm(n,  2, 1),
       dfun = function(x) dnorm(x,  2, 1))
))

# 2. 0.3*N(0,1) + 0.7*N(4,0.5) 
mix2 <- make_mixture(list(
  list(weight = 0.3,
       rfun = function(n) rnorm(n, 0, 1),
       dfun = function(x) dnorm(x, 0, 1)),
  list(weight = 0.7,
       rfun = function(n) rnorm(n, 4, 0.5),
       dfun = function(x) dnorm(x, 4, 0.5))
))

# 3. Three-component Normal: 0.3*N(-3,0.5) + 0.4*N(0,1) + 0.3*N(3,0.5) 
mix3 <- make_mixture(list(
  list(weight = 0.3,
       rfun = function(n) rnorm(n, -3, 0.5),
       dfun = function(x) dnorm(x, -3, 0.5)),
  list(weight = 0.4,
       rfun = function(n) rnorm(n,  0, 1.0),
       dfun = function(x) dnorm(x,  0, 1.0)),
  list(weight = 0.3,
       rfun = function(n) rnorm(n,  3, 0.5),
       dfun = function(x) dnorm(x,  3, 0.5))
))

# 4. Normal + Exponential: 0.6*N(0,1) + 0.4*Exp(1) 
mix4 <- make_mixture(list(
  list(weight = 0.6,
       rfun = function(n) rnorm(n, 0, 1),
       dfun = function(x) dnorm(x, 0, 1)),
  list(weight = 0.4,
       rfun = function(n) rexp(n, 1),
       dfun = function(x) dexp(x, 1))
))

# 5. Normal + Gamma: 0.5*N(5,1) + 0.5*Gamma(2,1) 
mix5 <- make_mixture(list(
  list(weight = 0.5,
       rfun = function(n) rnorm(n, 5, 1),
       dfun = function(x) dnorm(x, 5, 1)),
  list(weight = 0.5,
       rfun = function(n) rgamma(n, shape = 2, rate = 1),
       dfun = function(x) dgamma(x, shape = 2, rate = 1))
))

# 6. Beta + Uniform: 0.7*Beta(2,5) + 0.3*U(0,1) 
mix6 <- make_mixture(list(
  list(weight = 0.7,
       rfun = function(n) rbeta(n, 2, 5),
       dfun = function(x) dbeta(x, 2, 5)),
  list(weight = 0.3,
       rfun = function(n) runif(n, 0, 1),
       dfun = function(x) dunif(x, 0, 1))
))

# 7. Logistic + Cauchy: 0.8*Logistic(0,1) + 0.2*Cauchy(0,1) 
mix7 <- make_mixture(list(
  list(weight = 0.8,
       rfun = function(n) rlogis(n, 0, 1),
       dfun = function(x) dlogis(x, 0, 1)),
  list(weight = 0.2,
       rfun = function(n) rcauchy(n, 0, 1),
       dfun = function(x) dcauchy(x, 0, 1))
))

# 8. Weibull + Gamma: 0.5*Weibull(2,1) + 0.5*Gamma(3,2) 
mix8 <- make_mixture(list(
  list(weight = 0.5,
       rfun = function(n) rweibull(n, shape = 2, scale = 1),
       dfun = function(x) dweibull(x, shape = 2, scale = 1)),
  list(weight = 0.5,
       rfun = function(n) rgamma(n, shape = 3, rate = 2),
       dfun = function(x) dgamma(x, shape = 3, rate = 2))
))

# 9. Normal + Cauchy: 0.5*Normal(0,1) + 0.5*Cauchy(0,1) 
mix9 <- make_mixture(list(
  list(weight = 0.5,
       rfun = function(n) rnorm(n,0,1),
       dfun = function(x) dnorm(x,0,1)),
  list(weight = 0.5,
       rfun = function(n) rcauchy(n,0,1),
       dfun = function(x) dcauchy(x,0,1))
))

