##Q2
library(mice)
library(boot)

coverages_stochastic <- numeric(100)
coverages_bootstrap <- numeric(100)

# set the parameters
m <- 20
n_boot <- 100 
true_beta1 <- 3
seed <- 1

for (k in 1:100) {
  data <- dataex2[,,k]
  imputed_data <- mice(data, m=m, seed=seed, print=FALSE)
  stochastic_coverage <- bootstrap_coverage <- 0
  for (i in 1:m) {
    completed_data <- complete(imputed_data, action=i)
    fit <- lm(Y ~ X, data=completed_data)
    conf_int_stochastic <- confint(fit, level=0.95)["X", ]
    stochastic_coverage <- stochastic_coverage +
      as.numeric(conf_int_stochastic[1] <= true_beta1 && 
                   true_beta1 <= conf_int_stochastic[2])
    # Bootstrap Method
    bootstrap_results <- boot(completed_data, statistic=function(data, indices) {
      fit <- lm(Y ~ X, data=data[indices, ])
      coef(fit)[2]
    }, R=n_boot)
    # check if the confidence interval have
    if (!is.null(bootstrap_results$t0)) {
      boot_conf_int <- boot.ci(bootstrap_results, type="perc", conf=0.95)
      if ("percent" %in% names(boot_conf_int)) {
        lower_ci <- boot_conf_int$percent[4]
        upper_ci <- boot_conf_int$percent[5]
        bootstrap_coverage <- bootstrap_coverage +
          as.numeric(lower_ci <= true_beta1 && true_beta1 <= upper_ci)
      }
    }
  }
  
  # calculate the probability for each dataset
  coverages_stochastic[k] <- stochastic_coverage / m
  coverages_bootstrap[k] <- bootstrap_coverage / m
}

# calculate the mean probability of all dataset
empirical_coverage_stochastic <- mean(coverages_stochastic)
empirical_coverage_bootstrap <- mean(coverages_bootstrap)

cat("Stochastic Regression Imputation method coverage:", empirical_coverage_stochastic, "\n")
cat("Bootstrap method coverage:", empirical_coverage_bootstrap, "\n")



##Q3
load('dataex3.Rdata')
X <- dataex3$X
R <- dataex3$R
# Assume the value of sigma is 1.5, D is 4
sigma2 <- 1.5
D <- 4

# Define the likelihood function
log_likelihood <- function(mu, X, R, sigma2, D) {
  phi <- dnorm(X, mean = mu, sd = sqrt(sigma2))
  Phi <- pnorm(D, mean = mu, sd = sqrt(sigma2))
  ll <- ifelse(R == 1, log(phi), log(Phi))
  return(sum(ll))
}

neg_log_likelihood <- function(mu) {
  -log_likelihood(mu, X, R, sigma2, D)
}

# Assume the initial value is the mean of X
mu_initial <- mean(X)
# Use optim to find the MLE of mu
result <- optim(par = mu_initial, fn = neg_log_likelihood, method = "L-BFGS-B")
# MLE of mu
mu_mle <- result$par
mu_mle


##Q5
load('dataex5.Rdata')
# Load the data
data <- dataex5

# Logistic function, to make sure the value of p is between 0 and 1
logistic <- function(z) {
  return(1 / (1 + exp(-z)))
}

# Log-likelihood function for observed data
logLikelihood <- function(beta, data) {
  X <- data$X
  Y <- data$Y
  p <- logistic((beta[1] + beta[2] * X)/(1+(beta[1] + beta[2] * X)))
  logL <- sum(Y * log(p) + (1 - Y) * log(1 - p), na.rm = TRUE)
  return(logL)
}

# EM algorithm
em_algorithm <- function(data, max.iter = 1000, tol = 1e-8) {
  # Assume the initial value of beta is 0 and 0
  beta <- c(0, 0)
  # Repeat until convergence
  for (i in 1:max.iter) {
    result <- optim(beta, logLikelihood, data = data, control = list(fnscale = -1))
    new_beta <- result$par
    if (max(abs(new_beta - beta)) < tol) {
      break
    }
    beta <- new_beta
  }
  return(beta)
}

# Run the EM algorithm
beta_estimates <- em_algorithm(data)
beta_estimates