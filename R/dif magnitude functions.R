# Function to pull item parameters out of H5RESULTS object
#   for continuous dependent variables (pulls intercepts)
#
# Rich Jones
# Nov 25 2024
#
# Example usage:
# pull_itemParameters("mm06.h5","NHW","HISPSPAN","VDMDE5Z","GENMEM") 
#    where
#        "mm05.h5"    is the name of the h5 file 
#        "NHW"        is the label assigned (in the Mplus output) to 
#                     the reference group
#        "HISPAN"     is the label assigned (in the Mplus output) to
#                     the focal group 
#        "VDMDE5Z"    is the name of the CONTINUOUS dependent variable
#        "GENMEM"     is the label assigned to the underlying latent 
#                     variable
#
pull_itemParameters <- function(h5,group1,group2,item,factorname) {
  E <- mplus.print.model.results(h5) 
  a1 <- E |> filter(grepl(item, Statement) & grepl(group1,Group) & grepl("Intercepts",Section)) |> pull(Estimate)
  a2 <- E |> filter(grepl(item, Statement) & grepl(group2,Group) & grepl("Intercepts",Section)) |> pull(Estimate)
  b1 <- E |> filter(grepl(item, Statement) & grepl(group1,Group) & grepl(paste(factorname,"BY"),Statement)) |> pull(Estimate)
  b2 <- E |> filter(grepl(item, Statement) & grepl(group2,Group) & grepl(paste(factorname,"BY"),Statement)) |> pull(Estimate)
  c(a1,b1,a2,b2)
  
  
  pull_muSigma <- function(h5,group2,factorname) {
    E <- mplus.print.model.results(h5) 
    mu <- E |> filter(grepl(factorname, Statement) & grepl(group2,Group) & grepl("Means",Section)) |> pull(Estimate)
    sigma <- E |> filter(grepl(factorname, Statement) & grepl(group2,Group) & grepl("Variances",Section)) |> pull(Estimate)
    c(mu,sigma)
  }
  
  pooled_sd <- function(group_df, groupvar, item) {
    # Load necessary package
    library(dplyr)
    # force varnames and items to uppercase
    # because Mplus is case insensitivive and this is written to
    # work with H5RESULTS
    item <- toupper(item)
    groupvar <- toupper(groupvar)
    DF <- group_df %>% rename_all(toupper)
    DF <- DF %>% select(all_of(c(groupvar, item))) 
    # Ensure the grouping variable is a factor
    DF[[groupvar]] <- as.factor(DF[[groupvar]])
    # Calculate the number of levels in the grouping variable
    levels <- levels(DF[[groupvar]])
    k <- length(levels)
    # Initialize variables to store sums of squared deviations and total sample size
    sum_sq_dev <- 0
    total_n <- 0
    # Loop over each level to calculate the sum of squared deviations and total sample size
    for (level in levels) {
      subset_data <- DF %>% filter(DF[[groupvar]] == level)
      n <- sum(!is.na(subset_data[[item]])) # Count non-missing values
      sd_val <- sd(subset_data[[item]], na.rm = TRUE)
      sum_sq_dev <- sum_sq_dev + (n - 1) * sd_val^2
      total_n <- total_n + n
    }
    # Calculate the pooled standard deviation
    pooled_sd_val <- sqrt(sum_sq_dev / (total_n - k))
    print(paste("Pooled sd over",k,"levels of",groupvar," = ",pooled_sd_val))
    return(pooled_sd_val)
  }
}

computeAreas <- function(params, musigma, sd = NULL, cov_matrix = NULL) {
  # Extract parameters
  a1 <- params[1]
  b1 <- params[2]
  a2 <- params[3]
  b2 <- params[4]
  mu <- musigma[1]
  sigma <- musigma[2]
  
  # Define the squared difference function
  squared_f <- function(x) {
    diff <- (a1 - a2) + (b1 - b2) * x
    phi_x <- dnorm(x, mean = mu, sd = sigma)  # Latent trait density
    (diff^2) * phi_x
  }
  
  # Solve for crossover point
  if (b1 != b2) {
    xc <- (a2 - a1) / (b1 - b2)
  } else {
    xc <- NULL  # Parallel lines; no crossover point
  }
  
  # Compute unsigned area (full integral)
  total_result <- integrate(squared_f, lower = -Inf, upper = Inf)
  unsigned_area <- sqrt(total_result$value)
  
  # Compute signed area if slopes are different
  if (!is.null(xc)) {
    # Above region
    above_result <- integrate(squared_f, lower = xc, upper = Inf)
    D_above <- sqrt(above_result$value)
    
    # Below region
    below_result <- integrate(squared_f, lower = -Inf, upper = xc)
    D_below <- sqrt(below_result$value)
    
    # Signed area
    signed_area <- D_above - D_below
  } else {
    signed_area <- NA
  }
  
  # Compute standardized areas if pooled SD is supplied
  if (!is.null(sd)) {
    std_unsigned_area <- unsigned_area / sd
    std_signed_area <- if (!is.null(xc)) signed_area / sd else NA
  } else {
    std_unsigned_area <- NA
    std_signed_area <- NA
  }
  
  # Compute standard errors if covariance matrix is supplied
  if (!is.null(cov_matrix)) {
    # Gradient for unsigned area
    unsigned_gradient <- numDeriv::grad(function(p) {
      a1 <- p[1]
      b1 <- p[2]
      a2 <- p[3]
      b2 <- p[4]
      integrate(squared_f, lower = -Inf, upper = Inf)$value
    }, params)
    
    unsigned_variance <- t(unsigned_gradient) %*% cov_matrix %*% unsigned_gradient
    unsigned_se <- sqrt(unsigned_variance)
    
    # Gradient for signed area (only if slopes are different)
    if (!is.null(xc)) {
      signed_gradient <- numDeriv::grad(function(p) {
        a1 <- p[1]
        b1 <- p[2]
        a2 <- p[3]
        b2 <- p[4]
        if (b1 != b2) {
          xc <- (a2 - a1) / (b1 - b2)
          above <- integrate(squared_f, lower = xc, upper = Inf)$value
          below <- integrate(squared_f, lower = -Inf, upper = xc)$value
          sqrt(above) - sqrt(below)
        } else {
          NA
        }
      }, params)
      
      signed_variance <- t(signed_gradient) %*% cov_matrix %*% signed_gradient
      signed_se <- sqrt(signed_variance)
    } else {
      signed_se <- NA
    }
    
    # Standardized SEs
    if (!is.null(sd)) {
      std_unsigned_se <- unsigned_se / sd
      std_signed_se <- if (!is.null(xc)) signed_se / sd else NA
    } else {
      std_unsigned_se <- NA
      std_signed_se <- NA
    }
  } else {
    unsigned_se <- NA
    signed_se <- NA
    std_unsigned_se <- NA
    std_signed_se <- NA
  }
  
  # Return all results as a list
  return(list(
    unsigned_area = unsigned_area,
    signed_area = signed_area,
    std_unsigned_area = std_unsigned_area,
    std_signed_area = std_signed_area,
    unsigned_se = unsigned_se,
    signed_se = signed_se,
    std_unsigned_se = std_unsigned_se,
    std_signed_se = std_signed_se
  ))
}

expectedScore <- function(f, alpha, beta) {
  # ordinal response model assumes a0=-Inf, ak+1=+Inf
  alpha <- c(-Inf, alpha, Inf)
  k <- length(alpha)
  score <- 0
  for (j in 2:k) {  # Start from 2 to account for -Inf
    # The pnorm let's you know we are expecting 
    # parameters from the ordered probit regression
    # model. For example Mplus WLSMV theta.
    prob_j <- pnorm(alpha[j] - beta * f) - pnorm(alpha[j - 1] - beta * f)
    score <- score + (j - 2) * prob_j  # (j - 2) to adjust for the starting index
  }
  score <- score / (k - 2)
  return(score)
}

# Function: areaMeasures 
# The area measure is the weighted sum of Cohen's h effect
# size statistic over the range of ability (defined by the
# user, in fis) and weighted according to the distribution
# of ability in the focal group (as implied by the mean and 
# standard deviation of ability in the focal group provided
# by the user in `mu2` and `sigmasquared2`).
#
# Cohen's h effect size statistic is for the difference in
# proportions. It is computed as 
#
#    h = 2*asin(sqrt(p2)) - 2*asin(sqrt(p1))
#
# and is interpreted like Cohen's d statistic where .2, .5 
# and .8 demarcate small, medium and large effects (Cohen, 1988).
# In our framework, we have proportions for specific quadrature
# points on the ability dimension, which are expected proportions 
# of the total possible number of points on the item at that 
# ability level for members of the reference group (score1)
# and the focal group (score2). 
#
# We compute the weighted -- to the focal group ability 
# distribution -- mean h across the quadrature points, and 
# label this "SAh". This is analogous to the signed area 
# between the ICCs. 
#
# We also take the absolute value of h, and compute the weighted
# mean of that. We label this the unsigned area ("UAh").
#
# Both SAh and UAh are returned.
#
#
#  use like this:
#
# areaMeasures(score1, score2, mu2, sigmasquared2, fis)
# $SAh
# [1] 0.1202699
# 
# $UAh
# [1] 0.3825437
# where score1 and score2 come from the ExpectedScore function given
# a1, b1, and a2, b2, respectively, and mu2 is the mean of the latent
# trait in the focal group and sigmasquared2 is the variance of the
# latent trait in the focal group.
#
library(stats)
areaMeasures <- function(score1, score2, mu2, sigmasquared2, fis) {
  # Calculate the weights based on the normal distribution
  weight <- dnorm(fis, mean = mu2, sd = sqrt(sigmasquared2))  # Weights
  
  # calculate difference in scores
  SAp <- sum((score2-score1) * weight) / sum(weight)
  UAp <- sum(abs(score2-score1) * weight) / sum(weight)
  
  # Calculate h
  h <- 2 * (asin(sqrt(score2)) - asin(sqrt(score1)))
  # Calculate the absolute value of h
  absh <- abs(h)
  
  
  # Calculate the weighted mean of h
  SAh <- sum(h * weight) / sum(weight)
  
  # Calculate the weighted mean of abs(h)
  UAh <- sum(absh * weight) / sum(weight)
  
  # Return the calculated values as a list
  return(list(SAh = SAh, UAh = UAh, SAp = SAp, UAp = UAp))
}
#
# Function: plotERF
# This is a simple wrapper for plot that 
# will make an expected item response function (ERF)
# plot for the two groups.
#
#  example use:
#
#  plotERF(fis, score1, score2)
plotERF <- function(fis, score1, score2, 
                    xlab = "Latent Trait (f)", 
                    ylab = "Expected Score/Total Item Score", 
                    legend_labels = c("Focal", "Reference")) {
  # Plot score2 over fis as a simple line plot
  plot(fis, score2, type = "l", col = "red", ylim = c(0, 1), 
       xlab = xlab, ylab = ylab)
  
  # Add score1 to the plot in red
  lines(fis, score1, col = "blue")
  
  # Add a legend to differentiate the two lines
  legend("bottomright", legend = legend_labels, col = c("red", "blue"), lty = 1)
}
# plotERF(fis, score1, score2)

# Running through the three functions
# Define input parameters
# These are results from Mplus or other software
### This example is a toy example where the
### thresholds and loadings are different
### but the expected response functions are
### nearly identical. It's a mind-bender
### see below for something more expected
### a1 <- c(-1, 0, 1)
### b1 <- 0.85
### a2 <- c(-0.5, 0, 0.5)
### b2 <- 0.7
a1 <- c(-0.25, .25, 0.75)
b1 <- 0.6
a2 <- c(-0.5, 0, 0.5)
b2 <- 0.8
mu2 <- -1
sigmasquared2 <- 1
# Set quadrature points, the range of theta you'd like to 
# see the plots and define the ERFs. From -4 to +4 is 
# typical, so you might not want to edit the line below
fis <- seq(-4, 4, by = 0.01)  # Quadrature points
#
# Use the functions
score1 <- sapply(fis, function(f) expectedScore(f, a1, b1))
score2 <- sapply(fis, function(f) expectedScore(f, a2, b2))
areaMeasures(score1, score2, mu2, sigmasquared2, fis)
plotERF(fis, score1, score2)