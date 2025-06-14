# This is the first thing that should be run for this simulation study
# Here we load all the necessary packages and define all the necessary functions
# This will keep everything clean and organized, and that's what we want, right?

# Here are the packages that we need for R (our) life. Very useful things!
library(tidyverse)
library(ggplot2)
library(MASS)
library(purrr) #for looping
library(tictoc) #for checking runtime

# Here are some specialized packages for ease of life to extract visualizations from our simulations
library(gt)
library(gtsummary)
library(broom)

# Here are the packages we need for specific methods:
library(quantreg) #for logistic quantile regression
library(missMethods) #general missing methods
library(mice) #for PMM

# Specialized package to speed up calculations... may or may not use
library(Rcpp) #future implementation for faster calculations in CLQI

###########################################################################################

# takes in a probability and returns its logit value
logit <- function(prob) {
  value <- log(prob / (1 - prob))
  return(value)
}

# takes in a value and gives its expit value
expit <- function(x) {
  value <- 1 / (1 + exp(-x))
  return(value)
}

# transforms a value to its log quantile
log_quant_transform <- function(value, min, max) {
  if (is.na(value)) return(NA)  # short-circuit if missing
  if (value <= min | value >= max) return(NA)
  return(log((value - min) / (max - value)))
}

# untransforms
inv_log_quant_transform <- function(value, min, max) {
  new_value <- (exp(value)*max + min) / (1+exp(value))
  return(new_value)
}




# this is purely for the data generating mechanism
get_mixture_quantile <- function(p) {
  uniroot(
    function(q) 0.6 * pchisq(q, 5) + 0.4 * pchisq(q, 8) - p, #this is specific to our DGM
    c(0, 1000)
  )$root
}

# There will be more functions in later scripts, but that is to make life easier, not because they're essential

