#this qmd is purely to send to a cluster to perform the simulations

library(tidyverse)
library(ggplot2)
library(MASS)
library(purrr) #for looping
library(tictoc) #for checking runtime
library(LaplacesDemon)

library(future)
library(furrr) #for parallelization

# Here are some specialized packages for ease of life to extract visualizations from our simulations
library(gt)
library(gtsummary)
library(broom)

# Here are the packages we need for specific methods:
library(quantreg) #for logistic quantile regression
library(missMethods) #general missing methods
library(mice) #for PMM

cores <- 9 #for number of simulation settings: vectorize based on that
future::plan(multisession, workers = 1) #change to cores when this becomes implemented

############################################################################

num_sim <- 2 #will be 5000

#OUR SIMULATION SETTINGS
sample_sizes <- c(200, 500, 1000) # three settings
prop_missing <- c(0.30, 0.45, 0.60) # three settings


weak_effects <- c(logit(0.1), log(1.1), log(0.7), log(0.85))
strong_effects <- c(logit(0.1), log(1.5), log(0.7), log(0.85))

effect_lookup <- list(
  weak = weak_effects,
  strong = strong_effects
)

effect_strengths <- c("weak", "strong")

grid_design_MCAR <- expand_grid(
  sample_size = sample_sizes,
  prop_missing = prop_missing,
  missingness_type = "MCAR",
  effect_strength = effect_strengths
)

grid_design_MAR <- expand_grid(
  sample_size = sample_sizes,
  prop_missing = prop_missing,
  missingness_type = "MAR",
  effect_strength = effect_strengths
)

grid_design <- bindrows(grid_design_MCAR, grid_design_MAR)

############################################################################ All the functions here (blah)




############################################################################









############################################################################ DEFINE RUNNING EVERYTHING HERE

#furrr across all possible settings provided by design_grid
#before adding this, make sure you remember what the output of the last function is

