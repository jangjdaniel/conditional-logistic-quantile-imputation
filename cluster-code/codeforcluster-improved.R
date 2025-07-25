#this qmd is purely to send to a cluster to perform the simulations

library(tidyverse)
library(ggplot2)
library(MASS)
library(purrr) #for looping
library(tictoc) #for checking runtime
library(LaplacesDemon)

library(future)
library(furrr) #for parallelization
library(glue)

# Here are some specialized packages for ease of life to extract visualizations from our simulations
library(gt)
library(gtsummary)
library(broom)

# Here are the packages we need for specific methods:
library(quantreg) #for logistic quantile regression
library(missMethods) #general missing methods
library(mice) #for PMM

#install packages just in case the cluster doesn't have it
#install.packages(c("LaplacesDemon", "gt", "gtsummary", "broom", "quantreg", "missMethods", "mice"))

############################################################################ Set up all simulation settings
num_sim <- 5000

#sample size
sample_sizes <- c(200, 500, 1000) # three settings

#MCAR
prop_missing <- c(0.30, 0.45, 0.60) # three settings

#effects
no_effects <- c(logit(0.1), log(1), log(0.7), log(0.85))
weak_effects <- c(logit(0.1), log(1.1), log(0.7), log(0.85))
strong_effects <- c(logit(0.1), log(1.5), log(0.7), log(0.85))

effect_lookup <- list(
  none - no_effects,
  weak = weak_effects,
  strong = strong_effects
)

effect_strengths <- c("none", "weak", "strong")

grid_design_MCAR <- expand_grid(
  sample_size = sample_sizes,
  prop_missing = prop_missing,
  missingness_type = "MCAR",
  effect_strength = effect_strengths
)

grid_design_MAR <- expand_grid(
  sample_size = sample_sizes,
  prop_missing = prop_missing, #look at my documentation for corresponding values
  missingness_type = "MAR",
  effect_strength = effect_strengths
)

grid_design <- bind_rows(grid_design_MCAR, grid_design_MAR)

############################################################################ All fundamental functions here (00)

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

############################################################################ Data Generating Mechanism (01)

# establish the formula relationships for performing regression after imputation
formula_biomarker <- "outcome ~ biomarker + confounder + predictor"
formula_biomarker_MCAR <- "outcome ~ biomarker_MCAR + confounder + predictor"
formula_biomarker_MAR <- "outcome ~ biomarker_MAR + confounder + predictor"

# define the min and max for the biomarker variable 
min <- 0
max <- 10

# define the data generating mechanism
data_generator <- function(sample_size, coefficient_vector) {
  
  # Extract coefficients from the coefficient vector *b*
  b_0 <- coefficient_vector[1]
  b_1 <- coefficient_vector[2]
  b_2 <- coefficient_vector[3]
  b_3 <- coefficient_vector[4]
  
  # Generate values from distributions
  confounder <- rbinom(sample_size, size = 1, prob = 0.4) #prob of success is 0.4
  predictor <- rnorm(sample_size, mean = 0, sd = 1) #this is our continuous predictor
  biomarker <- rtrunc(sample_size, spec = "chisq", #use truncated distribution instead
                      a = 0, b = 20, #our min and max
                      df = 5 + 3*confounder) #if B = 1, X ~ chisq(8). else, X ~ chisq(5)
  
  # Generate the outcome
  outcome <- plogis(b_0 + b_1*biomarker + b_2*confounder + b_3*predictor) 
  outcome_binary <- rbinom(sample_size, size = 1, prob = outcome)
  
  # Into the dataframe
  my_data <- data.frame(
    confounder = confounder,
    predictor = predictor,
    biomarker = biomarker,
    outcome = outcome_binary
  )
  
  # Get our data
  return(my_data)
}

# define the missingness generator

missing_generator <- function(my_data, missing_prop_for_MCAR) {
  
  # Add biomarker variants to induce missingness
  my_data <- my_data |>
    mutate(biomarker_MCAR = biomarker,
           biomarker_MAR = biomarker)
  
  # Add a row index 
  my_data <- my_data |>
    mutate(row_index = row_number())
  
  
  # MCAR
  my_data <- missMethods::delete_MCAR(my_data,  #name of data
                                      missing_prop_for_MCAR,  #proportion of data you want missing
                                      "biomarker_MCAR") #name of variable you want to do this
  
  # MAR: bad fix
  if (missing_prop_for_MCAR == 0.3) {
    alpha <- log(0.2 / (1 - 0.2))        # logit(0.2)
    beta <- log(2)                       # doubles odds when C = 1
    
  } else if (missing_prop_for_MCAR == 0.45) {
    alpha <- log(0.3 / (1 - 0.3))        # logit(0.3) ≈ -0.847
    beta <- log(2.5)                     # increases odds ~2.5× when C = 1
    
  } else if (missing_prop_for_MCAR == 0.6) {
    alpha <- log(0.4 / (1 - 0.4))        # logit(0.4) ≈ -0.405
    beta <- log(3)                       # increases odds ~3× when C = 1
  }
  
  #calculate the missingness based on the probs given
  my_data <- my_data |>
    mutate(logit_p = alpha + beta * confounder,
           prob_missing = 1 / (1 + exp(-logit_p)), #calculate predicted prob of missing
           missing_MAR = rbinom(n(), size = 1, prob = prob_missing), #make indicator
           biomarker_MAR = ifelse(missing_MAR == 1, NA, biomarker_MAR) #done
    ) |>
    dplyr::select(-logit_p, prob_missing, missing_MAR)
  
  return(my_data)
}

data_generating_mechanism <- function(my_sample, beta_coefficients, prop_missing_MCAR,
                                      min_val = 0, max_val) {
  # Generate the data
  my_data <- data_generator(my_sample, beta_coefficients)
  
  # Induce MCAR and MAR missingness
  my_data <- missing_generator(my_data, prop_missing_MCAR)
  
  # Apply log transform to the biomarker data
  # As a caveat, I don't know if this is the best way to go about it
  my_data <- my_data |>
    mutate(biomarker_transformed = sapply(biomarker, log_quant_transform, #3 mutates for organization only
                                          min = min_val, max = max_val)) |> 
    mutate(biomarker_MCAR_transformed = sapply(biomarker_MCAR, log_quant_transform,
                                               min = min_val, max = max_val)) |>
    mutate(biomarker_MAR_transformed = sapply(biomarker_MAR, log_quant_transform,
                                              min = min_val, max = max_val)) 
  
  return(my_data)
}

############################################################################ Baseline Comparisons (CC, CD): 02

complete_case_or_data <- function(my_data, correct_formula) {
  
  #perform logistic regression, get everything into a tibble, and return that
  CCD <-   glm(as.formula(correct_formula),
               data = my_data,
               family = "binomial") |>
    summary() |>
    (\(x) as_tibble(x$coefficients, rownames = "term"))()
  
  return(CCD)
}

############################################################################ Predictive Mean Matching: 02

PMM <- function(my_data_filtered, correct_formula, num_imp) {
  
  # Generate imputed datasets
  capture.output({ #the mice package has weird output that isn't relevant for us. This is largely ignored
    PMM_df <- mice::mice(my_data_filtered, 
                         method = "pmm", 
                         m = num_imp) 
  })
  
  # Run regression!
  PMM_regression_results <- purrr::map_dfr(1:num_imp, ~ { #we repeat this num_imp times
    
    #extract the ith dataset (.x refers to the 1:num_imp)
    data_iteration <- complete(PMM_df, action = .x)
    
    #perform logistic regression, get everything into a tibble, and return that
    glm(as.formula(correct_formula),
        data = data_iteration,
        family = "binomial") |>
      summary() |>
      (\(x) {
        df <- as_tibble(x$coefficients, rownames = "term") 
        df[2, , drop = FALSE] #select only the second row: base R is fast
      })()
  }
  )
  
  #make sure to return it!
  return(PMM_regression_results)
}

############################################################################ Conditional Logistic Quantile Imp: 02

logistic_quantile_imputation <- function(my_data, transformed_imputation_relationship, row_index) {
  
  # Step 1) generate a random uniform value
  u <- runif(1, min = 0, max = 0.99) 
  
  # Step 2) perform quantile regression with the uth quantile
  reg_coeff <- quantreg::rq(as.formula(transformed_imputation_relationship), 
                            data = my_data, 
                            tau = u) #the uth quantile for regression
  
  # save all regression coefficients
  b_intercept <-  as.numeric(reg_coeff$coefficients["(Intercept)"])
  b_outcome <-    as.numeric(reg_coeff$coefficients["outcome"])
  b_confounder <- as.numeric(reg_coeff$coefficients["confounder"])
  b_predictor <-  as.numeric(reg_coeff$coefficients["predictor"])
  
  #predicted imputed value (which is in its transformed state)
  imputation_value_transformed <- b_intercept + (b_outcome * my_data[row_index,]$outcome) + (b_confounder * my_data[row_index,]$confounder) + (b_predictor * my_data[row_index,]$predictor)
  
  #return imputed value in its transformed state
  return(imputation_value_transformed)
}

implement_logistic_quantile_imputation <- function(my_data, transformed_imputation_relationship, var_for_imp,
                                                   min_val = 0, max_val) {
  # Convert var_for_imp to a symbol for tidy evaluation
  var_sym <- sym(var_for_imp)
  imputed_var_name <- paste0(var_for_imp, "_CLQI")
  
  # next, create that variable name to perform our operation  
  my_data <- my_data |>
    mutate(!!sym(imputed_var_name) := NA_real_)
  
  # Perform the algorithm
  for(row_index in 1:nrow(my_data)) { #for all indicies
    if(is.na(my_data[[var_for_imp]][row_index])) { #if the cell in the vector for the variable for imputation is empty
      
      imputed_value <- logistic_quantile_imputation( #implement algorithm
        my_data = my_data, 
        transformed_imputation_relationship = transformed_imputation_relationship,
        row_index = row_index)
      
      # then impute it!
      my_data[[imputed_var_name]][row_index] <- imputed_value
    }
    else { # or just assign the regular value
      my_data[[imputed_var_name]][row_index] <- my_data[[var_for_imp]][row_index]
    }
  }
  
  # Don't forget to UNTRANSFORM this new variable
  imputed_var_name_untransformed <- paste0(sub("_transformed", "", imputed_var_name), "_untransformed")
  
  my_data <- my_data |>
    mutate(!!sym(imputed_var_name_untransformed) := sapply(
      !!sym(imputed_var_name), #same name
      inv_log_quant_transform, #perform the inverse log transform
      min = min_val, max = max_val))
  
  return(my_data)
}

CLQI <- function(my_data, var_for_imp, transformed_imputation_relationship, correct_formula, 
                 num_imp, min_val = 0, max_val) {
  
  CLQI_regression_results <- purrr::map_dfr(1:num_imp, ~ { #we repeat this num_imp times
    
    # copy same missing data for each iteration
    copy_my_data <- my_data
    
    # implement algorithm
    data_iteration <- implement_logistic_quantile_imputation(
      my_data = copy_my_data,
      transformed_imputation_relationship = transformed_imputation_relationship,
      var_for_imp = var_for_imp,
      min_val = min_val,
      max_val = max_val
    )
    
    #print out the results
    glm(as.formula(correct_formula),
        data = data_iteration,
        family = "binomial") |>
      summary() |>
      (\(x) {
        df <- as_tibble(x$coefficients, rownames = "term") 
        df[2, , drop = FALSE] #select only the second row: base R is fast
      })()
  }
  )
  
  return(CLQI_regression_results)
}

############################################################################ Extracting Summary Statistics: 03

rubin_rule_estimate <- function(extracted_statistics) {
  return(mean(extracted_statistics$Estimate))
}

rubin_rule_SE <- function(extracted_statistics) {
  #extract necessary vectors
  estimate_vector <- as.vector(extracted_statistics$Estimate)
  SE_vector <- as.vector(extracted_statistics$`Std. Error`)
  m <- length(estimate_vector) #number of imputations
  
  #necessary calculations
  variance_within_MI <- mean((SE_vector)^2)
  variance_between_MI <- sum((estimate_vector - mean(estimate_vector))^2) / (m - 1)
  
  #actual value
  total_variance <- variance_within_MI + variance_between_MI + (variance_between_MI / m)
  
  return(sqrt(total_variance))
}

coverage <- function(RR_estimate, SE_estimate, orig_sample_size, true_value) {
  t_star <- qt(0.975, df = orig_sample_size - 1) #good enough estimate
  
  coverage_indicator <- ifelse((RR_estimate - t_star*SE_estimate) <= true_value & true_value <= (RR_estimate + t_star*SE_estimate),
                               1, 0) #if true value is covered, it's 1. 
  
  return(coverage_indicator)
}

power <- function(extracted_statistics) {
  
  # extract important information
  estimate_vector <- as.vector(extracted_statistics$Estimate)
  SE_vector <- as.vector(extracted_statistics$`Std. Error`)
  m <- length(estimate_vector) #number of imputations
  
  # df calculation
  variance_within_MI <- mean((SE_vector)^2)
  variance_between_MI <- sum((estimate_vector - mean(estimate_vector))^2) / (m - 1)
  my_frac <- variance_between_MI / (m * variance_within_MI)
  df <- ((m - 1) / (1 + my_frac)^2) 
  
  # wald t statistic: repeat rubin rules again for function cleaniness
  RR_estimate <- mean(estimate_vector)
  SE_estimate <- sqrt(variance_within_MI + variance_between_MI + (variance_between_MI / m))
  wald <- RR_estimate / SE_estimate
  
  # p-value calculation
  p_value <- 2 * (1 - pt(abs(wald), df))
  
  #return 1 if we reject H0 (from DGM, we know H1 to be true)
  return(ifelse(p_value < 0.05, 1, 0))
}

rel_efficiency <- function(simulation_data_A, simulation_data_B) {
  empirical_variance_A <- var(simulation_data_A$Estimate)
  empirical_variance_B <- var(simulation_data_B$Estimate)
  
  # Then return our ratio
  return(empirical_variance_B/empirical_variance_A)
}

MCSE <- function(vec_values_interest, num_sim) {
  return(sd(vec_values_interest) / sqrt(num_sim))
}


############################################################################ Functions for simulation study: 04

complete_performance_one_iter <- function(my_reg_coefficients, orig_sample_size, 
                                          true_effect, biomarker_name) {
  
  #generate t star value for CI evaluation
  t_star <- qt(0.975, df = orig_sample_size - 1)
  
  # calculate all the performance measures!
  sim_run <- my_reg_coefficients |>
    dplyr::filter(term == biomarker_name) |>
    dplyr::mutate(
      bias = Estimate - true_effect,
      rel_bias = (bias / true_effect) * 100,
      lower = Estimate - t_star * `Std. Error`,
      upper = Estimate + t_star * `Std. Error`,
      coverage = ifelse(lower <= true_effect & true_effect <= upper, 1, 0),
      power = ifelse(`Pr(>|z|)` < 0.05, 1, 0)
    ) |>
    dplyr::select(estimate = Estimate, bias, rel_bias, model_SE = `Std. Error`, coverage, power)
  
  return(sim_run)
}

MI_performance_one_iter <- function(my_list_of_reg, var_of_interest, true_value) {
  estimate <- rubin_rule_estimate(my_list_of_reg)
  model_SE <- rubin_rule_SE(my_list_of_reg)
  
  #now the dataframe
  sim_run <- data.frame(
    estimate = estimate,
    bias = estimate - true_value,
    model_SE = model_SE,
    coverage = coverage(estimate, model_SE,
                        orig_sample_size = 1000,
                        true_value = true_value) ,
    power =  power(my_list_of_reg)
  ) 
  
  return(sim_run)
}


full_performance <- function(simulation_results, rel_eff_comparison) {
  
  full_simulation_results <- simulation_results |>
    summarize(
      estimate_a = mean(estimate),
      MCSE_estimate = MCSE(simulation_results$estimate, num_sim = 30),
      
      bias = mean(bias),
      
      empirical_SE = sd(estimate),
      model_SE = mean(model_SE),
      MCSE_model_SE = MCSE(simulation_results$model_SE, num_sim = 30),
      
      coverage = sum(coverage) / n(),
      power = sum(power) / n(),
    ) |>
    mutate(RMSE = sqrt(bias^2 + empirical_SE^2))
  
  return(full_simulation_results)
  
}

############################################################################ Full implementation: 05

repeat_full_simulation_MCAR <- function(num_simulations = 1000, my_sample_size, num_imp,
                                        coefficient_effects, true_effect, prop_missing_MCAR) {
  
  # STEP 0) Define datasets for all four methods
  CC_total_result <- data.frame()
  CD_total_result <- data.frame()
  PMM_total_result <- data.frame()
  CLQI_total_result <- data.frame()
  
  #have an iteration indicator to see if function is still running
  iteration_indicator <- 0
  
  #now the loop
  for(i in 1:num_simulations) {
    
    # STEP 1) Generate data
    new_data <- data_generating_mechanism(
      my_sample = my_sample_size,
      beta_coefficients = coefficient_effects,
      prop_missing_MCAR = prop_missing_MCAR,
      min_val = 0,
      max_val = 20)
    
    # STEP 2) Do CC
    CC <- "outcome ~ biomarker_MCAR + confounder + predictor"
    CC_sim_result <- complete_case_or_data(new_data, CC)
    CC_sim_result <- complete_performance_one_iter(CC_sim_result, orig_sample_size = my_sample_size, 
                                                   true_effect = true_effect, biomarker_name = "biomarker_MCAR")
    
    CC_total_result <- dplyr::bind_rows(CC_total_result, CC_sim_result)
    
    # STEP 2.1) Do CD
    CD <- "outcome ~ biomarker + confounder + predictor"
    CD_sim_result <- complete_case_or_data(new_data, CD)
    CD_sim_result <- complete_performance_one_iter(CD_sim_result, orig_sample_size = my_sample_size, 
                                                   true_effect = true_effect, biomarker_name = "biomarker")
    
    CD_total_result <- dplyr::bind_rows(CD_total_result, CD_sim_result)
    
    # STEP 3) Do PMM
    PMM_formula <- "outcome ~ biomarker_MCAR + confounder + predictor"
    new_data_PMM <- new_data |>
      dplyr::select(outcome, biomarker_MCAR, confounder, predictor)
    
    PMM_sim_result <- PMM(new_data_PMM, PMM_formula, num_imp)
    PMM_sim_result <- MI_performance_one_iter(PMM_sim_result,
                                              "biomarker_MCAR",
                                              true_value = true_effect)
    
    # Append to PMM
    PMM_total_result <- dplyr::bind_rows(PMM_total_result, PMM_sim_result)
    
    # STEP 4) Do CLQI
    quant_reg_formula <- "biomarker_MCAR_transformed ~ outcome + confounder + predictor"
    estimation_formula <- "outcome ~ biomarker_MCAR_CLQI_untransformed + confounder + predictor"
    
    CLQI_sim_result <- CLQI(my_data = new_data, 
                            var_for_imp = "biomarker_MCAR_transformed",
                            transformed_imputation_relationship = quant_reg_formula, 
                            correct_formula = estimation_formula, 
                            num_imp = num_imp,
                            min_val = 0,
                            max_val = 20)
    CLQI_sim_result <- MI_performance_one_iter(CLQI_sim_result,
                                               "biomarker_MCAR_CLQI_untransformed",
                                               true_value = true_effect)
    
    # Append to CLQI
    CLQI_total_result <- dplyr::bind_rows(CLQI_total_result, CLQI_sim_result)
    
    #end this by showing imputation
    iteration_indicator <- iteration_indicator + 1
    print(iteration_indicator)
  }
  
  #gather all into a list
  sim_list <- list(CC_total_result, CD_total_result, PMM_total_result, CLQI_total_result)
  
  
  #obtain full performance measure calculations
  a <- full_performance(CC_total_result, 
                        rel_eff_comparison = CD_total_result)
  
  b <- full_performance(CD_total_result, 
                        rel_eff_comparison = CD_total_result)
  
  c <- full_performance(PMM_total_result, 
                        rel_eff_comparison = CD_total_result)
  
  d <- full_performance(CLQI_total_result, 
                        rel_eff_comparison = CD_total_result)
  
  #return the list
  return(list(CC_total_result, a, 
              CD_total_result, b,
              PMM_total_result, c,
              CLQI_total_result, d))
}

#*DEFINE MAR VERSION HERE!!*



############################################################################ TEST CASE
tictoc::tic()
tester <- repeat_full_simulation_MCAR(num_simulations = 5, #this is solely a test 
                                      my_sample_size = 1000, #again a test
                                      num_imp = 30, #again a test
                                      coefficient_effects = c(logit(0.1), log(1.1), log(0.7), log(0.85)),
                                      true_effect = log(1.1), 
                                      prop_missing_MCAR = 0.3)
tictoc::toc()

############################################################################ furrr across all possible settings from grid_design

run_subset_simulations_MCAR <- function(subset_design_grid, num_sim = 5000, num_imp = 30) {
  
  # Lookup for effects
  effect_lookup <- list(
    none = c(logit(0.1), log(1), log(0.7), log(0.85)),
    weak = c(logit(0.1), log(1.1), log(0.7), log(0.85)),
    strong = c(logit(0.1), log(1.5), log(0.7), log(0.85))
  )
  
  # Add num_sim and effect vectors to grid
  subset_design_grid <- subset_design_grid |>
    mutate(
      num_sim = num_sim,
      coefficient_effects = map(effect_strength, ~ effect_lookup[[.x]]),
      true_effect = map_dbl(effect_strength, ~ effect_lookup[[.x]][2])
    )
  
  # Parallel simulation run
  subset_sim_results <- furrr::future_pmap(
    .l = list(
      num_sim = subset_design_grid$num_sim,
      sample_size = subset_design_grid$sample_size,
      prop_missing_MCAR = subset_design_grid$prop_missing,
      coefficient_effects = subset_design_grid$coefficient_effects,
      true_effect = subset_design_grid$true_effect,
      effect_strength = subset_design_grid$effect_strength
    ),
    .f = function(num_sim, sample_size, prop_missing_MCAR, coefficient_effects, true_effect, effect_strength) {
      
      dir.create("individual_results_MCAR") #make sure directory exists
      file_name <- glue("n{sample_size}_prop{round(prop_missing_MCAR * 100)}_effect{effect_strength}.rds")
      file_path <- file.path("individual_results_MCAR", file_name)
      
      message("Running: ", file_name)
      
      # Run your simulation (returns a list)
      sim_result <- repeat_full_simulation_MCAR(
        num_simulations = num_sim,
        my_sample_size = sample_size,
        num_imp = num_imp,
        coefficient_effects = coefficient_effects,
        true_effect = true_effect,
        prop_missing_MCAR = prop_missing_MCAR
      )
      
      # Append metadata to list (no mutate)
      sim_result$sample_size <- sample_size
      sim_result$prop_missing_MCAR <- prop_missing_MCAR
      sim_result$effect_strength <- effect_strength
      
      # save results
      saveRDS(sim_result, file_path)
      
      return(sim_result)
    },
    .options = furrr_options(seed = TRUE)
  ) %>%
    set_names(paste0(
      "n", subset_design_grid$sample_size,
      "_prop", round(subset_design_grid$prop_missing * 100),
      "_effect", subset_design_grid$effect_strength
    ))
  
  return(subset_sim_results)
}

############################################################################ Get my results

set.seed(525) #set seed for reproducibility
cores <- 9 #for number of simulation settings: parallelize based on this value
future::plan(multisession, workers = cores) #change to cores when this becomes implemented

run_subset_simulations_MCAR(grid_design_MCAR, num_sim = 5000, num_imp = 30) #this saves all the results into a folder for loading