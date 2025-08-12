## ----------------------------------------------------------------------------------------------------------
knitr::purl("02_method_demonstration.qmd", output = "temp_scripts/02_temp_script.R")
source("temp_scripts/02_temp_script.R")


## ----------------------------------------------------------------------------------------------------------
set.seed(525)

my_sample <- 1000
weak_effects <- c(logit(0.1), log(1.1), log(0.7), log(0.85))
generated_data_for_testing <- data_generating_mechanism(my_sample = my_sample,
                                                        beta_coefficients = weak_effects,
                                                        prop_missing_MCAR = 0.3,
                                                        min_val = 0,
                                                        max_val = 20)


quant_reg_formula <- "biomarker_MCAR_transformed ~ outcome + confounder + predictor"
estimation_formula <- "outcome ~ biomarker_MCAR_CLQI_untransformed + confounder + predictor"

tictoc::tic()
CLQI_results <- CLQI(my_data = generated_data_for_testing, 
                     var_for_imp = "biomarker_MCAR_transformed",
                     transformed_imputation_relationship = quant_reg_formula, 
                     correct_formula = estimation_formula, 
                     num_imp = 10,
                     min_val = 0,
                     max_val = 20)
tictoc::toc()

CLQI_results

#comparison with complete data
glm(outcome ~ biomarker + confounder + predictor,
    data = generated_data_for_testing,
    family = "binomial")


## ----------------------------------------------------------------------------------------------------------
rubin_rule_estimate <- function(extracted_statistics) {
  return(mean(extracted_statistics$Estimate))
}


## ----eval=FALSE--------------------------------------------------------------------------------------------
# rubin_rule_estimate(CLQI_results)


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# rubin_rule_SE(CLQI_results)


## ----------------------------------------------------------------------------------------------------------
coverage <- function(RR_estimate, SE_estimate, orig_sample_size, true_value) {
  t_star <- qt(0.975, df = orig_sample_size - 1) #good enough estimate
  
  coverage_indicator <- ifelse((RR_estimate - t_star*SE_estimate) <= true_value & true_value <= (RR_estimate + t_star*SE_estimate),
                               1, 0) #if true value is covered, it's 1. 
  
  return(coverage_indicator)
}


## ----------------------------------------------------------------------------------------------------------
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


## ----------------------------------------------------------------------------------------------------------
rel_efficiency <- function(simulation_data_A, simulation_data_B) {
  empirical_variance_A <- var(simulation_data_A$Estimate)
  empirical_variance_B <- var(simulation_data_B$Estimate)
  
  # Then return our ratio
  return(empirical_variance_B/empirical_variance_A)
}


## ----------------------------------------------------------------------------------------------------------
MCSE <- function(vec_values_interest, num_sim) {
  return(sd(vec_values_interest) / sqrt(num_sim))
}

