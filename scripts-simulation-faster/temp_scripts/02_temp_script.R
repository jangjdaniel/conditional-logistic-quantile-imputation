## ----------------------------------------------------------------------------------------------------------
library(knitr)
knitr::purl("01_data_generating_mechanism.qmd", output = "temp_scripts/01_temp_script.R")
source("temp_scripts/01_temp_script.R")


## ----------------------------------------------------------------------------------------------------------
set.seed(525)

my_sample <- 1000
weak_effects <- c(logit(0.1), log(1.1), log(0.7), log(0.85))
strong_effects <- c(logit(0.1), log(1.5), log(0.7), log(0.85))

num_imp <- 10 #number of imputations
generated_data_for_testing <- data_generating_mechanism(my_sample = my_sample,
                                                        beta_coefficients = weak_effects,
                                                        prop_missing_MCAR = 0.3,
                                                        min_val = 0,
                                                        max_val = 20)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# glm(outcome ~ biomarker_MCAR + confounder + predictor,
#     data = generated_data_for_testing,
#     family = "binomial") |>
#   broom::tidy() |>
#   filter(term == "biomarker_MCAR")


## ----eval=FALSE--------------------------------------------------------------------------------------------
# glm(outcome ~ biomarker + confounder + predictor,
#     data = generated_data_for_testing,
#     family = "binomial") |>
#   broom::tidy() |>
#   filter(term == "biomarker")


## ----------------------------------------------------------------------------------------------------------
complete_case_or_data <- function(my_data, correct_formula) {

  #perform logistic regression, get everything into a tibble, and return that
  CCD <-   glm(as.formula(correct_formula),
               data = my_data,
               family = "binomial") |>
    summary() |>
    (\(x) as_tibble(x$coefficients, rownames = "term"))()
  
  return(CCD)
}


## ----eval=FALSE--------------------------------------------------------------------------------------------
# CC <- "outcome ~ biomarker_MCAR + confounder + predictor"
# complete_case_or_data(generated_data_for_testing, CC)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# CD <- "outcome ~ biomarker + confounder + predictor"
# complete_case_or_data(generated_data_for_testing, CD)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # mice requires the dataframe to only contain missing variable and its predictors
# generated_data_for_PMM <- generated_data_for_testing |>
#   dplyr::select(confounder, biomarker_MCAR, outcome, predictor)
# 
# #perform the mice function!
# PMM_test <- mice::mice(generated_data_for_PMM,
#                        method = "pmm",
#                        m = num_imp)
# 
# #check that a dataset has been made: there should be no NA's
# sum(is.na(complete(PMM_test, action = 1) |> dplyr::select(biomarker_MCAR)))


## ----eval=FALSE--------------------------------------------------------------------------------------------
# data_iteration <- complete(PMM_test, action = 1)
# log_reg_estimate <- glm(as.formula(formula_biomarker_MCAR),
#                           data = data_iteration,
#                           family = "binomial")
# 
# summary(log_reg_estimate)$coefficients |> as.tibble(rownames = "term") #this is faster
# log_reg_estimate |> broom::tidy()


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # purrr function to get coefficients: let me break it down for you
# PMM_coefficients <- purrr::map_dbl(1:num_imp, ~ { #we repeat this 10 times
# 
#   #extract the ith dataset (.x refers to the 1:num_imp)
#   data_iteration <- complete(PMM_test, action = .x)
# 
#   #perform logistic regression
#   log_reg_estimate <- glm(as.formula(formula_biomarker_MCAR),
#                           data = data_iteration,
#                           family = "binomial") |>
#     broom::tidy() |> #make it into a nice data frame
#     filter(term == "biomarker_MCAR") |> #extract biomarker_MCAR row
#     pull(estimate) #we only care about the estimate
# 
#   #this is not necessary, but for organization
#   log_reg_estimate
#   }
# )
# 
# # after this, we can perform rubin's rule for the estimate to check for bias
# mean(PMM_coefficients) - weak_effects[2]


## ----eval=FALSE--------------------------------------------------------------------------------------------
# tictoc::tic()
# PMM_regression_results <- purrr::map(1:num_imp, ~ { #we repeat this 10 times
# 
#   #extract the ith dataset (.x refers to the 1:num_imp)
#   data_iteration <- complete(PMM_test, action = .x)
# 
#   #perform logistic regression, extract results, and have that be the printout
#   glm(as.formula(formula_biomarker_MCAR),
#       data = data_iteration,
#       family = "binomial") |>
#     summary() |>
#     (\(x) as_tibble(x$coefficients, rownames = "term"))()
#   }
# )
# tictoc::toc()
# 
# PMM_regression_results


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# generated_data_MCAR <- generated_data_for_testing |>
#   dplyr::select(confounder, biomarker_MCAR, outcome, predictor)
# 
# PMM(my_data_filtered = generated_data_MCAR,
#     correct_formula = formula_biomarker_MCAR,
#     num_imp = 10)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# generated_data_MCAR <- generated_data_for_testing |>
#   dplyr::select(confounder, biomarker_MCAR, outcome, predictor)
# 
# tictoc::tic()
# PMM(my_data_filtered = generated_data_MCAR,
#     correct_formula = formula_biomarker_MCAR,
#     num_imp = 20)
# tictoc::toc()


## ----eval=FALSE--------------------------------------------------------------------------------------------
# generated_data_MAR <- generated_data_for_testing |>
#   dplyr::select(confounder, biomarker_MAR, outcome, predictor)
# 
# PMM(my_data_filtered = generated_data_MAR,
#     correct_formula = formula_biomarker_MAR,
#     num_imp = 10)


## ----------------------------------------------------------------------------------------------------------
# quantile regression!
reg_coeff <- quantreg::rq(biomarker_MCAR_transformed ~ outcome + confounder + predictor, 
                          data = generated_data_for_testing, 
                          tau = 1) #max quantile

as.numeric(reg_coeff$coefficients["(Intercept)"])

# or with broom:
reg_coeff <- reg_coeff |> broom::tidy()

# and here is a way to easily extract regression coefficients
# example of intercept estimate!
reg_coeff |> filter(term == "(Intercept)") |> pull(estimate) #we only care about the estimate


## ----------------------------------------------------------------------------------------------------------
logistic_quantile_imputation <- function(my_data, transformed_imputation_relationship, row_index) {
  
  # Step 1) generate a random uniform value
  u <- runif(1, min = 0, max = 1) 
  
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# tictoc::tic()
# 
# my_formula <- "biomarker_MCAR_transformed ~ outcome + confounder + predictor"
# 
# generated_data_for_testing <- generated_data_for_testing |>
#   mutate(biomarker_MCAR_transformed_CLQI = biomarker_MCAR_transformed) #create another copy to impute and compare
# 
# # Perform the algorithm
# for(row_index in 1:nrow(generated_data_for_testing)) { #for all indicies
#   if(is.na(generated_data_for_testing$biomarker_MCAR_transformed[row_index])) { #if the cell in the vector is empty
# 
#       imputed_value <- logistic_quantile_imputation(my_data = generated_data_for_testing, #implement algorithm
#                                             transformed_imputation_relationship = my_formula,
#                                             row_index = row_index)
# 
#       generated_data_for_testing$biomarker_MCAR_transformed_CLQI[row_index] <- imputed_value #and apply it
#   }
#   else {
#       generated_data_for_testing$biomarker_MCAR_transformed_CLQI[row_index] <- #if it has a value
#         generated_data_for_testing$biomarker_MCAR_transformed[row_index] # don't change anything
# 
#   }
# }
# 
# # Step 4) Sanity check
# generated_data_for_testing |>
#   dplyr::select(biomarker_MCAR_transformed, biomarker_MCAR_transformed_CLQI) |>
#   mutate(diff = biomarker_MCAR_transformed_CLQI - biomarker_MCAR_transformed)
# 
# hist(generated_data_for_testing$biomarker_MCAR_transformed_CLQI)
# 
# tictoc::toc()


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# my_formula <- "biomarker_MCAR_transformed ~ outcome + confounder + predictor"
# generated_data_for_CLQI <- implement_logistic_quantile_imputation(generated_data_for_testing,
#                                                                   my_formula,
#                                                                   "biomarker_MCAR_transformed",
#                                                                   min_val = 0,
#                                                                   max_val = 20)
# # mutating to see what happened
# generated_data_for_CLQI <- generated_data_for_CLQI |>
#   mutate(diff = biomarker_MCAR_CLQI_untransformed - biomarker_MCAR)
# 
# # this should be zero
# sum(generated_data_for_CLQI$diff, na.rm=TRUE)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # visualize what the imputed biomarker distributions look like
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker_MCAR_transformed_CLQI)) +
#   geom_density()
# 
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker_MCAR_transformed)) +
#   geom_density()
# 
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker_transformed)) +
#   geom_density()


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# quant_reg_formula <- "biomarker_MCAR_transformed ~ outcome + confounder + predictor"
# estimation_formula <- "outcome ~ biomarker_MCAR_CLQI_untransformed + confounder + predictor"
# 
# tictoc::tic()
# CLQI_results <- CLQI(my_data = generated_data_for_testing,
#                      var_for_imp = "biomarker_MCAR_transformed",
#                      transformed_imputation_relationship = quant_reg_formula,
#                      correct_formula = estimation_formula,
#                      num_imp = 10,
#                      min_val = 0,
#                      max_val = 20)
# tictoc::toc()
# 
# CLQI_results


## ----eval=FALSE--------------------------------------------------------------------------------------------
# CLQI_results[[1]]
# my_vector <- c()
# 
# for(i in 1:length(CLQI_results)) {
#   my_vector[i] <- CLQI_results[[i]] |> filter(term == "biomarker_MCAR_transformed_CLQI") |> pull(Estimate)
# }
# 
# mean(my_vector) - weak_effects[2]


## ----eval=FALSE--------------------------------------------------------------------------------------------
# estimation_formula <- "outcome ~ biomarker_MAR_transformed_CLQI + confounder + predictor"
# my_formula <- "biomarker_MAR_transformed ~ outcome + confounder + predictor"
# 
# tictoc::tic()
# CLQI_results <- CLQI(my_data = generated_data_for_testing,
#                      var_for_imp = "biomarker_MAR_transformed",
#                      correct_formula = estimation_formula,
#                      transformed_imputation_relationship = my_formula,
#                      num_imp = 20,
#                      min_val = 0,
#                      max_val = 20)
# tictoc::toc()
# 
# CLQI_results


## ----eval=FALSE--------------------------------------------------------------------------------------------
# my_formula <- "biomarker_MCAR ~ outcome + confounder + predictor"
# 
# generated_data_for_CLQI <- implement_logistic_quantile_imputation(generated_data_for_testing,
#                                                                   my_formula,
#                                                                   "biomarker_MCAR",
#                                                                   min_val = 0,
#                                                                   max_val = 20)
# # mutating to see what happened
# generated_data_for_CLQI <- generated_data_for_CLQI |>
#   dplyr::select(biomarker, biomarker_MCAR, biomarker_MCAR_CLQI) |>
#   mutate(diff = biomarker_MCAR_CLQI - biomarker_MCAR)
# 
# # this should be zero
# sum(generated_data_for_CLQI$diff, na.rm=TRUE)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # visualize what the imputed biomarker distributions look like
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker_MCAR_CLQI)) +
#   geom_density()
# 
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker_MCAR)) +
#   geom_density()
# 
# ggplot(data = generated_data_for_CLQI, aes(x = biomarker)) +
#   geom_density()

