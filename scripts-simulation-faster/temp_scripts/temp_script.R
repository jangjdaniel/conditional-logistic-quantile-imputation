## ----------------------------------------------------------------------------------------------------------
source("00_CLQI_functions_and_packages.R")
set.seed(525)


## ----------------------------------------------------------------------------------------------------------
my_sample <- 1000 
min <- 0
max <- 20 #change the data_generator function to include this
weak_effects <- c(logit(0.1), log(1.1), log(0.7), log(0.85))
strong_effects <- c(logit(0.1), log(1.5), log(0.7), log(0.85))


## ----------------------------------------------------------------------------------------------------------
#in general, to check the true logistic regression for the data before missingness, we need this
formula_biomarker <- "outcome ~ biomarker + confounder + predictor"
formula_biomarker_MCAR <- "outcome ~ biomarker_MCAR + confounder + predictor"
formula_biomarker_MAR <- "outcome ~ biomarker_MAR + confounder + predictor"


## ----------------------------------------------------------------------------------------------------------
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
                      a = 0, b = 20, 
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # Generate data
# weak_effects_data <- data_generator(my_sample, weak_effects)
# 
# # Run LOGISTIC regression: isn't good...
# glm(as.formula(formula_biomarker),
#     data = weak_effects_data,
#     family = "binomial") |>
#   tbl_regression(exponentiate = TRUE)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # Generate data
# strong_effects_data <- data_generator(my_sample, strong_effects)
# 
# # Run LOGISTIC regression
# glm(as.formula(formula_biomarker),
#     data = strong_effects_data,
#     family = "binomial") |>
#   tbl_regression(exponentiate = TRUE)


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# weak_effects_data <- data_generator(my_sample, weak_effects)
# 
# missing_test_case <- missing_generator(weak_effects_data, 0.30) |>
#   dplyr::select(biomarker_MCAR)
# 
# #should be 300
# sum(is.na(missing_test_case$biomarker_MCAR))
# rm(missing_test_case)


## ----eval=FALSE--------------------------------------------------------------------------------------------
# weak_effects_data <- data_generator(my_sample, weak_effects)
# 
# weak_effects_data |>
#   mutate(biomarker_transformed = sapply(biomarker, log_quant_transform,
#                                         min = 0, max = 20)) |>
#   mutate(biomarker_untransformed = sapply(biomarker_transformed, inv_log_quant_transform,
#                                           min = 0, max = 20)) |>
#   mutate(diff = biomarker - biomarker_untransformed) |> # see how close the original and untransformed vals are
#   dplyr::select(diff) #look into it.


## ----------------------------------------------------------------------------------------------------------
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


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # Generate data!
# DGM_weak <- data_generating_mechanism(my_sample = my_sample,
#                                       beta_coefficients = weak_effects,
#                                       prop_missing_MCAR = 0.3,
#                                       min_val = 0,
#                                       max_val = 20)
# 
# # With biomarker data
# summary(glm(as.formula(formula_biomarker),
#             data = DGM_weak,
#             family = "binomial"))
# 
# # With MCAR biomarker data: 95% CI will be wider
# summary(glm(as.formula(formula_biomarker),
#             data = DGM_weak,
#             family = "binomial"))


## ----eval=FALSE--------------------------------------------------------------------------------------------
# # Generate data!
# DGM_strong <- data_generating_mechanism(my_sample = my_sample,
#                                         beta_coefficients = weak_effects,
#                                         prop_missing_MCAR = 0.3,
#                                         min_val = 0,
#                                         max_val = 20)
# 
# # With biomarker data
# glm(as.formula(formula_biomarker),
#     data = DGM_strong,
#     family = "binomial") |>
#   broom::tidy()
# 
# # With MCAR biomarker data: 95% CI will be wider
# glm(as.formula(formula_biomarker),
#     data = DGM_strong,
#     family = "binomial") |>
#   broom::tidy()


## ----------------------------------------------------------------------------------------------------------
DGM_weak <- data_generating_mechanism(my_sample = my_sample,
                                      beta_coefficients = weak_effects,
                                      prop_missing_MCAR = 0.3,
                                      min_val = 0,
                                      max_val = 20)

#this is just some chatgpt code for a quick check
DGM_weak$confounder <- factor(DGM_weak$confounder, levels = c(0, 1), labels = c("Confounder = 0", "Confounder = 1"))

# Create the plot with overlaid density curves
ggplot(data = DGM_weak, aes(x = biomarker, fill = confounder, color = confounder)) +
  geom_density(alpha = 0.5) + 
  labs(title = "Biomarker Distribution by Confounder",
       x = "Biomarker",
       y = "Density") +
  scale_fill_manual(values = c("skyblue", "orange")) +  # Custom colors for filling
  scale_color_manual(values = c("blue", "red"))        # Custom border colors


## ----------------------------------------------------------------------------------------------------------
rm(my_sample)

# weak
rm(weak_effects)
rm(weak_effects_data)
rm(DGM_weak)

# strong
rm(strong_effects)
rm(strong_effects_data)
rm(DGM_strong)

