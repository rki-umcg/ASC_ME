### General information ----
# Title: Establishing a random forest model to predict the maturity of adult, human plasma cells using various cell surface marker
# Author: Tobit D. Steinmetz
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: d.t.steinmetz@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Date created: 19-10-2023
# Date last adjustment: 20-06-2024
# RStudio version: 2023.06.1 Build 524
# References: 

### introduction ----
# this scripts uses bootstrapping as internal validation method for the ASC-ME model
# 

### This R-script requires the following packages:
library(randomForest)

## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

set.seed(555)

### use bootstrapping to evaluate the performance of the ASC-ME model
# pool all the imputed datasets into one
combined_df <- do.call(rbind, data.full)

# set the number of bootstrapping
num_bootstrap_samples <- 100

# Function to perform bootstrapping
bootstrap_random_forest <- function(data, response_variable, predictors, num_bootstrap_samples) {
  n <- nrow(data)
  
  # Initialize vectors to store error metrics
  mse_values <- numeric(num_bootstrap_samples)
  rmse_values <- numeric(num_bootstrap_samples)
  mae_values <- numeric(num_bootstrap_samples)
  
  # Initialize list to store models and their hyperparameters
  models <- list()
  hyperparameters <- list()
  
  # Create formula for the random forest model
  formula <- as.formula(paste(response_variable, "~", paste(predictors, collapse = " + ")))
  
  for (i in 1:num_bootstrap_samples) {
    # Generate bootstrap sample
    bootstrap_sample <- data[sample(1:n, replace = TRUE), ]
    
    # Identify OOB samples
    oob_indices <- setdiff(1:n, unique(as.numeric(rownames(bootstrap_sample))))
    oob_sample <- data[oob_indices, ]
    
    # Train random forest on bootstrap sample
    rf_model <- randomForest(formula, data = bootstrap_sample, ntree = 100, mtry = 6)
    
    # Predict on OOB sample
    oob_predictions <- predict(rf_model, oob_sample)
    
    # Calculate error metrics
    actual_values <- oob_sample[[response_variable]]
    mse_values[i] <- mean((actual_values - oob_predictions)^2)
    rmse_values[i] <- sqrt(mse_values[i])
    mae_values[i] <- mean(abs(actual_values - oob_predictions))
    
    # Save the model and hyperparameters
    models[[i]] <- rf_model
    hyperparameters[[i]] <- list(ntree = rf_model$ntree, mtry = rf_model$mtry)
  }
  
  # Identify the best model based on the lowest MSE
  best_model_index <- which.min(mse_values)
  best_model <- models[[best_model_index]]
  best_hyperparameters <- hyperparameters[[best_model_index]]
  
  # Return the best model, its hyperparameters, and error metrics
  return(list(
    best_model = best_model,
    best_hyperparameters = best_hyperparameters,
    mse = mse_values,
    rmse = rmse_values,
    mae = mae_values
  ))
}

# Define response variable and predictors
response_variable <- "dpi"
predictors <- c("CD19", "CD20", "CD28", "CD45", "CD56", "CD138", "HLA.DR", "Ki67")

# Perform bootstrapping
bootstrap_results <- bootstrap_random_forest(combined_df, response_variable, predictors, num_bootstrap_samples)

# Summarize error metrics
mean_mse <- mean(bootstrap_results$mse)
mean_rmse <- mean(bootstrap_results$rmse)
mean_mae <- mean(bootstrap_results$mae)
std_mse <- sd(bootstrap_results$mse)
std_rmse <- sd(bootstrap_results$rmse)
std_mae <- sd(bootstrap_results$mae)

# Print results
print(paste("Mean MSE:", mean_mse))
print(paste("Mean RMSE:", mean_rmse))
print(paste("Mean MAE:", mean_mae))
print(paste("Standard Deviation of MSE:", std_mse))
print(paste("Standard Deviation of RMSE:", std_rmse))
print(paste("Standard Deviation of MAE:", std_mae))

# Print the best model and its hyperparameters
print("Best Model:")
print(bootstrap_results$best_model)
#print("Best Hyperparameters:")
#print(bootstrap_results$best_hyperparameters)

