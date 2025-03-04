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

### introduction & Read-ME ###
# this scripts contains all necessary steps to apply the ASC-ME (antibody-secreting cell maturity index)
# to your own data, given that your data is adequately processed and formatted
# if necessary, use the R script "example_cytometry_code_Z6YW.R" to extrate the raw data from your FCS files 
### 

### This R-script requires the following packages:
library(randomForest)
library(mice)
library(readxl)
library(dplyr)
library(tidyverse)
library(writexl)

## continue with R environment from training of the ASC-ME model
# the R environment "env_train_model_final.RData" contains required standardization data and the random forest model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

set.seed(555)

## import your own data (adjust file name if necessary)
# !IMPORTANT!:
# your excel sheet requires to following layout with 17 columns (A-Q):
# A: "sample" (name of your samples / sample ID; required information)
# B: "dpi" (days since last vaccination or infection; required information - enter "0" if unknown)
# C-J: "maturity markers" (CD19, CD20, CD28, CD45, CD56, CD138, HLA.DR, Ki67 - as normalized expression value)
# !IMPORTANT!: a minimum of one marker is required! 
# !IMPORTANT!: Cytometry-acquired markers presumably give a more reliable result than imputed values

# K: "StudyID" (name of your Study)
# L: "organ" (e.g. "blood" or a specific tissue 
# M: "Gender" (sex/gender of sampled individual)
# N: "Age" (age of individual even sample was taken)
# O: "File" (name of fcs file)
# P: "group" (e.g. healthy or disease groups)
# Q: "datatype" ("flow" or "CyTOF", dependent on your data from flow or mass cytometry)
## if data is not present leave the respective entry blank

data_new <- read_excel("my_own_data.xlsx")

# trim input for missing data imputation
new_merging <- data_new %>% select(sample,File)
new_imp <- data_new %>% select(-c(File))
# calculate percentage of missing data
percent_miss<-(sum(is.na(new_imp)) / (nrow(new_imp) * ncol(new_imp))) *100

# align your data with standardization data necessary for reliable data imputation
iter=20
data.norm<-cbind(data.imputed_mean[,c(1,2,4:11)],dataforimp[,c(11:16)])
data.total<-rbind(data.norm,new_imp)
imp.new <- mice(data.total, method = "rf", m = iter, seed = 555)

# extract your imputed data
new.imputed <- list()
for(m in 1:20){
  new.imputed[[m]] <- complete(imp.new, action = m)
} 
for(m in 1:20){  
  new.imputed[[m]] <- cbind(new.imputed[[m]][c((nrow(data.norm)+1):nrow(new.imputed[[m]])),], new_merging)
}

#use the ASC-ME model on your imputed data
prediction_new <- list()
for (m in 1:length(new.imputed)) {
  prediction_new[[m]] <- average_predictions(new.imputed[[m]], rf_ascmi)
}

# average out predictions from the RF bundle
new_predictions <- matrix(0, nrow(prediction_new[[1]]), ncol(prediction_new[[1]]))
for (m in 1:length(prediction_new)) {
  new_predictions <- new_predictions + prediction_new[[m]][,3]
}
new_predictions <- new_predictions / length(prediction_new)
predictions_new <- cbind(new_merging,new_imp[,c(11:12,15)], prediction_new[[1]][, 1:2], new_predictions)

# trim result table
new_output<-predictions_new %>%
  select(sample, StudyID, organ, group, dpi, p) %>%
  mutate(Age=new_imp$Age) %>%
  mutate(Gender=new_imp$Gender)
new_output$Age<-as.numeric(new_output$Age)

## generate your data output and export results
# the predicted DPI/ ASC maturity value is added in column R to your input data 
my_own_data_ASCed<-cbind(data_new,ASC_maturity=new_output$p)
write_xlsx(my_own_data_ASCed, "my_own_data_ASCed.xlsx")
