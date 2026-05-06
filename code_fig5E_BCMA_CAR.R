### General information ----
# Title:            Development, validation and implementation of the antibody-secreting cell maturity index (ASC-ME): 
#                   universal prediction of human plasma cell maturity
# Author:           Tobit D. Steinmetz
# Department:       Rheumatology and Clinical Immunology
# Affiliation:      University Medical Center Groningen
# Email:            d.t.steinmetz@umcg.nl
# Collaboration:    please ask permission from the author before using this script
# Last adjustment:  05-05-2026
# Remark:           part of this code was generated or troubleshooted using ChatGPT

### Reference:
# iScience 2026

### This R-script requires the following packages:
library(randomForest)
library(mice)
library(readxl)
library(dplyr)
library(tidyverse)
library(dunn.test)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

# import required data from model development
set.seed(557)
rf_ascmi<-readRDS("final_model.rds")
joined_data<-read_xlsx("final_joined_data.xlsx")
data.norm<-read_xlsx("final_data.norm.xlsx")
model_output<-read_xlsx("final_model_output.xlsx")
prediction_marker<-c("CD19","CD20","CD28","CD45","CD56","CD138","HLA.DR","Ki67")

# Define a function to make and average predictions
average_predictions <- function(data, models) {
  predictions_list <- map(models, ~ predict(.x, newdata = data))
  predictions_df <- bind_cols(predictions_list)
  predictions_df <- predictions_df %>%
    mutate(p = rowMeans(across(everything()))) %>%
    mutate(ID = data$sample,dpi = data$dpi) %>%
    select(ID,dpi,p)
  return(predictions_df)
}

data_new<- joined_data %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86,subfolder,V2,Race,Race.Specify,Ethnicity,ARM.Accession,Expsample.Accession,Original.File.Name)) %>% 
  filter(!(sample=="NA")) %>%
  rename(dpi = Study.Time.Collected) %>%
  rename(group = ARM.Name) %>%
  rename(datatype = File.Detail) %>%
  rename(Age = Subject.Age) %>%
  rename(StudyID = Study.Accession) %>%
  rename(SubjectID = Subject.Accession) %>%
  rename(organ = Biosample.Type) %>%
  mutate(across(3:10, ~ ifelse(.x < -2 | .x > 5, NA, .x))) %>%
  filter(count>9) %>%
  filter(StudyID=="Z4KB")

new_merging <- data_new %>% 
  select(sample,count,source_file,SubjectID,Age,dpi)
new_imp <- data_new %>% 
  select(-c(count,source_file,SubjectID,Age,dpi))

iter=20
data.total<-rbind(data.norm[,-c(11,15)], cbind(new_imp[,c(1:14)]))
imp.new <- mice(data.total, method = "rf", m = iter, seed = 557)

new.imputed <- list()
for(m in 1:20){
  new.imputed[[m]] <- complete(imp.new, action = m)
} 
for(m in 1:20){  
  new.imputed[[m]] <- cbind(new.imputed[[m]][c((nrow(data.norm)+1):nrow(new.imputed[[m]])),], new_merging)
}

#get average imputed values for prediction markers

new.imputed_mean <- new.imputed[[1]][,c(1,15,11)]
for (o in 1:length(prediction_marker)) {
  new_average_marker <- matrix(0, nrow(new.imputed[[1]]))
  for (m in 1:length(new.imputed)) {
    new_average_marker <- new_average_marker + new.imputed[[m]][,o+1]
  }
  new_average_marker <- new_average_marker / length(new.imputed)
  new.imputed_mean <- cbind(new.imputed_mean,new_average_marker)
  colnames(new.imputed_mean)[o+3] <- paste(prediction_marker[o])
}

#
#use model on all imputed data
prediction_new <- list()
for (m in 1:length(new.imputed)) {
  prediction_new[[m]] <- average_predictions(new.imputed[[m]], rf_ascmi)
}

#
new_predictions <- matrix(0, nrow(prediction_new[[1]]), ncol(prediction_new[[1]]))
for (m in 1:length(prediction_new)) {
  new_predictions <- new_predictions + prediction_new[[m]][,3]
}
new_predictions <- new_predictions / length(prediction_new)
predictions_new <- cbind(new_merging,new_imp[,c(13,14,12)], prediction_new[[1]][, 1:2], new_predictions)

# trim output table
new_output<-predictions_new[,-11]

ASC_per<-read_xlsx("BCMA_Z4KB_ASCper.xlsx") %>%  filter(pop=="g.ASC")

plot_data<-left_join(new_output, ASC_per[,c(1,4)], by = "sample")

# scaling factor so both variables fit in same plot space
scale_factor <- max(plot_data$p, na.rm = TRUE) / max(plot_data$ASCper, na.rm = TRUE)

ggplot(plot_data, aes(x = dpi)) +
  # left y-axis: p
  geom_point(aes(y = p, fill = "p"), shape = 21, size = 3, color = "black") +
  # right y-axis: ASCper (rescaled)
  geom_point(aes(y = ASCper * scale_factor, fill = "ASCper"), shape = 21, size = 3, color = "black") +
  
  scale_y_continuous(
    name = "predicted DPI",
    sec.axis = sec_axis(~ . / scale_factor, name = "ASC (% of lymphocytes)")
  ) +
  scale_fill_manual(
    name = "Variable",
    values = c("p" = "blue", "ASCper" = "red")
  ) +
  labs(x = "days post BCMA CAR-T cell treatment") +
  theme_bw() +
  theme(legend.position = "none",
      axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
      axis.title.y.left = element_text(size = 20,face = "bold", colour = "blue"),
      axis.title.y.right = element_text(size = 20,face = "bold", colour = "red"),
      axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
      axis.title.x = element_text(size = 20,face = "bold"),
      plot.title = element_text(size = 15,face = "bold"))
ggsave("output/Fig5E_plot_BCMA_CAR.png", width = 3000, height = 1920, units = "px")
