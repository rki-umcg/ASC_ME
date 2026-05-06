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
  
### import output of cytometry analysis

# Path to the parent folder containing the 68 subfolders
parent_folder <- paste(data_path, "AG_spec", sep = "")

# Get the list of subfolders
subfolders <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)

# Loop through subfolders and read in the excel files
spec_data <- lapply(subfolders, function(folder) {
  file_path <- file.path(folder, "final_results.xlsx")
  
  if (file.exists(file_path)) {
    df <- read_excel(file_path)
    df[,c(3:15)] <- lapply(df[,c(3:15)], as.numeric)
    df$subfolder <- basename(folder)  # add subfolder name as a new column
    return(df)
  } else {
    message("No final_results.xlsx found in: ", folder)
    return(NULL)
  }
}) %>% 
  bind_rows()  # combine into one big data frame


### import sample meta data
{
  data_folder <- paste(data_path, ".sample_data", sep = "")
  txt_files <- list.files(data_folder, pattern = "\\.txt$", full.names = TRUE)
  xlsx_files<- list.files(data_folder, pattern = "\\.xlsx$", full.names = TRUE)
  
  # Loop through and read each txt file
  sample_annotations <- lapply(txt_files, function(file) {
    df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    df$Planned.Visit.Name <- lapply(df$Planned.Visit.Name, as.character)
    # Add a column with the filename (without extension)
    df$source_file <- tools::file_path_sans_ext(basename(file))
    
    return(df)
  }) %>% 
    bind_rows()
  
  colnames(sample_annotations)[colnames(sample_annotations)=="File.Name"] <- "sample"
  
  xlsx_annotations <- lapply(xlsx_files, function(file) {
    df <- read_excel(file)
    df$Subject.Age <- as.numeric(df$Subject.Age)
    df$Study.Time.Collected <- as.numeric(df$Study.Time.Collected)
    df$Subject.Accession <- as.character(df$Subject.Accession)
    df$source_file <- tools::file_path_sans_ext(basename(file))
    
    return(df)
  }) %>% 
    bind_rows()
  
  all_annotations<-rbind(sample_annotations[,c(1,3:5,7,10,13:15,20,22,30,38:41)],xlsx_annotations)
  }

joined_spec<-left_join(spec_data, all_annotations, by = "sample") %>%
  mutate(across(Gender, ~ gsub("Not Specified", "NA",.))) %>%
  mutate(across(Gender, ~ gsub("Female", "female",.))) %>%
  mutate(across(Gender, ~ gsub("Male", "male",.))) %>%
  mutate(across(Race, ~ gsub("Other", "NA/other",.))) %>%
  mutate(across(Race, ~ gsub("American Indian or Alaska Native", "NA/other",.))) %>%
  mutate(across(Race, ~ gsub("Not Specified", "NA/other",.))) %>%
  mutate(across(Race, ~ gsub("Native Hawaiian or NA/other Pacific Islander", "NA/other",.))) %>%
  mutate(across(Race, ~ gsub("Unknown", "NA/other",.))) %>%
  mutate(Race   = ifelse(is.na(Race), "NA/other", Race), Gender = ifelse(is.na(Gender), "NA/other", Gender)) %>%
  mutate(across(File.Detail, ~ gsub("flow", "flow_cyto",.))) %>%
  mutate(across(File.Detail, ~ gsub("Flow cytometry result", "flow_cyto",.))) %>%
  mutate(across(File.Detail, ~ gsub("CyTOF result", "mass_cyto",.))) %>%
  mutate(across(File.Detail, ~ gsub("CyTOF", "mass_cyto",.))) %>%
  mutate(source_file = subfolder)

data_new<- joined_spec %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86,subfolder,Race,Race.Specify,Ethnicity,ARM.Accession,Expsample.Accession,Original.File.Name)) %>% 
  filter(!(sample=="NA")) %>%
  rename(dpi = Study.Time.Collected) %>%
  rename(group = ARM.Name) %>%
  rename(datatype = File.Detail) %>%
  rename(Age = Subject.Age) %>%
  rename(StudyID = Study.Accession) %>%
  rename(SubjectID = Subject.Accession) %>%
  rename(organ = Biosample.Type) %>%
  mutate(across(3:10, ~ ifelse(.x < -2 | .x > 5, NA, .x)))

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
new_output<-predictions_new[,-11] #%>%  select(sample, StudyID,SubjectID, organ, group, dpi, p)

plot_data<-new_output %>% filter((source_file=="Z68Q_RBD")|(source_file=="Z68Q_total"))

plot_data_RBD<-new_output %>% filter((source_file=="Z7A5_RBD_spec")|(source_file=="Z7A5_RBD_total"))
plot_data_TT<-new_output %>% filter((source_file=="Z7A5_TT_spec")|(source_file=="Z7A5_TT_total"))
test_RBD<-wilcox.test(plot_data_RBD[c(1:20),11],plot_data_RBD[c(21:40),11], paired = TRUE)
test_TT<-wilcox.test(plot_data_TT[c(1:20),11],plot_data_TT[c(21:40),11], paired = TRUE)
test_Z68Q_RBD<-wilcox.test(plot_data[c(1:22),11],plot_data[c(23:44),11], paired = TRUE)

### plot TT-specific ASC

plot_TT<-ggplot(plot_data_TT, aes(x = source_file, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "Z7A5_TT",       x = NULL,       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("Z7A5_TT_total","Z7A5_TT_spec") +
  annotate("text", x=1, y=max(plot_data_TT$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(plot_data_TT$p)*1.1, xend = 1.6, yend = max(plot_data_TT$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x=2, y=max(plot_data_TT$p)*1.1, label=paste("p =",format(test_TT[["p.value"]], digits=3)), size=5)
ggsave(filename = "output/Fig5B_plot_Z7A5_TT.png", plot = plot_TT,
       width = 3000, height = 1920, units = "px")

###
# plot combined data of Z68Q and Z7A5 RBD-specific ASC together

plot_data_RBD<-new_output %>% filter((source_file=="Z68Q_RBD")|(source_file=="Z68Q_total")|(source_file=="Z7A5_RBD_spec")|(source_file=="Z7A5_RBD_total"))
plot_RBD<-ggplot(plot_data_RBD, aes(x = source_file, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "Z68Q and Z7A5",       x = NULL,       y = "predicted DPI",
       color = "StudyID") +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("Z68Q_total", "Z68Q_RBD","Z7A5_RBD_total","Z7A5_RBD_spec") +
  annotate("text", x=1, y=max(plot_data_RBD$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(plot_data_RBD$p)*1.1, xend = 1.6, yend = max(plot_data_RBD$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x=2, y=max(plot_data_RBD$p)*1.1, label=paste("p =",format(test_Z68Q_RBD[["p.value"]], digits=3)), size=5)+
  annotate("text", x=3, y=max(plot_data_RBD$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 3.4, y = max(plot_data_RBD$p)*1.1, xend = 3.6, yend = max(plot_data_RBD$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x=4, y=max(plot_data_RBD$p)*1.1, label=paste("p =",format(test_RBD[["p.value"]], digits=3)), size=5)
ggsave(filename = "output/Fig5B_plot_RBD_spec_comp.png",
       width = 3000, height = 1920, units = "px")
