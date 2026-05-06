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
library(caret)
library(mice)
library(readxl)
library(writexl)
library(gtools)
library(VIM)
library(dplyr)
library(randomForestExplainer)
library(tidyverse)
library(skimr)
library(ggplot2)
library(viridis)
library(dunn.test)
library(gridExtra)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

### import output of cytometry analysis

# Path to the parent folder containing the 68 subfolders
parent_folder <- paste(data_path,".new_output", sep = "")

# Get the list of subfolders
subfolders <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)

# Loop through subfolders and read in the excel files
all_data <- lapply(subfolders, function(folder) {
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
  data_folder <- paste(data_path,".sample_data", sep = "")
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

# joining both dataframes togehter
# also streamline Gender/Sex, Race/Ethnicity and cytometry type
joined_data<-left_join(all_data, all_annotations, by = "sample") %>%
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
  mutate(across(File.Detail, ~ gsub("CyTOF", "mass_cyto",.)))


# filter kinetic data for model training and validation
filtered_kinetic_data<-joined_data %>% filter(count>49) %>% filter(Subject.Age>9.9999) %>%
  filter((source_file=="AID_UMCG")|(source_file=="CER_UMCG")|(source_file=="SDY1086_sample")|(source_file=="SDY1288_sample")|
           (source_file=="SDY144_sample")|(source_file=="SDY1669_sample")|(source_file=="SDY1734_sample")|
           (source_file=="SDY224_sample")|(source_file=="SDY364_sample")|(source_file=="SDY387_sample")|
           (source_file=="SDY522_sample")|(source_file=="SDY272_sample")|(source_file=="SDY984_sample")|
           (source_file=="SDY180_sample")|(source_file=="SDY296_sample")|(source_file=="SDY301_sample")|
           (source_file=="SDY819_sample")|(source_file=="ZY77")|(source_file=="SDY1397_sample")|
           (source_file=="SDY368_sample")|(source_file=="SDY788_sample")|(source_file=="SDY34_sample")|
           (source_file=="SDY648_sample")|(source_file=="SDY739_sample"))

# streamline sample annotations for kinetic information
filtered_kinetic_data<-filtered_kinetic_data %>%
  mutate(across(ARM.Name, ~ gsub("Vax002_group_A", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Vax002_group_B", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("CHIKV patient", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Infected individual", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("TIV Vaccine", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Emory cohort", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Hong Kong cohort", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Africans", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Europeans", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Study group 1 Pneunomax23", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Study group 1 Saline", "excl",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Study group 1 2009-2010 Fluzone", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("TIV 2010", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Young", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Aged", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("AIRFV 2011-12", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("AIRFV 2012-13", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Healthy Controls", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("6 month post-transplant subjects receiving trivalent influenza vaccine", "TP_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("NCH-2012-13", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("NCH-2013-14", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("NCH-2010-11", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("LAIV Vaccine", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Immune response to Influenza vaccination in aged populations - Year 3", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Immune response to Influenza vaccination in aged populations - Year 4", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("control", "HC",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Non-responder", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Responder", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Cohort2", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("Immune response to Influenza vaccination in aged populations - Year 5", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("elderly", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("young", "HC_Vac",.))) %>%
  mutate(across(ARM.Name, ~ gsub("uncomplicated", "HC",.))) %>%
  mutate(across(ARM.Name, ~ gsub("infections", "HC_Inf",.))) %>%
  mutate(across(ARM.Name, ~ gsub("inf_relapse", "HC_Inf",.)))

filtered_kinetic_data<-filtered_kinetic_data %>%
  filter((ARM.Name=="HC_Vac")|(ARM.Name=="HC")|(ARM.Name=="HC_Inf")) %>%
  filter(Study.Time.Collected<181)

split_by_sourcefile <- function(filtered_kinetic_data,
                                source_col = "source_file",
                                cat_cols = c("Race", "Gender"),
                                num_cols = c("Subject.Age"),
                                balance_col = "Study.Time.Collected",  # NEW: main variable to balance
                                train_frac = 0.75,
                                iterations = 5000,
                                tol = 0.02,
                                seed = NULL) {
  library(dplyr)
  if (!is.null(seed)) set.seed(seed)
  
  df <- filtered_kinetic_data
  total_n <- nrow(df)
  target <- train_frac * total_n
  
  # summarize counts per source file
  file_summary <- df %>%
    group_by_at(source_col) %>%
    summarise(n = n(), .groups = "drop")
  files <- file_summary[[source_col]]
  
  best_score <- Inf
  best_train_files <- NULL
  valid_candidates <- 0L
  
  for (i in seq_len(iterations)) {
    perm <- sample(files)                                 # random order of files
    perm_counts <- file_summary$n[match(perm, files)]    # counts in that order
    cum_counts <- cumsum(perm_counts)
    k <- which.min(abs(cum_counts - target))             # choose cut closest to target
    n_train <- cum_counts[k]
    
    if (abs(n_train - target) > tol * total_n) next
    
    train_files <- perm[1:k]
    train_df <- df[df[[source_col]] %in% train_files, , drop = FALSE]
    val_df   <- df[!df[[source_col]] %in% train_files, , drop = FALSE]
    
    ### ----------------------------
    ### Main balance: Study.Time.Collected
    ### ----------------------------
    train_bal <- suppressWarnings(as.numeric(train_df[[balance_col]]))
    val_bal   <- suppressWarnings(as.numeric(val_df[[balance_col]]))
    train_bal <- train_bal[!is.na(train_bal)]
    val_bal   <- val_bal[!is.na(val_bal)]
    
    if (length(train_bal) < 2 || length(val_bal) < 2) {
      main_score <- 0
    } else {
      pooled_sd <- sqrt(((length(train_bal)-1)*var(train_bal) + (length(val_bal)-1)*var(val_bal)) /
                          (length(train_bal) + length(val_bal) - 2))
      if (is.na(pooled_sd) || pooled_sd == 0) pooled_sd <- 1
      main_score <- abs(mean(train_bal) - mean(val_bal)) / pooled_sd +
        abs(sd(train_bal) - sd(val_bal)) / pooled_sd
    }
    
    ### ----------------------------
    ### Secondary balance: Race, Gender, Age
    ### ----------------------------
    # categorical imbalance
    total_cat_score <- 0
    for (cc in cat_cols) {
      levs <- union(unique(train_df[[cc]]), unique(val_df[[cc]]))
      p_train <- prop.table(table(factor(train_df[[cc]], levels = levs)))
      p_val   <- prop.table(table(factor(val_df[[cc]],   levels = levs)))
      total_cat_score <- total_cat_score + sum(abs(as.numeric(p_train) - as.numeric(p_val)))
    }
    
    # numeric imbalance for Age
    total_num_score <- 0
    for (nc in num_cols) {
      train_num <- suppressWarnings(as.numeric(train_df[[nc]]))
      val_num   <- suppressWarnings(as.numeric(val_df[[nc]]))
      train_num <- train_num[!is.na(train_num)]
      val_num   <- val_num[!is.na(val_num)]
      if (length(train_num) >= 2 && length(val_num) >= 2) {
        pooled_sd <- sqrt(((length(train_num)-1)*var(train_num) + (length(val_num)-1)*var(val_num)) /
                            (length(train_num) + length(val_num) - 2))
        if (is.na(pooled_sd) || pooled_sd == 0) pooled_sd <- 1
        total_num_score <- total_num_score +
          abs(mean(train_num) - mean(val_num)) / pooled_sd
      }
    }
    
    # penalty for deviating from exact train fraction
    sample_penalty <- abs(n_train/total_n - train_frac)
    
    ### ----------------------------
    ### Weighted total score
    ### ----------------------------
    score <- (main_score * 5) +          # strong weight on Study.Time.Collected
      (total_cat_score * 1) +     # lighter weight on Race & Gender
      (total_num_score * 1) +     # lighter weight on Age
      (sample_penalty * 10)
    
    if (score < best_score) {
      best_score <- score
      best_train_files <- train_files
    }
    valid_candidates <- valid_candidates + 1L
  }
  
  if (is.null(best_train_files)) {
    stop("No valid split found within tolerance. Try increasing `tol` or `iterations`.")
  }
  
  train_df <- df[df[[source_col]] %in% best_train_files, , drop = FALSE]
  val_df   <- df[!df[[source_col]] %in% best_train_files, , drop = FALSE]
  assignment <- data.frame(
    source_file = files,
    assigned = ifelse(files %in% best_train_files, "train", "val"),
    n = file_summary$n,
    stringsAsFactors = FALSE
  )
  
  cat("Split chosen:\n")
  cat(" - train samples:", nrow(train_df), "\n")
  cat(" - val   samples:", nrow(val_df), "\n")
  cat(" - achieved train fraction:", round(nrow(train_df)/total_n, 4), "\n")
  cat(" - best score:", round(best_score, 4), "\n")
  cat(" - valid candidates evaluated:", valid_candidates, " (of", iterations, ")\n\n")
  
  return(list(train = train_df, val = val_df,
              assignment = assignment, score = best_score))
}


### apply the splitting

res <- split_by_sourcefile(filtered_kinetic_data,
                           source_col = "source_file",
                           cat_cols = c("Race","Gender"),
                           num_cols = c("Subject.Age"),
                           balance_col = "Study.Time.Collected",
                           train_frac = 0.75,
                           iterations = 10000,
                           tol = 0.05,
                           seed = 937)

train_df <- res$train
val_df   <- res$val

set.seed(937)
data <- train_df %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86,subfolder,V2,Race,Race.Specify,Ethnicity,ARM.Accession,Expsample.Accession,Original.File.Name)) %>% 
  filter(!(sample=="NA")) %>%
  rename(dpi = Study.Time.Collected) %>%
  rename(group = ARM.Name) %>%
  rename(datatype = File.Detail) %>%
  rename(Age = Subject.Age) %>%
  rename(StudyID = Study.Accession) %>%
  rename(SubjectID = Subject.Accession) %>%
  rename(organ = Biosample.Type) %>%
  mutate(across(3:10, ~ ifelse(.x < -2 | .x > 5, NA, .x)))

#use multiple imputation (with random forest) to complete missing training data
dataformerging <- data %>% 
  select(sample,count,source_file,SubjectID, Ki67)
dataforimp <- data %>% 
  select(-c(count,source_file,SubjectID, Ki67))

iter=20
imp2025 <- mice(dataforimp, method = "rf", m = iter, seed = 937)

#get all imputed data
data.imputed <- list()
for(m in 1:20){
  data.imputed[[m]] <- complete(imp2025, action = m)
} 
for(m in 1:20){  
  data.imputed[[m]] <- cbind(data.imputed[[m]], dataformerging)
}
#get average imputed values for prediction markers
prediction_marker<-c(colnames(dataforimp[,c(2:8)]))
data.imputed_mean <- data.imputed[[1]][,c(1,15,11)]
for (o in 1:length(prediction_marker)) {
  average_marker <- matrix(0, nrow(data.imputed[[1]]))
  for (m in 1:length(data.imputed)) {
    average_marker <- average_marker + data.imputed[[m]][,o+1]
  }
  average_marker <- average_marker / length(data.imputed)
  data.imputed_mean <- cbind(data.imputed_mean,average_marker)
  colnames(data.imputed_mean)[o+3] <- paste(prediction_marker[o])
}
prediction_marker_age<-c(colnames(data.imputed[[1]][,c(2:8)]),"Age")

#turn the MIDS into a list to bundle imputed data for training the prediction model
data.full <- list()
for(m in 1:iter){
  data.full[[m]] <- complete(imp2025, action = m) %>%
    left_join(dataformerging, by = c("sample")) %>% 
    mutate(group = as.factor(group)) %>% 
    filter(dpi>4)
}

# Define a function to make and average predictions
average_predictions <- function(data, models) {
  # Apply predictions using each model in the list
  predictions_list <- map(models, ~ predict(.x, newdata = data))
  
  # Bind all predictions into a dataframe
  predictions_df <- bind_cols(predictions_list)
  
  # Calculate the row-wise average of all predictions
  predictions_df <- predictions_df %>%
    mutate(p = rowMeans(across(everything()))) %>%
    mutate(ID = data$sample,dpi = data$dpi) %>%
    select(ID,dpi,p)
  
  return(predictions_df)
}

#train the final model
rf_ascmi <- data.full %>% 
  map(~
        randomForest(dpi ~ CD19+CD20+CD28+CD45+CD56+CD138+HLA.DR, 
                     data = ., 
                     ntree = 100,
                     mtry = 6,
                     importance = TRUE, 
                     proximity =TRUE
        ))

#extract importance of predictors
df_list <- rf_ascmi %>%  map(~.$importance)
df_list <- lapply(df_list, as.data.frame)
add_id_variable <- function(df) {
  df <- df %>% 
    mutate(ID = row.names(.))  # Create "ID" variable with row names
  return(df)
}
df_list <- lapply(df_list, add_id_variable)
variable_importance <- df_list %>% 
  bind_rows() %>% 
  group_by(ID) %>%
  summarise_all(mean, na.rm = TRUE)

#use model on imputed data and average out predictions
prediction_all <- list()
for (m in 1:length(data.imputed)) {
  prediction_all[[m]] <- average_predictions(data.imputed[[m]], rf_ascmi)
}
mean_predictions <- matrix(0, nrow(prediction_all[[1]]), ncol(prediction_all[[1]]))
for (m in 1:length(prediction_all)) {
  mean_predictions <- mean_predictions + prediction_all[[m]][,3]
}
mean_predictions <- mean_predictions / length(prediction_all)
predictions_mean <- cbind(dataformerging,dataforimp[,c(12,13,11)], prediction_all[[1]][, 1:2], mean_predictions)

# trim output table
model_output<-predictions_mean %>%
  select(sample, StudyID,SubjectID, organ, group, dpi, p) %>%
  mutate(Age=dataforimp$Age) %>%
  mutate(Gender=dataforimp$Gender)
model_output$Age<-as.numeric(model_output$Age)

###
### check performance of model and biological validation

#add average imputed values for prediction markers
data.imputed_mean<-data.imputed_mean %>%
  mutate(predict=model_output$p)

# correlate actual and predicted dpi
model_HC_kinetic<-model_output %>%
  filter(dpi>4)

r_HC_kinetic<-cor.test(model_HC_kinetic$dpi, model_HC_kinetic$p)

plot_HC_kinetic<-ggplot(model_HC_kinetic, aes(x=dpi, y=p)) +
  geom_point(shape=1) + #none-filled circles
  geom_smooth(method = "lm", color = "red", se = TRUE , level = 0.99) + #add regression line /w confidence interval
  labs(title = paste("pearson correlation, r =",format(r_HC_kinetic$estimate, digits=3)," p =",format(r_HC_kinetic$p.value, digits=3)),
       x = "actual DPI",
       y = "predicted DPI") +
  theme_bw() +
  #xlim(5,250)+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))

ggsave(filename = "output/train_noKi67.png")

# continue with validation

data_new <- val_df %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86,subfolder,V2,Race,Race.Specify,Ethnicity,ARM.Accession,Expsample.Accession,Original.File.Name)) %>% 
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
  select(sample,count,source_file,SubjectID, Ki67)
new_imp <- data_new %>% 
  select(-c(count,source_file,SubjectID, Ki67))

percent_miss<-(sum(is.na(new_imp)) / (nrow(new_imp) * ncol(new_imp))) *100

iter=20
data.norm<-cbind(data.imputed_mean[,c(1,4:10)],dataforimp[,c(9:15)])
data.total<-rbind(data.norm,new_imp)
imp.new <- mice(data.total, method = "rf", m = iter, seed = 937)

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
predictions_new <- cbind(new_merging,new_imp[,c(12,13,11)], prediction_new[[1]][, 1:2], new_predictions)

# trim output table
new_output<-predictions_new %>%
  select(sample, SubjectID, group, dpi, p) %>%
  mutate(Age=new_imp$Age) %>%
  mutate(Gender=new_imp$Gender)
new_output$Age<-as.numeric(new_output$Age)

val_data<-new_output %>%
  filter(dpi>4)

r_ex_kinetic<-cor.test(val_data$dpi, val_data$p, method = "pearson")
val_plot_scatter<-ggplot(val_data, aes(x=dpi, y=p)) +
  geom_point(shape=1) + #none-filled circles
  geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
  labs(title = paste("pearson r =",format(r_ex_kinetic$estimate, digits=3)," p =",format(r_ex_kinetic$p.value, digits=3)),
       x = "actual DPI",
       y = "predicted DPI") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))

ggsave(filename = "output/val_plot_noKi67.png")

###: evaluate the trained model 
#
##extract model characteristics
tree_depth_min<-list()
for (f in 1:iter) {
  tree_depth_min[[f]]<-min_depth_distribution(rf_ascmi[[f]])
}
all_tree_depth <- do.call(rbind, tree_depth_min)
marker_depth<-cbind(data.frame(CD138 = all_tree_depth[seq(1, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD19 = all_tree_depth[seq(2, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD20 = all_tree_depth[seq(3, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD28 = all_tree_depth[seq(4, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD45 = all_tree_depth[seq(5, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD56 = all_tree_depth[seq(6, nrow(all_tree_depth), by = 7), 3]))
marker_depth<-cbind(marker_depth,data.frame(HLA = all_tree_depth[seq(7, nrow(all_tree_depth), by = 7), 3]))
#marker_depth<-cbind(marker_depth,data.frame(Ki67 = all_tree_depth[seq(8, nrow(all_tree_depth), by = 8), 3]))
#
forest_importance<-list()
for (f in 1:iter) {
  forest_importance[[f]]<-measure_importance(rf_ascmi[[f]])
}
all_forest_importance <- do.call(rbind, forest_importance)
write_xlsx(marker_depth, "output/marker_depth_noKi67.xlsx")
write_xlsx(all_forest_importance, "output/all_forest_importance_noKi67.xlsx")

#generate min_tree_depth distribution heatmap for Figure 4A

md_long <- marker_depth %>%
  mutate(sample_id = row_number()) %>%        # keep original row id if desired
  pivot_longer(cols = -sample_id, names_to = "Marker", values_to = "Depth") %>%
  group_by(Marker) %>%
  arrange(Depth, .by_group = TRUE) %>%        # sort values low -> high within each marker
  mutate(rank_within_marker = row_number()) %>%
  ungroup()

md_long$Marker <- factor(md_long$Marker,
                         levels = c("CD20", "CD28", "HLA", "CD19", "CD138", "CD45", "CD56"))

# Make Depth a factor so we can use a discrete palette (0..8)
md_long$Depth <- factor(md_long$Depth, levels = 0:8)

# Basic plot
{
  p <- ggplot(md_long, aes(x = rank_within_marker, y = Marker, fill = Depth)) +
    geom_tile(width = 1, height = 1, color = NA) +
    geom_vline(xintercept = seq(200, 1800, by = 200),
               linetype = "dotted", size = 1.5, color = "black", alpha = 0.8) +
    scale_fill_viridis_d(option = "turbo", direction = -1, begin = 0, end = 1,
                         name = "minimal tree depth of marker:") +
    scale_x_continuous(
      name = "number of trees within random forest model",
      expand = c(0, 0),
      breaks = seq(0, 2000, by = 200),
      limits = c(-6, 2004),
      position = "top"
    ) +
    scale_y_discrete(name = NULL, expand = c(0, 0), position = "right") +
    theme_bw(base_rect_size = 4, base_line_size = 1.5) +
    theme(
      title = element_text(size=15, face = "bold", colour = "black",hjust = 1.06),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0, "pt"),
      axis.text.x = element_text(size = 15, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 14, face = "bold", colour = "white"),
      axis.ticks.length = grid::unit(10, "pt"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.key.size = grid::unit(0.9, "lines"),
      legend.text = element_text(size = 14,face = "bold", colour = "black"),
      legend.title = element_text(size = 14,face = "bold", colour = "black"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"))
  print(p)
}

ggsave("output/tree_depth_noKi67.png", plot = p)

#generate importance overview for markers in Fig 4B

{
  mean_importance<-all_forest_importance[c(1:8),-c(6,8)]
  colnames(mean_importance)[colnames(mean_importance)=="mean_min_depth"] <- "mean depth"
  colnames(mean_importance)[colnames(mean_importance)=="no_of_nodes"] <- "tree nodes (%)"
  colnames(mean_importance)[colnames(mean_importance)=="mse_increase"] <- "MSE increase"
  colnames(mean_importance)[colnames(mean_importance)=="node_purity_increase"] <- "node purity"
  colnames(mean_importance)[colnames(mean_importance)=="times_a_root"] <- "tree root (%)"
  for (a in 1:8) {
    average<-all_forest_importance %>% filter(variable==paste(mean_importance[a,1]))
    mean_importance[a,2]<-mean(average$mean_min_depth)
    mean_importance[a,3]<-(sum(average$no_of_nodes)/sum(all_forest_importance$no_of_nodes)*100)
    mean_importance[a,4]<-mean(average$mse_increase)
    mean_importance[a,5]<-mean(average$node_purity_increase)
    mean_importance[a,6]<-mean(average$times_a_root)
  }
  mean_importance[,c(2:6)]<-round(mean_importance[,c(2:6)],2)
  mean_importance[,c(5)]<-round(mean_importance[,c(5)],0)
}

mean_importance_scale<-mean_importance
for (a in 1:8) {
  mean_importance_scale[a,2]<-(1-(mean_importance[a,2]-min(mean_importance$`mean depth`))/(max(mean_importance$`mean depth`)-min(mean_importance$`mean depth`)))
  mean_importance_scale[a,3]<-(mean_importance[a,3]-min(mean_importance$`tree nodes (%)`))/(max(mean_importance$`tree nodes (%)`)-min(mean_importance$`tree nodes (%)`))
  mean_importance_scale[a,4]<-(mean_importance[a,4]-min(mean_importance$`MSE increase`))/(max(mean_importance$`MSE increase`)-min(mean_importance$`MSE increase`))
  mean_importance_scale[a,5]<-(mean_importance[a,5]-min(mean_importance$`node purity`))/(max(mean_importance$`node purity`)-min(mean_importance$`node purity`))
  mean_importance_scale[a,6]<-(mean_importance[a,6]-min(mean_importance$`tree root (%)`))/(max(mean_importance$`tree root (%)`)-min(mean_importance$`tree root (%)`))
}

mis_long <- mean_importance_scale %>%
  pivot_longer(
    cols = -variable,
    names_to = "Metric",
    values_to = "Value"
  )

mis_long2 <- mean_importance %>%
  pivot_longer(
    cols = -variable,
    names_to = "Metric",
    values_to = "Value"
  )

mis_join<-cbind(mis_long, raw=mis_long2$Value)

mis_join$variable <- factor(mis_join$variable, 
                            levels = c("CD20", "CD28", "HLA.DR", "CD19", "CD138", "CD45", "CD56"))

p2 <- ggplot(mis_join, aes(x = Metric, y = variable, fill = Value)) +
  geom_tile(width = 1, height = 1, color = NA) +  # add subtle borders
  geom_text(aes(x = Metric,y = variable, label = raw), color = "black", size = 5, inherit.aes = TRUE, fontface = "bold") +
  scale_fill_gradient(low = "#FAD6D6", high = "red3",name = "low importance                                      high importance") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw(base_rect_size = 4, base_line_size = 1.5) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = grid::unit(0, "pt"),
    axis.ticks.length = grid::unit(10, "pt"),
    axis.text.x = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black",hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.align = 0.5,
    legend.key.size = grid::unit(5, "lines"),
    legend.text = element_text(size = 0,face = "bold", colour = "white"),
    legend.title = element_text(size = 14,face = "bold", colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  guides(fill = guide_colorbar(title.position = "top",barwidth = 25,barheight = 1))

print(p2)
ggsave("output/importance_noKi67.png", plot = p2)

###
# organ application

data_new<- joined_data %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86,subfolder,V2,Race,Race.Specify,Ethnicity,ARM.Accession,Expsample.Accession,Original.File.Name)) %>% 
  select(-Ki67) %>%
  filter(!(sample=="NA")) %>%
  rename(dpi = Study.Time.Collected) %>%
  rename(group = ARM.Name) %>%
  rename(datatype = File.Detail) %>%
  rename(Age = Subject.Age) %>%
  rename(StudyID = Study.Accession) %>%
  rename(SubjectID = Subject.Accession) %>%
  rename(organ = Biosample.Type) %>%
  mutate(across(3:9, ~ ifelse(.x < -2 | .x > 5, NA, .x))) %>%
  filter(count>49) %>%
  filter((StudyID=="NA")|(StudyID=="Festen")|(StudyID=="JJS")|(StudyID=="SDY1389")|(StudyID=="SDY2107")
         |(StudyID=="Z24H")|(StudyID=="Z3EP")|(StudyID=="Z4JJ")|(StudyID=="Z4JK")|(StudyID=="Z68Q")
         |(StudyID=="Z7A5")|(StudyID=="ZYQ9"))

data_organ<-data_new %>%
  filter((organ=="gut_ileum")|(organ=="gut_rectum")|(organ=="gut_sigmoid")|(organ=="gut_cdescen")
         |(organ=="Jejunum")|(organ=="BM")|(organ=="Bone Marrow")) %>%
  mutate(across(organ, ~ gsub("gut_ileum", "gut",.))) %>%
  mutate(across(organ, ~ gsub("gut_rectum", "gut",.))) %>%
  mutate(across(organ, ~ gsub("gut_sigmoid", "gut",.))) %>%
  mutate(across(organ, ~ gsub("gut_cdescen", "gut",.))) %>%
  mutate(across(organ, ~ gsub("Jejunum", "gut",.))) %>%
  mutate(across(organ, ~ gsub("Bone Marrow", "BM",.)))

new_merging <- data_organ %>% 
  select(sample,count,source_file,SubjectID,Age,dpi)
new_imp <- data_organ %>% 
  select(-c(count,source_file,SubjectID,Age,dpi))

iter=20
data.total<-rbind(data.norm[,-c(10,14)], new_imp)
imp.new <- mice(data.total, method = "rf", m = iter, seed = 937)

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
predictions_new <- cbind(new_merging,new_imp[,c(12,13,11)], prediction_new[[1]][, 1:2], new_predictions)

# trim output table
new_output<-predictions_new[,-11] #%>%  select(sample, StudyID,SubjectID, organ, group, dpi, p)

comp_BLD_BM_gut<-rbind(model_output[,-c(5,9)] %>% 
                         filter(dpi<2) %>%
                         mutate(across(organ, ~ gsub("PBMC", "blood",.))) %>%
                         mutate(across(organ, ~ gsub("Other", "blood",.))) %>%
                         mutate(across(organ, ~ gsub("Whole blood", "blood",.))) %>%
                         mutate(across(organ, ~ gsub("B cell", "blood",.)))
                       ,new_output[,c(1,9,4,7,6,11,5)])

ANOVA_imp<- lm(comp_BLD_BM_gut[[6]] ~ comp_BLD_BM_gut[[4]])  |> aov() |> TukeyHSD()
plot_HC_organ<-ggplot(comp_BLD_BM_gut, aes(x = organ, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(comp_BLD_BM_gut[[2]])), width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  annotate("text", x="BM", y=max(comp_BLD_BM_gut$p)*1.05, label=paste("p=",formatC(ANOVA_imp[["comp_BLD_BM_gut[[4]]"]][1,4], format = "e", digits=2)), size=5)+
  annotate("text", x="gut", y=max(comp_BLD_BM_gut$p)*1.05, label=paste("p=",formatC(ANOVA_imp[["comp_BLD_BM_gut[[4]]"]][2,4], format = "e", digits=2)), size=5)
print(plot_HC_organ) 

ggsave(filename = paste("output/organ_comp_noKi67.png", sep = ""))

