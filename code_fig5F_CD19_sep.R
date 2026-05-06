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

  ###
  parent_folder <- paste(data_path, "CD19_separation", sep = "")
  subfolders <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)
  
  # Loop through subfolders and read in the excel files
  all_data <- lapply(subfolders, function(folder) {
    file_path <- file.path(folder, "final_results.xlsx")
    
    if (file.exists(file_path)) {
      df <- read_excel(file_path)
      df[,c(3:15)] <- lapply(df[,c(3:15)], as.numeric)
      df$subfolder <- basename(folder)  # add subfolder name as a new column
      df$type <- "total.ASCs"  # add ASC type as a new column
      return(df)
    } else {
      message("No final_results.xlsx found in: ", folder)
      return(NULL)
    }
  }) %>% 
    bind_rows()  # combine into one big data frame
  
  all_data_CD19pos <- lapply(subfolders, function(folder) {
    file_path <- file.path(folder, "CD19pos_final_results.xlsx")
    
    if (file.exists(file_path)) {
      df <- read_excel(file_path)
      df[,c(3:15)] <- lapply(df[,c(3:15)], as.numeric)
      df$subfolder <- basename(folder)  # add subfolder name as a new column
      df$type <- "CD19pos.ASCs"  # add ASC type as a new column
      return(df)
    } else {
      message("No final_results.xlsx found in: ", folder)
      return(NULL)
    }
  }) %>% 
    bind_rows()  # combine into one big data frame
  
  all_data_CD19neg <- lapply(subfolders, function(folder) {
    file_path <- file.path(folder, "CD19neg_final_results.xlsx")
    
    if (file.exists(file_path)) {
      df <- read_excel(file_path)
      df[,c(3:15)] <- lapply(df[,c(3:15)], as.numeric)
      df$subfolder <- basename(folder)  # add subfolder name as a new column
      df$type <- "CD19neg.ASCs"  # add ASC type as a new column
      return(df)
    } else {
      message("No final_results.xlsx found in: ", folder)
      return(NULL)
    }
  }) %>% 
    bind_rows()  # combine into one big data frame
  
  data_merge<- rbind(all_data, all_data_CD19neg, all_data_CD19pos) %>%
    select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,CD86)) %>% 
    mutate(across(3:10, ~ ifelse(.x < -2 | .x > 5, NA, .x))) %>%
    filter(!(CD19=="NA")) %>%
    #rename(dpi = Study.Time.Collected) %>%
    #rename(group = ARM.Name) %>%
    #rename(datatype = File.Detail) %>%
    #rename(Age = Subject.Age) %>%
    #rename(StudyID = Study.Accession) %>%
    #rename(SubjectID = Subject.Accession) %>%
    #rename(organ = Biosample.Type) 
    filter(count>49)
  
  new_merging <- data_merge %>% 
    select(sample,count,subfolder,type)
  new_imp <- data_merge %>% 
    select(-c(count,subfolder,type))
  
  iter=20
  data.total<-rbind(data.norm[,c(1:9)], new_imp)
  imp.new <- mice(data.total, method = "rf", m = iter, seed = 557)
  
  new.imputed <- list()
  for(m in 1:20){
    new.imputed[[m]] <- complete(imp.new, action = m)
    new.imputed[[m]] <- cbind(new.imputed[[m]][c((nrow(data.norm)+1):nrow(new.imputed[[m]])),], new_merging[,-1])
    new.imputed[[m]] <- left_join(new.imputed[[m]], joined_data[,c(1,32)], by = "sample")
    colnames(new.imputed[[m]])[13]<-"dpi"
  }
  
  #get average imputed values for prediction markers
  
  new.imputed_mean <- new.imputed[[1]][,c(1,12,10,11)]
  for (o in 1:length(prediction_marker)) {
    new_average_marker <- matrix(0, nrow(new.imputed[[1]]))
    for (m in 1:length(new.imputed)) {
      new_average_marker <- new_average_marker + new.imputed[[m]][,o+1]
    }
    new_average_marker <- new_average_marker / length(new.imputed)
    new.imputed_mean <- cbind(new.imputed_mean,new_average_marker)
    colnames(new.imputed_mean)[o+4] <- paste(prediction_marker[o])
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
  predictions_new <- cbind(new_merging,prediction_new[[1]][, 1], new_predictions)
  
  # trim output table
  new_output<-predictions_new[,-5] 
  
  #
  # select data for plotting
  #
  
  plot_data<-new_output
  kt<- lm(plot_data[[5]] ~ plot_data[[4]])  |> aov() |> TukeyHSD()
  plot<-ggplot(plot_data, aes(x = type, y = p)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25, alpha = 1, size = 2) +
    theme_bw() +
    labs(x = NULL,
         y = "predicted DPI",
         color = "StudyID") +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold")) +
    xlim("total.ASCs", "CD19pos.ASCs", "CD19neg.ASCs")+
    annotate("text", x=1, y=max(plot_data$p)*1.2, label=paste("comparison"), size=5)+
    geom_segment(aes(x = 1.4, y = max(plot_data$p)*1.2, xend = 1.6, yend = max(plot_data$p)*1.2),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
    annotate("text", x="CD19pos.ASCs", y=max(plot_data$p)*1.2, label=paste("p=",format(kt[["plot_data[[4]]"]][3,4], digits=3)), size=5)+
    annotate("text", x="CD19neg.ASCs", y=max(plot_data$p)*1.2, label=paste("p=",format(kt[["plot_data[[4]]"]][2,4], digits=3)), size=5)+
    annotate("text", x=3, y=max(plot_data$p)*1.1, label=paste("comparison"), size=5)+
    geom_segment(aes(x = 2.6, y = max(plot_data$p)*1.1, xend = 2.4, yend = max(plot_data$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
    annotate("text", x=2, y=max(plot_data$p)*1.1, label=paste("p=",format(kt[["plot_data[[4]]"]][1,4], digits=3)), size=5)
  ggsave(filename = "output/Fig5F_plot_all.data_CD19.png", plot = plot,
         width = 3000, height = 1920, units = "px")

## plot CD19 + and - by dpi group
merge_dpi<-left_join(plot_data, joined_data[,c(1,32)], by = "sample")
plot_sep<-merge_dpi %>% 
    filter(type=="CD19pos.ASCs") %>%
    filter(Study.Time.Collected<181) %>%
    mutate(across(Study.Time.Collected, ~ ifelse(.x < 0, 0, .x)))
comp_CD19pos_group<-plot_sep %>%
  mutate(dpi_group = cut(Study.Time.Collected,  
                         breaks = c(-Inf, 1, 4.9, 10, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("BL_pos", "excl", "VE_pos", "E_pos", "I_pos", "L_pos", "VL_pos"))) %>%# Define the labels for the age groups
  filter(!(dpi_group=="excl"))
plot_sep2<-merge_dpi %>% 
  filter(type=="CD19neg.ASCs") %>%
  filter(Study.Time.Collected<181) %>%
  mutate(across(Study.Time.Collected, ~ ifelse(.x < 0, 0, .x)))
comp_CD19neg_group<-plot_sep2 %>%
  mutate(dpi_group = cut(Study.Time.Collected,  
                         breaks = c(-Inf, 1, 4.9, 10, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("BL_neg", "excl", "VE_neg", "E_neg", "I_neg", "L_neg", "VL_neg"))) %>%# Define the labels for the age groups
  filter(!(dpi_group=="excl"))

plot_comp.data<-rbind(comp_CD19pos_group, comp_CD19neg_group)
plot_comp.data$group_type <- ifelse(grepl("_pos$", plot_comp.data$dpi_group),
                                    "CD19+ ASCs", "CD19- ASCs")
plot_comp <- ggplot(plot_comp.data, aes(x = dpi_group, y = p, color = group_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2, shape = 1) +
  theme_bw() +
  labs(
    x = "immune response stage",
    y = "predicted DPI",
    color = "Group type"
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_color_manual(values = c("CD19+ ASCs" = "blue", "CD19- ASCs" = "red")) + 
  ylim(5, 95) +
  xlim("BL_pos", "BL_neg", "VE_pos", "VE_neg", "E_pos", "E_neg",
       "I_pos", "I_neg", "L_pos", "L_neg", "VL_pos", "VL_neg")
ggsave("output/Fig5F_plot_CD19_kinetic.png",
       width = 3000, height = 1920, units = "px")

kt<- kruskal.test(plot_comp.data[[5]] ~ plot_comp.data[[7]])
dt <- dunn.test(plot_comp.data[[5]], plot_comp.data[[7]], method = "bonferroni", list = TRUE)

