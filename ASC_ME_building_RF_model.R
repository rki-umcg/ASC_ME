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
# this scripts contains the setup strategies and codes for the random forest ASC maturity index model
# raw data from flow cytometry files were previously analyzed, processed and imported below 
# 

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

###Part 1: adjusting data and train the random forest prediction modedl
#
#set work direction and load training datasets
setwd(".../ASC_ME")

###
set.seed(555)
data <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
  filter(!(sample=="NA")) %>%
  filter(((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>%
  filter((StudyID=="AID")|(StudyID=="CER")|(StudyID=="SDY144")|(StudyID=="SDY180")|(StudyID=="SDY224")
         |(StudyID=="SDY272")|(StudyID=="SDY296")|(StudyID=="SDY301")|(StudyID=="SDY364")|(StudyID=="SDY387")
         |(StudyID=="SDY522")|(StudyID=="SDY80")|(StudyID=="SDY819")|(StudyID=="SDY984")|(StudyID=="SDY1086")
         |(StudyID=="SDY1288")|(StudyID=="SDY1669")|(StudyID=="SDY1734")|(StudyID=="ZY77")
  )

#use multiple imputation (with random forest) to complete missing training data
dataformerging <- data %>% 
  select(sample,check,SubjectID,File,condition)
dataforimp <- data %>% 
  select(-c(check,SubjectID,File,condition))
print(percent_miss<-(sum(is.na(dataforimp)) / (nrow(dataforimp) * ncol(dataforimp))) *100)
iter=20
imp0911 <- mice(dataforimp, method = "rf", m = iter, seed = 555)
#get all imputed data
data.imputed <- list()
for(m in 1:20){
  data.imputed[[m]] <- complete(imp0911, action = m)
} 
for(m in 1:20){  
  data.imputed[[m]] <- cbind(data.imputed[[m]], dataformerging)
}
#get average imputed values for prediction markers
prediction_marker<-c(colnames(dataforimp[,c(3:10)]))
data.imputed_mean <- data.imputed[[1]][,c(1,2,14)]
for (o in 1:length(prediction_marker)) {
  average_marker <- matrix(0, nrow(data.imputed[[1]]))
  for (m in 1:length(data.imputed)) {
    average_marker <- average_marker + data.imputed[[m]][,o+2]
  }
  average_marker <- average_marker / length(data.imputed)
  data.imputed_mean <- cbind(data.imputed_mean,average_marker)
  colnames(data.imputed_mean)[o+3] <- paste(prediction_marker[o])
}
prediction_marker_age<-c(colnames(data.imputed[[1]][,c(3:10)]),"Age")

#turn the MIDS into a list to bundle imputed data for training the prediction model
data.full <- list()
for(m in 1:iter){
  data.full[[m]] <- complete(imp0911, action = m) %>%
    left_join(dataformerging, by = c("sample")) %>% 
    mutate(group = as.factor(group)) %>% 
    filter(dpi>4 & ((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>%
    filter(!(StudyID=="SDY648")) %>%
    filter(!(StudyID=="SDY739"))%>%
    filter(!(StudyID=="SDY622a"))  %>%
    filter(!(StudyID=="SDY622b"))  %>%
    filter(!(StudyID=="SDY368"))  %>%
    filter(!(StudyID=="SDY34"))  %>%
    filter(!(StudyID=="SDY1397")) 
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
set.seed(555)
rf_ascmi <- data.full %>% 
  map(~
        randomForest(dpi ~ CD19+CD20+CD28+CD45+CD56+CD138+HLA.DR+Ki67, 
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
predictions_mean <- cbind(dataformerging,dataforimp[,c(11:12,15)], prediction_all[[1]][, 1:2], mean_predictions)
# trim output table
model_output<-predictions_mean %>%
  select(sample, StudyID,SubjectID, organ, group, dpi, p) %>%
  mutate(Age=dataforimp$Age) %>%
  mutate(Gender=dataforimp$Gender)
model_output$Age<-as.numeric(model_output$Age)

###
###
### check performance of model and biological validation
library(ggplot2)
library(viridis)
library(dunn.test)
library(gridExtra)

#add average imputed values for prediction markers
data.imputed_mean<-data.imputed_mean %>%
    mutate(predict=model_output$p)

# correlate actual and predicted dpi
model_HC_kinetic<-model_output %>%
  filter(dpi>4 & ((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>% 
  #filter(dpi<200) %>%
  filter(!(StudyID=="SDY648")) %>%
  filter(!(StudyID=="SDY739"))%>%
  filter(!(StudyID=="SDY622a"))  %>%
  filter(!(StudyID=="SDY622b"))  %>%
  filter(!(StudyID=="SDY368"))  %>%
  filter(!(StudyID=="SDY34"))  %>%
  filter(!(StudyID=="SDY1397"))  %>%
  filter(!(StudyID=="SDY788"))
  
r_HC_kinetic<-cor.test(model_HC_kinetic$dpi, model_HC_kinetic$p)

plot_HC_kinetic<-ggplot(model_HC_kinetic, aes(x=dpi, y=p)) +
  geom_point(shape=1) + #none-filled circles
  #geom_smooth(method = "gam", n=10, color = "red", se = TRUE , level = 0.99) + #add regression line /w confidence interval
  geom_smooth(method = "lm", color = "red", se = TRUE , level = 0.99) + #add regression line /w confidence interval
  labs(title = paste("pearson correlation, r =",format(r_HC_kinetic$estimate, digits=3)," p =",format(r_HC_kinetic$p.value, digits=3)),
       x = "actual DPI",
       y = "predicted DPI") +
  theme_bw() +
  xlim(5,250)+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))

print(plot_HC_kinetic)
ggsave("plot_HC_kinetic.png", plot = plot_HC_kinetic)
## -> Figure 3A

# check longitudinal/kinetic samples
model_HC_long<-model_output %>%
  filter((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT")) %>%
  filter(!(StudyID=="SDY648")) %>%
  filter(!(StudyID=="SDY739"))%>%
  filter(!(StudyID=="SDY622a"))  %>%
  filter(!(StudyID=="SDY622b"))  %>%
  filter(!(StudyID=="SDY368"))  %>%
  filter(!(StudyID=="SDY34"))  %>%
  filter(!(StudyID=="SDY1397")) %>%
  filter(dpi==0|dpi>4)

#Identify non-unique SubjectID values in model_HC_long
non_unique_subjects <- model_HC_long %>%
  count(SubjectID) %>%  # Count occurrences of each SubjectID
  filter(n > 1) %>%  # Keep only SubjectID with more than 1 occurrence
  pull(SubjectID)  # Extract the SubjectID values
#filter out non-unique samples
model_HC_long <- model_HC_long %>%
  filter(SubjectID %in% non_unique_subjects)

plot_HC_long<-ggplot(model_HC_long, aes(x=dpi, y=p, group=SubjectID, color=SubjectID)) +
  geom_point(shape=1) + #none-filled circles
  geom_line(size=0.25) + #, aes(color=SubjectID)) +
  scale_color_viridis(discrete = TRUE, option = "H") +
  labs(title = "linked longitudinal data points",
       x = "actual DPI",
       y = "predicted DPI") +
  xlim(0,250)+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))
print(plot_HC_long)
ggsave("plot_HC_long.png", plot = plot_HC_long)
## -> Figure 3B

#by dpi groups
model_HC_dpi_gp<-model_output %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl", "early", "intermediate", "late", "very late"))) %>%# Define the labels for the age groups
  filter(!(dpi_group=="excl"))

ANOVA_dev<- lm(model_HC_dpi_gp[[7]] ~ model_HC_dpi_gp[[10]])  |> aov() |> TukeyHSD()

plot_HC_dpi_gp<-ggplot(model_HC_dpi_gp, aes(x = dpi_group, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(model_HC_dpi_gp[[2]])), width = 0.25, alpha = 1, size = 2,shape=1) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "immune response stage",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  #scale_y_continuous(breaks = c(0,30,45,60,90,120), limits = c(4,116)) +
  ylim(5,145) +
  xlim("baseline", "early", "intermediate", "late", "very late") +
  annotate("text", x="baseline", y=145, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = 144, xend = 1.6, yend = 144),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x="early", y=145, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][1,4], digits=3)), size=5)+
  annotate("text", x="intermediate", y=145, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][2,4], digits=3)), size=5)+
  annotate("text", x="late", y=145, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][3,4], digits=3)), size=5)+
  annotate("text", x="very late", y=145, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][4,4], digits=3)), size=5)+
  annotate("text", x="baseline", y=135, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][4,4], digits=3)), size=5)+
  annotate("text", x="early", y=135, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][7,4], digits=3)), size=5)+
  annotate("text", x="intermediate", y=135, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][9,4], digits=3)), size=5)+
  annotate("text", x="late", y=135, label=paste("p=",format(ANOVA_dev[["model_HC_dpi_gp[[10]]"]][10,4], digits=3)), size=5)+
  geom_segment(aes(x = 4.6, y = 134, xend = 4.4, yend = 134),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x="very late", y=135, label=paste("comparison"), size=5)
print(plot_HC_dpi_gp) 
ggsave("plot_HC_dpi_gp.png", plot = plot_HC_dpi_gp)
## -> Figure 3C

###Part 2: evaluate the trained model 
#
##extract model characteristics
tree_depth_min<-list()
for (f in 1:iter) {
  tree_depth_min[[f]]<-min_depth_distribution(rf_ascmi[[f]])
}
all_tree_depth <- do.call(rbind, tree_depth_min)
marker_depth<-cbind(data.frame(CD138 = all_tree_depth[seq(1, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD19 = all_tree_depth[seq(2, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD20 = all_tree_depth[seq(3, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD28 = all_tree_depth[seq(4, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD45 = all_tree_depth[seq(5, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(CD56 = all_tree_depth[seq(6, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(HLA = all_tree_depth[seq(7, nrow(all_tree_depth), by = 8), 3]))
marker_depth<-cbind(marker_depth,data.frame(Ki67 = all_tree_depth[seq(8, nrow(all_tree_depth), by = 8), 3]))
#
forest_importance<-list()
for (f in 1:iter) {
  forest_importance[[f]]<-measure_importance(rf_ascmi[[f]])
}
all_forest_importance <- do.call(rbind, forest_importance)
write_xlsx(marker_depth, "marker_depth.xlsx")
write_xlsx(all_forest_importance, "all_forest_importance.xlsx")

## -> rearrange data in Excel and plot heatmaps for Figure 4 with GraphPad Prism 9
