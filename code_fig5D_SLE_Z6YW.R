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
  filter((StudyID=="Z6YW")|(StudyID=="AID"))

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

plot_data<-new_output %>%
  filter(!(SubjectID=="anchorstim")) %>%
  filter((str_detect(sample, "_A_"))) %>% #|(str_detect(sample, "PSL"))|(str_detect(sample, "CTL"))|(str_detect(sample, "VC"))) %>%
  filter(dpi==0) %>%
  mutate(ID = SubjectID) %>%
  mutate(across(ID, ~ gsub("525","HC",.))) %>%
  mutate(across(ID, ~ gsub("323","LN",.))) %>%
  mutate(across(ID, ~ gsub("526","HC",.))) %>%
  mutate(across(ID, ~ gsub("328","SLE",.))) %>%
  mutate(across(ID, ~ gsub("326","SLE",.))) %>%
  mutate(across(ID, ~ gsub("523","HC",.))) %>%
  mutate(across(ID, ~ gsub("330","LN",.)))  %>%
  mutate(across(ID, ~ gsub("528","HC",.))) %>%
  mutate(across(ID, ~ gsub("336","LN",.)))  %>%
  mutate(across(ID, ~ gsub("504","HC",.))) %>%
  mutate(across(ID, ~ gsub("303","SLE",.))) %>%
  mutate(across(ID, ~ gsub("501","HC",.))) %>%
  mutate(across(ID, ~ gsub("508","HC",.))) %>%
  mutate(across(ID, ~ gsub("516","HC",.))) %>%
  mutate(across(ID, ~ gsub("311","SLE",.))) %>%
  mutate(across(ID, ~ gsub("506","HC",.))) %>%
  mutate(across(ID, ~ gsub("514","HC",.))) %>%
  mutate(across(ID, ~ gsub("307","SLE",.))) %>%
  mutate(across(ID, ~ gsub("512","HC",.))) %>%
  mutate(across(ID, ~ gsub("306","SLE",.))) %>%
  mutate(across(ID, ~ gsub("308","SLE",.))) %>%
  mutate(across(ID, ~ gsub("507","HC",.))) %>%
  mutate(across(ID, ~ gsub("515","HC",.))) %>%
  mutate(across(ID, ~ gsub("505","HC",.))) %>%
  mutate(across(ID, ~ gsub("513","HC",.))) %>%
  mutate(across(ID, ~ gsub("316","SLE",.))) %>%
  mutate(across(ID, ~ gsub("517","HC",.))) %>%
  mutate(across(ID, ~ gsub("317","LN",.)))  %>%
  mutate(across(ID, ~ gsub("518","HC",.))) %>%
  mutate(across(ID, ~ gsub("319","SLE",.))) %>%
  mutate(across(ID, ~ gsub("321","SLE",.))) %>%
  mutate(across(ID, ~ gsub("322","SLE",.))) %>%
  mutate(across(ID, ~ gsub("327","LN",.)))  %>%
  mutate(across(ID, ~ gsub("329","LN",.)))  %>%
  mutate(across(ID, ~ gsub("524","HC",.)))   %>%
  mutate(across(ID, ~ gsub("527","HC",.)))  %>%
  mutate(across(ID, ~ gsub("509","HC",.)))  %>%
  mutate(across(ID, ~ gsub("519","HC",.)))  %>%
  mutate(across(ID, ~ gsub("521","HC",.)))  %>%
  mutate(across(ID, ~ gsub("522","HC",.)))  %>%
  mutate(across(ID, ~ gsub("302","SLE",.)))  %>%
  mutate(across(ID, ~ gsub("305","SLE",.)))  %>%
  mutate(across(ID, ~ gsub("309","SLE",.)))  %>%
  mutate(across(ID, ~ gsub("310","SLE",.)))  %>%
  mutate(across(ID, ~ gsub("325","LN",.)))

plot_data2<-plot_data %>% filter(!(ID=="LN"))
dt2<-t.test((plot_data2 %>% filter(ID=="HC"))$p, (plot_data %>% filter(ID=="SLE"))$p, paired = FALSE)

plot_SLE2<-ggplot(plot_data2, aes(x = ID, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = NULL,
       x = NULL,
       y = "predicted DPI") +
  #ylim(2,52) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  xlim("HC", "SLE") +
  annotate("text", x=1, y=max(plot_data2$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(plot_data2$p)*1.1, xend = 1.6, yend = max(plot_data2$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x=2, y=max(plot_data2$p)*1.1, label=paste("p=",format(dt2$p.value, digits=3)), size=5)
ggsave("output/Fig5D_SLE_Z6YW.png", width = 3000, height = 1920, units = "px")
