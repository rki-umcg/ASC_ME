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
  filter(count>49) %>%
  filter((StudyID=="No49_Mitsialis_LP")|(StudyID=="Festen")|(StudyID=="JJS")|(StudyID=="SDY1389")|(StudyID=="SDY2107")
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
  labs(title = NULL,
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  annotate("text", x="BM", y=max(comp_BLD_BM_gut$p)*1.1, label=paste("p=",formatC(ANOVA_imp[["comp_BLD_BM_gut[[4]]"]][1,4], format = "e", digits=2)), size=5)+
  annotate("text", x=1, y=max(comp_BLD_BM_gut$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(comp_BLD_BM_gut$p)*1.1, xend = 1.6, yend = max(comp_BLD_BM_gut$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x="gut", y=max(comp_BLD_BM_gut$p)*1.1, label=paste("p=",formatC(ANOVA_imp[["comp_BLD_BM_gut[[4]]"]][2,4], format = "e", digits=2)), size=5)
print(plot_HC_organ) 
ggsave("output/Fig5A_organ_comp.png", plot = plot_HC_organ, width = 3000, height = 1920, units = "px")

#gut section
Festen_data<-new_output %>%
  filter(StudyID == "Festen") %>%
  mutate(organ = sample) %>%
  mutate(across(organ, ~ gsub("V010 T0 ileum NI 3360 MC A.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V010 T0 ileum NI 3360 MCint_1 B.fcs","ileum",.))) %>%
  #mutate(across(organ, ~ gsub("V011 T0 rectum I 3394 MCint_1 B.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V011 T4 sigmoid NI 3445 MCint_1 B.fcs","colon",.))) %>%
  #mutate(across(organ, ~ gsub("V013 T0 rectum I 3397 MCint_1 B 1.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V014 T4 sigmoid NI 3458 MCint_1 B 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V016 T4 ileum NI 3606 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V018 T0 ileum NI 3402 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V018 T0 sigmoid I 3402 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V019 T0 ileum I 3417 MC A 1.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V019 T0 ileum I 3417 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V019 T0 sigmoid I 3417 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V019 T4 sigmoid I 3497 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V020 T0 cdescen NI 3453 blanco_1 B 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V020 T0 cdescen NI 3453 MCint_1 B 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V020 T0 ileum NI 3453 MCint_1 B 1.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V020 T0 sigmoid I 3453 MCint_1 B 1.fcs","colon",.))) %>%
  #mutate(across(organ, ~ gsub("V020 T4 rectum NI 3537 MC A 1.fcs","rectum",.))) %>%
  #mutate(across(organ, ~ gsub("V020 T4 rectum NI 3537 MCint_1 B 1.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V020 T4 sigmoid NI 3537 MC A 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V020 T4 sigmoid NI 3537 MCint_1 B 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V022 T0 ileum NI 3452 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V022 T4 sigmoid I 3582 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V023 T0 ileum I 3494 MC A.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V023 T0 ileum I 3494 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V023 T0 sigmoid I 3494 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V023 T4 sigmoid NI 3606 blanco_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V023 T4 sigmoid NI 3606 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V024 T0 sigmoid I 9B MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V031 T0 ileum NI 3534 MC A.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V031 T0 ileum NI 3534 MCint_1 B 1.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V031 T0 sigmoid NI 3534 MCint_1 B 1.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V031 T4 sigmoid NI 3608 MCint_1 B 1.fcs","colon",.))) %>%
  #mutate(across(organ, ~ gsub("V032 T0 rectum I 3390 MCint_1 B.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V032 T0 sigmoid NI 3390 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V032 T4 sigmoid NI 3638 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V033 T0 ileum NI 3549 MCint_1 B.fcs","ileum",.))) %>%
  #mutate(across(organ, ~ gsub("V033 T0 rectum NI 3549 MCint_1 B.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V033 T0 sigmoid NI 3549 MCint_1 B.fcs","colon",.))) %>%
  #mutate(across(organ, ~ gsub("V033 T4 rectum NI 28B MCint_1 B.fcs","rectum",.))) %>%
  mutate(across(organ, ~ gsub("V033 T4 sigmoid NI 28A MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V035 T0 cdescen I 19A MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V035 T4 sigmoid NI 3684 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V039 T0 ileum NI 3640 MCint_1 B.fcs","ileum",.))) %>%
  mutate(across(organ, ~ gsub("V039 T4 sigmoid I 3711 MCint_1 B.fcs","colon",.))) %>%
  mutate(across(organ, ~ gsub("V039 T4 sigmoid NI 3711 MCint_1 B.fcs","colon",.))) %>%
  filter((organ=="ileum")|(organ=="colon"))

dt<-wilcox.test((Festen_data %>% filter(organ=="ileum"))[,11],(Festen_data %>% filter(organ=="colon"))[,11], paired = FALSE)

plot_Festen_organ<-ggplot(Festen_data, aes(x = organ, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 1, size = 2) +
  theme_bw() +
  labs(title = NULL,       x = NULL,       y = "predicted DPI") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("ileum", "colon")+
  annotate("text", x=1, y=max(Festen_data$p)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(Festen_data$p)*1.1, xend = 1.6, yend = max(Festen_data$p)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x=2, y=max(Festen_data$p)*1.1, label=paste("p=",format(dt$p.value, digits=3)), size=5)
print(plot_Festen_organ) 
ggsave(filename = "output/Fig5C_plot_Festen_organ.png", plot = plot_Festen_organ, width = 3000, height = 1920, units = "px")

