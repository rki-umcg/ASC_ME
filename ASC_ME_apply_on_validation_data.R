## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

library(dplyr)
library(readxl)
library(mice)
library(randomForest)
library(tidyverse)

data_new <- read_excel("all_extended_data.xlsx") %>% 
                    select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
                    filter(!(sample=="NA")) %>%
                    filter(((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>%
                    filter((StudyID=="SDY648")
                           |(StudyID=="SDY739")
                           |(StudyID=="SDY368")
                           |(StudyID=="SDY34")
                           |(StudyID=="SDY1397")
                           |(StudyID=="SDY788"))

new_merging <- data_new %>% 
  select(sample,check,SubjectID,File,condition)

new_imp <- data_new %>% 
  select(-c(check,SubjectID,File,condition))

percent_miss<-(sum(is.na(new_imp)) / (nrow(new_imp) * ncol(new_imp))) *100

iter=20
data.norm<-cbind(data.imputed_mean[,c(1,2,4:11)],dataforimp[,c(11:16)])
data.total<-rbind(data.norm,new_imp)
imp.new <- mice(data.total, method = "rf", m = iter, seed = 555)

new.imputed <- list()
for(m in 1:20){
  new.imputed[[m]] <- complete(imp.new, action = m)
} 
for(m in 1:20){  
  new.imputed[[m]] <- cbind(new.imputed[[m]][c((nrow(data.norm)+1):nrow(new.imputed[[m]])),], new_merging)
}

#get average imputed values for prediction markers
#prediction_marker<-c(colnames(dataforimp[,c(3:10)]))

new.imputed_mean <- new.imputed[[1]][,c(1,2,14)]
for (o in 1:length(prediction_marker)) {
  new_average_marker <- matrix(0, nrow(new.imputed[[1]]))
  for (m in 1:length(new.imputed)) {
    new_average_marker <- new_average_marker + new.imputed[[m]][,o+2]
  }
  new_average_marker <- new_average_marker / length(new.imputed)
  new.imputed_mean <- cbind(new.imputed_mean,new_average_marker)
  colnames(new.imputed_mean)[o+3] <- paste(prediction_marker[o])
}

### TDS-code:
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
predictions_new <- cbind(new_merging,new_imp[,c(11:12,15)], prediction_new[[1]][, 1:2], new_predictions)

# trim output table
new_output<-predictions_new %>%
  select(sample, StudyID,SubjectID, organ, group, dpi, p) %>%
  mutate(Age=new_imp$Age) %>%
  mutate(Gender=new_imp$Gender)
new_output$Age<-as.numeric(new_output$Age)

val_data<-new_output %>%
  filter(dpi>4 & ((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>% 
  filter((StudyID=="SDY648")|(StudyID=="SDY739")|(StudyID=="SDY368")|(StudyID=="SDY34")|(StudyID=="SDY1397")|(StudyID=="SDY788"))
  
r_ex_kinetic<-cor.test(val_data$dpi, val_data$p, method = "pearson")
val_plot_scatter<-ggplot(val_data, aes(x=dpi, y=p)) +
        geom_point(shape=1) + #none-filled circles
        geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
        labs(title = paste("pearson r =",format(r_ex_kinetic$estimate, digits=3)," p =",format(r_ex_kinetic$p.value, digits=3)),
             x = "actual DPI",
             y = "predicted DPI") +
        theme_bw()+
        xlim(5,250) +
        ylim(10,128) +
        theme(legend.position = "none",
              axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
              axis.title.y = element_text(size = 20,face = "bold"),
              axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
              axis.title.x = element_text(size = 20,face = "bold"),
              plot.title = element_text(size = 15,face = "bold"))
        
print(val_plot_scatter)
ggsave("plot_val_kinetic.png", plot = val_plot_scatter)
## -> Figure 5A

# check longitudinal/kinetic samples
val_data_long<-new_output %>%
  filter((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT")) %>% 
  filter((StudyID=="SDY648")|(StudyID=="SDY739")|(StudyID=="SDY368")|(StudyID=="SDY34")|(StudyID=="SDY1397")|(StudyID=="SDY788")) %>%
  filter(dpi==0|dpi>4)

#Identify non-unique SubjectID values in model_HC_long
non_unique_subjects <- val_data_long %>%
  count(SubjectID) %>%  # Count occurrences of each SubjectID
  filter(n > 1) %>%  # Keep only SubjectID with more than 1 occurrence
  pull(SubjectID)  # Extract the SubjectID values
#filter out non-unique samples
val_data_long <- val_data_long %>%
  filter(SubjectID %in% non_unique_subjects)

plot_val_long<-ggplot(val_data_long, aes(x=dpi, y=p, group=SubjectID, color=SubjectID)) +
  geom_point(shape=1) + #none-filled circles
  geom_line(size=0.25) + #, aes(color=SubjectID)) +
  scale_color_viridis(discrete = TRUE, option = "H") +
  labs(title = "linked longitudinal data points",
       x = "actual DPI",
       y = "predicted DPI") +
  xlim(0,250) +
  ylim(10,128) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))
print(plot_val_long)
ggsave("plot_val_long.png", plot = plot_val_long)
## -> Figure 5B

#by dpi groups
model_val_gp<-new_output %>%
  filter((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")) %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl", "early", "intermediate", "late", "very late"))) %>%  # Define the labels for the age groups
  filter(dpi==0|dpi<300) #%>%
#filter(p<47)

ANOVA_dev<- lm(model_val_gp[[7]] ~ model_val_gp[[10]])  |> aov() |> TukeyHSD()

plot_val_dpi_gp<-ggplot(model_val_gp, aes(x = dpi_group, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(model_val_gp[[2]])), width = 0.25, alpha = 1, size = 2,shape=1) +
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
  annotate("text", x="early", y=145, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][1,4], format = "e", digits=2)), size=5)+
  annotate("text", x="intermediate", y=145, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][2,4], format = "e", digits=2)), size=5)+
  annotate("text", x="late", y=145, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][3,4], format = "e", digits=2)), size=5)+
  annotate("text", x="very late", y=145, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][4,4], format = "e", digits=2)), size=5)+
  annotate("text", x="baseline", y=135, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][4,4], format = "e", digits=2)), size=5)+
  annotate("text", x="early", y=135, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][7,4], format = "e", digits=2)), size=5)+
  annotate("text", x="intermediate", y=135, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][9,4], format = "e", digits=2)), size=5)+
  annotate("text", x="late", y=135, label=paste("p=",formatC(ANOVA_dev[["model_val_gp[[10]]"]][10,4], format = "e", digits=2)), size=5)+
  geom_segment(aes(x = 4.6, y = 134, xend = 4.4, yend = 134),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x="very late", y=135, label=paste("comparison"), size=5)
print(plot_val_dpi_gp) 
ggsave("plot_val_dpi_gp.png", plot = plot_val_dpi_gp)
## -> Figure 5C

