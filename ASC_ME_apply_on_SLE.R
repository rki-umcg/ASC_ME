## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

library(randomForest)
library(mice)
library(readxl)
library(dplyr)
library(tidyverse)

set.seed(555)
data_new <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(CD98,Bcl.2,BAFFR,BCMA,TACI,CD44,CD95,CD69,No,CD86)) %>% 
  filter((StudyID=="Z6YW")|(StudyID=="AID")|(StudyID=="SDY997")) %>% filter(!(group=="anchor")) %>%
  filter(!(condition=="pSLE_T6")) %>%
  filter(!(condition=="pSLE_T6R")) %>%
  filter(!(group=="HC_Vac")) %>% filter(!(group=="SSc")) %>% filter(!(group=="SjD"))

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

#plot SLE data
AID_data<-new_output %>%
  filter(!(group=="anchor")) %>%
  filter((StudyID=="Z6YW")) %>%
  filter(!(str_detect(sample, "_B_"))) %>% filter(!(str_detect(sample, "_C_"))) %>% filter(!(str_detect(sample, "_D_"))) %>% filter(!(str_detect(sample, "_E_"))) %>% filter(!(str_detect(sample, "_F_"))) %>% filter(!(str_detect(sample, "_G_")))
AID_data<-AID_data %>% mutate(across(group, ~ gsub("SLE_LN", "lupus nephritis", .)))
                                                                                                                                                                                                               
ANOVA_AID<- dunn.test(AID_data[[7]], AID_data[[5]], method = "bonferroni", list = TRUE)

plot_AID_own<-ggplot(AID_data, aes(x = group, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(AID_data[[2]])), width = 0.25, alpha = 1, size = 2, shape=16) +
  theme_bw() +
  labs(title = "Kruskal-Wallis test with Bonferroni correction",
       x = "",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  annotate("text", x="HC", y=60, label=paste("p=",format(ANOVA_AID[["P.adjusted"]][2], digits=3)), size=5)+
  annotate("text", x="lupus nephritis", y=60, label=paste("p=",format(ANOVA_AID[["P.adjusted"]][3], digits=3)), size=5)+
  xlim("HC","SLE","lupus nephritis")
print(plot_AID_own) 
ggsave(filename = "plot_SLE_baseline.png", plot = plot_AID_own)
## -> Figure 6C

SLE_scores <- read_excel("C:/Users/SteinmetzDT/OneDrive - UMCG/own publications/Steinmetz_prep_ASC MI study/SLE_Z6YW_stats_sorted.xlsx")

SLE_merged<-merge(new_output%>%filter(StudyID=="Z6YW"), SLE_scores, by="sample", all.x = TRUE)# %>% filter(MMF_time>0)
r_AID_plot<-cor.test(SLE_merged$p, SLE_merged$Timepoint, method = "pearson", exact = FALSE)
plot_AID_example<-ggplot(SLE_merged, aes(x = Timepoint, y = p)) +
  geom_jitter(width = 0.1, alpha = 1, size = 3) +
  theme_bw() +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = paste("pearson r =",format(r_AID_plot$estimate, digits=3)," p =",format(r_AID_plot$p.value, digits=3)),
       x = "Timepoint",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  scale_x_continuous(breaks = seq(1,7, by=1))
  #xlim(0.5,7.5)
print(plot_AID_example) 
ggsave(filename = "plot_SLE_timeline.png", plot = plot_AID_example)
## -> Figure 6C

r_AID_plot<-cor.test(SLE_merged[SLE_merged$Steroid>0,7], SLE_merged[SLE_merged$Steroid>0,25], method = "spearman", exact = FALSE)
plot_AID_example<-ggplot(SLE_merged, aes(x = Steroid, y = p)) +
  geom_jitter(width = 0.1, alpha = 1, size = 3) +
  theme_bw() +
  geom_smooth(method = "lm", color = "red", se = TRUE, 
              data = SLE_merged %>% filter(Steroid > 0), 
              aes(x = Steroid, y = p)) +
  labs(title = paste("spearman r =",format(r_AID_plot$estimate, digits=3)," p =",format(r_AID_plot$p.value, digits=3)),
       x = "steroid usage intensity",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  xlim(-0.2,3.2)
print(plot_AID_example) 
ggsave(filename = "plot_SLE_steroids.png", plot = plot_AID_example)
## -> Figure 6C

SLE_scores <- read_excel("C:/Users/SteinmetzDT/OneDrive - UMCG/own publications/Steinmetz_prep_ASC MI study/- manuscript prep 15Okt/.final/codes/backup/ASC_ME_prediction_model/SLE_Z6YW_stats_sorted.xlsx")
SLE_merged<-merge(new_output, SLE_scores, by="sample", all.x = TRUE) 
SLE_merged<-SLE_merged %>% filter(StudyID=="Z6YW") %>% filter(Timepoint==1)

#LDA, Arthritis, Complement, ds_DNA, Leukopenia, Steroid_use, MMF_MPA, RTX_BLY, Cytoxan, AZA_MTX
{r_test<-wilcox.test(p ~ LDA, data = SLE_merged)
print(ggplot(SLE_merged, aes(x = LDA, y = p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(SLE_merged[[2]])), width = 0.25, alpha = 1, size = 2, shape=16) +
  theme_bw() +
  #geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Mann Whitney U-test",
       x = "LDA",
       y = "predicted DPI",
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))+
  annotate("text", x = 1.5, y = 60, label = paste("p=", format(r_test$p.value, digits = 3)), size = 5)
)}
ggsave(filename = "plot_SLE_Proteinuria.png")

#C3, C4, ESR, dsDNA, SLEDAI
{r_test<-cor.test(SLE_merged$p, SLE_merged$C3, method = "spearman", exact = FALSE)
  print(ggplot(SLE_merged, aes(x = C3, y = p)) +
          geom_point(shape=16, size = 3) + 
          geom_smooth(method = "lm", color = "red", se = TRUE , level = 0.95) + #add regression line /w confidence interval
          theme_bw() +
          #geom_smooth(method = "lm", color = "red", se = TRUE) +
          labs(title = paste("spearman correlation, r =",format(r_test$estimate, digits=3)," p =",format(r_test$p.value, digits=3)),
               x = "serum C3 [mg/dl]",
               y = "predicted DPI",
               color = "StudyID") +
          theme(legend.position = "none",
                axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.y = element_text(size = 20,face = "bold"),
                axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.x = element_text(size = 20,face = "bold"),
                plot.title = element_text(size = 15,face = "bold"))
    )}
ggsave(filename = "plot_SLE_BL_C3.png")

