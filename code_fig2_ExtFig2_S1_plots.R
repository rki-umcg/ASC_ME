## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

library(dplyr)
library(ggplot2)
library(readxl)

#dpi grouping for imputed data
imp_markers<-data.imputed_mean %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl","early", "intermediate", "late", "very late"))) # Define the labels for the age groups
imp_markers<-cbind(imp_markers,StudyID=dataforimp$StudyID)

for (a in 4:11) {
ANOVA_plot<- lm(imp_markers[[a]] ~ imp_markers[[13]])  |> aov() |> TukeyHSD()
#set.seed(555)
plot_imp_marker<-ggplot(imp_markers, aes(x = dpi_group, y = imp_markers[,a])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(imp_markers[[14]])), width = 0.25, alpha = 1, size = 2,shape=1) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "immune response stage",
       y = paste(colnames(imp_markers[a]), "(normalized mean expression)"),
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  #scale_y_continuous(breaks = c(0,30,45,60,90,120), limits = c(4,116)) +
  #ylim(-1,5.2) +
  xlim("baseline", "early", "intermediate", "late", "very late") +
  annotate("text", x="early", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers[[13]]"]][2,4], format = "e", digits=2)), size=5) +
  annotate("text", x="intermediate", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers[[13]]"]][3,4], format = "e", digits=2)), size=5) +
  annotate("text", x="late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers[[13]]"]][4,4], format = "e", digits=2)), size=5) +
  annotate("text", x="very late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers[[13]]"]][5,4], format = "e", digits=2)), size=5)
print(plot_imp_marker) 
ggsave(filename = paste("plot_imp_",colnames(imp_markers[a]),".png"), plot = plot_imp_marker)
}
## -> Figure 2

#age correlation for imputed data
imp_markers$Age<-as.numeric(imp_markers$Age)
imp_markers_BL<-imp_markers %>% filter(dpi_group=="baseline")
for (a in 4:11) {
  r_plot<-cor.test(imp_markers_BL$Age, imp_markers_BL[[a]])
  print(ggplot(imp_markers_BL, aes(x=Age, y=imp_markers_BL[,a])) +
    geom_point(shape=1) + #none-filled circles
    geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
    labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
         x = "sample age",
         y = paste(colnames(imp_markers_BL[a]), "(normalized mean expression)")) +
    theme_bw() +
    ylim(-1,5)+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("plot_BL_",colnames(imp_markers_BL[a]),"_age.png", sep = ""))
}
## -> Extended Figure 2


#dpi grouping for unimputed data
unimp_markers <- read_excel("all_extended_data.xlsx") %>% 
  select(-c(No)) %>% 
  filter(!(sample=="NA")) %>%
  filter(((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT"))) %>%
  filter((StudyID=="AID")|(StudyID=="CER")|(StudyID=="SDY144")|(StudyID=="SDY180")|(StudyID=="SDY224")
         |(StudyID=="SDY272")|(StudyID=="SDY296")|(StudyID=="SDY301")|(StudyID=="SDY364")|(StudyID=="SDY387")
         |(StudyID=="SDY522")|(StudyID=="SDY80")|(StudyID=="SDY819")|(StudyID=="SDY984")|(StudyID=="SDY1086")
         |(StudyID=="SDY1288")|(StudyID=="SDY1669")|(StudyID=="SDY1734")|(StudyID=="ZY77")
  )
unimp_markers<-unimp_markers %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "early", "intermediate", "late", "very late"))) # Define the labels for the age groups

for (b in c(3,4,7,8,10,13,14,15)) {
#col 5,6,9,16,18,19
#b<-6
ANOVA_plot<- lm(unimp_markers[[b]] ~ unimp_markers[[30]])  |> aov() |> TukeyHSD()
#set.seed(555)
plot_unimp_marker<-ggplot(unimp_markers, aes(x = dpi_group, y = unimp_markers[[b]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(unimp_markers[[20]])), width = 0.25, alpha = 1, size = 2,shape=1) +
  theme_bw() +
  labs(title = "one-way ANOVA with Tukey correction",
       x = "immune response stage",
       y = paste(colnames(unimp_markers[b]), "(normalized mean expression)"),
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  #scale_y_continuous(breaks = c(0,30,45,60,90,120), limits = c(4,116)) +
  #ylim(-0.0,3.0) +
  xlim("baseline", "early", "intermediate", "late", "very late") +
  annotate("text", x="early", y=max(unimp_markers[,b],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["unimp_markers[[30]]"]][1,4], format = "e", digits=2)), size=5)+
  annotate("text", x="intermediate", y=max(unimp_markers[,b],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["unimp_markers[[30]]"]][2,4], format = "e", digits=2)), size=5)+
  annotate("text", x="late", y=max(unimp_markers[,b],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["unimp_markers[[30]]"]][3,4], format = "e", digits=3)), size=5)+
  annotate("text", x="very late", y=max(unimp_markers[,b],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["unimp_markers[[30]]"]][4,4], format = "e", digits=2)), size=5)
print(plot_unimp_marker) 
ggsave(filename = paste("plot_unimp_",colnames(unimp_markers[b]),".png", sep = ""), plot = plot_unimp_marker)   
}
## -> Figure Supp.1A

#age correlation for unimputed data
unimp_markers$Age<-as.numeric(unimp_markers$Age)
unimp_markers_BL<-unimp_markers %>% filter(dpi_group=="baseline")
for (b in c(3:10,13:16,18,19)) {
  r_plot<-cor.test(unimp_markers_BL$Age, unimp_markers_BL[[b]])
  print(ggplot(unimp_markers_BL, aes(x=Age, y=unimp_markers_BL[[b]])) +
          geom_point(shape=1) + #none-filled circles
          geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
          labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
               x = "sample age",
               y = paste(colnames(unimp_markers_BL[b]), "(normalized mean expression)")) +
          theme_bw() +
          ylim(-1,5)+
          theme(legend.position = "none",
                axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.y = element_text(size = 20,face = "bold"),
                axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.x = element_text(size = 20,face = "bold"),
                plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("plot_BL.unimp_",colnames(unimp_markers_BL[b]),"_age.png", sep = ""))
}
#b<-13, 18
