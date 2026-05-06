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
library(dplyr)
library(ggplot2)
library(readxl)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

# import data from model development
data.imputed_mean<-read_xlsx("data.imputed_mean.xlsx")
data.unimputed<-read_xlsx("data.unimputed.xlsx")

###dpi grouping for imputed data
imp_markers<-data.imputed_mean %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 10, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl","very early", "early", "intermediate", "late", "very late"))) # Define the labels for the age groups
imp_markers<-cbind(imp_markers,StudyID=data.unimputed$StudyID)

#
for (a in 4:11) {
ANOVA_plot<- lm(imp_markers[[a]] ~ imp_markers$dpi_group)  |> aov() |> TukeyHSD()
plot_imp_marker<-ggplot(imp_markers, aes(x = dpi_group, y = imp_markers[,a])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(imp_markers$StudyID)), width = 0.25, alpha = 1, size = 2,shape=1) +
  theme_bw() +
  labs(title = NULL,
       x = "immune response stage",
       y = paste(colnames(imp_markers[a]), "(normalized median expression)"),
       color = "StudyID") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold")) +
  xlim("baseline", "very early", "early", "intermediate", "late", "very late") +
  annotate("text", x="baseline", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("comparison"), size=5)+
  geom_segment(aes(x = 1.4, y = max(imp_markers[,a],na.rm = TRUE)*1.1, xend = 1.6, yend = max(imp_markers[,a],na.rm = TRUE)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text", x="very early", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][2,4], format = "e", digits=2)), size=5) +
  annotate("text", x="early", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][3,4], format = "e", digits=2)), size=5) +
  annotate("text", x="intermediate", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][4,4], format = "e", digits=2)), size=5) +
  annotate("text", x="late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][5,4], format = "e", digits=2)), size=5) +
  annotate("text", x="very late", y=max(imp_markers[,a],na.rm = TRUE)*1.1, label=paste("p=",formatC(ANOVA_plot[["imp_markers$dpi_group"]][6,4], format = "e", digits=2)), size=5)
print(plot_imp_marker) 
ggsave(filename = paste("output/Fig2B_plot_imp_marker_",paste(colnames(imp_markers[a])),".png", sep = ""),
       width = 3000, height = 1920, units = "px")
}

###age correlation for imputed data
imp_markers$Age<-as.numeric(imp_markers$Age)
imp_markers_BL<-imp_markers %>% filter(dpi_group=="baseline")
for (a in 4:11) {
  r_plot<-cor.test(imp_markers_BL$Age, imp_markers_BL[[a]])
  print(ggplot(imp_markers_BL, aes(x=Age, y=imp_markers_BL[,a])) +
    geom_point(shape=1) + #none-filled circles
    geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
    labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
         x = "sample age",
         y = paste(colnames(imp_markers_BL[a]), "(normalized median expression)")) +
    theme_bw() +
    #ylim(-1,5)+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("output/FigS4_plot_BL_",colnames(imp_markers_BL[a]),"_age.png", sep = ""),
         width = 3000, height = 1920, units = "px")
}

###
###dpi grouping for unimputed data
unimp_markers<-data.unimputed %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 10, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "very early", "early", "intermediate", "late", "very late"))) # Define the labels for the age groups
#
for (b in c(3:10)) {
  ANOVA_plot<- lm(unimp_markers[[b]] ~ unimp_markers[[20]])  |> aov() |> TukeyHSD()
  annotations <- list()
  annotations2 <- list()
  groups <- c("very early", "early", "intermediate", "late", "very late")
  groups2 <- c("baseline", "very early", "early", "intermediate", "late")
  for (group in groups) {
    if (sum(!is.na(unimp_markers[unimp_markers$dpi_group == group, b])) > 0) {
      annotations <- append(annotations,annotate("text", x = group, y = max(unimp_markers[, b], na.rm = TRUE) * 1.2,
                                                 label = paste("p=", formatC(ANOVA_plot[["unimp_markers[[20]]"]][paste0(group, "-baseline"), 4],format = "e",digits = 2)
                                                 ),size = 5))}}
  for (group in groups2) {
    if (sum(!is.na(unimp_markers[unimp_markers$dpi_group == group, b])) > 0) {
      annotations2 <- append(annotations2,annotate("text", x = group, y = max(unimp_markers[, b], na.rm = TRUE) * 1.1,
                                                 label = paste("p=", formatC(ANOVA_plot[["unimp_markers[[20]]"]][paste0("very late-", group), 4],format = "e",digits = 2)
                                                 ),size = 5))}}
  plot_unimp_marker<-ggplot(unimp_markers, aes(x = dpi_group, y = unimp_markers[[b]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = as.factor(unimp_markers[[15]])), width = 0.25, alpha = 1, size = 2,shape=1) +
    theme_bw() +
    labs(title = NULL,
         x = "immune response stage",
         y = paste(colnames(unimp_markers[b]), "(normalized median expression)"),
         color = "StudyID") +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold")) +
    xlim("baseline", "very early", "early", "intermediate", "late", "very late") +
    annotations +
    annotations2 +
    annotate("text", x=1, y=max(unimp_markers[,b],na.rm = TRUE)*1.2, label=paste("comparison"), size=5)+
    geom_segment(aes(x = 1.4, y = max(unimp_markers[,b],na.rm = TRUE)*1.2, xend = 1.6, yend = max(unimp_markers[,b],na.rm = TRUE)*1.2),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
    annotate("text", x=6, y=max(unimp_markers[,b],na.rm = TRUE)*1.1, label=paste("comparison"), size=5)+
    geom_segment(aes(x = 5.6, y = max(unimp_markers[,b],na.rm = TRUE)*1.1, xend = 5.4, yend = max(unimp_markers[,b],na.rm = TRUE)*1.1),arrow = arrow(length=unit(0.30,"cm"), type = "closed"))
  print(plot_unimp_marker) 
ggsave(filename = paste("output/FigS5_plot_unimp_",colnames(unimp_markers[b]),".png", sep = ""), plot = plot_unimp_marker,
       width = 3000, height = 1920, units = "px")  
}

#age correlation for unimputed data
unimp_markers$Age<-as.numeric(unimp_markers$Age)
unimp_markers_BL<-unimp_markers %>% filter(dpi_group=="baseline")
for (b in c(3:10)) {
  r_plot<-cor.test(unimp_markers_BL$Age, unimp_markers_BL[[b]])
  print(ggplot(unimp_markers_BL, aes(x=Age, y=unimp_markers_BL[[b]])) +
          geom_point(shape=1) + #none-filled circles
          geom_smooth(method = "lm", color = "red", se = TRUE ) + #add regression line /w confidence interval
          labs(title = paste("pearson correlation, r =",format(r_plot$estimate, digits=3)," p =",format(r_plot$p.value, digits=3)),
               x = "sample age",
               y = paste(colnames(unimp_markers_BL[b]), "(normalized median expression)")) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.y = element_text(size = 20,face = "bold"),
                axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.x = element_text(size = 20,face = "bold"),
                plot.title = element_text(size = 15,face = "bold"))
  )
  ggsave(filename = paste("output/FigS5_plot_BL.unimp_",colnames(unimp_markers_BL[b]),"_age.png", sep = ""),
         width = 3000, height = 1920, units = "px")
}
