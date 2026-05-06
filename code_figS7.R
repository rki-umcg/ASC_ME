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
library(readxl)
library(viridis)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

## import data
model_output <- read_excel("final_model_output.xlsx")

# check longitudinal/kinetic samples
model_HC_long<-model_output %>%
  filter((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT")) %>%
  mutate(across(dpi, ~ ifelse(.x < 0, 0, .x)))

studies<-unique(model_HC_long$StudyID)

for (a in 1:length(studies)) {
  plot_data<-model_HC_long %>% filter(StudyID==paste(studies[[a]]))
  
plot_HC_long<-ggplot(plot_data, aes(x=dpi, y=p, group=SubjectID, color=SubjectID)) +
  geom_point(shape=1) + #none-filled circles
  geom_line(size=0.25) + #, aes(color=SubjectID)) +
  scale_color_viridis(discrete = TRUE, option = "H") +
  labs(title = paste("linked longitudinal data points for study", studies[[a]]),
       x = "actual DPI",
       y = "predicted DPI") +
  xlim(0,185)+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
        axis.title.x = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 15,face = "bold"))
print(plot_HC_long)
ggsave(paste("output/plot_",studies[[a]],".png", sep = ""))
}

###
# check longitudinal/kinetic samples
new_output<-read_excel("val_output.xlsx")
val_df <- read_excel("val_df.xlsx")

model_HC_val_long<-left_join(new_output, val_df[,c(1,30)], by = "sample") %>%
  filter((group == "HC")|(group == "HC_CoV")|(group == "HC_Inf")|(group == "HC_Vac")|(group == "HC_SCT")) %>%
  mutate(across(dpi, ~ ifelse(.x < 0, 0, .x)))

studies_val<-unique(model_HC_val_long$Study.Accession)

for (a in 1:length(studies_val)) {
  plot_data<-model_HC_val_long %>% filter(Study.Accession==paste(studies_val[[a]]))
  
  plot_HC_long<-ggplot(plot_data, aes(x=dpi, y=p, group=SubjectID, color=SubjectID)) +
    geom_point(shape=1) + #none-filled circles
    geom_line(size=0.25) + #, aes(color=SubjectID)) +
    scale_color_viridis(discrete = TRUE, option = "H") +
    labs(title = paste("linked longitudinal data points for study", studies_val[[a]]),
         x = "actual DPI",
         y = "predicted DPI") +
    xlim(0,185)+
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.y = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold", colour = "black"),
          axis.title.x = element_text(size = 20,face = "bold"),
          plot.title = element_text(size = 15,face = "bold"))
  print(plot_HC_long)
  ggsave(paste("output/plot_",studies_val[[a]],".png", sep = ""))
}

