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
library(tidyr)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

## contine with R environment from training the RF model
all_data.imputed_mean<-read_xlsx("train_val_data.imputed_mean.xlsx")
data.unimputed<-read_xlsx("data.unimputed.xlsx")

###dpi grouping for imputed data
imp_markers<-all_data.imputed_mean %>%
  mutate(dpi_group = cut(dpi,  
                         breaks = c(-Inf, 1, 4.9, 10, 19, 50, 99, 500),  # Define the breaks for the age groups
                         labels = c("baseline", "excl","very early", "early", "intermediate", "late", "very late"))) # Define the labels for the age groups

#generate dot plots with median+IQR range for markers/studies/timepoints

summary_stats <- imp_markers %>%
  pivot_longer(cols = 4:11, 
               names_to = "marker", 
               values_to = "value") %>%
  
  # subgroup-level summaries
  group_by(general = !!sym(names(imp_markers)[12]),
           subgroup = !!sym(names(imp_markers)[13]),
           marker) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(level = "subgroup") %>%
  
  bind_rows(
    # general-level summaries
    imp_markers %>%
      pivot_longer(cols = 4:11, 
                   names_to = "marker", 
                   values_to = "value") %>%
      group_by(general = !!sym(names(imp_markers)[12]), marker) %>%
      summarise(
        median = median(value, na.rm = TRUE),
        IQR = IQR(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(subgroup = "all data",   # use "ALL" instead of NA for plotting
             level = "general")
  ) %>%
  arrange(marker, general, level, subgroup)

plots <- lapply(unique(summary_stats$marker), function(m) {
  
  df <- summary_stats %>% filter(marker == m) %>% filter(!(subgroup=="excl"))
  
  ggplot(df, aes(x = subgroup, y = general)) +
    geom_point(aes(color = median, size = IQR), shape = 19, alpha = 0.9) +
    
    # fixed color scale for median
    scale_color_viridis_c(option = "plasma", limits = c(-1.25, 5)) +
    
    # fixed size scale for IQR
    scale_size_continuous(limits = c(0, 3), range = c(2, 8)) +
    geom_vline(xintercept = 1.5) +
    theme_bw() +
    geom_hline(yintercept = 6.5) +
    ylim("AID","SDY301","SDY34","SDY387","SDY819","SDY984",
         "CER","SDY1086","SDY1288","SDY144","SDY1669","SDY1734","SDY180","SDY224","SDY272","SDY296","SDY364","SDY368","SDY522","SDY648","SDY739","SDY788","ZY77") +
    xlim("all data", "baseline","very early","early","intermediate","late","very late") + 
    theme(
      axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
      plot.title = element_text(size = 15,face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = paste("Marker:", m),
      x = NULL,      y = NULL,
      color = "Median",
      size = "IQR"
    )
})

# export plots
for (a in 1:8) {
  print(plots[[a]])
  ggsave(filename = paste("output/marker_",paste(a),".png", sep = ""))
}
