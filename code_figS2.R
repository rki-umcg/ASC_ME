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
library(ggplot2)
library(dplyr)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

####
train_df <- readxl::read_excel("train_df.xlsx")
val_df <- readxl::read_excel("val_df.xlsx")

train_df$V2<-paste(train_df$Race, train_df$Race.Specify, train_df$Ethnicity)

train_df<-train_df %>%
  mutate(across(V2, ~ gsub("NA/other NA NA", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("White  Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("NA/other  Not Hispanic or Latino", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("White  Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("Black or African American  Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("Asian  Not Hispanic or Latino", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other W/AA Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Not Specified", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("Black or African American NA Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("White NA Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("Asian NA Not Hispanic or Latino", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other More than one race Not Hispanic or Latino", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("NA/other American Indian/White Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("White NA Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("NA/other White/Black or African American Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("Asian  Other", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Not Hispanic or Latino", "NA/other",.)))


### --- Prepare outer data: Sex ---
sex_df_train <- train_df %>%
  mutate(Gender = ifelse(is.na(Gender), "NA", Gender)) %>%
  count(Gender) %>%
  mutate(prop = n / sum(n),
         ring = "outer")

### --- Prepare inner data: Race ---
race_df_train <- train_df %>%
  mutate(V2 = ifelse(is.na(V2), "NA", V2)) %>%
  count(V2) %>%
  mutate(prop = n / sum(n),
         ring = "inner")

### --- Plot ---
ggplot() +
  # OUTER DONUT (Sex)
  geom_col(
    data = sex_df_train,
    aes(x = 2, y = prop, fill = Gender),
    width = 1,
    color = "black",     # black outline
    size = 1.2
  ) +
  # INNER PIE (Race)
  geom_col(
    data = race_df_train,
    aes(x = 1, y = prop, fill = V2),
    width = 1,
    color = "black",     # black outline
    size = 1.2
  ) +
  coord_polar(theta = "y") +
  # Adjust hole size (inner ring)
  xlim(0.5, 2.5) +
  theme_void() +
  theme(
    legend.position = "right"
  ) +
  guides(fill = guide_legend(title = "Categories"))
ggsave("output/pie_train_data.png")


####
val_df$V2<-paste(val_df$Race, val_df$Race.Specify, val_df$Ethnicity)

val_df<-val_df %>%
  mutate(across(V2, ~ gsub("NA/other NA NA", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("White  Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("NA/other  Not Hispanic or Latino", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("White  Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("Black or African American  Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("Asian  Not Hispanic or Latino", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other W/AA Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Not Specified", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("Black or African American NA Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("White NA Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("Asian NA Not Hispanic or Latino", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other More than one race Not Hispanic or Latino", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("NA/other American Indian/White Not Hispanic or Latino", "White",.))) %>%
  mutate(across(V2, ~ gsub("White NA Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("NA/other White/Black or African American Not Hispanic or Latino", "Black",.))) %>%
  mutate(across(V2, ~ gsub("Asian  Other", "Asian",.))) %>%
  mutate(across(V2, ~ gsub("NA/other NA Not Hispanic or Latino", "NA/other",.))) %>%
  mutate(across(V2, ~ gsub("NA/other Hispanic Hispanic or Latino", "Hispanic",.))) %>%
  mutate(across(V2, ~ gsub("White NA Other", "White",.))) %>%
  mutate(across(V2, ~ gsub("Black or African American NA Other", "Black",.))) %>%
  mutate(across(V2, ~ gsub("Asian NA Other", "Asian",.)))


### --- Prepare outer data: Sex ---
sex_df_val <- val_df %>%
  mutate(Gender = ifelse(is.na(Gender), "NA", Gender)) %>%
  count(Gender) %>%
  mutate(prop = n / sum(n),
         ring = "outer")

### --- Prepare inner data: Race ---
race_df_val <- val_df %>%
  mutate(V2 = ifelse(is.na(V2), "NA", V2)) %>%
  count(V2) %>%
  mutate(prop = n / sum(n),
         ring = "inner")

### --- Plot ---
ggplot() +
  # OUTER DONUT (Sex)
  geom_col(
    data = sex_df_val,
    aes(x = 2, y = prop, fill = Gender),
    width = 1,
    color = "black",     # black outline
    size = 1.2
  ) +
  # INNER PIE (Race)
  geom_col(
    data = race_df_val,
    aes(x = 1, y = prop, fill = V2),
    width = 1,
    color = "black",     # black outline
    size = 1.2
  ) +
  coord_polar(theta = "y") +
  # Adjust hole size (inner ring)
  xlim(0.5, 2.5) +
  theme_void() +
  theme(
    legend.position = "right"
  ) +
  guides(fill = guide_legend(title = "Categories"))

ggsave("output/pie_val_data.png")

