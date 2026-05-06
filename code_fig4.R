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
library(randomForestExplainer)
library(tidyr)
library(dplyr)
library(ggplot2)

### 
# insert here the location of the folder containing the files from https://github.com/rki-umcg/ASC_ME

data_path<-"..."  #for example: data_path<-"C:/ASC_ME/"
setwd(data_path)
if (file.exists("output")){} else {dir.create(file.path(data_path, "output"))}

#import the RF model (rf_ascmi object)
rf_ascmi<-readRDS("final_model.rds")

###: evaluate the trained model 
#
##extract model characteristics
tree_depth_min<-list()
for (f in 1:20) {
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
for (f in 1:20) {
  forest_importance[[f]]<-measure_importance(rf_ascmi[[f]])
}
all_forest_importance <- do.call(rbind, forest_importance)

write_xlsx(marker_depth, "output/marker_depth.xlsx")
write_xlsx(all_forest_importance, "output/all_forest_importance.xlsx")

#generate min_tree_depth distribution heatmap for Figure 4A

md_long <- marker_depth %>%
  mutate(sample_id = row_number()) %>%        # keep original row id if desired
  pivot_longer(cols = -sample_id, names_to = "Marker", values_to = "Depth") %>%
  group_by(Marker) %>%
  arrange(Depth, .by_group = TRUE) %>%        # sort values low -> high within each marker
  mutate(rank_within_marker = row_number()) %>%
  ungroup()

md_long$Marker <- factor(md_long$Marker,
                         levels = c("CD28", "HLA", "CD20", "CD19", "CD138", "CD45", "Ki67", "CD56"))

# Make Depth a factor so we can use a discrete palette (0..8)
md_long$Depth <- factor(md_long$Depth, levels = 0:8)

# Basic plot
{
p <- ggplot(md_long, aes(x = rank_within_marker, y = Marker, fill = Depth)) +
  geom_tile(width = 1, height = 1, color = NA) +
  geom_vline(xintercept = seq(200, 1800, by = 200),
             linetype = "dotted", size = 1.5, color = "black", alpha = 0.8) +
  scale_fill_viridis_d(option = "turbo", direction = -1, begin = 0, end = 1,
                       name = "minimal tree depth of marker:") +
  scale_x_continuous(
    name = "number of trees within random forest model",
    expand = c(0, 0),
    breaks = seq(0, 2000, by = 200),
    limits = c(-6, 2004),
    position = "top"
  ) +
  scale_y_discrete(name = NULL, expand = c(0, 0), position = "right") +
  theme_bw(base_rect_size = 4, base_line_size = 1.5) +
  theme(
    title = element_text(size=15, face = "bold", colour = "black",hjust = 1.06),
    panel.grid = element_blank(),
    panel.spacing = grid::unit(0, "pt"),
    axis.text.x = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 14, face = "bold", colour = "white"),
    axis.ticks.length = grid::unit(10, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.align = 0.5,
    legend.key.size = grid::unit(0.9, "lines"),
    legend.text = element_text(size = 14,face = "bold", colour = "black"),
    legend.title = element_text(size = 14,face = "bold", colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"))
}

ggsave("output/Fig4A_tree_depth.png", plot = p)

#generate importance overview for markers in Fig 4B

{
mean_importance<-all_forest_importance[c(1:8),-c(6,8)]
colnames(mean_importance)[colnames(mean_importance)=="mean_min_depth"] <- "mean depth"
colnames(mean_importance)[colnames(mean_importance)=="no_of_nodes"] <- "tree nodes (%)"
colnames(mean_importance)[colnames(mean_importance)=="mse_increase"] <- "MSE increase"
colnames(mean_importance)[colnames(mean_importance)=="node_purity_increase"] <- "node purity"
colnames(mean_importance)[colnames(mean_importance)=="times_a_root"] <- "tree root (%)"
for (a in 1:8) {
  average<-all_forest_importance %>% filter(variable==paste(mean_importance[a,1]))
  mean_importance[a,2]<-mean(average$mean_min_depth)
  mean_importance[a,3]<-(sum(average$no_of_nodes)/sum(all_forest_importance$no_of_nodes)*100)
  mean_importance[a,4]<-mean(average$mse_increase)
  mean_importance[a,5]<-mean(average$node_purity_increase)
  mean_importance[a,6]<-mean(average$times_a_root)
}
mean_importance[,c(2:6)]<-round(mean_importance[,c(2:6)],2)
mean_importance[,c(5)]<-round(mean_importance[,c(5)],0)
}

mean_importance_scale<-mean_importance
for (a in 1:8) {
  mean_importance_scale[a,2]<-(1-(mean_importance[a,2]-min(mean_importance$`mean depth`))/(max(mean_importance$`mean depth`)-min(mean_importance$`mean depth`)))
  mean_importance_scale[a,3]<-(mean_importance[a,3]-min(mean_importance$`tree nodes (%)`))/(max(mean_importance$`tree nodes (%)`)-min(mean_importance$`tree nodes (%)`))
  mean_importance_scale[a,4]<-(mean_importance[a,4]-min(mean_importance$`MSE increase`))/(max(mean_importance$`MSE increase`)-min(mean_importance$`MSE increase`))
  mean_importance_scale[a,5]<-(mean_importance[a,5]-min(mean_importance$`node purity`))/(max(mean_importance$`node purity`)-min(mean_importance$`node purity`))
  mean_importance_scale[a,6]<-(mean_importance[a,6]-min(mean_importance$`tree root (%)`))/(max(mean_importance$`tree root (%)`)-min(mean_importance$`tree root (%)`))
}

mis_long <- mean_importance_scale %>%
  pivot_longer(
    cols = -variable,
    names_to = "Metric",
    values_to = "Value"
  )

mis_long2 <- mean_importance %>%
  pivot_longer(
    cols = -variable,
    names_to = "Metric",
    values_to = "Value"
  )

mis_join<-cbind(mis_long, raw=mis_long2$Value)

mis_join$variable <- factor(mis_join$variable, 
                            levels = c("CD28", "HLA.DR", "CD20", "CD19", "CD138", "CD45", "Ki67", "CD56"))

p2 <- ggplot(mis_join, aes(x = Metric, y = variable, fill = Value)) +
  geom_tile(width = 1, height = 1, color = NA) +  # add subtle borders
  geom_text(aes(x = Metric,y = variable, label = raw), color = "black", size = 5, inherit.aes = TRUE, fontface = "bold") +
  scale_fill_gradient(low = "#FAD6D6", high = "red3",name = "low importance                                      high importance") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw(base_rect_size = 4, base_line_size = 1.5) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = grid::unit(0, "pt"),
    axis.ticks.length = grid::unit(10, "pt"),
    axis.text.x = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black",hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.align = 0.5,
    legend.key.size = grid::unit(5, "lines"),
    legend.text = element_text(size = 0,face = "bold", colour = "white"),
    legend.title = element_text(size = 14,face = "bold", colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  guides(fill = guide_colorbar(title.position = "top",barwidth = 25,barheight = 1))

ggsave("output/Fig4B_importance.png", plot = p2)

