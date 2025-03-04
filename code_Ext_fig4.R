# code for Extended Data Figure 3
## contine with R environment from training the RF model
load(".../ASC_ME/ASC_ME_environment.RData")
setwd(".../ASC_ME")

# Bin the 'predict' values into sections of width 5 and calculate the midpoint of each bin
binned_data <- data.imputed_mean %>%
  mutate(predict_bin = cut(predict, breaks=seq(5, 125, by=5), right=FALSE),
         predict_mid = as.numeric(as.character(cut(predict, breaks=seq(5, 125, by=5), labels=seq(7.5, 122.5, by=5)))))
binned_data$Age<-as.numeric(binned_data$Age)

# Plot boxplots with the actual data points plotted at their true x-values
for (a in 3:3) {
  print(ggplot(binned_data, aes(x=predict, y=binned_data[,a])) +
          geom_boxplot(aes(x=predict_mid, group=predict_bin), outlier.shape = NA, width=4.5, color="red") +  # Box plots at midpoints
          geom_point(alpha=0.5, shape=1, size=1) +  # Points at actual x-values
          theme_bw() +
          labs(x="predict DPI",
               y="sample Age") +
          theme(axis.title.y = element_text(size = 20,face = "bold"),
                axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
                axis.title.x = element_text(size = 20,face = "bold"),
                axis.text.x = element_text(size = 20,face = "bold", colour = "black", angle = 0)))  # Rotate x-axis labels if necessary
  ggsave(filename = paste("plot_",colnames(binned_data[a]),".by.predict_low.png", sep = ""))
}

for (a in 4:11) {
print(ggplot(binned_data, aes(x=predict, y=binned_data[,a])) +
        geom_boxplot(aes(x=predict_mid, group=predict_bin), outlier.shape = NA, width=4.5, color="red") +  # Box plots at midpoints
        geom_point(alpha=0.5, shape=1, size=1) +  # Points at actual x-values
        theme_bw() +
        labs(x="predict DPI",
             y=paste(colnames(binned_data[a]), "(normalized mean expression)")) +
        ylim(-1,5)+
        theme(axis.title.y = element_text(size = 20,face = "bold"),
              axis.text.y = element_text(size = 20,face = "bold", colour = "black"),
              axis.title.x = element_text(size = 20,face = "bold"),
              axis.text.x = element_text(size = 20,face = "bold", colour = "black", angle = 0)))  # Rotate x-axis labels if necessary
ggsave(filename = paste("plot_",colnames(binned_data[a]),".by.predict_low.png", sep = ""))
}
