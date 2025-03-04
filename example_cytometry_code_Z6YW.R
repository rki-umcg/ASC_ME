### General information ----
# Title: Establishing a random forest model to predict the maturity of adult, human plasma cells using various cell surface marker
# Author: Tobit D. Steinmetz
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: d.t.steinmetz@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Date created: 19-10-2023
# Date last adjustment: 02-12-2024
# RStudio version: 2023.06.1 Build 524
# References: 

### introduction ----
# this scripts contains the setup strategies and codes for the random forest ASC maturity index (ASC-ME) model
# raw data from flow cytometry files were previously analyzed, processed and imported below 
# 

### This R-script requires the following packages:
library(flowCore)
library(flowAI)
library(ggplot2)
library(ggcyto)
library(flowViz)
library(openCyto)
library(flowWorkspace)
library(gridExtra)
library(writexl)

###Part 1: importing cytometry data, adjusting marker/channel names and gate required cell populations
#
# load file-list of dataset (select folder containing your fcs-files)
files <- list.files(path="D:/Z6YW_set2", pattern = ".FCS", ignore.case = TRUE)
# write into flowSet
fs_all <- read.flowSet(files, path="D:/Z6YW_set2", truncate_max_range = FALSE)
# use autoplot to double check which marker is in which channel
autoplot(fs_all[[1]])
# define gating set (containing all further gating hierarchy)
gs_all<-GatingSet(fs_all)

### rename channels
# adjust/link marker to respective channel/metal
# check/uncheck using "#" according to what is present in the respective dataset
# adjust marker to actual channel
colnames(gs_all)[colnames(gs_all)=="Nd144Di"] <- "CD3"
colnames(gs_all)[colnames(gs_all)=="Sm154Di"] <- "CD14"
colnames(gs_all)[colnames(gs_all)=="Ho165Di"] <- "CD16"
colnames(gs_all)[colnames(gs_all)=="Gd155Di"] <- "CD27"
colnames(gs_all)[colnames(gs_all)=="Yb172Di"] <- "CD38"
colnames(gs_all)[colnames(gs_all)=="Ir191Di"] <- "DNA1"
colnames(fs_all)[colnames(fs_all)=="Ir191Di"] <- "DNA1"
colnames(gs_all)[colnames(gs_all)=="Ir193Di"] <- "DNA2"
colnames(fs_all)[colnames(fs_all)=="Ir193Di"] <- "DNA2"
#colnames(gs_all)[colnames(gs_all)=="Pt198Di"] <- "viability"
#colnames(fs_all)[colnames(fs_all)=="Pt198Di"] <- "viability"
colnames(gs_all)[colnames(gs_all)=="Ce140Di"] <- "Beads"
# maturity marker
colnames(gs_all)[colnames(gs_all)=="Dy163Di"] <- "CD19"
#colnames(gs_all)[colnames(gs_all)=="Nd146Di"] <- "CD20"
#colnames(gs_all)[colnames(gs_all)=="Er167Di"] <- "CD28"
#colnames(gs_all)[colnames(gs_all)=="Nd144Di"] <- "CD44"
colnames(gs_all)[colnames(gs_all)=="Eu153Di"] <- "CD45"
colnames(gs_all)[colnames(gs_all)=="Er168Di"] <- "CD56"
#colnames(gs_all)[colnames(gs_all)=="Sm154Di"] <- "CD69"
colnames(gs_all)[colnames(gs_all)=="Nd150Di"] <- "CD86"
#colnames(gs_all)[colnames(gs_all)=="Sm152Di"] <- "CD95"
#colnames(gs_all)[colnames(gs_all)=="Eu151Di"] <- "CD98"
#colnames(gs_all)[colnames(gs_all)=="Nd150Di"] <- "CD138"
colnames(gs_all)[colnames(gs_all)=="Er167Di"] <- "HLA.DR"
#colnames(gs_all)[colnames(gs_all)=="Yb176Di"] <- "Ki67"
#colnames(gs_all)[colnames(gs_all)=="Er166Di"] <- "Bcl.2"
#colnames(gs_all)[colnames(gs_all)=="Dy162Di"] <- "BAFFR"
#colnames(gs_all)[colnames(gs_all)=="Gd158Di"] <- "BCMA"
#colnames(gs_all)[colnames(gs_all)=="Pr141Di"] <- "TACI"

### set output path used to export scatter plots and end results
output_path="D:/Routput"

### set up gating strategy
# general gating prodecure
# name gate and define boundaries
# check gate in a representative sample, adjust gate if necessary or go further
# chack gate in all samples, adjust gate if necessary or go further
# save scatter plots as quality check
# add cell population to the gating set object
# recompute gating hierarchy (if cell population is required for further sub-gating)

# 1. clean/single gate
g.clean<-polygonGate("DNA1"=c(400,10000,6000,300,60,80),"DNA2"=c(400,10000,14000,1000,150,80), filterId = "g.clean") 
ggcyto(fs_all[[1]],aes(x="DNA1",y="DNA2"))+geom_hex(bins = 256)+geom_gate(g.clean)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh() 
ggcyto(fs_all,aes(x="DNA1",y="DNA2"))+geom_hex(bins = 256)+geom_gate(g.clean)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh() 
ggsave("lymph_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.clean, parent="root") 
recompute(gs_all) 

# 2.live gate
g.live<-rectangleGate("Beads"=c(-10,20),"Event_length"=c(-10,100), filterId = "g.live")
ggcyto(gs_all[[1]],aes(x="Beads",y="Event_length"), subset = "g.clean")+geom_hex(bins = 128)+geom_gate(g.live)+geom_stats(adjust = 1)+geom_hex(bins = 128)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+ggcyto_par_set(limits = "instrument") #check gate
ggcyto(gs_all,aes(x="Beads",y="Event_length"), subset = "g.clean")+geom_hex(bins = 128)+geom_gate(g.live)+geom_stats(adjust = 1)+geom_hex(bins = 128)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+ggcyto_par_set(limits = "instrument") #check gate
ggsave("via_length_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.live, parent="g.clean") #add to gating set
recompute(gs_all) #apply changes

# 3. CD3/CD14 plot
#T cell gate
g.Tcells<-polygonGate("CD3"=c(15,2000,2000,15),"CD14"=c(-10,-10,10,10), filterId = "g.Tcells")
ggcyto(gs_all[[1]],aes(x="CD3",y="CD14"), subset = "g.live")+geom_hex(bins = 512)+geom_stats(adjust = 1)+ggcyto_par_set(limits = "instrument")+geom_gate(g.Tcells)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
# none T cell gate
g.noT<-rectangleGate("CD3"=c(-2,8),"CD14"=c(-2,2e3), filterId = "g.noT")
ggcyto(gs_all[[1]],aes(x="CD3",y="CD14"), subset = "g.live")+geom_hex(bins = 512)+ggcyto_par_set(limits = "instrument")+geom_gate(g.Tcells)+geom_gate(g.noT)+geom_stats()+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggcyto(gs_all,aes(x="CD3",y="CD14"), subset = "g.live")+geom_hex(bins = 512)+ggcyto_par_set(limits = "instrument")+geom_gate(g.Tcells)+geom_gate(g.noT)+geom_stats()+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggsave("CD3_CD19_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.Tcells, parent="g.live") #add to gating set
gs_pop_add(gs_all, g.noT, parent="g.live") #add to gating set
recompute(gs_all) #apply changes

# 4. gating CD16+ cells
g.CD16cells<-polygonGate("CD16"=c(15,15,6000,6000),"CD3"=c(-2,70,200,-2), filterId = "g.CD16cells")
ggcyto(gs_all[[1]],aes(x="CD16",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+geom_stats(adjust = 1)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD16cells)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggcyto(gs_all,aes(x="CD16",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+geom_stats(adjust = 1)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD16cells)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggsave("CD16_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.CD16cells, parent="g.live") #add to gating set
recompute(gs_all) #apply changes

# 5. plasma cell / ASC gate
# ASC gate
g.ASC<-polygonGate("CD38"=c(900,50000,1e6,1e6,50,300),"CD27"=c(30,10,10,7e5,7e5,100), filterId = "g.ASC")
ggcyto(gs_all[[3]],aes(x="CD38",y="CD27"), subset = "g.noT")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.ASC)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
# not-ASC gate
g.not.ASC<-polygonGate("CD38"=c(-10,1000,1000,50,10,-10),"CD27"=c(-10,-10,0.5,50,400,2000), filterId = "g.not.ASC")
ggcyto(gs_all[[1]],aes(x="CD38",y="CD27"), subset = "g.noT")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.ASC)+geom_gate(g.not.ASC)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggcyto(gs_all,aes(x="CD38",y="CD27"), subset = "g.noT")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.ASC)+geom_gate(g.not.ASC)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggsave("CD38_CD27_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.ASC, parent="g.noT") #add to gating set
gs_pop_add(gs_all, g.not.ASC, parent="g.noT") #add to gating set
recompute(gs_all) #apply changes

# 6. B cell gate
g.Bcells<-polygonGate("CD19"=c(7,3000,3000,7),"CD45"=c(40,40,800,800), filterId = "g.Bcells")
ggcyto(gs_all[[1]],aes(x="CD19",y="CD45"), subset = "g.noT")+geom_hex(bins = 256)+geom_stats(adjust = 1)+ggcyto_par_set(limits = "instrument")+geom_gate(g.Bcells)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggcyto(gs_all,aes(x="CD19",y="CD45"), subset = "g.noT")+geom_hex(bins = 256)+geom_stats(adjust = 1)+ggcyto_par_set(limits = "instrument")+geom_gate(g.Bcells)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggsave("CD19_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.Bcells, parent="g.noT") #add to gating set
recompute(gs_all) #apply changes

### gating for control populations
colnames(gs_all)
# Markers present in this dataset: CD19, CD45, CD56, CD86, HLA
# check/uncheck using "#" according to what control populations are necessary for the markers above

#live+CD14+ cells (for CD45- cells)
g.live.CD14pos<-rectangleGate("CD14"=c(5,500),"CD3"=c(-1.5,15), filterId = "g.live.CD14pos")
ggcyto(gs_all[[1]],aes(x="CD14",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD14pos)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggcyto(gs_all,aes(x="CD14",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD14pos)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggsave("CD14_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.live.CD14pos, parent="g.live")

#live+CD56- cells
g.live.CD56neg<-rectangleGate("CD56"=c(-1.5,1),"CD3"=c(-1.5,8000), filterId = "g.live.CD56neg")
ggcyto(gs_all[[1]],aes(x="CD56",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD56neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggcyto(gs_all,aes(x="CD56",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD56neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
ggsave("CD56_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.live.CD56neg, parent="g.live") 

#CD16+CD56+ cells
g.live.CD16CD56<-rectangleGate("CD56"=c(15,20000),"CD16"=c(15,3000), filterId = "g.live.CD16CD56")
ggcyto(gs_all[[1]],aes(x="CD56",y="CD16"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD16CD56)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggcyto(gs_all,aes(x="CD56",y="CD16"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD16CD56)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggsave("CD56_CD16_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.live.CD16CD56, parent="g.live") 

#CD3+HLA- cells
g.T.HLAneg<-rectangleGate("CD3"=c(15,3000),"HLA.DR"=c(-1.5,15), filterId = "g.T.HLAneg")
ggcyto(gs_all[[1]],aes(x="CD3",y="HLA.DR"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.T.HLAneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggcyto(gs_all,aes(x="CD3",y="HLA.DR"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.T.HLAneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggsave("CD3_HLA.DR_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.T.HLAneg, parent="g.live") 

#live+CD138- cells
#g.live.CD138neg<-rectangleGate("CD138"=c(-10,1),"Cell_length"=c(1,300), filterId = "g.live.CD138neg")
#ggcyto(gs_all[[1]],aes(x="CD138",y="Cell_length"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD138neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#gs_pop_add(gs_all, g.live.CD138neg, parent="g.live")
#live+CD138+ cells
#g.live.CD138pos<-rectangleGate("CD138"=c(10,30000),"Cell_length"=c(1,300), filterId = "g.live.CD138pos")
#ggcyto(gs_all[[1]],aes(x="CD138",y="Cell_length"), subset = "g.not.ASC")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD138neg)+geom_gate(g.live.CD138pos)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
#ggcyto(gs_all,aes(x="CD138",y="Cell_length"), subset = "g.not.ASC")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD138neg)+geom_gate(g.live.CD138pos)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
#ggsave("CD138_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.CD138pos, parent="g.live") 

#live+CD86- cells
g.live.CD86neg<-rectangleGate("CD86"=c(-1,1),"CD19"=c(-0.1,10000), filterId = "g.live.CD86neg")
ggcyto(gs_all[[1]],aes(x="CD86",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD86neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggcyto(gs_all,aes(x="CD86",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD86neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
ggsave("CD86_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
gs_pop_add(gs_all, g.live.CD86neg, parent="g.live") 

#live+CD44+ cells
#g.live.CD44pos<-rectangleGate("CD44"=c(0.8,10),"CD3"=c(1.5,100), filterId = "g.live.CD44pos")
#ggcyto(gs_all[[1]],aes(x="CD44",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD44pos)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#gs_pop_add(gs_all, g.live.CD44pos, parent="g.live") 

#live+CD44- cells
#g.live.CD44neg<-rectangleGate("CD44"=c(-10,0.7),"CD3"=c(1.5,100), filterId="g.live.CD44neg")
#ggcyto(gs_all[[1]],aes(x="CD44",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD44pos)+geom_gate(g.live.CD44neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="CD44",y="CD3"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.CD44pos)+geom_gate(g.live.CD44neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("CD44_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.CD44neg, parent="g.live") 

#CD14+CD69+ cells
#g.CD14CD69<-polygonGate("CD69"=c(5,20,5000,5000,5),"CD14"=c(14,-1,-1,2500,2500), filterId = "g.CD14CD69")
#ggcyto(gs_all[[1]],aes(x="CD69",y="CD14"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD14CD69)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="CD69",y="CD14"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD14CD69)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#gs_pop_add(gs_all, g.CD14CD69, parent="g.live") 
#live+CD69- cells
#g.live.CD69neg<-rectangleGate("CD69"=c(-2,0.5),"CD14"=c(-2,2500), filterId = "g.live.CD69neg")
#ggcyto(gs_all[[1]],aes(x="CD69",y="CD14"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD14CD69)+geom_gate(g.live.CD69neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="CD69",y="CD14"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.CD14CD69)+geom_gate(g.live.CD69neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("CD69_CD14_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.CD69neg, parent="g.live") 

#live+Ki67+ cells
#g.live.Ki67pos<-rectangleGate("Ki67"=c(0.5,10),"CD3"=c(1,10000), filterId = "g.live.Ki67pos")
#ggcyto(gs_all[[1]],aes(x="Ki67",y="CD3"), subset = "g.live")+geom_gate(g.live.Ki67pos)+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="Ki67",y="CD3"), subset = "g.live")+geom_gate(g.live.Ki67pos)+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#gs_pop_add(gs_all, g.live.Ki67pos, parent="g.live") 

#live+Ki67- cells
#g.live.Ki67neg<-rectangleGate("Ki67"=c(-1,0.01),"CD3"=c(1,10000), filterId = "g.live.Ki67neg")
#ggcyto(gs_all[[1]],aes(x="Ki67",y="CD3"), subset = "g.live")+geom_gate(g.live.Ki67neg)+geom_gate(g.live.Ki67pos)+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="Ki67",y="CD3"), subset = "g.live")+geom_gate(g.live.Ki67neg)+geom_gate(g.live.Ki67pos)+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("Ki67_CD3_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.Ki67neg, parent="g.live") 

#live+Bcl2- cells
#g.live.BCL2neg<-rectangleGate("Bcl.2"=c(-1,0.01),"CD19"=c(-0.1,10000), filterId = "g.live.BCL2neg")
#ggcyto(gs_all[[1]],aes(x="Bcl.2",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.BCL2neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="Bcl.2",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.BCL2neg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("Bcl.2_CD19_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.BCL2neg, parent="g.live") 

#live+BCMA- cells
#g.live.BCMAneg<-rectangleGate("BCMA"=c(-1,0.01),"CD19"=c(-0.1,10000), filterId = "g.live.BCMAneg")
#ggcyto(gs_all[[1]],aes(x="BCMA",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.BCMAneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="BCMA",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.BCMAneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("BCMA_CD19_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.BCMAneg, parent="g.live") 

#live+TACI- cells
#g.live.TACIneg<-rectangleGate("TACI"=c(-1,0.01),"CD19"=c(-0.1,10000), filterId = "g.live.TACIneg")
#ggcyto(gs_all[[1]],aes(x="TACI",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.TACIneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggcyto(gs_all,aes(x="TACI",y="CD19"), subset = "g.live")+geom_hex(bins = 256)+ggcyto_par_set(limits = "instrument")+geom_gate(g.live.TACIneg)+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()
#ggsave("TACI_CD19_plot.png",height = 4000, width = 6000,units = "px", path = output_path)
#gs_pop_add(gs_all, g.live.TACIneg, parent="g.live") 

recompute(gs_all)

### Part 2: extract Median fluorescence intensity (MedFI) for each sample and cell population
###         calculate relative expression by ASC according to respective control populations
###         quality control output

# defining function to read out MEAN FIs
pop.means <- function(gs_all){
  chnls <- colnames(gs_all)
  res <- colMeans(exprs(gs_all))
  names(res) <- chnls
  res
}
# display all gated populations within this dataset
plot(gs_all)
# summarize all gated populations within this dataset
all.pop<-c("g.live",
           "g.Tcells",
           "g.ASC",
           "g.Bcells",
           "g.live.CD14pos",
           "g.live.CD16CD56",
           "g.live.CD56neg",
        #  "g.live.CD138neg",
        #  "g.live.CD138pos",
           "g.live.CD86neg",
        #  "g.live.CD44neg",
        #  "g.live.CD44pos",
        #  "g.CD14CD69",
        #  "g.live.CD69neg",    
        #  "g.live.Ki67pos",
        #  "g.live.Ki67neg",
        #  "g.live.BCL2neg",
        #  "g.live.BCMAneg",
        #  "g.live.TACIneg",
           "g.T.HLAneg"
          )
# write all MedFI results (all markers, all cell populations) into a table and add ASC counts
results_mean<-gs_pop_get_stats(gs_all, all.pop, type = pop.means)
counts<-gs_pop_get_stats(gs_all, all.pop)
N.col.res_mean=ncol(results_mean)
results<-cbind(counts,results_mean[,3:N.col.res_mean])

### Set up output table with adjust relative expression values
# write new table with samples + ASC counts
N.smpl=length(files)
N.pops=length(all.pop)
line.ASC=3
sampleID<-results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),c(1,3)]

# generate dummy table to fill in for markers not present in dataset
dummy_tab = data.frame(
  CD19 = matrix("NA",N.smpl),
  CD20 = matrix("NA",N.smpl),
  CD28 = matrix("NA",N.smpl),
  CD44 = matrix("NA",N.smpl),
  CD45 = matrix("NA",N.smpl),
  CD56 = matrix("NA",N.smpl),
  CD69 = matrix("NA",N.smpl),
  CD86 = matrix("NA",N.smpl),
  CD95 = matrix("NA",N.smpl),
  CD98 = matrix("NA",N.smpl),
  CD138 = matrix("NA",N.smpl),
  HLA.DR = matrix("NA",N.smpl),
  Ki67 = matrix("NA",N.smpl),
  Bcl.2 = matrix("NA",N.smpl),
  BAFFR = matrix("NA",N.smpl),
  BCMA = matrix("NA",N.smpl),
  TACI = matrix("NA",N.smpl)
)
# define lines of gated pops within all.pop object
print(all.pop)
line.ASC=3
line.Bcells=4
line.Tcells=2
line.CD14=5
line.CD56neg=7
line.CD16CD56=6
line.T.HLAneg=9
#line.Ki67neg=6
#line.Ki67pos=5
line.CD86neg=8
#line.CD138pos=5
#line.CD138neg=6
#line.CD14CD69=10
#line.CD69neg=11
#line.BCL2neg=14
#line.TACIneg=15
#line.BCMAneg=16
#line.CD44neg=9
#line.CD44pos=10

# calculate normalized ASC expression values and write into table
# in general: relative value = (ASC-neg.pop)/(pos.pop-neg.pop)
# for each marker:
#     choose between calculation or dummy based on availability of marker within this dataset

{
# CD19
sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD19"]-results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"CD19"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD19"]-results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"CD19"]))
#sampleID<-cbind(sampleID,dummy_tab[,1])
#names(sampleID)[3] <- "CD19"
# CD20
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD20"]-results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"CD20"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD20"]-results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"CD20"]))
sampleID<-cbind(sampleID,dummy_tab[,2])
names(sampleID)[4] <- "CD20"
# CD28 
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD28"]-results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD28"])/(results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"CD28"]-results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD28"]))
sampleID<-cbind(sampleID,dummy_tab[,3])
names(sampleID)[5] <- "CD28"
# CD44 
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD44"]-results[as.numeric(paste(line.CD44neg+N.pops*0:(N.smpl-1))),"CD44"])/(results[as.numeric(paste(line.CD44pos+N.pops*0:(N.smpl-1))),"CD44"]-results[as.numeric(paste(line.CD44neg+N.pops*0:(N.smpl-1))),"CD44"]))
sampleID<-cbind(sampleID,dummy_tab[,4])
names(sampleID)[6] <- "CD44"
# CD45 
sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD45"]-results[as.numeric(paste(line.CD14+N.pops*0:(N.smpl-1))),"CD45"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD45"]-results[as.numeric(paste(line.CD14+N.pops*0:(N.smpl-1))),"CD45"]))
#sampleID<-cbind(sampleID,dummy_tab[,5])
#names(sampleID)[7] <- "CD45"
# CD56
sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD56"]-results[as.numeric(paste(line.CD56neg+N.pops*0:(N.smpl-1))),"CD56"])/(results[as.numeric(paste(line.CD16CD56+N.pops*0:(N.smpl-1))),"CD56"]-results[as.numeric(paste(line.CD56neg+N.pops*0:(N.smpl-1))),"CD56"]))
#sampleID<-cbind(sampleID,dummy_tab[,6])
#names(sampleID)[8] <- "CD56"
# CD69 
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD69"]-results[as.numeric(paste(line.CD69neg+N.pops*0:(N.smpl-1))),"CD69"])/(results[as.numeric(paste(line.CD14CD69+N.pops*0:(N.smpl-1))),"CD69"]-results[as.numeric(paste(line.CD69neg+N.pops*0:(N.smpl-1))),"CD69"]))
sampleID<-cbind(sampleID,dummy_tab[,7])
names(sampleID)[9] <- "CD69"
# CD86
sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD86"]-results[as.numeric(paste(line.CD86neg+N.pops*0:(N.smpl-1))),"CD86"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"CD86"]-results[as.numeric(paste(line.CD86neg+N.pops*0:(N.smpl-1))),"CD86"]))
#sampleID<-cbind(sampleID,dummy_tab[,8])
#names(sampleID)[10] <- "CD86"
# CD138
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"CD138"]-results[as.numeric(paste(line.CD138neg+N.pops*0:(N.smpl-1))),"CD138"])/(results[as.numeric(paste(line.CD138pos+N.pops*0:(N.smpl-1))),"CD138"]-results[as.numeric(paste(line.CD138neg+N.pops*0:(N.smpl-1))),"CD138"]))
sampleID<-cbind(sampleID,dummy_tab[,11])
names(sampleID)[11] <- "CD138"
# HLA.DR
sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"HLA.DR"]-results[as.numeric(paste(line.T.HLAneg+N.pops*0:(N.smpl-1))),"HLA.DR"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"HLA.DR"]-results[as.numeric(paste(line.T.HLAneg+N.pops*0:(N.smpl-1))),"HLA.DR"]))
#sampleID<-cbind(sampleID,dummy_tab[,12])
#names(sampleID)[12] <- "HLA.DR"
# Ki67
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"Ki67"]-results[as.numeric(paste(line.Ki67neg+N.pops*0:(N.smpl-1))),"Ki67"])/(results[as.numeric(paste(line.Ki67pos+N.pops*0:(N.smpl-1))),"Ki67"]-results[as.numeric(paste(line.Ki67neg+N.pops*0:(N.smpl-1))),"Ki67"]))
sampleID<-cbind(sampleID,dummy_tab[,13])
names(sampleID)[13] <- "Ki67"
# Bcl.2
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"Bcl.2"]-results[as.numeric(paste(line.BCL2neg+N.pops*0:(N.smpl-1))),"Bcl.2"])/(results[as.numeric(paste(line.Tcells+N.pops*0:(N.smpl-1))),"Bcl.2"]-results[as.numeric(paste(line.BCL2neg+N.pops*0:(N.smpl-1))),"Bcl.2"]))
sampleID<-cbind(sampleID,dummy_tab[,14])
names(sampleID)[14] <- "Bcl.2"
# BCMA
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"BCMA"]-results[as.numeric(paste(line.BCMAneg+N.pops*0:(N.smpl-1))),"BCMA"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"BCMA"]-results[as.numeric(paste(line.BCMAneg+N.pops*0:(N.smpl-1))),"BCMA"]))
sampleID<-cbind(sampleID,dummy_tab[,16])
names(sampleID)[15] <- "BCMA"
# TACI
#sampleID<-cbind(sampleID,(results[as.numeric(paste(line.ASC+N.pops*0:(N.smpl-1))),"TACI"]-results[as.numeric(paste(line.TACIneg+N.pops*0:(N.smpl-1))),"TACI"])/(results[as.numeric(paste(line.Bcells+N.pops*0:(N.smpl-1))),"TACI"]-results[as.numeric(paste(line.TACIneg+N.pops*0:(N.smpl-1))),"TACI"]))
sampleID<-cbind(sampleID,dummy_tab[,17])
names(sampleID)[16] <- "TACI"
}

### calculate mean, min and max of normalized values as QC of the results
# check/uncheck markers that apply to respective data set
{
QC_sampleID = data.frame(QC=c("mean","max","min"))
QC_sampleID<-cbind(QC_sampleID, CD19=c(colMeans(sampleID[,3],na.rm = TRUE),max(sampleID[,3],na.rm = TRUE),min(sampleID[,3],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, CD20=c(colMeans(sampleID[,4],na.rm = TRUE),max(sampleID[,4],na.rm = TRUE),min(sampleID[,4],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, CD28=c(colMeans(sampleID[,5],na.rm = TRUE),max(sampleID[,5],na.rm = TRUE),min(sampleID[,5],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, CD44=c(colMeans(sampleID[,6],na.rm = TRUE),max(sampleID[,6],na.rm = TRUE),min(sampleID[,6],na.rm = TRUE)))
QC_sampleID<-cbind(QC_sampleID, CD45=c(colMeans(sampleID[,7],na.rm = TRUE),max(sampleID[,7],na.rm = TRUE),min(sampleID[,7],na.rm = TRUE)))
QC_sampleID<-cbind(QC_sampleID, CD56=c(colMeans(sampleID[,8],na.rm = TRUE),max(sampleID[,8],na.rm = TRUE),min(sampleID[,8],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, CD69=c(colMeans(sampleID[,9],na.rm = TRUE),max(sampleID[,9],na.rm = TRUE),min(sampleID[,9],na.rm = TRUE)))
QC_sampleID<-cbind(QC_sampleID, CD86=c(colMeans(sampleID[,10],na.rm = TRUE),max(sampleID[,10],na.rm = TRUE),min(sampleID[,10],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, CD138=c(colMeans(sampleID[,11],na.rm = TRUE),max(sampleID[,11],na.rm = TRUE),min(sampleID[,11],na.rm = TRUE)))
QC_sampleID<-cbind(QC_sampleID, HLA.DR=c(colMeans(sampleID[,12],na.rm = TRUE),max(sampleID[,12],na.rm = TRUE),min(sampleID[,12],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, Ki67=c(colMeans(sampleID[,13],na.rm = TRUE),max(sampleID[,13],na.rm = TRUE),min(sampleID[,13],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, Bcl.2=c(colMeans(sampleID[,14],na.rm = TRUE),max(sampleID[,14],na.rm = TRUE),min(sampleID[,14],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, BAFFR=c(colMeans(sampleID[,15],na.rm = TRUE),max(sampleID[,15],na.rm = TRUE),min(sampleID[,15],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, BCMA=c(colMeans(sampleID[,16],na.rm = TRUE),max(sampleID[,16],na.rm = TRUE),min(sampleID[,16],na.rm = TRUE)))
#QC_sampleID<-cbind(QC_sampleID, TACI=c(colMeans(sampleID[,17],na.rm = TRUE),max(sampleID[,17],na.rm = TRUE),min(sampleID[,17],na.rm = TRUE)))
}
print(QC_sampleID)

###Part 3: export results and generate overlay plots of control populations for quality control
#
# export results
write_xlsx(QC_sampleID, ".../Routput/QC_results_Z6YW.xlsx")
write_xlsx(results,".../Routput/results_Z6YW.xlsx")
write_xlsx(sampleID,".../Routput/final_results_Z6YW.xlsx")

# plot overlays for each sample and marker
dev.off(dev.list()["RStudioGD"])  #empty plot display
gc()    #empty unused memory
for(n in 1:length(files)) {
  ggcyto(gs_all[[n]],aes(x="CD19",y="CD3"), subset = "g.ASC")+scale_fill_gradient(trans = "sqrt", high = "black")+geom_overlay(data="g.Tcells", size=0.01, color="yellowgreen")+geom_overlay(data="g.Bcells", size=0.01,color="tomato3")+geom_hex(bins=64)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  ggsave(filename = paste("CD19_overlay_",n,".png", sep = ""), height = 800, width = 1000,units = "px", path = output_path)
  ggcyto(gs_all[[n]],aes(x="CD45",y="CD3"), subset = "g.ASC")+scale_fill_gradient(trans = "sqrt", high = "black")+geom_overlay(data="g.live.CD14pos", size=0.01, color="yellowgreen")+geom_overlay(data="g.Bcells", size=0.01,color="tomato3")+geom_hex(bins=64)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  ggsave(filename = paste("CD45_overlay_",n,".png", sep = ""), height = 800, width = 1000,units = "px", path = output_path)
  ggcyto(gs_all[[n]],aes(x="CD56",y="CD19"), subset = "g.ASC")+scale_fill_gradient(trans = "sqrt", high = "black")+geom_overlay(data="g.live.CD56neg", size=0.01, color="yellowgreen")+geom_overlay(data="g.live.CD16CD56", size=0.01,color="tomato3")+geom_hex(binwidth=c(3,3))+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  ggsave(filename = paste("CD56_overlay_",n,".png", sep = ""), height = 800, width = 1000,units = "px", path = output_path)
  ggcyto(gs_all[[n]],aes(x="CD86",y="CD19"), subset = "g.ASC")+scale_fill_gradient(trans = "sqrt", high = "black")+geom_overlay(data="g.live.CD86neg", size=0.01, color="yellowgreen")+geom_overlay(data="g.Bcells", size=0.01,color="tomato3")+geom_hex(bins=64)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  ggsave(filename = paste("CD86_overlay_",n,".png", sep = ""), height = 800, width = 1000,units = "px", path = output_path)
  ggcyto(gs_all[[n]],aes(x="HLA.DR",y="CD19"), subset = "g.ASC")+scale_fill_gradient(trans = "sqrt", high = "black")+geom_overlay(data="g.T.HLAneg", size=0.01, color="yellowgreen")+geom_overlay(data="g.Bcells", size=0.01,color="tomato3")+geom_hex(binwidth=c(3,3))+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_fasinh()+scale_y_flowjo_fasinh()+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  ggsave(filename = paste("HLA_overlay_",n,".png", sep = ""), height = 800, width = 1000,units = "px", path = output_path)
}
