## Script for plot Figure 2
##Yulan Deng, last updated 2024-2-23
##my e-mail:kndeajs@163.com

#############################
#Figure 2A (R version 4.1.1)#
#############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${clinicalFile} is file of clinical information
#${EpiFile} is the seurat object of all the Epithelium

#Load required packages
library(reticulate)
use_python(pythonDir,required=T)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library('ggrastr')
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the clinical data
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
#Load the seurat object of epithelium
epi.integrated <- readRDS(file = EpiFile)

#label epithlial cells with copy number
label <- as.character(epi.integrated@meta.data[,"orig.ident"])
names(label) <- rownames(epi.integrated@meta.data)
label2 <- clinical[match(label,clinical[,"rawDataID"]),"class"]
names(label2) <- rownames(epi.integrated@meta.data)
label2[label2!="Normal"] <- "precancerous"
label2[grep("c",as.character(epi.integrated@meta.data[,"combAnno"]))] <- "cancer"
epi.integrated$clinical <- label2

#plot for benign patients
nor <- subset(epi.integrated,cells=rownames(epi.integrated@meta.data)[as.character(epi.integrated@meta.data[,"orig.ident"])%in%clinical[clinical[,"class"]=="Normal","rawDataID"]])
b_embed1 <- nor@"reductions"[["tsne"]]@"cell.embeddings"
b_df1 <- data.frame(tSNE_1=b_embed1[,"tSNE_1"],tSNE_2=b_embed1[,"tSNE_2"],stringsAsFactors=F)
p1 <- ggplot(b_df1, aes(tSNE_1, tSNE_2))+
geom_point( size = 0.01,col = brewer.pal(6, "Set3")[5]) +
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p1,scale = 0.25,dpi=300)

#plot for cancer patients
ca <- subset(epi.integrated,cells=rownames(epi.integrated@meta.data)[as.character(epi.integrated@meta.data[,"orig.ident"])%in%clinical[clinical[,"class"]!="Normal","rawDataID"]])
b_embed2 <- ca@"reductions"[["tsne"]]@"cell.embeddings"
b_df2 <- data.frame(tSNE_1=b_embed2[,"tSNE_1"],tSNE_2=b_embed2[,"tSNE_2"],colo=ca@meta.data[,"clinical"],stringsAsFactors=F)
p2 <- ggplot(b_df2, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point(size = 0.01) +scale_colour_manual(values=brewer.pal(7, "Set3")[c(4,7)]) +
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p2,scale = 0.25,dpi=300) 