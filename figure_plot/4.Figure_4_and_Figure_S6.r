## Script for plot Figure 4 and Figure S6
##Yulan Deng, last updated 2022-6-28

#############################
#Figure 4B (R version 4.1.1)#
#############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${moFile} is the seurat object of monocyte and macrophage
#${BFile} is the seurat object of B cells
#${CD4Tfile} is the seurat object of CD4+ T cell
#${CD8Tfile} is the seurat object of CD8+ T cell

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

#monocyte and macrophage
B.integrated <- readRDS(file = moFile)
DefaultAssay(B.integrated)<-"RNA"
b_embed <- B.integrated@"reductions"[["umap"]]@"cell.embeddings"
colo <- as.character(B.integrated@meta.data[,"assigned_cell_type_byNC"])
b_df <- data.frame(UMAP_1=b_embed[,"UMAP_1"],UMAP_2=b_embed[,"UMAP_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(UMAP_1, UMAP_2, colour = colo))+
geom_point() +scale_colour_manual(values=c(brewer.pal(8, "Set2")[1:3],"gray",brewer.pal(8, "Set2")[5:8]))+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#B cells
B.integrated <- readRDS(file = BFile)
DefaultAssay(B.integrated)<-"RNA"
b_embed <- B.integrated@"reductions"[["umap"]]@"cell.embeddings"
colo <- as.character(B.integrated@meta.data[,"assigned_cell_type_byNC"])
colo[colo=="MALT B cells"] <- "Follicular B cells"
b_df <- data.frame(UMAP_1=b_embed[,"UMAP_1"],UMAP_2=b_embed[,"UMAP_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(UMAP_1, UMAP_2, colour = colo))+
geom_point() +scale_colour_manual(values=c(brewer.pal(3, "Set2")[1:2],"gray"))+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#CD4+ T cells
B.integrated <- readRDS(file = CD4Tfile)
DefaultAssay(B.integrated)<-"RNA"
b_embed <- B.integrated@"reductions"[["umap"]]@"cell.embeddings"
colo <- as.character(B.integrated@meta.data[,"assigned_cell_type_byNC"])
b_df <- data.frame(UMAP_1=b_embed[,"UMAP_1"],UMAP_2=b_embed[,"UMAP_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(UMAP_1, UMAP_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())

#CD8+ T cells
B.integrated <- readRDS(file = CD8Tfile)
DefaultAssay(B.integrated)<-"RNA"
b_embed <- B.integrated@"reductions"[["umap"]]@"cell.embeddings"
colo <- as.character(B.integrated@meta.data[,"assigned_cell_type_byNC"])
b_df <- data.frame(UMAP_1=b_embed[,"UMAP_1"],UMAP_2=b_embed[,"UMAP_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(UMAP_1, UMAP_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#############################
#Figure 4C (R version 3.6.1)#
#############################
#${workDir} is the working directory
#${moFile} is the seurat object of monocyte and macrophage
#${BFile} is the seurat object of B cells
#${CD4Tfile} is the seurat object of CD4+ T cell
#${CD8Tfile} is the seurat object of CD8+ T cell
#${EcoMoFile} is the Ecotyper result of monocyte and macrophage
#${EcoBFile} is the Ecotyper result of B cells
#${EcoCD4Tfile} is the Ecotyper result of CD4+ T cells
#${EcoCD8Tfile} is the Ecotyper result of CD8+ T cells

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(RColorBrewer)
library(riverplot)

#Set working directory
setwd(workDir)

#monocyte and macrophage
immune.integrated <- readRDS(file = moFile)
label_new <- as.character(immune.integrated@meta.data[,"assigned_cell_type_byNC"])
label_new[label_new%in%c("CD14_monocyte","CD16_monocyte")] <- "monocyte"
label_new[label_new%in%c("TAM_anti_inflammatory")] <- "SPP1_macrophage"
label_new[label_new%in%c("cDC")] <- "DC_like_monocyte"
immune.integrated$label_new<- label_new
ecotyper_anno <- read.table(file =EcoMoFile,sep="\t",stringsAsFactors=F, header=T)
mono.ecotyper.integrated <- subset(immune.integrated,cells=rownames(immune.integrated@meta.data)[rownames(immune.integrated@meta.data)%in%ecotyper_anno[,1]]) 
label2 <- ecotyper_anno[,2]
names(label2) <- ecotyper_anno[,1]
mono.ecotyper.integrated$ecotyper <- label2[rownames(mono.ecotyper.integrated@meta.data)]
stat <- table(mono.ecotyper.integrated@meta.data[,c("label_new","ecotyper")])
edge_cn <- matrix(unlist(strsplit(outer(c("DC_like_monocyte","early_stage_macrophages","monocyte","SPP1_macrophage",
"tissue_resident_Macrophage"),paste0("S0",1:5),paste)," ")),ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(5,5)), y=rep(1:5,2))
cols <- c(brewer.pal(5, "Pastel1")[c(4,5,5,3,1)],brewer.pal(5, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) 
plot(r, plot_area=1, nodewidth=10, default_style=d)

#B cells
immune.integrated <- readRDS(file = BFile)
ecotyper_anno <- read.table(file =EcoBFile,sep="\t",stringsAsFactors=F, header=T)
B.ecotyper.integrated <- subset(immune.integrated,cells=rownames(immune.integrated@meta.data)[rownames(immune.integrated@meta.data)%in%ecotyper_anno[,1]]) 
labelF <- as.character(B.ecotyper.integrated@meta.data[,"assigned_cell_type_byNC"])
labelF[labelF%in%c("MALT B cells")] <- "Plasma cells"
labelF[labelF%in%c("Plasma cells")] <- "Plasma_cells"
labelF[labelF%in%c("Follicular B cells")] <- "Follicular_B_cells"
B.ecotyper.integrated$"new_label" <- labelF
label2 <- ecotyper_anno[,2]
names(label2) <- ecotyper_anno[,1]
B.ecotyper.integrated$ecotyper <- label2[rownames(B.ecotyper.integrated@meta.data)]
stat <- table(B.ecotyper.integrated@meta.data[,c("new_label","ecotyper")])
edge_cn <- matrix(unlist(strsplit(outer(c("Follicular_B_cells","Plasma_cells"),paste0("S0",1:3),paste)," ")),
ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(2,3)), y=c(1:2,1:3))
cols <- c(brewer.pal(5, "Pastel1")[c(3,1)],brewer.pal(5, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) 
plot(r, plot_area=1, nodewidth=10, default_style=d)

#CD4+ T
immune.integrated <- readRDS(file = CD4Tfile)
DefaultAssay(immune.integrated) <- "RNA"
ecotyper_anno <- read.table(file =EcoCD4Tfile,sep="\t",stringsAsFactors=F, header=T)
CD4T.ecotyper.integrated <- subset(immune.integrated,cells=rownames(immune.integrated@meta.data)[rownames(immune.integrated@meta.data)%in%ecotyper_anno[,1]]) 
label2 <- ecotyper_anno[,2]
names(label2) <- ecotyper_anno[,1]
CD4T.ecotyper.integrated$ecotyper <- label2[rownames(CD4T.ecotyper.integrated@meta.data)]
stat <- table(CD4T.ecotyper.integrated@meta.data[,c("assigned_cell_type_byNC","ecotyper")])
edge_cn <- matrix(unlist(strsplit(outer(c("CD4_Th","CXCL13_CD4T","Naive_CD4T","Treg"),paste0("S0",1:5),paste)," ")),ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(4,5)), y=c(1:4,1:5))
cols <- c(brewer.pal(7, "Pastel1")[c(4,7,2,1)],brewer.pal(5, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) 
plot(r, plot_area=1, nodewidth=10, default_style=d)

#CD8+ T cell
CD8T.integrated <- readRDS(file=CD8Tfile)
ecotyper_anno <- read.table(file =EcoCD8Tfile,sep="\t",stringsAsFactors=F, header=T)
label_new <- as.character(CD8T.integrated@meta.data[,"assigned_cell_type_byNC"])
label_new[label_new=="pre-effector_CD8"] <- "transitional_effector_memory_T"
CD8T.integrated$label_new<- label_new
label <- rep("unassigned",nrow(CD8T.integrated@meta.data))
names(label) = rownames(CD8T.integrated@meta.data)
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S01",1],names(label))] = "S01"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S02",1],names(label))] = "S02"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S03",1],names(label))] = "S03"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S04",1],names(label))] = "S04"
CD8T.integrated$ecoRecovery = label
stat <- table(CD8T.integrated@meta.data[,c("label_new","ecoRecovery")])
edge_cn <- matrix(unlist(strsplit(outer(c("effector_CD8","exhausted_CD8",
"transitional_effector_memory_T"),paste0("S0",1:4),paste)," ")),
ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(3,4)), y=c(1:3,1:4))
cols <- c(brewer.pal(9, "Pastel1")[c(1,2,3)],brewer.pal(4, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) # default style
plot(r, plot_area=1, nodewidth=10, default_style=d)

#############################
#Figure 4D (R version 4.0.3)#
#############################
#${workDir} is the working directory
#${moFile} is the seurat object of monocyte and macrophage
#${BFile} is the seurat object of B cells
#${CD4Tfile} is the seurat object of CD4+ T cell
#${CD8Tfile} is the seurat object of CD8+ T cell
#${clinicalFile} is the file of clinical information 

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

#Set working directory
setwd(workDir)

#SPP1+macrophage
CD8T.integrated <-readRDS(file =  moFile)
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical[clinical[,"class"]=="s25GGO","class"] <- "GGO25"
clinical[clinical[,"class"]=="s50GGO","class"] <- "GGO50"
clinical[clinical[,"class"]=="s75GGO","class"] <- "GGO75"
clinical[clinical[,"class"]=="s100GGO","class"] <- "GGO100"
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NoAden"
CD8T.integrated$"stage" <- factor(clinical[match(as.character(CD8T.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
stage_list <- split(CD8T.integrated@meta.data[,c("assigned_cell_type_byNC","orig.ident")],
f=factor(CD8T.integrated@meta.data[,"stage"])) 
stage_fraction_list <- lapply(stage_list,function(x) {
	mt_tmp <- table(x)
	mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"doublets"),]
	fraction_tmp <- apply(mt_tmp,2,function(y) y/sum(y))
	return(fraction_tmp)
})
names(stage_fraction_list) <- names(stage_list)
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN")
sampleName <- unlist(lapply(stage_fraction_list[names(colo)],function(x) colnames(x)))
mono_S3_df <- data.frame(sample=sampleName,
stage=clinical[match(sampleName,clinical[,"rawDataID"]),"class"],
fraction=unlist(lapply(stage_fraction_list[names(colo)],function(x) unlist(x[5,]))),stringsAsFactors=F)
mono_S3_df[,"stage"] <- factor(mono_S3_df[,"stage"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
my_comparisons <- list( c("Normal", "SolidN"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("the fraction of moMac SPP1 Macrophage")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#B cells
CD8T.integrated <- readRDS(file = BFile)
CD8T.integrated$"stage" <- factor(clinical[match(as.character(CD8T.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
labelF <- as.character(CD8T.integrated@meta.data[,"assigned_cell_type_byNC"])
labelF[labelF%in%c("MALT B cells")] <- "Plasma cells"
labelF[labelF%in%c("Plasma cells")] <- "Plasma_cells"
labelF[labelF%in%c("Follicular B cells")] <- "Follicular_B_cells"
CD8T.integrated$"labelF" <- labelF
stage_list <- split(CD8T.integrated@meta.data[,c("labelF","orig.ident")],
f=factor(CD8T.integrated@meta.data[,"stage"])) 
stage_fraction_list <- lapply(stage_list,function(x) {
	mt_tmp <- table(x)
	mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"Undetermined"),]
	fraction_tmp <- apply(mt_tmp,2,function(y) y/sum(y))
	return(fraction_tmp)
})
names(stage_fraction_list) <- names(stage_list)
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN")
sampleName <- unlist(lapply(stage_fraction_list[names(colo)],function(x) colnames(x)))
mono_S3_df <- data.frame(sample=sampleName,
stage=clinical[match(sampleName,clinical[,"rawDataID"]),"class"],
fraction=unlist(lapply(stage_fraction_list[names(colo)],function(x) unlist(x[2,]))),stringsAsFactors=F)
mono_S3_df[,"stage"] <- factor(mono_S3_df[,"stage"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
my_comparisons <- list( c("Normal", "SolidN"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("the fraction of Plasma cells")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#CD4+ T cell
CD8T.integrated <- readRDS(file = CD4Tfile)
CD8T.integrated$"stage" <- factor(clinical[match(as.character(CD8T.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
stage_list <- split(CD8T.integrated@meta.data[,c("assigned_cell_type_byNC","orig.ident")],
f=factor(CD8T.integrated@meta.data[,"stage"])) 
stage_fraction_list <- lapply(stage_list,function(x) {
	mt_tmp <- table(x)
	mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"unassigned"),]
	fraction_tmp <- apply(mt_tmp,2,function(y) y/sum(y))
	return(fraction_tmp)
})
names(stage_fraction_list) <- names(stage_list)
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN")
sampleName <- unlist(lapply(stage_fraction_list[names(colo)],function(x) colnames(x)))
mono_S3_df <- data.frame(sample=sampleName,
stage=clinical[match(sampleName,clinical[,"rawDataID"]),"class"],
fraction=unlist(lapply(stage_fraction_list[names(colo)],function(x) unlist(x[4,]))),stringsAsFactors=F)
mono_S3_df[,"stage"] <- factor(mono_S3_df[,"stage"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN"))
my_comparisons <- list( c("pGGO", "SolidN"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("the fraction of Treg")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#CD8+ T cell
CD8T.integrated <-readRDS(file =  CD8Tfile)
CD8T.integrated$"stage" <- factor(clinical[match(as.character(CD8T.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
stage_list <- split(CD8T.integrated@meta.data[,c("assigned_cell_type_byNC","orig.ident")],
f=factor(CD8T.integrated@meta.data[,"stage"])) 
stage_fraction_list <- lapply(stage_list,function(x) {
	mt_tmp <- table(x)
	mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"doublets"),]
	fraction_tmp <- apply(mt_tmp,2,function(y) y/sum(y))
	return(fraction_tmp)
})
names(stage_fraction_list) <- names(stage_list)
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN")
sampleName <- unlist(lapply(stage_fraction_list[names(colo)],function(x) colnames(x)))
mono_S3_df <- data.frame(sample=sampleName,
stage=clinical[match(sampleName,clinical[,"rawDataID"]),"class"],
fraction=unlist(lapply(stage_fraction_list[names(colo)],function(x) unlist(x[2,]))),stringsAsFactors=F)
mono_S3_df[,"stage"] <- factor(mono_S3_df[,"stage"],levels=c("Normal","pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN"))
my_comparisons <- list( c("Normal", "SolidN"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("exhausted CD8T S02")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#############################
#Figure 4E (R version 4.0.3)#
#############################
#${workDir} is the working directory
#${EcoMoFile} is the Ecotyper result of monocyte and macrophage
#${EcoBFile} is the Ecotyper result of B cells
#${EcoCD4Tfile} is the Ecotyper result of CD4+ T cells
#${EcoCD8Tfile} is the Ecotyper result of CD8+ T cells
#${clinicalFile} is the file of clinical information 

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

#Set working directory
setwd(workDir)


fraction_CD8T <- read.table(file=EcoCD8Tfile,sep="\t",stringsAsFactors=F,header=T)
fraction_CD8T <- fraction_CD8T[,!(colnames(fraction_CD8T)%in%c("X10C","X10F","X16C","X16F"))]
fraction_B <- read.table(file=EcoBFile,sep="\t",stringsAsFactors=F,header=T)
fraction_B <- fraction_B[,!(colnames(fraction_B)%in%c("X10C","X10F","X16C","X16F"))]
fraction_Macro <- read.table(file=EcoMoFile,sep="\t",stringsAsFactors=F,header=T)
fraction_Macro <- fraction_Macro[,!(colnames(fraction_Macro)%in%c("X10C","X10F","X16C","X16F"))]
fraction_CD4T <- read.table(file=EcoCD4Tfile,sep="\t",stringsAsFactors=F,header=T)
fraction_CD4T <- fraction_CD4T[,!(colnames(fraction_CD4T)%in%c("X10C","X10F","X16C","X16F"))]
clinical_info <- read.table(file=clinicalFile,
sep="\t",stringsAsFactors=F,header=T)
clinical_info[,1] <- paste0("X",clinical_info[,1],"C")
clinical_info[clinical_info[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical_info[clinical_info[,"class"]=="s25GGO","class"] <- "GGO25"
clinical_info[clinical_info[,"class"]=="s50GGO","class"] <- "GGO50"
clinical_info[clinical_info[,"class"]=="s75GGO","class"] <- "GGO75"
clinical_label <- rep("Normal",ncol(fraction_CD8T))
clinical_label[clinical_label=="HiDenGGO"] <- "dGGO"
clinical_label[clinical_label=="s25GGO"] <- "GGO25"
clinical_label[clinical_label=="s50GGO"] <- "GGO50"
clinical_label[clinical_label=="s75GGO"] <- "GGO75"
names(clinical_label) <- colnames(fraction_CD8T)
clinical_label[clinical_info[,1]] <- clinical_info[,2]
n_clinical_label <- table(clinical_label)
clinical_labelN <- paste0(clinical_label,"(n=",n_clinical_label[clinical_label],")")
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(3,4,5,6,7)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")

#monocyte and macrophage 
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_Macro[3,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of Macro_S03")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#B cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_B[1,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of B_S01")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

##CD4+ T cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_CD4T[1,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of CD4T_S01")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#CD8+ T cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_CD8T[2,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of CD8T S02")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

##############################
#Figure S6E (R version 4.1.1)#
##############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${EndoFile} is the seurat object of endothelium

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

cancer.integrated <- readRDS(file =EndoFile)
b_embed <- cancer.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(cancer.integrated@meta.data[,"assigned_cell_type"])
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

##############################
#Figure S6G (R version 4.1.1)#
##############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${FibroFile} is the seurat object of fibroblast

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

cancer.integrated <- readRDS(file =FibroFile)
b_embed <- cancer.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(cancer.integrated@meta.data[,"assigned_cell_type"])
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

##############################
#Figure S6I (R version 3.6.1)#
##############################
#${workDir} is the working directory
#${endoFile} is the seurat object of endothelium
#${EcoEndoFile} is the Ecotyper result of endothelium

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(RColorBrewer)
library(riverplot)

#Set working directory
setwd(workDir)

DC.integrated <- readRDS(file=endoFile)
ecotyper_anno <- read.table(file =EcoEndoFile,sep="\t",stringsAsFactors=F, header=T)
label <- rep("unassigned",nrow(DC.integrated@meta.data))
names(label) = rownames(DC.integrated@meta.data)
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S01",1],names(label))] = "S01"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S02",1],names(label))] = "S02"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S03",1],names(label))] = "S03"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S04",1],names(label))] = "S04"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S05",1],names(label))] = "S05"
DC.integrated$ecoRecovery = label
stat <- table(DC.integrated@meta.data[,c("assigned_cell_type","ecoRecovery")])
edge_cn <- matrix(unlist(strsplit(outer(c("alveolar_type_I","arteries","lymphatic_EC","scavenging",
"stalk_like_EC","tip_like_EC","tumor_EC"),paste0("S0",1:5),paste)," ")),
ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(7,5)), y=c(1:7,1:5))
cols <- c(brewer.pal(9, "Pastel1")[c(4,5,1,9,3,9,2)],brewer.pal(5, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) 
plot(r, plot_area=1, nodewidth=10, default_style=d)

##############################
#Figure S6J (R version 3.6.1)#
##############################
#${workDir} is the working directory
#${FibroFile} is the seurat object of fibroblast
#${EcoFibroFile} is the Ecotyper result of fibroblast

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(RColorBrewer)
library(riverplot)

#Set working directory
setwd(workDir)

DC.integrated <- readRDS(file=FibroFile)
ecotyper_anno <- read.table(file =EcoFibroFile,sep="\t",stringsAsFactors=F, header=T)
label <- rep("unassigned",nrow(DC.integrated@meta.data))
names(label) = rownames(DC.integrated@meta.data)
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S01",1],names(label))] = "S01"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S02",1],names(label))] = "S02"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S03",1],names(label))] = "S03"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S04",1],names(label))] = "S04"
DC.integrated$ecoRecovery = label
stat <- table(DC.integrated@meta.data[,c("assigned_cell_type","ecoRecovery")])
edge_cn <- matrix(unlist(strsplit(outer(c("COL10A1+_FBs","COL13A1+_matrix_FBs","COL14A1+_matrix_FBs",
"Myofibroblasts","Pericytes","PLA2G2A_FBs"),paste0("S0",1:4),paste)," ")),
ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(6,4)), y=c(1:6,1:4))
cols <- c(brewer.pal(9, "Pastel1")[c(2,3,4,9,1,9)],brewer.pal(4, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1) 
plot(r, plot_area=1, nodewidth=10, default_style=d)
