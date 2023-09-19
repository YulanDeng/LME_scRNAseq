## Script for plot Figure 3 and Figure S4
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

###############################
#Figure 3a (R version R.4.0.3)#
###############################
#${ggoFile} is the seurat object of all the file
#${workDir} is the working directory
#${rankFile} is the rank date from Ecotyper
#${EcoyterDir} is the directory of Ecotyper result

#Load required packages
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(gplots)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load the seurat object
ggo.integrated <- readRDS(file=ggoFile)

#gene sigature
gene_list <- list()
gene_list[[1]] <- c("ZNF683","EOMES","IL7R","TCF7","SELL","CCR7","LEF1")
names(gene_list)[1] <- c("Memory markers")
gene_list[[2]] <- c("ENTPD1","ID2","TOX","HAVCR2","LAG3","PDCD1","TIGIT")
names(gene_list)[2] <- c("Exhaustion markers")
gene_list[[3]] <- c("CD28","PRF1","CD226","MKI67","GZMB","IFNG")
names(gene_list)[3] <- c("Other markers")
gene_list[[4]] <- c("GZMA","INPP4B","IL32","PLEKHF1","SIT1","MYL12A",
"ITM2B","IL16","TRAF3IP3","MTRNR2L12","PTGDR","FCRL3",
"RNPC3","ITGB2","LAMP1","KLF12","ZNF302","PLEK",
"FCMR","KLRG1","NKG7","GZMK","CD27","GIMAP7",
"CMC1","LIMD2","IFITM1","CHI3L2","PLAAT4","DENND2D")
names(gene_list)[4] <- c("Progenitor exhausted cluster markers")


rank_data <- read.table(file=rankFile,sep="\t",stringsAsFactors=F,header=T)

	x=3
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	ecotyper_df=read.table(file =paste0(EcoyterDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	
	label <- rep("unassigned",nrow(ggo.integrated@meta.data))
	names(label)=rownames(ggo.integrated@meta.data)
	label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
	ggo.integrated$Ecotype = label

	mono.ecotyper.integrated <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[ggo.integrated@meta.data[,"Ecotype"]!="unassigned"]) 
	
	DefaultAssay(mono.ecotyper.integrated)<-"RNA"
	mono.ecotyper.integrated <- ScaleData(mono.ecotyper.integrated)
	mt <- mono.ecotyper.integrated@assays$RNA@"scale.data"
	#quantile(mt)
	
#heatmap 1
	mt1 <- mt[gene_list[[1]],]
	mt1sum <- apply(mt1,1,function(x){
		tmp_list <- split(x,f=factor(mono.ecotyper.integrated@meta.data[,"Ecotype"]))
		res <- sapply(tmp_list,mean)
		return(res)
	})
	legend_col <- c(-1,-1,1,1)
	mt1sum <- cbind(mt1sum,legend_col)
	
	pheatmap(mt1sum, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), border_color = "grey60",
  cellwidth = 10, cellheight = 10, scale = "none", cluster_rows = FALSE,
  cluster_cols = FALSE,legend = TRUE, show_rownames = T, show_colnames = T,
  filename = "CD8_Memory_markers.pdf")
  
##heatmap2
  mt2 <- mt[gene_list[[2]],]
	mt2sum <- apply(mt2,1,function(x){
		tmp_list <- split(x,f=factor(mono.ecotyper.integrated@meta.data[,"Ecotype"]))
		res <- sapply(tmp_list,mean)
		return(res)
	})
	legend_col <- c(-1,-1,1,1)
	mt2sum <- cbind(mt2sum,legend_col)
	
	pheatmap(mt2sum, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), border_color = "grey60",
  cellwidth = 10, cellheight = 10, scale = "none", cluster_rows = FALSE,
  cluster_cols = FALSE,legend = TRUE, show_rownames = T, show_colnames = T,
  filename = "CD8_Exhaustion_markers.pdf")
  
##heatmap3
  mt3 <- mt[gene_list[[3]],]
	mt3sum <- apply(mt3,1,function(x){
		tmp_list <- split(x,f=factor(mono.ecotyper.integrated@meta.data[,"Ecotype"]))
		res <- sapply(tmp_list,mean)
		return(res)
	})
	legend_col <- c(-1,-1,1,1)
	mt3sum <- cbind(mt3sum,legend_col)
	
	pheatmap(mt3sum, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), border_color = "grey60",
  cellwidth = 10, cellheight = 10, scale = "none", cluster_rows = FALSE,
  cluster_cols = FALSE,legend = TRUE, show_rownames = T, show_colnames = T,
  filename = "CD8_Other_markers.pdf")
  
##heatmap4
  mt4 <- mt[gene_list[[4]],]
	mt4sum <- apply(mt4,1,function(x){
		tmp_list <- split(x,f=factor(mono.ecotyper.integrated@meta.data[,"Ecotype"]))
		res <- sapply(tmp_list,mean)
		return(res)
	})
	legend_col <- c(-1,-1,1,1)
	mt4sum <- cbind(mt4sum,legend_col)
	
	pheatmap(mt4sum, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), border_color = "grey60",
  cellwidth = 10, cellheight = 10, scale = "none", cluster_rows = FALSE,
  cluster_cols = FALSE,legend = TRUE, show_rownames = T, show_colnames = T,
  filename = "CD8_Progenitor_exhausted_cluster_markers.pdf")


###############################
#Figure 3e (R version R.4.0.3)#
###############################
#${ggoFile} is the seurat object of all the file
#${workDir} is the working directory
#${rankFile} is the rank date from Ecotyper
#${EcoyterDir} is the directory of Ecotyper result

#Load required packages
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(gplots)

#Set working directory
setwd(workDir)

#Load the clinical information and rename
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)

#Load the seurat object of all the cell
ggo.integrated <- readRDS(file=ggoFile)

#Load T cell signatures
load("Tsig_list.RData")

rank_data <- read.table(file=rankFile,sep="\t",stringsAsFactors=F,header=T)

	x=3
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	ecotyper_df=read.table(file =paste0(EcoyterDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	
	label <- rep("unassigned",nrow(ggo.integrated@meta.data))
	names(label)=rownames(ggo.integrated@meta.data)
	label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
	ggo.integrated$Ecotype = label

	mono.ecotyper.integrated <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[ggo.integrated@meta.data[,"Ecotype"]!="unassigned"]) 

	DefaultAssay(mono.ecotyper.integrated)<-"RNA"
	
	mono.ecotyper.integrated <- AddModuleScore(object = mono.ecotyper.integrated, 
	features = Tsig_list, name = "Tsig")

#violin plot
	VlnPlot(mono.ecotyper.integrated, features = paste0("Tsig",1:9), ncol=3,   pt.size = 0,group.by="Ecotype",cols=brewer.pal(8,"Pastel1"))

#--------------------------------------------coding-------------------------------------------------
################################
#Figure S4e (R version R.4.1.3)#
################################
#${ggoFile} is the seurat object of all the file
#${workDir} is the working directory
#${CD4TsigFile} is the signatures of CD4 T cells
#${CD4TmonocleDir} is the directory of monocle3 results for CD4 T cells

#Load required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(magrittr)

#Set working directory
setwd(workDir)

#Load monocle3 result and signature of CD4 T cells
CD4T_P1 <- readRDS(file=CD4TmonocleDir)
load(CD4TsigFile)

DefaultAssay(CD4T_P1)<-"RNA"	
CD4T_P1 <- AddModuleScore(object = CD4T_P1, features =CD4Tsig_list, name = "CD4Tsig")


CD4T_P1_path1 <- subset(CD4T_P1,cells=rownames(CD4T_P1@meta.data)[as.character(CD4T_P1@meta.data[,"monocle3_clusters"])%in%c("23","6","18","13","11","24","5","14","13","17","19","21",
"4","12","10","9","7")]) 
CD4T_P1_path2 <- subset(CD4T_P1,cells=rownames(CD4T_P1@meta.data)[as.character(CD4T_P1@meta.data[,"monocle3_clusters"])%in%c("23","6","18","13","24","5","14",
"16","1")]) 
CD4T_P1_path3 <- subset(CD4T_P1,cells=rownames(CD4T_P1@meta.data)[as.character(CD4T_P1@meta.data[,"monocle3_clusters"])%in%c("23","6","18","13","11","24","5","14","22",
"20","22","8","25","15","3","2")]) 


Tex_cyto_df_path1 <- data.frame(pseudotime=CD4T_P1_path1@meta.data[,"monocle3_pseudotime"],
	Th1=CD4T_P1_path1@meta.data[,"CD4Tsig3"],TFH=CD4T_P1_path1@meta.data[,"CD4Tsig2"],
	Treg=CD4T_P1_path1@meta.data[,"CD4Tsig1"],Th2=CD4T_P1_path1@meta.data[,"CD4Tsig4"],colour="path1")

Tex_cyto_df_path2 <- data.frame(pseudotime=CD4T_P1_path2@meta.data[,"monocle3_pseudotime"],
	Th1=CD4T_P1_path2@meta.data[,"CD4Tsig3"],TFH=CD4T_P1_path2@meta.data[,"CD4Tsig2"],
	Treg=CD4T_P1_path2@meta.data[,"CD4Tsig1"],Th2=CD4T_P1_path2@meta.data[,"CD4Tsig4"],colour="path2")
	
Tex_cyto_df_path3 <- data.frame(pseudotime=CD4T_P1_path3@meta.data[,"monocle3_pseudotime"],
	Th1=CD4T_P1_path3@meta.data[,"CD4Tsig3"],TFH=CD4T_P1_path3@meta.data[,"CD4Tsig2"],
	Treg=CD4T_P1_path3@meta.data[,"CD4Tsig1"],Th2=CD4T_P1_path3@meta.data[,"CD4Tsig4"],colour="path3")

	
Tex_cyto_df <- rbind(Tex_cyto_df_path1,rbind(Tex_cyto_df_path2,Tex_cyto_df_path3))

#plot the results for CD4T_Th1
ggplot(Tex_cyto_df, aes(pseudotime, Th1)) +
	geom_smooth(aes(x=pseudotime,y=Th1, colour=colour),method='lm', formula=y ~ splines::bs(x, 3),se=F)+
	theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18))+ ylab("Signature score")

#plot the results for CD4T_Th2
ggplot(Tex_cyto_df, aes(pseudotime, Th2)) +
	geom_smooth(aes(x=pseudotime,y=Th2, colour=colour),method='lm', formula=y ~ splines::bs(x, 3),se=F)+
	theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18))+ ylab("Signature score")

#plot the results for Treg
ggplot(Tex_cyto_df, aes(pseudotime, Treg)) +
	geom_smooth(aes(x=pseudotime,y=Treg, colour=colour),method='lm', formula=y ~ splines::bs(x, 3),se=F)+
	theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18))+ ylab("Signature score")
