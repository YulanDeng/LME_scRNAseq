## Script for plot Figure 7 and Figure S9
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

###############################
#Figure 7F-H (R version 3.6.1)#
###############################
#${workDir} is the working directory
#${EcoFileGSE12604} is the results of Ecotyper for GSE12604
#${clinicalFileGSE126044} is the clinical information of GSE126044
#${EcoFileGSE166449} is the results of Ecotyper for GSE166449
#${clinicalFileGSE166449} is the clinical information of GSE166449
#${smpId2_clin} is the sample information for GSE166449
#${EcoFileGSE135222} is the results of Ecotyper for GSE135222
#${clinicalFileGSE135222} is the clinical information of GSE135222
#${smpId1_clin} is the sample information for GSE135222
#${survivalFileGSE135222} is the GSE135222 data used for KM plot

#Load required packages
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library("survminer")

#Set working directory
setwd(workDir)

#Load the results of Ecotyper for GSE12604
fraction <- read.table(file=EcoFileGSE12604,sep="\t",stringsAsFactors=F,header=T)

#Load the clinical information of GSE126044
load(clinicalFileGSE126044)

#Calculate the LME04 score
score11_list <- split(unlist(fraction[4,rownames(Cho_clinic)]),f=factor(Cho_clinic[,"Clinicalbenefit"]))
clinical_label1 <- Cho_clinic[,"Clinicalbenefit"]
score11_df <- data.frame(fraction=unlist(fraction[4,rownames(Cho_clinic)]),
clinical_label=factor(clinical_label1),stringsAsFactors=F)
my_comparisons <- list( c("non-responder","responder"))

#plot the result for GSE12604
p<-ggboxplot(score11_df, x = "clinical_label", y = "fraction", color="clinical_label",outlier.shape = NA,
palette = brewer.pal(3,"Pastel2"),legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("LME4 score")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#Load the results of Ecotyper for GSE166449
fraction <- read.table(file=EcoFileGSE166449,sep="\t",stringsAsFactors=F,header=T)

#Load the clinical data for GSE166449
clini2 <- read.table(file=clinicalFileGSE166449,sep="\t",stringsAsFactors=F,header=F) 
clini2[,2] <- gsub("[0-9]","",clini2[,2])
clinical_label2 <- smpId2_clin[,"response"]
names(clinical_label2) <- smpId2_clin[,"SRR"]
#Calculate the LME04 score
score12_list <- split(unlist(fraction[4,names(clinical_label2)]),f=factor(clinical_label2))
score12_df <- data.frame(fraction=unlist(fraction[4,names(clinical_label2)]),
clinical_label=factor(clinical_label2),stringsAsFactors=F)
my_comparisons <- list( c("Immunotherapy_nonResponder","Immunotherapy_Responder"))
#plot the result for GSE166449
p<-ggboxplot(score12_df, x = "clinical_label", y = "fraction", color="clinical_label",outlier.shape = NA,
palette = brewer.pal(3,"Pastel2"),legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("LME4 score")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#Load the results of Ecotyper for GSE135222
fraction <- read.table(file=EcoFileGSE135222,sep="\t",stringsAsFactors=F,header=T)
#Load the clinical information of GSE135222
clini1 <- read.table(file=smpId1_clin,sep=",",stringsAsFactors=F,header=T) 
clinical_label1 <- smpId1_clin[,"response"]
names(clinical_label1) <- smpId1_clin[,"SRR"]
#Calculate the LME04 score
score11_list <- split(unlist(fraction[4,names(clinical_label1)]),f=factor(clinical_label1))
score11_df <- data.frame(fraction=unlist(fraction[4,names(clinical_label1)]),
clinical_label=factor(clinical_label1),stringsAsFactors=F)

#plot the result for GSE135222
my_comparisons <- list( c("Response", "Non_response"))
p<-ggboxplot(score11_df, x = "clinical_label", y = "fraction", color="clinical_label",outlier.shape = NA,
palette = brewer.pal(3,"Pastel2"),legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("LME4 score")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#Load GSE135222 data used for KM plot
load(survivalFileGSE135222)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,ylab="E4_signature_GSE135222")

###############################
#Figure S10C (R version 3.6.1)#
###############################
#${workDir} is the working directory
#${CD4Tfile} is the seraut object of CD4+ T cells
#${CD4TcellState} is the cell state of CD4+ T cells

#Load required packages
library(monocle)
library(tidyverse)
library(Seurat)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the seurat object and cell state of CD4+ T cells
CD4T.integrated <- readRDS(file = CD4Tfile)
ecotyper_CD4T <- read.table(file =CD4TcellState,sep="\t",stringsAsFactors=F, header=T)

#transfer the seurat obeject to monocle object
CD4T.state.integrated <- subset(CD4T.integrated,cells=intersect(rownames(CD4T.integrated@meta.data),
ecotyper_CD4T[,1])) 
label <- ecotyper_CD4T[match(rownames(CD4T.state.integrated@meta.data),ecotyper_CD4T[,1]),2]
CD4T.state.integrated$"cell_state" <- label
RA_matrix <-as(as.matrix(CD4T.state.integrated@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-CD4T.state.integrated@meta.data
sample_ann<- cbind(sample_ann,label)
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
RA.cds<- newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())

#trajectories analysis
RA.cds <- estimateSizeFactors(RA.cds) 
RA.cds <- estimateDispersions(RA.cds)
RA.cds=detectGenes(RA.cds,min_expr = 1)

CD4_marker_S01 = FindMarkers(CD4T.state.integrated,ident.1 = "S01", group.by = "cell_state" )
CD4_marker_S02 = FindMarkers(CD4T.state.integrated,ident.1 = "S02", group.by = "cell_state" )
CD4_marker_S03 = FindMarkers(CD4T.state.integrated,ident.1 = "S03", group.by = "cell_state" )
CD4_marker_S04 = FindMarkers(CD4T.state.integrated,ident.1 = "S04", group.by = "cell_state" )
CD4_marker_S05 = FindMarkers(CD4T.state.integrated,ident.1 = "S05", group.by = "cell_state" )

marker_gene=rbind(CD4_marker_S01,rbind(CD4_marker_S02,rbind(CD4_marker_S03,rbind(CD4_marker_S04,CD4_marker_S05))))
test_ordering_genes=unique(rownames(marker_gene[(marker_gene[,2]>log(1.6))&(marker_gene[,5]<0.05),]))
RA.cds=setOrderingFilter(RA.cds,ordering_genes = test_ordering_genes) 

RA.cds=reduceDimension(RA.cds,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 

#residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次
RA.cds=orderCells(RA.cds)
RA.cds=orderCells(RA.cds,root_state=4)
plot_cell_trajectory(RA.cds,color_by = "cell_state")
label2 <- RA.cds@"phenoData"@"data"[match(rownames(CD4T.state.integrated@meta.data),
rownames(RA.cds@"phenoData"@"data")),"State"]
CD4T.state.integrated$"monocle" <- label2
CD4_monocle_1 = FindMarkers(CD4T.state.integrated,ident.1 = "1", group.by = "monocle" )
CD4_monocle_2 = FindMarkers(CD4T.state.integrated,ident.1 = "2", group.by = "monocle" )
CD4_monocle_3 = FindMarkers(CD4T.state.integrated,ident.1 = "3", group.by = "monocle" )
CD4_monocle_4 = FindMarkers(CD4T.state.integrated,ident.1 = "4", group.by = "monocle" )
CD4_monocle_5 = FindMarkers(CD4T.state.integrated,ident.1 = "5", group.by = "monocle" )
plot_cell_trajectory(RA.cds,color_by = "Pseudotime")