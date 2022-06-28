## Script for plot Figure 7 and Figure S9
##Yulan Deng, last updated 2022-6-28

#############################
#Figure 7B (R version 3.6.1)#
#############################
#${workDir} is the working directory
#${ggoFile} is the seurat object of all the file
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state
#${clinicalFile} is the file of clinical information
#${EcoFile} is the results of Ecotyper
#${FindMarkerDir} is the result directory of findmarkers

#Load required packages
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(gplots)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load the clinical data
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
#Load the seurat object of all cells
ggo.integrated <- readRDS(file=ggoFile)
#Load the the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load the results of Ecotyper
ecotype <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)

#For LME04, plot heatmap for feature genes of cell state
k=4
	#extract the cell state information for the LME
	E1_cell_state <- ecotype[ecotype[,"Ecotype"]==paste0("E",k),]
	E1_cell <- c()
	for(i in seq(nrow(E1_cell_state)))
	{
		x <- which(rank_data[,1]==E1_cell_state[i,"CellType"])
		cell_type=rank_data[x,1]
		choose_rank=rank_data[x,2]
		ecotyper_df=read.table(file =paste0(cellStateDir,
		cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
		E1_cell_tmp <- rep(E1_cell_state[i,"ID"],sum(ecotyper_df[,"State"]==E1_cell_state[i,"State"]))
		names(E1_cell_tmp) <- ecotyper_df[ecotyper_df[,"State"]==E1_cell_state[i,"State"],"ID"]
		E1_cell <- c(E1_cell,E1_cell_tmp)
	}
	E1.ecotyper.integrated <- subset(ggo.integrated,cells=names(E1_cell))
	label <- rep("0",length(E1_cell))
	names(label) <- rownames(E1.ecotyper.integrated@meta.data)
	label[names(E1_cell)] <- E1_cell
	E1.ecotyper.integrated$Ecotype = label
	
	#Extract matrix of the ecotype
	log2TPMf <-  as.matrix(E1.ecotyper.integrated@assays$RNA@data)
	
	#Extract marker genes  of cell states
	ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
		print(x)
		cell_type=rank_data[x,1]
		choose_rank=rank_data[x,2]
		pathD=cellStateDir
		ecotyper_df=read.table(file =paste0(pathD,
		cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
		geneTopFC <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/gene_info.txt"),
		sep="\t",stringsAsFactors=F,header=T)
		markers1_list <- lapply(seq(length(unique(geneTopFC[,"State"]))),function(y) { 
			res <- geneTopFC[(geneTopFC[,"State"]==paste0("S0",y))&(geneTopFC[,"MaxFC"]>0.9),"Gene"]
			return(res)
		})		
		markers2_list <- lapply(seq(length(unique(geneTopFC[,"State"]))),function(y) { 
			res <- readRDS(file=paste0(FindMarkerDir,cell_type,".S0",
			y,".markerGenes.rds"))
			resT <- rownames(res)[(res[,"avg_logFC"]>log(1.5))&(res[,"p_val_adj"]<0.01)]
			return(resT)
		})
		markers3_list <- sapply(seq(length(unique(geneTopFC[,"State"]))),function(y){
			res <- unique(c(markers1_list[[y]],markers2_list[[y]]))
			return(res)
		})
		names(markers3_list) <- paste0("S0",seq(length(unique(geneTopFC[,"State"]))))
		return(markers3_list)
	})
	names(ecotyper_list) <- rank_data[,1]
	cell_type_list <- unique(rank_data[,1])
	cell_type_gene_list <- lapply(cell_type_list,function(x){
		pathD2=cellStateDir
		cell_type_gene_df=read.table(file =paste0(pathD2,x,"_cell_type_specific_genes.txt"),
		sep="\t",stringsAsFactors=F, header=T)
		gene <- cell_type_gene_df[(cell_type_gene_df[,"FC"]>1)&(cell_type_gene_df[,"Q"]<0.05),"Gene"]
		if(length(gene)>60)
		{
			gene <- gene[1:60]
		}
		return(gene)
	})
	names(cell_type_gene_list) <- cell_type_list
	gene_list <- lapply(seq(nrow(E1_cell_state)),function(x){
			return(c(cell_type_gene_list[E1_cell_state[x,1]][[1]],
			ecotyper_list[E1_cell_state[x,1]][[1]][E1_cell_state[x,2]][[1]])
			)
		})
	
	#order samples of expression matrix by the top gene of cell state
	log2TPMf <- log2TPMf[apply(log2TPMf,1,sum)>0,]
	gene_list <- lapply(gene_list,function(x) intersect(x,rownames(log2TPMf)))
	cellType_log2TPM <- log2TPMf[unlist(gene_list),names(label)]
	order_list <- c()
	n_cluster <- 0
	for(j in seq(nrow(E1_cell_state))){
		mt_tmp <- cellType_log2TPM[,label==E1_cell_state[j,"ID"]]
		gene <- intersect(gene_list[[j]],rownames(log2TPMf))[1]
		order_tmp <- which(label==E1_cell_state[j,"ID"])[order(unlist(mt_tmp[gene,]),decreasing=T)]
		order_list <- c(order_list,order_tmp)
		n_cluster <- c(n_cluster,n_cluster[length(n_cluster)]+sum(label==E1_cell_state[j,"ID"]))
	}
	cellType_log2TPM_order <- cellType_log2TPM[,order_list]
	
	#Set the parameter for heatmap
	annotation_col = data.frame(
		cell_state = factor(label[order_list]))
		rownames(annotation_col) = names(label)[order_list]
	cell_state_colo <- c(brewer.pal(8, "Set1"),brewer.pal(5, "Set2"))
		names(cell_state_colo) <- E1_cell_state[,"ID"]
	ann_colors = list(
			 cell_state = cell_state_colo[unique(label)])	
	ngene <- c(0,cumsum(sapply(gene_list,length)))	
	cellType_log2TPM_order <- t(apply(cellType_log2TPM_order,1,function(x) (x-min(x))/(max(x)-min(x))))
		
	#plot heatmap
		pheatmap(cellType_log2TPM_order,border_color=NA, 
		color=colorRampPalette(c(brewer.pal(6, "Set1")[2],
		"white",
		brewer.pal(6, "Set1")[1]))(50),
		scale="none",cluster_rows=F,cluster_cols=F,
		gaps_col=n_cluster[-c(1,length(n_cluster))],
		gaps_row=ngene[-length(ngene)],
		annotation_col = annotation_col,
		annotation_colors = ann_colors)

#############################
#Figure 7H (R version 3.6.1)#
#############################
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

##############################
#Figure S9E (R version 3.6.1)#
##############################
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