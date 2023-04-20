## Script for plot Figure 4 and Figure S6
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

#############################
#Figure 4I (R version 4.0.4)#
#############################
#${workDir} is the working directory
#${survivalFileTCGA} is the TCGA data used for KM plot
#${survivalFileGSE31210} is the GSE31210 data used for KM plot
#${survivalFileGSE140343} is the GSE140343 data used for KM plot

#Load required packages
library("survminer")

#Set working directory
setwd(workDir)

#Load TCGA data used for KM plot
load(survivalFileTCGA)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="E1_TCGA_DFS_stage_I_II")

#Load GSE31210 data used for KM plot
load(survivalFileGSE31210)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="E1_GSE31210_DFS_stage_I_II")

#Load GSE140343 data used for KM plot
load(survivalFileGSE140343)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="E1_GSE140343_DFS")

#################################
#Figure 4A S6C (R version 4.0.4)#
#################################
#${workDir} is the working directory
#${ecoFile} is the result of Ecotyper

#Load required packages
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

ecotype <- read.table(file=ecoFile,sep="\t",stringsAsFactors=F,header=T)
E1_df <- ecotype[ecotype[,"Ecotype"]=="E1",]
cell_list <- E1_df[,"ID"]

E1_combine_df <- matrix(unlist(strsplit(outer(cell_list,cell_list,paste)," ")),byrow=T,ncol=2)
nichenet_cellchat_res_list <- lapply(seq(nrow(E1_combine_df)),function(x){
	fl <- paste0(E1_combine_df[x,1],"@",E1_combine_df[x,2],"@nicheNet_cellChat.RData")
	if(file.exists(fl))
	{
		load(fl)
		ligand_cell <- rep(E1_combine_df[x,1],nrow(nicheNet_cellChat_con_df))
		receptor_cell <- rep(E1_combine_df[x,2],nrow(nicheNet_cellChat_con_df))
		res_df <- cbind(ligand_cell,cbind(receptor_cell,nicheNet_cellChat_con_df))
		res_df <- res_df[order(res_df[,"cellChat"],decreasing=T),]
		if(nrow(res_df)>10)
		{
			res_df <- res_df[1:10,]
		}
		return(res_df)
	}
})

nichenet_cellchat_res_df <- do.call(rbind,nichenet_cellchat_res_list)

label <- paste(nichenet_cellchat_res_df[,3],nichenet_cellchat_res_df[,4],sep="_")
label_value <- nichenet_cellchat_res_df[,5]*nichenet_cellchat_res_df[,6]
label_tab <- table(label)
label_U <- names(label_tab)[label_tab==1]
label_M <- names(label_tab)[label_tab>1]

nichenet_cellchat_res_dfU <- nichenet_cellchat_res_df[label%in%label_U,]

nichenet_cellchat_res_df_de <- c()
for(i in 1:length(label_M))
{
	index_M <- which(label==label_M[i])
	value_M <- label_value[index_M]
	order_sel <- index_M[which(value_M==max(value_M))[1]]
	nichenet_cellchat_res_df_de <- rbind(nichenet_cellchat_res_df_de,nichenet_cellchat_res_df[order_sel,])
}

nichenet_cellchat_res_df <- rbind(nichenet_cellchat_res_dfU,nichenet_cellchat_res_df_de)

mt1 <- matrix(0,nrow=length(unique(paste(nichenet_cellchat_res_df[,1],nichenet_cellchat_res_df[,3],sep="_"))),
ncol=length(unique(paste(nichenet_cellchat_res_df[,2],nichenet_cellchat_res_df[,4],sep="_"))))
rownames(mt1) <- unique(paste(nichenet_cellchat_res_df[,1],nichenet_cellchat_res_df[,3],sep="_"))
colnames(mt1) <- unique(paste(nichenet_cellchat_res_df[,2],nichenet_cellchat_res_df[,4],sep="_"))

mt2 <- mt1

for(i in seq(nrow(nichenet_cellchat_res_df)))
{
	mt1[paste(nichenet_cellchat_res_df[i,1],nichenet_cellchat_res_df[i,3],sep="_"),
	paste(nichenet_cellchat_res_df[i,2],nichenet_cellchat_res_df[i,4],sep="_")] <- nichenet_cellchat_res_df[i,5]
	mt2[paste(nichenet_cellchat_res_df[i,1],nichenet_cellchat_res_df[i,3],sep="_"),
	paste(nichenet_cellchat_res_df[i,2],nichenet_cellchat_res_df[i,4],sep="_")] <- nichenet_cellchat_res_df[i,6]
}

col_fun2 = colorRamp2(c(0,max(mt1)), c( "white", "purple"))
col_fun = colorRamp2(c( 0,max(mt2)), c("white","red" ))

mt1_r_label <- sapply(strsplit(rownames(mt1),"_"),function(x) paste(x[1],x[2],sep="_"))
mt1_c_label <- sapply(strsplit(colnames(mt1),"_"),function(x) paste(x[1],x[2],sep="_"))

mt1_r_order <- order(mt1_r_label)
mt1_c_order <- order(mt1_c_label)

mt1_order <- mt1[mt1_r_order,mt1_c_order]
mt2_order <- mt1[mt1_r_order,mt1_c_order]

row_ha = rowAnnotation(sender = mt1_r_label[mt1_r_order],
col = list(B.cells_S01 = "#00111111",
CD4.T.cells_S01 = "#00333333",
Endothelial.cells_S02 = "#00555555",
Fibroblasts_S02 = "#00777777",
Mast.cells_S04 = "#00999999",
Monocytes.and.Macrophages_S03 = "#00AAAAAA",
NK.cells_S03 = "#00CCCCCC"))

column_ha = HeatmapAnnotation(receiver = mt1_c_label[mt1_c_order],
col = list(B.cells_S01 = "#00111111",
CD4.T.cells_S01 = "#00333333",
Endothelial.cells_S02 = "#00555555",
Fibroblasts_S02 = "#00777777",
Mast.cells_S04 = "#00999999",
Monocytes.and.Macrophages_S03 = "#00AAAAAA",
NK.cells_S03 = "#00CCCCCC"))

rownames(mt1_order) <- sapply(strsplit(rownames(mt1_order),"_"),function(x) x[3])
colnames(mt1_order) <- sapply(strsplit(colnames(mt1_order),"_"),function(x) x[3])

rownames(mt2_order) <- sapply(strsplit(rownames(mt2_order),"_"),function(x) x[3])
colnames(mt2_order) <- sapply(strsplit(colnames(mt2_order),"_"),function(x) x[3])

#Figure S6C
Heatmap(mt1_order, width = unit(35, "cm"), height = unit(35, "cm"),top_annotation = column_ha, left_annotation = row_ha,
name = "NicheNet", col = col_fun2,cluster_rows = FALSE,cluster_row_slices = FALSE,
 cluster_columns = FALSE,cluster_column_slices = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
            grid.polygon(
                unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
				 gp=gpar(fill=col_fun2(mt1_order[i, j]),col="white"))
			grid.polygon(
                unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
				 gp=gpar(fill=col_fun(mt2_order[i,j]),col="white"))
	}
)

#Figure 4A
Heatmap(mt2_order, width = unit(35, "cm"), height = unit(35, "cm"),top_annotation = column_ha, left_annotation = row_ha,
name = "cellChat", col = col_fun,cluster_rows = FALSE,cluster_row_slices = FALSE,
 cluster_columns = FALSE,cluster_column_slices = FALSE)

#############################
#Figure S6A (R version 3.6.1)#
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

#For LME01, plot heatmap for feature genes of cell state
k=1
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