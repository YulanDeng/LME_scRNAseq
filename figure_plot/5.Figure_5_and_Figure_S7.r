## Script for plot Figure 5 and Figure S7
##Yulan Deng, last updated 2022-6-28

#############################
#Figure 5A (R version 3.6.1)#
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

#############################
#Figure 5I (R version 4.0.4)#
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

##############################
#Figure S7A (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${ggoFile} is the seurat object of all the cells
#${EcoFile} is the result of Ecotyper
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#Load required packages
library(CellChat)
library(patchwork)
library(Seurat)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the seurat object of all the cells
ggo.integrated <- readRDS(file=ggoFile)
#Load the result of Ecotyper
ecotye <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)
#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load the information of cell state
ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	res=read.table(file =paste0(cellStateDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	res[,"State"]=paste(cell_type,res[,"State"],sep="_")
	return(res)
})
ecotyper_df <- do.call(rbind,ecotyper_list)
ecotyper_df <- ecotyper_df[ecotyper_df[,"State"]%in%ecotye[,"ID"],]
label <- rep("unassigned",nrow(ggo.integrated@meta.data))
names(label)=rownames(ggo.integrated@meta.data)
label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
ggo.integrated$Ecotype = label
ecotype.integrate <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"Ecotype"])!="unassigned"])

#run cellchat
cellchat <- createCellChat(object = ecotype.integrate, group.by = "Ecotype")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 16) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat) 
cellchat <- aggregateNet(cellchat)

#select cell interaction from monocyte and macrophage to NK cells at different cell state
LR <- colnames(cellchat@net[[2]][1, ,])
m3nk3sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S03","NK.cells_S03" ,]<0.05]
m4nk1sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S04","NK.cells_S01" ,]<0.05]
m4nk2sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S04","NK.cells_S02" ,]<0.05]
m5nk4sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S05","NK.cells_S04" ,]<0.05]
m3sigu <- setdiff(m3nk3sig,c(m4nk1sig,m4nk2sig,m5nk4sig))
m4sigu <- setdiff(c(m4nk1sig,m4nk2sig),c(m3nk3sig,m5nk4sig))
m5sigu <- setdiff(m5nk4sig,c(m3nk3sig,m4nk1sig,m4nk2sig))
LR_sel <- c(m3sigu,m4sigu,m5sigu)
p_raw <- c(cellchat@net[[2]]["Monocytes.and.Macrophages_S03","NK.cells_S03" ,LR_sel],
cellchat@net[[2]]["Monocytes.and.Macrophages_S04","NK.cells_S01" ,LR_sel],
cellchat@net[[2]]["Monocytes.and.Macrophages_S04","NK.cells_S02"  ,LR_sel],
cellchat@net[[2]]["Monocytes.and.Macrophages_S05","NK.cells_S04" ,LR_sel])
p_pseu <- p_raw
p_pseu[p_pseu<=0.001] <- 0.001
logp <- -log10(p_pseu)
prob <- c(cellchat@net[[1]]["Monocytes.and.Macrophages_S03","NK.cells_S03" ,LR_sel],
cellchat@net[[1]]["Monocytes.and.Macrophages_S04","NK.cells_S01" ,LR_sel],
cellchat@net[[1]]["Monocytes.and.Macrophages_S04","NK.cells_S02"  ,LR_sel],
cellchat@net[[1]]["Monocytes.and.Macrophages_S05","NK.cells_S04" ,LR_sel])
cluster <- c(rep("Monocytes.and.Macrophages_S03 -> NK.cells_S03",length(LR_sel)),
rep("Monocytes.and.Macrophages_S04 -> NK.cells_S01",length(LR_sel)),
rep("Monocytes.and.Macrophages_S04 -> NK.cells_S02",length(LR_sel)),
rep("Monocytes.and.Macrophages_S05 -> NK.cells_S04",length(LR_sel)))
Term <- factor(c(LR_sel,LR_sel,LR_sel,LR_sel),levels=c("BMP4_BMPR1B_ACVR2B",
"CCL13_CCR1","CCL5_CCR5" ,"CXCL1_CXCR2", "HGF_MET","IGF1_IGF1R","IGF1_ITGA6_ITGB4",  
"KLK3_NTRK1","PDGFC_PDGFRA","PROS1_AXL", "VEGFC_VEGFR3",
 "CD40LG_ITGAM_ITGB2","C3_CR2","C3_ITGAM_ITGB2" ,"TGFB1_TGFBR1_TGFBR2", "TGFB3_TGFBR1_TGFBR2","TGFB1_ACVR1_TGFBR1","TGFB2_ACVR1_TGFBR1",
"TGFB3_ACVR1_TGFBR1","CX3CL1_CX3CR1","IL21_IL21R_IL2RG",
"LTA_TNFRSF14", "F2_F2R","OSM_LIFR_IL6ST","IL1B_IL1R1_IL1RAP","IL1B_IL1R2","CCL20_CCR6","CXCL16_CXCR6",
"TNFSF15_TNFRSF25","RETN_TLR4","RETN_CAP1"))
flres_df <- data.frame(cluster=cluster,
Term=Term,
log10P=logp,Probability=prob/max(prob),stringsAsFactors=F)
pp <- ggplot(flres_df,aes(cluster,Term)) 
#plot
pp + geom_point(aes(size = log10P,colour=Probability))+ylab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
panel.grid.major=element_line(size=0.1,linetype =1,color=brewer.pal(6, "Greys")[2]),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=90))+scale_color_gradient(low = brewer.pal(6, "Reds")[2],high = brewer.pal(6, "Reds")[5])

##############################
#Figure S7B (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${ggoFile} is the seurat object of all the cells
#${EcoFile} is the result of Ecotyper
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#Load required packages
library(CellChat)
library(patchwork)
library(Seurat)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the seurat object of all the cells
ggo.integrated <- readRDS(file=ggoFile)
#Load the result of Ecotyper
ecotye <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)
#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load the information of cell state
ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	res=read.table(file =paste0(cellStateDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	res[,"State"]=paste(cell_type,res[,"State"],sep="_")
	return(res)
})
ecotyper_df <- do.call(rbind,ecotyper_list)
ecotyper_df <- ecotyper_df[ecotyper_df[,"State"]%in%ecotye[,"ID"],]
label <- rep("unassigned",nrow(ggo.integrated@meta.data))
names(label)=rownames(ggo.integrated@meta.data)
label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
ggo.integrated$Ecotype = label
ecotype.integrate <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"Ecotype"])!="unassigned"])

#run cellchat
cellchat <- createCellChat(object = ecotype.integrate, group.by = "Ecotype")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 16) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat) 
cellchat <- aggregateNet(cellchat)

#select cell interatcion among macrophage S03,fibroblast S02  and Endothelial.cells S02 
#within LME01
LR <- colnames(cellchat@net[[2]][1, ,])
m3f2sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S03" ,"Fibroblasts_S02",]<0.05]
f2m3sig <- LR[cellchat@net[[2]]["Fibroblasts_S02" ,"Monocytes.and.Macrophages_S03",]<0.05]
m3e2sig <- LR[cellchat@net[[2]]["Monocytes.and.Macrophages_S03" ,"Endothelial.cells_S02",]<0.05]
e2m3sig <- LR[cellchat@net[[2]]["Endothelial.cells_S02" ,"Monocytes.and.Macrophages_S03",]<0.05]
e2f2sig <- LR[cellchat@net[[2]]["Endothelial.cells_S02" ,"Fibroblasts_S02",]<0.05]
f2e2sig <- LR[cellchat@net[[2]]["Fibroblasts_S02" ,"Endothelial.cells_S02",]<0.05]
m3f2sigu <- setdiff(m3f2sig,c(f2m3sig,m3e2sig,e2m3sig,e2f2sig,f2e2sig))
f2m3sigu <- setdiff(f2m3sig,c(m3f2sig,m3e2sig,e2m3sig,e2f2sig,f2e2sig))
m3e2sigu <- setdiff(m3e2sig,c(f2m3sig,m3f2sig,e2m3sig,e2f2sig,f2e2sig))
e2m3sigu <- setdiff(e2m3sig,c(f2m3sig,m3e2sig,m3f2sig,e2f2sig,f2e2sig))
e2f2sigu <- setdiff(e2f2sig,c(f2m3sig,m3e2sig,e2m3sig,m3f2sig,f2e2sig))
f2e2sigu <- setdiff(f2e2sig,c(f2m3sig,m3e2sig,e2m3sig,e2f2sig,m3f2sig))
LR_sel <- c(m3f2sigu,f2m3sigu,m3e2sigu,e2m3sigu,e2f2sigu,f2e2sigu )
p_raw <- c(cellchat@net[[2]]["Monocytes.and.Macrophages_S03" ,"Fibroblasts_S02",LR_sel],
cellchat@net[[2]]["Fibroblasts_S02" ,"Monocytes.and.Macrophages_S03",LR_sel],
cellchat@net[[2]]["Monocytes.and.Macrophages_S03" ,"Endothelial.cells_S02",LR_sel],
cellchat@net[[2]]["Endothelial.cells_S02" ,"Monocytes.and.Macrophages_S03",LR_sel],
cellchat@net[[2]]["Endothelial.cells_S02" ,"Fibroblasts_S02",LR_sel],
cellchat@net[[2]]["Fibroblasts_S02" ,"Endothelial.cells_S02",LR_sel])
p_pseu <- p_raw
p_pseu[p_pseu<=0.001] <- 0.001
logp <- -log10(p_pseu)
prob <- c(cellchat@net[[1]]["Monocytes.and.Macrophages_S03" ,"Fibroblasts_S02",LR_sel],
cellchat@net[[1]]["Fibroblasts_S02" ,"Monocytes.and.Macrophages_S03",LR_sel],
cellchat@net[[1]]["Monocytes.and.Macrophages_S03" ,"Endothelial.cells_S02",LR_sel],
cellchat@net[[1]]["Endothelial.cells_S02" ,"Monocytes.and.Macrophages_S03",LR_sel],
cellchat@net[[1]]["Endothelial.cells_S02" ,"Fibroblasts_S02",LR_sel],
cellchat@net[[1]]["Fibroblasts_S02" ,"Endothelial.cells_S02",LR_sel])
cluster <- c(rep("Monocytes.and.Macrophages_S03 -> Fibroblasts_S02",length(LR_sel)),
rep("Fibroblasts_S02 -> Monocytes.and.Macrophages_S03",length(LR_sel)),
rep("Monocytes.and.Macrophages_S03 -> Endothelial.cells_S02",length(LR_sel)),
rep("Endothelial.cells_S02 -> Monocytes.and.Macrophages_S03",length(LR_sel)),
rep("Endothelial.cells_S02 -> Fibroblasts_S02",length(LR_sel)),
rep("Fibroblasts_S02 ->Endothelial.cells_S02",length(LR_sel)))
Term <- factor(c(LR_sel,LR_sel),levels=c("MDK_SDC1","PTN_SDC2","PTN_SDC3",
"PTN_SDC4","RARRES2_CMKLR1",
"WNT2_FZD4_LRP5","WNT2_FZD6_LRP5","WNT2_FZD4_LRP6","WNT2_FZD6_LRP6",
"WNT5A_FZD4","WNT5A_FZD6","WNT5A_MCAM","WNT5B_FZD4",
"WNT5B_FZD6","MDK_PTPRZ1","PTN_PTPRZ1","EDN1_EDNRB",
"WNT2B_FZD4_LRP5","WNT2B_FZD4_LRP6",
"IL4_IL4R_IL13RA2","TNF_TNFRSF1B","ADM_CALCRL",
"AREG_EGFR","AREG_EGFR_ERBB2"))
flres_df <- data.frame(cluster=cluster,
Term=Term,
log10P=logp,Probability=prob/max(prob),stringsAsFactors=F)
pp <- ggplot(flres_df,aes(cluster,Term)) 
#plot
pp + geom_point(aes(size = log10P,colour=Probability))+ylab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
panel.grid.major=element_line(size=0.1,linetype =1,color=brewer.pal(6, "Greys")[2]),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=90))+scale_color_gradient(low = brewer.pal(6, "Reds")[2],high = brewer.pal(6, "Reds")[5])
