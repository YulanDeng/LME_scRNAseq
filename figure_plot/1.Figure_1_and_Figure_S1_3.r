## Script for plot Figure 1 and Figure S1-3
##Yulan Deng, last updated 2022-6-27

#############################
#Figure 1C (R version 4.1.1)#
#############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${ggoFile} is the seurat object of all the cell
#${lymFile} is the seurat object of lymphocyte
#${lymAnnoFile} is the annotation of lymphocyte
#${MyeFile} is the seurat object of myeloid
#${MyeAnnoFile} is the annotation of myeloid

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

#major cluster
ggo.integrated <- readRDS(file=ggoFile)
b_embed <- ggo.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(ggo.integrated@meta.data[,"assigned_cell_type"])
names(colo) <- rownames(ggo.integrated@meta.data)
colo[colo%in%c("T","NK","B")] <- "lymphocyte"
colo[colo%in%c("MAST")] <- "Myeloid"
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#non-immune cluster
ggo.integrated$"major0701" <- colo
ggo.non.integrated <- subset(ggo.integrated,
	cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"major0701"])%in%c("Endothelial",
		"Epithelial","Fibroblasts")])
b2_embed <- ggo.non.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo2 <- as.character(ggo.non.integrated@meta.data[,"assigned_cell_type"])
b2_df <- data.frame(tSNE_1=b2_embed[,"tSNE_1"],tSNE_2=b2_embed[,"tSNE_2"],colo=colo2,stringsAsFactors=F)
p <- ggplot(b2_df, aes(tSNE_1, tSNE_2, colour = colo2))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#cluster of lymphocyte
lym.integrated <- readRDS(file=lymFile)
cluster <- read.table(file = lymAnnoFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
lym.integrated$major<- cluster[match(as.character(lym.integrated@meta.data[,"seurat_clusters"]),cluster[,1]),4]
b_embed <- lym.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(lym.integrated@meta.data[,"major"])
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set3")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#cluster of myeloid
mye.integrated <- readRDS(file=MyeFile)
cluster <- read.table(file = MyeAnnoFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
mye.integrated$major<- cluster[match(as.character(mye.integrated@meta.data[,"seurat_clusters"]),cluster[,1]),4]
b_embed <- mye.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(mye.integrated@meta.data[,"major"])
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set3")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#############################
#Figure 1D (R version 3.6.1)#
#############################
#${workDir} is the working directory
#${clinicalFile} is the clinical information
#${ggoFile} is the seurat object of all the cell

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the clinical data
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)

#Load the seurat object of all the cells
ggo.integrated <- readRDS(file=ggoFile)

#basic processing
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NoAden"
ggo.integrated$"stage" <- factor(clinical[match(as.character(ggo.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","pGGO","HiDenGGO","s25GGO",
 "s50GGO","s75GGO","s100GGO","Solid1","Solid3","SolidN","NoAden"))
label <- as.character(ggo.integrated@meta.data[,"assigned_cell_type"])
names(label) <- rownames(ggo.integrated@meta.data)
label[label%in%c("T","NK","B")] <- "lymphocyte"
label[label%in%c("MAST")] <- "Myeloid"
ggo.integrated$"major1217" <- factor(label,levels=unique(label))

#bar plot
major_mt <- ggo.integrated@meta.data[,c("major1217","stage")]
major_mt_num <- table(major_mt)
major_mt_fra <- apply(major_mt_num,2,function(x) x/sum(x))
barplot(major_mt_fra,col=brewer.pal(8, "Set2"),las=2,legend.text = T)

#############################
#Figure 1F (R version 3.6.1)#
#############################
#${workDir} is the working directory
#${ecoFile} is the result of Ecotyper
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#Load required packages
library(RColorBrewer)
library(gplots)

#Set working directory
setwd(workDir)

#Load result of Ecotyper
ecotyper <- read.table(file=ecoFile,sep="\t",stringsAsFactors=F,header=T)
#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load cell state information
ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
	print(x)
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	ecotyper_df=read.table(file =paste0(cellStateDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)	
	return(ecotyper_df)
})
names(ecotyper_list) <- rank_data[,1]

#plot heatmap
mt <- matrix(0,nrow=nrow(rank_data),ncol=length(unique(ecotyper[,"Ecotype"])))
rownames(mt) <- rank_data[,1]
colnames(mt) <- unique(ecotyper[,"Ecotype"])
for(i in seq(nrow(ecotyper)))
{
	mt[ecotyper[i,"CellType"],ecotyper[i,"Ecotype"]] <- mt[ecotyper[i,"CellType"],ecotyper[i,"Ecotype"]] + sum(ecotyper_list[ecotyper[i,"CellType"]][[1]][,"State"]==ecotyper[i,"State"])
}
mtN <- apply(mt,2,function(x) x/sum(x))
heatmap.2(mtN,Rowv = F,Colv=F,dendrogram = "none",scale = "none",trace="none",
col=c("white",colorRampPalette(brewer.pal(5, "Greens")[1:4])(20)))

#############################
#Figure 1G (R version 3.6.1)#
#############################
#${workDir} is the working directory
#${clinicalFile} is the file of clinical information
#${ggoFile} is the seurat object of all the file
#${ecoFile} is the result of Ecotyper
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#Load required packages
library(patchwork)
library(Seurat)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
#Load the seurat object of all the cell
ggo.integrated <- readRDS(file=ggoFile)
#Load the result of Ecotype
ecotye <- read.table(file=ecoFile,sep="\t",stringsAsFactors=F,header=T)
#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load the result of cell states
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

#match LME to clinical stage
label <- rep("unassigned",nrow(ggo.integrated@meta.data))
names(label)=rownames(ggo.integrated@meta.data)
label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
ggo.integrated$Ecotype = label
ecotype.integrate <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"Ecotype"])!="unassigned"])
clin_meta <- ecotype.integrate@meta.data
clini_class <- clinical[match(clin_meta[,"orig.ident"],clinical[,"rawDataID"]),"class"]
eco_class <- ecotye[match(clin_meta[,"Ecotype"],ecotye[,"ID"]),"Ecotype"]
clin_meta <- cbind(cbind(clin_meta,clini_class,stringsAsFactors=F),eco_class,stringsAsFactors=F)

#fraction of LME in GGO group
GGO_mt <- clin_meta[clin_meta[,"clini_class"]%in%c("pGGO","s100GGO","s25GGO","HiDenGGO","s50GGO","s75GGO"),
c("clini_class","eco_class")]
GGO_mt[,"clini_class"] <- factor(GGO_mt[,"clini_class"],levels=c("pGGO","HiDenGGO","s25GGO","s50GGO","s75GGO","s100GGO"))
GGO_mt_num <- table(GGO_mt)
GGO_mt_fra <- apply(GGO_mt_num,1,function(x) x/sum(x))
barplot(GGO_mt_fra,col=brewer.pal(12,"Paired")[c(2,4,6,8,10,12)],legend.text = T)

#fraction of LME in major group
major_mt <- clin_meta[,c("clini_class","eco_class")]
major_mt[major_mt[,"clini_class"]%in%c("pGGO","HiDenGGO","s25GGO","s50GGO","s75GGO","s100GGO"),"clini_class"] <- "GGO"
major_mt[major_mt[,"clini_class"]%in%c("Solid1","Solid3","SolidN"),"clini_class"] <- "Solid"
major_mt[,"clini_class"] <- factor(major_mt[,"clini_class"],levels=c("Normal","GGO","Solid"))
major_mt_num <- table(major_mt)
major_mt_fra <- apply(major_mt_num,1,function(x) x/sum(x))
barplot(major_mt_fra,col=brewer.pal(12,"Paired")[c(2,4,6,8,10,12)],legend.text = T)

#fraction of LME in solid group
solid_mt <- clin_meta[clin_meta[,"clini_class"]%in%c("Solid1","Solid3","SolidN"),
c("clini_class","eco_class")]
solid_mt[,"clini_class"] <- factor(solid_mt[,"clini_class"],levels=c("Solid1","Solid3","SolidN"))
solid_mt_num <- table(solid_mt)
solid_mt_fra <- apply(solid_mt_num,1,function(x) x/sum(x))
barplot(solid_mt_fra,col=brewer.pal(12,"Paired")[c(2,4,6,8,10,12)],legend.text = T)

##############################
#Figure S1A (R version 3.6.1)#
##############################
#${clinicalFile} is the file of clinical information
#${ggoFile} is the seurat object of all the file
#${workDir} is the working directory

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)

#Set working directory
setwd(workDir)

#Load the clinical information and rename
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical[clinical[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical[clinical[,"class"]=="s25GGO","class"] <- "GGO25"
clinical[clinical[,"class"]=="s50GGO","class"] <- "GGO50"
clinical[clinical[,"class"]=="s75GGO","class"] <- "GGO75"
clinical[clinical[,"class"]=="s100GGO","class"] <- "GGO100"

#Load the seurat object of all the cell
ggo.integrated <- readRDS(file=ggoFile)

#Calculate statistics and qunality control
ggo.integrated$"stage" <- clinical[match(as.character(ggo.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"]
clinical_sub <- clinical[clinical[,"rawDataID"]%in%unique(as.character(ggo.integrated@meta.data[,"orig.ident"])),]
patient_list <- split(clinical_sub[,"rawDataID"],f=factor(clinical_sub[,"class"]))
patient_list <- sapply(patient_list,length)
patient_list <- patient_list[c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")]
Ncell_list <- split(table(as.character(ggo.integrated@meta.data[,"orig.ident"])),
f=factor(clinical[match(names(table(as.character(ggo.integrated@meta.data[,"orig.ident"]))),as.character(clinical[,"rawDataID"])),"class"]))
Ncell_list <- Ncell_list[c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")]
Ncell_list <- lapply(Ncell_list,as.numeric)
Ngene <- split(ggo.integrated@meta.data[,"nFeature_RNA"],f=ggo.integrated@meta.data[,"orig.ident"])
NgeneM <- sapply(Ngene,median)
Ngene_list <- split(NgeneM,
f=factor(clinical[match(names(NgeneM),as.character(clinical[,"rawDataID"])),"class"]))
Ngene_list <- Ngene_list[c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")]
Nmit <- split(ggo.integrated@meta.data[,"percent.mt"],f=ggo.integrated@meta.data[,"orig.ident"])
Nmit <- sapply(Nmit,median)
Nmit_list <- split(Nmit,
f=factor(clinical[match(names(Nmit),as.character(clinical[,"rawDataID"])),"class"]))
Nmit_list <- Nmit_list[c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")]
Ncell_df <- data.frame(stage=factor(clinical[match(names(table(as.character(ggo.integrated@meta.data[,"orig.ident"]))),as.character(clinical[,"rawDataID"])),"class"],
levels=c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")),
fraction=as.numeric(table(as.character(ggo.integrated@meta.data[,"orig.ident"]))),stringsAsFactors=F)
Ngene_df <- data.frame(stage=factor(clinical[match(names(NgeneM),
as.character(clinical[,"rawDataID"])),"class"],
levels=c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")),
fraction=NgeneM,stringsAsFactors=F)
Nmit_df <- data.frame(stage=factor(clinical[match(names(Nmit),
as.character(clinical[,"rawDataID"])),"class"],
levels=c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")),
fraction=Nmit,stringsAsFactors=F)
p_df <- data.frame(Npt=factor(clinical_sub[,"class"],
levels=c("Normal", "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN","NonAden")))

#plot statistics and qunality control
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)],brewer.pal(6, "Greys")[4])
p1 <- ggplot(p_df,aes(Npt))+geom_bar(fill=colo)+coord_flip()+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.x=element_text(size=70),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("Number of patients")
p2<-ggboxplot(Ncell_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.x=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("Number of cells")+coord_flip()		  
p3<-ggboxplot(Ngene_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.x=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("Number of genes")+coord_flip()		  
p4<-ggboxplot(Ncell_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+
xlab(NULL)+theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.x=element_text(size=19),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("% of mitochondrial")+coord_flip()
grid.arrange(p1, p2,p3, p4,nrow=1)


##############################
#Figure S2C (R version 3.6.1)#
##############################
#${workDir} is the working directory
#${ggoFile} is the seurat object of all the file
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state
#${clinicalFile} is the file of clinical information
#${FindMarkerDir} is the result directory of findmarkers

#Load required packages
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(gplots)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load the seurat object of all cell
ggo.integrated <- readRDS(file=ggoFile)
#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
pathD <- cellStateDir

#For each cell type, plot heatmap for feature genes of cell state
for(x in 1:11)
{
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	
	#Load cell state information
	ecotyper_df=read.table(file =paste0(pathD,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	abundance <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/state_abundances.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	
	#extract annotated cells for the cell type
	label <- rep("unassigned",nrow(ggo.integrated@meta.data))
	names(label)=rownames(ggo.integrated@meta.data)
	label[ecotyper_df[,"ID"]]=ecotyper_df[,"State"]
	ggo.integrated$Ecotype = label
	mono.ecotyper.integrated <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[ggo.integrated@meta.data[,"Ecotype"]!="unassigned"]) 
	
	#exact expression matrix
	log2TPMf <-  as.matrix(mono.ecotyper.integrated@assays$RNA@data)
	sample_inter <- intersect(colnames(log2TPMf),colnames(abundance))
	cell_state <- sapply(sample_inter,function(x) rownames(abundance)[which(abundance[,x]==max(abundance[,x]))[1]])
	log2TPMf <- log2TPMf[,sample_inter]
	log2TPMf <- log2TPMf[apply(log2TPMf,1,sum)!=0,]
	
	#extract marker genes
	geneTopFC <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/gene_info.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	markers1_list <- lapply(seq(length(unique(geneTopFC[,"State"]))),function(x) { 
		res <- geneTopFC[(geneTopFC[,"State"]==paste0("S0",x))&(geneTopFC[,"MaxFC"]>0.5),"Gene"]
		return(res)
	})
	markers2_list <- lapply(seq(length(unique(geneTopFC[,"State"]))),function(x) { 
		res <- readRDS(file=paste0(FindMarkerDir,cell_type,".S0",
		x,".markerGenes.rds"))
		resT <- rownames(res)[(res[,"avg_logFC"]>log(1.5))&(res[,"p_val_adj"]<0.01)]
		return(resT)
	})
	markers3_list <- sapply(seq(length(unique(geneTopFC[,"State"]))),function(x){
		res <- unique(c(markers1_list[[x]],markers2_list[[x]]))
		return(res)
	})
	markers3_list <- lapply(markers3_list,function(x) intersect(x,rownames(log2TPMf)))
	gene_inter <- intersect(unlist(markers3_list),rownames(log2TPMf))
	cellType_log2TPM <- log2TPMf[gene_inter,]
	
	#order the sample of expression matrix by the top marker genes of each cell state
	order_list <- c()
	n_cluster <- 0
	for(j in seq(nrow(abundance))){
		mt_tmp <- cellType_log2TPM[,cell_state==rownames(abundance)[j]]
		gene <- intersect(markers3_list[[j]],rownames(log2TPMf))[1]
		order_tmp <- which(cell_state==rownames(abundance)[j])[order(unlist(mt_tmp[gene,]),decreasing=T)]
		order_list <- c(order_list,order_tmp)
		n_cluster <- c(n_cluster,n_cluster[length(n_cluster)]+sum(cell_state==rownames(abundance)[j]))
	}
	cellType_log2TPM_order <- cellType_log2TPM[,order_list]
	
	annotation_col = data.frame(
	cell_state = factor(cell_state[order_list]))
	rownames(annotation_col) = sample_inter[order_list]
	
	#Set the parameter for heatmap
	cell_state_colo <- brewer.pal(8, "Pastel1")
	names(cell_state_colo) <- paste0("S0",1:8)
	ann_colors = list(
         cell_state = cell_state_colo[unique(cell_state)])
	ngene <- c(0,cumsum(sapply(markers3_list,length)))
	cellType_log2TPM_order <- t(apply(cellType_log2TPM_order,1,function(x) (x-min(x))/(max(x)-min(x))))
	
	#plot heatmap
	pheatmap(cellType_log2TPM_order,border_color=NA, 
	color=colorRampPalette(c(brewer.pal(6, "Set1")[2],
	"white",brewer.pal(6, "Set1")[1]))(30),
	scale="none",cluster_rows=F,cluster_cols=F,
	annotation_col = annotation_col,
	gaps_col=n_cluster[-c(1,length(n_cluster))],
	gaps_row=ngene[-length(ngene)],
	annotation_colors = ann_colors)
}

##############################
#Figure S2D (R version 3.6.1)#
##############################
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

#For each LME, plot heatmap for feature genes of cell state
for(k in 1:6)
{
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
}

##############################
#Figure S3A (R version 3.6.1)#
##############################
#${workDir} is the working directory
#${EcoFile} is the results of Ecotyper
#${clinicalFile} is the file of clinical information of bulk RNaseq
#${ggoFile} is the seurat object of all the file
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#Load required packages
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the result of Ecotyper
fraction <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)
fraction <- fraction[,!(colnames(fraction)%in%c("X10C","X10F","X16C","X16F"))]

#Load the clinical information and basic processing
clinical_info <- read.table(file=clinicalFile,
sep="\t",stringsAsFactors=F,header=T)
clinical_info[,1] <- paste0("X",clinical_info[,1],"C")
clinical_info[clinical_info[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical_info[clinical_info[,"class"]=="s25GGO","class"] <- "GGO25"
clinical_info[clinical_info[,"class"]=="s50GGO","class"] <- "GGO50"
clinical_info[clinical_info[,"class"]=="s75GGO","class"] <- "GGO75"
clinical_label <- rep("Normal",ncol(fraction))
names(clinical_label) <- colnames(fraction)
clinical_label[clinical_info[,1]] <- clinical_info[,2]
n_clinical_label <- table(clinical_label)
clinical_labelN <- paste0(clinical_label,"(n=",n_clinical_label[clinical_label],")")

#calculate the fraction of LME for each clinical stage
fraction_list <- lapply(seq(nrow(fraction)),function(x){
	res <- split(unlist(fraction[x,]),
	f=factor(clinical_labelN,
	levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")))
	return(res)
})

#Set colors for boxplot
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(3,4,5,6,7)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")

#LME01
E1_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
E1=unlist(fraction[1,]),stringsAsFactors=F)
p<-ggboxplot(E1_df, x = "stage", y = "E1", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of LME1")
p

#LME03
E3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
E3=unlist(fraction[3,]),stringsAsFactors=F)
p<-ggboxplot(E3_df, x = "stage", y = "E3", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of LME3")
p

#LME05
E5_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
E5=unlist(fraction[5,]),stringsAsFactors=F)
p<-ggboxplot(E5_df, x = "stage", y = "E5", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of LME5")
p

#LME06
E6_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
E6=unlist(fraction[6,]),stringsAsFactors=F)
p<-ggboxplot(E6_df, x = "stage", y = "E6", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of LME6")
p

##############################
#Figure S3B (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${pathR} is the result directory of recovery cell state
#${pathD} is the result directory of discovery cell state
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${expressionFile} is the expression matrix of bulk RNAseq,scaled to TPM
#${clinicalFile} is the file of clinical information of bulk RNAseq

#Load required packages
library(RColorBrewer)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load the parameter for the appropriate number of cell states
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
#Load expression matrix of bulk RNAseq and log2 transformed
TPMMat <- read.table(file=expressionFile,sep="\t",stringsAsFactors=F,header=T,row.names=1)
log2TPM <- log2(TPMMat+1)
#Load clinical information
clinical_info <- read.table(file=clinicalFile,
sep="\t",stringsAsFactors=F,header=T) 

#For each cell type, plot heatmap for common feature genes of cell state
for(i in seq(nrow(rank_data)))
{
	cellType <- rank_data[i,1]
	print(cellType)
	ncluster <- rank_data[i,2]
	abundance <- read.table(file=paste0(pathR,cellType,"/",ncluster,"/state_abundances.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	sample_inter <- intersect(c(paste0("X",clinical_info[,1],"C"),paste0("X",clinical_info[,1],"F")),colnames(abundance))
	cell_state <- sapply(sample_inter,function(x) rownames(abundance)[which(abundance[,x]==max(abundance[,x]))[1]])
	clinical_label <- rep("Normal",length(sample_inter))
	names(clinical_label) <- sample_inter
	sample_interT <- grep("C$",sample_inter,value=T)
	clinical_label[sample_interT] <- clinical_info[match(sample_interT,
	paste0("X",clinical_info[,1],"C")),2]
	#Load the scaled expression data generated by Ecotyper
	log2TPMf <- read.table(file=paste0(pathR,cellType,"/",ncluster,"/expression_matrix_scaled.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	log2TPMf <- log2TPMf[,sample_inter]
	
	#select the common feature genes
	geneTopFC <- read.table(file=paste0(pathD,cellType,"/",ncluster,"/gene_info.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	matchingInitial <- read.table(file=paste0(pathD,cellType,"/",ncluster,"/mapping_to_initial_states.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	gene_top <- lapply(seq(nrow(matchingInitial)),function(x){
		geneFCtmp <- geneTopFC[geneTopFC[,"InitialState"]==matchingInitial[x,2],c("Gene","MaxFC")]
		gene_inter_tmp <- intersect(geneFCtmp[,"Gene"],rownames(log2TPMf))[1:50]
		meanstate <- apply(log2TPMf[gene_inter_tmp,cell_state==matchingInitial[x,1]],1,median)
		meanothers <- apply(log2TPMf[gene_inter_tmp,cell_state!=matchingInitial[x,1]],1,median)
		geneFCtmpOrder <- meanstate- meanothers
		names(geneFCtmpOrder) <- gene_inter_tmp
		return(names(geneFCtmpOrder)[geneFCtmpOrder>0.5])
	})
	gene_inter <- intersect(unlist(gene_top),rownames(log2TPMf))
	cellType_log2TPM <- log2TPMf[gene_inter,]
	
	#order the sample of expression matrix by top gene
	order_list <- c()
	n_cluster <- 0
	for(j in seq(nrow(matchingInitial))){
		mt_tmp <- cellType_log2TPM[,cell_state==matchingInitial[j,1]]
		gene <- intersect(gene_top[[j]],rownames(log2TPMf))[1]
		order_tmp <- which(cell_state==matchingInitial[j,1])[order(unlist(mt_tmp[gene,]),decreasing=T)]
		order_list <- c(order_list,order_tmp)
		n_cluster <- c(n_cluster,n_cluster[length(n_cluster)]+sum(cell_state==matchingInitial[j,1]))
	}
	cellType_log2TPM_order <- cellType_log2TPM[,order_list]
	annotation_col = data.frame(
	cell_state = factor(cell_state[order_list]),
	clinical_label=factor(clinical_label[order_list],levels=c("Normal","pGGO","HiDenGGO",
	"s25GGO","s50GGO","s75GGO","Solid1","Solid3","SolidN")))
	rownames(annotation_col) = sample_inter[order_list]
	
	#Ser the parameter for heatmap
	colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(3,4,5,6,7)],
	brewer.pal(6, "OrRd")[c(2,4,6)])
	names(colo) <- c("Normal","pGGO","HiDenGGO","s25GGO","s50GGO",
	"s75GGO","Solid1","Solid3","SolidN")
	cell_state_colo <- brewer.pal(8, "Pastel1")
	names(cell_state_colo) <- paste0("S0",1:8)
	ann_colors = list(
         cell_state = cell_state_colo[unique(cell_state)],
         clinical_label=colo)
	cellType_log2TPM_order[cellType_log2TPM_order>4] <- 4
	cellType_log2TPM_order[cellType_log2TPM_order < -4] = -4
	ngene <- cumsum(sapply(gene_top,length))
	
	#plot heatmap
	pheatmap(cellType_log2TPM_order,border_color=NA, 
	color=colorRampPalette(c(brewer.pal(6, "Set1")[2],
	"black",brewer.pal(6, "Set1")[6]))(50),
	scale="none",cluster_rows=F,cluster_cols=F,
	gaps_col=n_cluster[-c(1,length(n_cluster))],
	gaps_row=ngene[-length(ngene)],
	annotation_col = annotation_col,
	annotation_colors = ann_colors)
}
