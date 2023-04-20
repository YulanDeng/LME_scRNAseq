## Script for the annotation of epithlium from scRNAseq
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

###########################
#copyKat (R version 4.0.3)#
###########################
#${smp} is sample ID
#${annoFile} is the data from cell annotation
#${workDir} is the working directory

#Load required packages
library(copykat)
library(Seurat)

#Set working directory
setwd(workDir)
dir.create(smp)
setwd(paste0(workDir,smp))

#Load annotation file 
ggo.integrated <- readRDS(file = annoFile)

#Run copyKAT for each sample
pbmc <- subset(ggo.integrated, cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"orig.ident"])==smp])
exp.rawdata <- as.matrix(pbmc@assays$RNA@counts)
copykat.ZG <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=smp, distance="euclidean", 
norm.cell.names=colnames(exp.rawdata)[as.character(pbmc@meta.data[,"assigned_cell_type"])!="Epithelial"],n.cores=32)
saveRDS(copykat.ZG, file = paste0(smp,".copyKAT.rds"))

############################################################
#extract aneuploid epithlium from copyKAT (R version 3.6.1)#
############################################################
#${workDir} is the working directory
#${clinicalFile} is the path for clinical data,where the forth column represents sample ID.
#${integratedFile} is seurat object with annotation of major cell type
#${copykatDir} is the path of copyKAT result

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)

#Set working directory
setwd(workDir)

#Load clinical data and annotation file
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
sample_list <- clinical[,4]
ggo.integrated <- readRDS(file=integratedFile)

#Extract epithlium 
epi.integrated <- subset(ggo.integrated,idents="Epithelial")
epi.integrated <- RunPCA(epi.integrated, features = VariableFeatures(object = epi.integrated))

#Extract result from copyKAT
anno_label_list <- c()
for(smp in sample_list) 
{
	pbmc <- readRDS(paste0(copykatDir,smp,"/",smp,".copyKAT.rds"))
	tumor.cells <- pbmc["prediction"][[1]][which(pbmc["prediction"][[1]][,"copykat.pred"]=="aneuploid"),"cell.names"]
	tumor.cells <- substr(tumor.cells,1,16)
	##remove suffix from barcode ID
	anno_label <- rownames(epi.integrated@meta.data)[(as.character(epi.integrated@meta.data[,"orig.ident"])==smp)]
	names(anno_label) <- substr(anno_label,1,16)
	##remove suffix from barcode ID
	intersmp <- anno_label[names(anno_label)%in%tumor.cells]
	names(intersmp) <- NULL
	anno_label_list <- c(anno_label_list,intersmp)
}
label <- rep("normal",nrow(epi.integrated@meta.data))
names(label) <- rownames(epi.integrated@meta.data)
label[names(label)%in%anno_label_list] <- "tumor"
epi.integrated$cnv <- label
cancer.integrated <- subset(epi.integrated,cells=rownames(epi.integrated@meta.data)[as.character(epi.integrated@meta.data[,"cnv"])=="tumor"])

#basic processing of epithelium
cancer.integrated <- RunUMAP(cancer.integrated, dims = 1:30, verbose = FALSE)
cancer.integrated <- RunTSNE(cancer.integrated, dims = 1:30, nthreads = 4, max_iter = 2000)
cancer.integrated <- FindNeighbors(cancer.integrated, dims = 1:30, verbose = FALSE)
saveRDS(cancer.integrated, file = "cancer.integrated.rds")

############################
#inferCNV (R version 4.0.3)#
############################
#${smp} is sample ID
#${annoFile} is the data from all cell annotation
#${epiAnnoFile} is the data from epithelium cell annotation
#${workDir} is the working directory

#Load required packages and set working directory
setwd(workDir)
dir.create(smp)
setwd(paste0(workDir,smp))
library(Seurat)
library(dplyr)
library(patchwork)
library(infercnv)
library(RColorBrewer)

#Load annotation file
ggo.integrated <- readRDS(file = annoFile)
cancer.integrated <- readRDS(file = epiAnnoFile)
ggo.RNA <- subset(ggo.integrated, cells=rownames(ggo.integrated@meta.data)[as.character(ggo.integrated@meta.data[,"orig.ident"])==smp])
cellanno <- data.frame(cell=rownames(ggo.RNA@meta.data),
cluster="others",stringsAsFactors=F)
inter_normal <- rownames(ggo.RNA@meta.data)[as.character(ggo.RNA@meta.data[,"assigned_cell_type"])!="Epithelial"]
cellanno[cellanno[,1]%in%inter_normal,2] <- "normal"
inter_epi <- intersect(rownames(ggo.RNA@meta.data)[as.character(ggo.RNA@meta.data[,"assigned_cell_type"])=="Epithelial"], rownames(cancer.integrated@meta.data))
cellanno[match(inter_epi,cellanno[,1]),2] <- as.character(cancer.integrated@meta.data[inter_epi,"integrated_snn_res.0.8"])

#Exclude cell clusters with less than 10 cells
ncluster <- table(cellanno[,2])
cluster_sel <- names(ncluster)[ncluster>=10]
cellanno <- cellanno[cellanno[,2]%in%cluster_sel,]

#Generate input data for inferCNV
counts <- ggo.RNA@"assays"[[1]]@"counts"
counts <- counts[,cellanno[,1]]
counts <- cbind(rownames(counts),as.matrix(counts))
colnames(counts)[1] <- ""
write.table(counts, file=paste0(smp,".singleCell.counts.matrix"),sep="\t" , quote=FALSE, row.names=F, col.names=T)
write.table(cellanno, file=paste0(smp,".cellAnnotations.txt"),
 sep="\t" , quote=FALSE, row.names=F, col.names=F )

#Create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(smp,
".singleCell.counts.matrix"),
annotations_file=paste0(smp,".cellAnnotations.txt"),delim="\t",
gene_order_file="/NAS/dyl/project/singlecell/GGO/10.infercnv/gene_ordering_fileM.txt",
ref_group_names=c("normal"))

#Perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1,
out_dir=paste0(smp,".infercnv"),  
cluster_by_groups=T,denoise=T,HMM=T,output_format = "pdf")
