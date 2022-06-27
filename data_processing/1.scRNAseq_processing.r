## Script for the basis processing of scRNAseq
##Yulan Deng, last updated 2022-6-23

###################
#cellranger (bash)#
###################
#${smp} is sample ID
#${fastqsDir} is the file directory for fastq files of scRNAseq dataset
#${smp_list} is the file name for fastq files
#${transcriptomeDir} is the directory for reference genome(hg38)
mkdir smp
cd smp
cellranger count --id=run_count_smp \
--fastqs=fastqsDir \
--localmem 2048 --localcores 80  --sample smp_list \
--include-introns --transcriptome=transcriptomeDir 

##########################
#seurat (R version 3.6.1)#
##########################
#${workDir} is the working directory
#${clinicalFile} is the path for clinical data,where the first column represents clinical stage,and
##the forth column represents sample ID.

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)

#Set working directory
setwd(workDir)

#Load clinical data and results from cellranger
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
fl_list <- sample_list
data_list <- lapply(seq(length(sample_list)),function(x){ 
	pbmc.data <- Read10X(data.dir = paste0("/NAS/dyl/project/singlecell/GGO/1.cellranger/2.cellranger_count/",
	fl_list[x],"/run_count_",sample_list[x],"/outs/filtered_feature_bc_matrix"))
	pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample_list[x], min.cells = 3, min.features = 200)
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
	pbmc <-subset(pbmc,features=setdiff(rownames(pbmc@assays[[1]]),
	c(grep("MT-",rownames(pbmc@assays[[1]]), value = T),
	grep("^RP[SL]",rownames(pbmc@assays[[1]]),value = T)))) 
	pbmc <- NormalizeData(pbmc, verbose = FALSE)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)
	return(pbmc)
})
names(data_list) <- sample_list

#scRNAseq integration
features <- SelectIntegrationFeatures(object.list = data_list)
data.reduce.list <- lapply(X = data_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
t1 <- which(!duplicated(clinical[,1]))
#one sample of each clinical stage was selected for refrence of anchors
anchors <- FindIntegrationAnchors(object.list = data.reduce.list, reference = t1), reduction = "rpca", 
dims = 1:50)
ggo.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

#clustering of integrated datasets
DefaultAssay(ggo.integrated) <- "integrated"
ggo.integrated <- ScaleData(ggo.integrated, verbose = FALSE)
ggo.integrated <- RunPCA(ggo.integrated,npcs=50,verbose = FALSE)
print(ElbowPlot(ggo.integrated,ndims=50))
#dim of PCA was chosen bt Elbow Plot
ggo.integrated <- RunUMAP(ggo.integrated, dims = 1:40, verbose = FALSE)
ggo.integrated <- RunTSNE(ggo.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
ggo.integrated <- FindNeighbors(ggo.integrated, dims = 1:40, verbose = FALSE)
ggo.integrated <- FindClusters(ggo.integrated, verbose = FALSE) 
saveRDS(ggo.integrated, file = "ggo.integrated.rds")

#################################################
#annotation of major cell type (R version 3.6.1)#
#################################################
#${annotationFile} is the annotation result for clustering, where the third column represents cell annotation.
markers.to.plot <- c("EPCAM","KRT19","KRT18","CDH1","CD3D","CD3E","CD3G","TRAC","LYZ","MARCO",
"CD68","FCGR3A","CD79A","DCN","THY1","COL1A1","COL1A2","PECAM1","FLT1",
"KIT","MS4A2","GATA2","NKG7","NCAM1","KLRD1")
DotPlot(ggo.integrated, features = markers.to.plot,  dot.scale = 8,assay="RNA") + 
    RotatedAxis()

#Load file for cell annotation
cluster <- read.table(file = annotationFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
new.cluster.ids <- cluster[,3]
names(new.cluster.ids) <- levels(ggo.integrated)
ggo.integrated <- RenameIdents(ggo.integrated, new.cluster.ids)
ggo.integrated$assigned_cell_type<-Idents(ggo.integrated)
saveRDS(ggo.integrated,file="ggo.integrated.anno.rds")

###############################################
#reclustering of immune cell (R version 3.6.1)#
###############################################
#${workDir} is the working directory
#${integratedFile} is the file for annotation of major cell type
#${immuneClusterFile} is the file for annotation of immune cell,
##where the second column respresents cell annotation.

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load results from annotation of major cell type
ggo.integrated <- readRDS(file=integratedFile)

#extract immune cells from annotation results
immu.integrated <- subset(ggo.integrated,idents=c("T","Myeloid","NK","B","MAST")) 

#clustering processing of immune cells
immu.integrated <- RunPCA(immu.integrated, features = VariableFeatures(object = immu.integrated))
print(ElbowPlot(immu.integrated,ndims=50))
#dim of PCA was chosen bt Elbow Plot
immu.integrated <- RunUMAP(immu.integrated, dims = 1:40, verbose = FALSE)
immu.integrated <- RunTSNE(immu.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
immu.integrated <- FindNeighbors(immu.integrated, dims = 1:40, verbose = FALSE)
immu.integrated <- FindClusters(immu.integrated, verbose = FALSE,resolution=0.5) 

#annotation of clustering results
markers.to.plot <- c( "CD3D","CD3E","CD3G","TRAC","CD4","COTL1", "LDHB","CCR7","LEF1","FOXP3",
"CD8A","GZMK", "DUSP2","EOMES","GZMH", "GZMB","SLC4A10","ZBTB16","TRDC","KLRB1","KLRG1","TRGC1","TRGC2","NKG7",
"NCAM1","TYROBP","KLRD1","FCER1G",
"RORA","ICOS","IL7R","NCR1","IL23R",
"CD79A","CD24", "MS4A1", "CD19","CD27", "SLAMF7","JCHAIN","MZB1",
"MARCO","MSR1", "MRC1","CD68","CD14","S100A8","FCGR3A","LILRB4", "IRF8", "LILRA4","CLEC9A", "LAMP3","CD1C",
"FCER1A","PLD4",
"S100A9","FCGR3B","IFITM2","NAMPT","KIT","MS4A2","GATA2","TPSAB1","TUBB3","TH","CALCA","PTPRC")
DotPlot(immu.integrated, features = markers.to.plot,  dot.scale = 8,assay="RNA") + 
    RotatedAxis()

#Load file for the annotation of immune cell
cluster <- read.table(file = immuneClusterFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
new.cluster.ids <- cluster[,2]
names(new.cluster.ids) <- levels(immu.integrated)
immu.integrated<- RenameIdents(immu.integrated, new.cluster.ids)
immu.integrated$assigned_cell_type<-Idents(immu.integrated)
saveRDS(immu.integrated,file="immu.integrated.anno.rds")

##############################################
#reclustering of lymphocyte (R version 3.6.1)#
##############################################
#${lymClusterFile} is the file for annotation of lymphocyte,
##where the forth column respresents cell annotation.
lymph <- subset(immu.integrated,cells=rownames(immu.integrated@meta.data)[as.character(immu.integrated@meta.data[,"major"])=="L"])
lymph <- ScaleData(lymph,  verbose = FALSE)
lymph <- FindVariableFeatures(lymph, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)		
lymph <- RunPCA(lymph,npcs=50,verbose = FALSE, features = VariableFeatures(object = lymph))
print(ElbowPlot(lymph,ndims=50))
#dim of PCA was chosen bt Elbow Plot
lymph.integrated <- RunUMAP(lymph, dims = 1:40, nthreads = 32,verbose = FALSE)
lymph.integrated <- RunTSNE(lymph.integrated, dims = 1:40, nthreads = 32, max_iter = 2000)
lymph.integrated  <- FindNeighbors(lymph.integrated , dims = 1:40,nthreads = 32, verbose = FALSE)
lymph.integrated  <- FindClusters(lymph.integrated , nthreads = 32,verbose = FALSE,resolution=0.2) 

#Load file for the annotation of lymphocyte
cluster <- read.table(file = lymClusterFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
lym.integrated$major<- cluster[match(as.character(lym.integrated@meta.data[,"seurat_clusters"]),cluster[,1]),4]

#Plot canonical markers
markers.to.plot <- c( "CD79A", "CD24", "MS4A1", "CD19","CD27", "SLAMF7","CD3E", "CD8A","CD4",
"KLRD1", "NKG7", "TYROBP","FCER1G")
DotPlot(lym.integrated, features = markers.to.plot,  group.by="major",dot.scale = 8,assay="RNA") + 
    RotatedAxis()
saveRDS(lym.integrated , file = "lym.integrated.rds")

##########################################
#reclustering of myloid (R version 3.6.1)#
##########################################
#${myeClusterFile} is the file for annotation of lymphocyte,
##where the forth column respresents cell annotation.
mye <- subset(immu.integrated,cells=rownames(immu.integrated@meta.data)[as.character(immu.integrated@meta.data[,"major"])=="M"])
mye <- ScaleData(mye,  verbose = FALSE)		
mye <- RunPCA(mye,npcs=50,verbose = FALSE, features = VariableFeatures(object = mye))
print(ElbowPlot(mye,ndims=50))
#dim of PCA was chosen bt Elbow Plot
mye.integrated <- RunUMAP(mye, dims = 1:40, verbose = FALSE)
mye.integrated <- RunTSNE(mye.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
mye.integrated  <- FindNeighbors(mye.integrated , dims = 1:40, verbose = FALSE)
mye.integrated  <- FindClusters(mye.integrated , verbose = FALSE,resolution=0.3) 

#Load file for myloid
cluster <- read.table(file = myeClusterFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
mye.integrated$major<- cluster[match(as.character(mye.integrated@meta.data[,"seurat_clusters"]),cluster[,1]),4]

#Plot canonical markers
markers.to.plot <- c( "S100A8", "S100A9", "IFITM2", "FCGR3B",
"MS4A2", "CPA3", "TPSAB1","SIGLEC8",
"NRGN", "PPBP", "PF4", "OST4",
"MARCO", "MSR1", "MRC1",
"LILRB4", "IRF8", "LILRA4","CLEC9A", "LAMP3",
 "CD1C", "PLD4","CD14","FCGR3A")
DotPlot(mye.integrated, features = markers.to.plot, group.by="major", dot.scale = 8,assay="RNA") + 
    RotatedAxis()
saveRDS(mye.integrated , file = "mye.integrated.rds")
