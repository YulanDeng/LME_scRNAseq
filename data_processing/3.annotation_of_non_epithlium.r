## Script for the annotation of non-epithlium from scRNAseq
##Yulan Deng, last updated 2022-6-23

#####################################
#B cell annotation (R version 4.0.3)#
#####################################
#${workDir} is the working directory
#${pythonDir} is the directory of python (version 3.6.5)
#${BFile} is the seurat object of B cell
#${GSE131907_File} is raw count data of GSE131907 
#${GSE131907_phenotype_File} is phenotype file of GSE131907 

#Load required packages
library(reticulate)
use_python(pythonDir)
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of B cell
immu.integrated <- readRDS(file =BFile)
#Set the default assay for the seurat object of B cell
DefaultAssay(immu.integrated)<-"RNA"

#Load raw count data of GSE131907 
ref.data <- readRDS(file = GSE131907_File)
#Load phenotype file of GSE131907 
pheno <- read.table(file = GSE131907_phenotype_File,sep="\t",quote = "",stringsAsFactors=F,header=T)
#Basic process for raw count data of GSE131907  
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
c(grep("MT-",rownames(ref@assays[[1]]), value = T),
grep("^RP[SL]",rownames(ref@assays[[1]]),value = T)))) 
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", 
        nfeatures = 3000, verbose = FALSE)

#Integration B cell from our scRNAseq and  GSE131907
data_list <- list()
data_list[[1]] <- immu.integrated
data_list[[2]] <- ref
names(data_list) <- c(label,"ref")
sub.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:40) 
sub.integrated <- IntegrateData(anchorset = sub.anchors, dims = 1:40)
DefaultAssay(sub.integrated) <- "integrated"
sub.integrated <- ScaleData(sub.integrated, verbose = FALSE)
sub.integrated <- RunPCA(sub.integrated, verbose = FALSE)
sub.integrated <- RunUMAP(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- RunTSNE(sub.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
sub.integrated <- FindNeighbors(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- FindClusters(sub.integrated, verbose = FALSE)

#annotation the clusters of integrated data by the phenotype of GSE131907
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(rownames(immu.integrated@meta.data),rownames(sub.integrated@meta.data))
interlabel[intersmp] <- immu.integrated@meta.data[intersmp,"combAnno"]
sub.integrated$"combAnno" <- interlabel
interlabel2 <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel2) <- rownames(sub.integrated@meta.data)
intersmp2 <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel2[intersmp2] <- pheno[match(intersmp2,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel2
msi3 <- table(sub.integrated@meta.data[,c("seurat_clusters","combAnno")])
msi3 <- cbind(rownames(msi3),msi3)
colnames(msi3)[1] <- "clusters"
msi4 <- table(sub.integrated@meta.data[,c("seurat_clusters","combAnno")])
msi4 <- msi4[,setdiff(colnames(msi4),"NA")]
msiM4 <- apply(msi4,1,function(x) {
	if(sum(x)<10)
	{
		return("NA")
	}else{
		return(colnames(msi4)[which(x==max(x))[1]])
	}
})
msi4res <- cbind(rownames(msi4),msiM4)
colnames(msi4res)[1] <- "clusters"
msi5 <- table(sub.integrated@meta.data[,c("seurat_clusters","NCanno")])
msi5 <- msi5[,setdiff(colnames(msi5),"NA")]
msiM5 <- apply(msi5,1,function(x) {
	if(sum(x)<10)
	{
		return("NA")
	}else{
		return(colnames(msi5)[which(x==max(x))[1]])
	}
})
msi5res <- cbind(rownames(msi5),msiM5)
colnames(msi5res)[1] <- "clusters"

########################################
#CD4T cell annotation (R version 4.0.3)#
########################################
#${workDir} is the working directory
#${pythonDir} is the directory of python (version 3.6.5)
#${CD4T_File} is the seurat object of CD4+ T cell
#${GSE131907_File} is raw count data of GSE131907 
#${GSE131907_phenotype_File} is phenotype file of GSE131907 

#Load required packages
library(reticulate)
use_python(pythonDir)
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of CD4 T cell
immu.integrated <- readRDS(file =CD4T_File)
#Set the default assay for the seurat object of CD4 T cell
DefaultAssay(immu.integrated)<-"RNA"

#Load raw count data of GSE131907 
ref.data <- readRDS(file = GSE131907_File)
#Load phenotype file of GSE131907 
pheno <- read.table(file = GSE131907_phenotype_File,sep="\t",quote = "",stringsAsFactors=F,header=T)
#Basic process for raw count data of GSE131907 
smp <- pheno[pheno[,"Cell_subtype"]%in%c("Naive CD4+ T","Treg","CD4+ Th"),1]
ref.data <- ref.data[,colnames(ref.data)%in%smp]
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
c(grep("MT-",rownames(ref@assays[[1]]), value = T),
grep("^RP[SL]",rownames(ref@assays[[1]]),value = T)))) 
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", 
        nfeatures = 3000, verbose = FALSE)

#Integration CD4T cell from our scRNAseq and  GSE131907
data_list <- list()
data_list[[1]] <- immu.integrated
data_list[[2]] <- ref
names(data_list) <- c(label,"ref")
sub.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:40) 
sub.integrated <- IntegrateData(anchorset = sub.anchors, dims = 1:40)
DefaultAssay(sub.integrated) <- "integrated"
sub.integrated <- ScaleData(sub.integrated, verbose = FALSE)
sub.integrated <- RunPCA(sub.integrated, verbose = FALSE)
sub.integrated <- RunUMAP(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- RunTSNE(sub.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
sub.integrated <- FindNeighbors(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- FindClusters(sub.integrated, verbose = FALSE) 

#annotation the clusters of integrated data by the phenotype of GSE131907
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
msi3 <- table(sub.integrated@meta.data[,c("seurat_clusters","RNA_snn_res.0.2")])
msi3 <- cbind(rownames(msi3),msi3)
colnames(msi3)[1] <- "clusters"
msi4 <- table(sub.integrated@meta.data[,c("seurat_clusters","RNA_snn_res.0.2")])
msi4 <- msi4[,setdiff(colnames(msi4),"NA")]
msiM4 <- apply(msi4,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi4)[which(x==max(x))[1]])
}
})
msi4res <- cbind(rownames(msi4),msiM4)
colnames(msi4res)[1] <- "clusters"
msi5 <- table(sub.integrated@meta.data[,c("seurat_clusters","NCanno")])
msi5 <- msi5[,setdiff(colnames(msi5),"NA")]
msiM5 <- apply(msi5,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi5)[which(x==max(x))[1]])
}
})
msi5res <- cbind(rownames(msi5),msiM5)
colnames(msi5res)[1] <- "clusters"
labelF <- as.character(immu.integrated@meta.data[,"RNA_snn_res.0.2"])
labelF[labelF%in%c("0","3","4","6")] <- "CD4_Th"
labelF[labelF=="1"] <- "Treg"
labelF[labelF=="2"] <- "Naive_CD4T"
labelF[labelF=="5"] <- "CXCL13_CD4T"
immu.integrated$"assigned_cell_type_byNC" <- labelF

########################################
#CD8T cell annotation (R version 4.0.3)#
########################################
#${workDir} is the working directory
#${pythonDir} is the directory of python (version 3.6.5)
#${CD8T_File} is the seurat object of CD8+ T cell
#${GSE131907_File} is raw count data of GSE131907 
#${GSE131907_phenotype_File} is phenotype file of GSE131907 

#Load required packages
library(reticulate)
use_python(pythonDir)
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of CD8 T cell
immu.integrated <- readRDS(file =CD8T_File)
#Set the default assay for the seurat object of CD8 T cell
DefaultAssay(immu.integrated)<-"RNA"

#Load raw count data of GSE131907 
ref.data <- readRDS(file = GSE131907_File)
#Load phenotype file of GSE131907 
pheno <- read.table(file = GSE131907_phenotype_File,sep="\t",quote = "",stringsAsFactors=F,header=T)
#Basic process for raw count data of GSE131907
smp <- pheno[pheno[,"Cell_subtype"]%in%c("Cytotoxic CD8+ T","Exhausted CD8+ T","Naive CD8+ T","CD8 low T"),1]
ref.data <- ref.data[,colnames(ref.data)%in%smp]
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
c(grep("MT-",rownames(ref@assays[[1]]), value = T),
grep("^RP[SL]",rownames(ref@assays[[1]]),value = T)))) 
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", 
        nfeatures = 3000, verbose = FALSE)

#Integration CD8T cell from our scRNAseq and  GSE131907
data_list <- list()
data_list[[1]] <- immu.integrated
data_list[[2]] <- ref
names(data_list) <- c(label,"ref")
sub.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:40) 
sub.integrated <- IntegrateData(anchorset = sub.anchors, dims = 1:40)
DefaultAssay(sub.integrated) <- "integrated"
sub.integrated <- ScaleData(sub.integrated, verbose = FALSE)
sub.integrated <- RunPCA(sub.integrated, verbose = FALSE)
sub.integrated <- RunUMAP(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- RunTSNE(sub.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
sub.integrated <- FindNeighbors(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- FindClusters(sub.integrated, verbose = FALSE) 

#annotation the clusters of integrated data by the phenotype of GSE131907
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
msi3 <- table(sub.integrated@meta.data[,c("seurat_clusters","RNA_snn_res.0.1")])
msi3 <- cbind(rownames(msi3),msi3)
colnames(msi3)[1] <- "clusters"
msi4 <- table(sub.integrated@meta.data[,c("seurat_clusters","RNA_snn_res.0.1")])
msi4 <- msi4[,setdiff(colnames(msi4),"NA")]
msiM4 <- apply(msi4,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi4)[which(x==max(x))[1]])
}
})
msi4res <- cbind(rownames(msi4),msiM4)
colnames(msi4res)[1] <- "clusters"
msi5 <- table(sub.integrated@meta.data[,c("seurat_clusters","NCanno")])
msi5 <- msi5[,setdiff(colnames(msi5),"NA")]
msiM5 <- apply(msi5,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi5)[which(x==max(x))[1]])
}
})
msi5res <- cbind(rownames(msi5),msiM5)
colnames(msi5res)[1] <- "clusters"
labelF <- as.character(immu.integrated@meta.data[,"RNA_snn_res.0.1"])
labelF[labelF%in%c("0","2","4")] <- "pre-effector_CD8"
#pre-effector_CD8 is also called transitory CD8 T cell in the manuscript
labelF[labelF%in%c("1","6","7","8")] <- "effector_CD8"
labelF[labelF%in%c("3","5")] <- "exhausted_CD8"
immu.integrated$"assigned_cell_type_byNC" <- labelF

######################################################
#monocyte and macrophage annotation (R version 4.0.3)#
######################################################
#${workDir} is the working directory
#${pythonDir} is the directory of python (version 3.6.5)
#${moMac_File} is the seurat object of monocyte and macrophage
#${GSE131907_File} is raw count data of GSE131907 
#${GSE131907_phenotype_File} is phenotype file of GSE131907 

#Load required packages
library(reticulate)
use_python(pythonDir)
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of  monocyte and macrophage cell
immu.integrated <- readRDS(file =moMac_File)
#Set the default assay for the seurat object of monocyte and macrophage
DefaultAssay(immu.integrated)<-"RNA"

#Load raw count data of GSE131907 
ref.data <- readRDS(file = GSE131907_File)
#Load phenotype file of GSE131907 
pheno <- read.table(file = GSE131907_phenotype_File,sep="\t",quote = "",stringsAsFactors=F,header=T)
#Basic process for raw count data of GSE131907
smp <- pheno[pheno[,"Cell_subtype"]%in%c("mo-Mac","Monocytes","Alveolar Mac",
"Pleural Mac","CD1c+ DCs","CD141+ DCs","pDCs","CD163+CD14+ DCs","Activated DCs"),1]
ref.data <- ref.data[,colnames(ref.data)%in%smp]
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
c(grep("MT-",rownames(ref@assays[[1]]), value = T),
grep("^RP[SL]",rownames(ref@assays[[1]]),value = T)))) 
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", 
        nfeatures = 3000, verbose = FALSE)

#Integration CD8T cell from our scRNAseq and  GSE131907
data_list <- list()
data_list[[1]] <- immu.integrated
data_list[[2]] <- ref
names(data_list) <- c(label,"ref")
sub.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:40) 
sub.integrated <- IntegrateData(anchorset = sub.anchors, dims = 1:40)
DefaultAssay(sub.integrated) <- "integrated"
sub.integrated <- ScaleData(sub.integrated, verbose = FALSE)
sub.integrated <- RunPCA(sub.integrated, verbose = FALSE)
sub.integrated <- RunUMAP(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- RunTSNE(sub.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
sub.integrated <- FindNeighbors(sub.integrated, dims = 1:40, verbose = FALSE)
sub.integrated <- FindClusters(sub.integrated, verbose = FALSE)

#annotation the clusters of integrated data by the phenotype of GSE131907
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
interlabel3 <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel3) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(rownames(immu.integrated@meta.data),rownames(sub.integrated@meta.data))
interlabel3[intersmp] <- immu.integrated@meta.data[intersmp,"combAnno"]
sub.integrated$"combAnno" <- interlabel3
interlabel <- rep("NA",nrow(sub.integrated@meta.data))
names(interlabel) <- rownames(sub.integrated@meta.data)
intersmp <- intersect(pheno[,1],rownames(sub.integrated@meta.data))
interlabel[intersmp] <- pheno[match(intersmp,pheno[,1]),"Cell_subtype"]
sub.integrated$"NCanno" <- interlabel
msi3 <- table(sub.integrated@meta.data[,c("seurat_clusters","combAnno")])
msi3 <- cbind(rownames(msi3),msi3)
colnames(msi3)[1] <- "clusters"
msi4 <- table(inte@meta.data[,c("seurat_clusters","combAnno")])
msi4 <- msi4[,setdiff(colnames(msi4),"NA")]
msiM4 <- apply(msi4,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi4)[which(x==max(x))[1]])
}
})
msi4res <- cbind(rownames(msi4),msiM4)
colnames(msi4res)[1] <- "clusters"
msi5 <- table(inte@meta.data[,c("seurat_clusters","NCanno")])
msi5 <- msi5[,setdiff(colnames(msi5),"NA")]
msiM5 <- apply(msi5,1,function(x) {
if(sum(x)<10)
{
return("NA")
}else{
return(colnames(msi5)[which(x==max(x))[1]])
}
})
msi5res <- cbind(rownames(msi5),msiM5)
colnames(msi5res)[1] <- "clusters"
labelF <- as.character(immu.integrated@meta.data[,"combAnno"])
labelF[labelF%in%c("1","2","5","6","8")] <- "tissue_resident_Macrophage"
labelF[labelF%in%c("3","9")] <- "doublets"
labelF[labelF%in%c("7","c0")] <- "TAM_anti_inflammatory"
labelF[labelF%in%c("c1")] <- "cDC"
labelF[labelF%in%c("c2")] <- "early_stage_macrophages"
labelF[labelF%in%c("s0")] <- "CD14_monocyte"
labelF[labelF%in%c("s1")] <- "CD16_monocyte"
#TAM_anti_inflammatory is also called SPP1+ macrophage in the manuscript
immu.integrated$"assigned_cell_type_byNC" <- labelF

##########################################
#endothelium annotation (R version 3.6.1)#
##########################################
#${workDir} is the working directory
#${allCellFile} is the seurat object of all cells

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of  all cells
ggo.integrated <- readRDS(file=allCellFile)
#select endothelium
endo.integrated <- subset(ggo.integrated,idents="Endothelial")
#Basic process for the object of emdothelium
DefaultAssay(endo.integrated)<-"RNA"
endo.integrated <- SCTransform(endo.integrated,return.only.var.genes = F,verbose = T,conserve.memory = T)
endo.integrated <- RunPCA(endo.integrated, features = VariableFeatures(object = endo.integrated))
endo.integrated <- RunUMAP(endo.integrated, dims = 1:30, verbose = FALSE)
endo.integrated <- RunTSNE(endo.integrated, dims = 1:30, nthreads = 4, max_iter = 2000)
endo.integrated <- FindNeighbors(endo.integrated, dims = 1:30, verbose = FALSE)
endo.integrated <- FindClusters(endo.integrated, resolution=seq(0.1,1,by=0.1), verbose = FALSE,algorithm=1)

#plot marker genes of endothelial subtype
markers.to.plot <- c( "INSR", "HSPG2",  "VWA1", "RGCC",
"RAMP3",   "ACKR1","SELP", "TYROBP", "C1QB", "LYVE1","PROX1",
"CXCR4",  "GJA5","CA4", "IL33","CPE","POSTN", "CCL14","PRCP","ENG","PGF","LXN",
"IGFBP3","SPARC","IGFBP7","EDN1","EPAS1","EGR1","HPGD","IL1RL1","EDNRB",
"FCN3","BTNL9","NOSTRIN","HLA-DRA","C1QA","MARCO","SFTPC","CCL21","LYVE1","PDPN")
DefaultAssay(endo.integrated)<-"RNA"
FeaturePlot(endo.integrated, features = markers.to.plot,
reduction="umap",ncol=3) + RotatedAxis()

#annotation of endothelium
labelF <- as.character(endo.integrated@meta.data[,"SCT_snn_res.0.5"])
labelF[labelF%in%c("0","4","12")] <- "tumor_EC"
labelF[labelF%in%c("6")] <- "tip_like_EC"
labelF[labelF%in%c("2","7")] <- "alveolar_type_I"
labelF[labelF%in%c("3","9","10")] <- "stalk_like_EC"
labelF[labelF=="5"] <- "lymphatic_EC"
labelF[labelF%in%c("1","8")] <- "arteries"
labelF[labelF=="11"] <- "scavenging"
endo.integrated$"assigned_cell_type" <- labelF

#########################################
#fibroblast annotation (R version 3.6.1)#
#########################################
#${workDir} is the working directory
#${allCellFile} is the seurat object of all cells

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load seurat object of  all cells
ggo.integrated <- readRDS(file=allCellFile)
#select fibroblast
fibro.integrated <- subset(ggo.integrated,idents="Fibroblasts")
#Basic process for the object of fibroblast
DefaultAssay(fibro.integrated)<-"RNA"
fibro.integrated <- SCTransform(fibro.integrated,return.only.var.genes = F,verbose = T,conserve.memory = T)
fibro.integrated <- RunPCA(fibro.integrated, features = VariableFeatures(object = fibro.integrated))
fibro.integrated <- RunUMAP(fibro.integrated, dims = 1:30, verbose = FALSE)
fibro.integrated <- RunTSNE(fibro.integrated, dims = 1:30, nthreads = 4, max_iter = 2000)
fibro.integrated <- FindNeighbors(fibro.integrated, dims = 1:30, verbose = FALSE)
fibro.integrated <- FindClusters(fibro.integrated, resolution=seq(0.1,1,by=0.1), verbose = FALSE,algorithm=1)

#plot marker genes offibroblast subtype
markers.to.plot <- c( "COL14A1","GSN","PI16","CYGB","PRRX1",
"COL13A1","TCF21","ITGA8","CXCL14","NPNT",
"TAGLN","ACTA2","ACTG2","MYH11","MYLK",
"SYNPO2","CRYAB","CNN1","DES",
"UPK3B","MSLN","CALB2","WT1",
"CYP1B1","APOD",
"RGS5","CSPG4","ABCC9","KCNJ8",
"COL10A1","COL4A1","PLA2G2A","MMP3","FIGF","CCL2")

#annotation of fibroblast
labelF <- as.character(fibro.integrated@meta.data[,"SCT_snn_res.0.5"])
labelF[labelF%in%c("4")] <- "COL14A1+_matrix_FBs"
labelF[labelF%in%c("0","3","6","10")] <- "COL13A1+_matrix_FBs"
labelF[labelF%in%c("8","9","11","13")] <- "Myofibroblasts"
labelF[labelF%in%c("1")] <- "Pericytes"
labelF[labelF%in%c("2","5")] <- "COL10A1+_FBs"
labelF[labelF%in%c("12")] <- "PLA2G2A_FBs"
labelF[labelF=="7"] <- "mixture"
fibro.integrated$"assigned_cell_type" <- labelF
