## Script for the annotation of epithlium from scRNAseq
##Yulan Deng, last updated 2022-6-23

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

##################################################################################
#determine the best resolution of epithelium clusterings by sihouette coefficient#
#and balance of number of cells within each cluster (R version 3.6.1)            #
##################################################################################
#${workDir} is the working directory
#${pythonDir} is the directory of python (version 3.6.5)
#${cancerFile} is seurat object with cancer cell from copyKAT

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(cluster)
library("reticulate")
use_python(pythonDir)

#Set working directory
setwd(workDir)

#Load seurat object for cancer cell
cancer.integrated <- readRDS(file =cancerFile)

#function of running sihouette coefficient
sil_subsample_v3=function(mat,clust){
  out=as.data.frame(matrix(rep(NA,20*ncol(clust)),nrow=20))
  for(x in 1:20){
    for(j in 1:ncol(clust)){
      i=c()
      a=clust[,j]
      for(lab in unique(a)){
        i=c(i,sample((1:ncol(mat))[which(a==lab)],length(which(a==lab))/length(a)*min(1000,ncol(mat))))
      }
      d=as.dist(1 - cor(mat[,i], method = "pearson"))
      if(length(table(clust[i,j]))==1){out[x,j]=0}
      else{
        sil=silhouette(as.numeric(clust[i,j]),d)
        out[x,j]=mean(sil[, "sil_width"])}}
  }
  means=apply(out,2,median)
  sd=apply(out,2,mad)
  return(list(means,sd))
}

#enumerate the resolution for findclusters, from 0.05 to 1,step by 0.05.
cancer.integrated <- FindClusters(cancer.integrated, resolution=seq(0.05,1,by=0.05), verbose = T,algorithm=1) 

#calculate sihouette coefficient
clust=cancer.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(cancer.integrated@meta.data)))]
mat=as.data.frame(t(cancer.integrated$pca@cell.embeddings))
out=sil_subsample_v3(mat,clust)

#plot result for sihouette coefficient
means=out[[1]]
std=out[[2]]
x=seq(0.05,1,by=0.05)
barplot(means,ylab="mean silhouette score",xlab="resolution parameter")

#plot the number of cells within each cluster at different candidate resolution, that is, 0.05,0.1,0.2,0.25
barplot(table(cancer.integrated@meta.data[,"integrated_snn_res.0.05"]))
barplot(table(cancer.integrated@meta.data[,"integrated_snn_res.0.1"]))
barplot(table(cancer.integrated@meta.data[,"integrated_snn_res.0.2"]))
barplot(table(cancer.integrated@meta.data[,"integrated_snn_res.0.25"]))

#Finally, resolution=0.25 is selected
#find marker genes of each cluster
cancer.integrated <- FindClusters(cancer.integrated, resolution=0.25, verbose = T,algorithm=1) 
DefaultAssay(cancer.integrated)<-"RNA"
cancer.markers<-FindAllMarkers(cancer.integrated,min.pct = 0.25, logfc.threshold = 0.25)
cancer.markers[,6] <- as.character(cancer.markers[,6])
saveRDS(cancer.markers, file = "cancer.markers.0.25.rds")
