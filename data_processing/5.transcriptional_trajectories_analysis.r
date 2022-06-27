## Script for the transcriptional trajectories analysis
##Yulan Deng, last updated 2022-6-27

###############################################################
#monocle analysis of monocyte and macrophage (R version 3.6.1)#
###############################################################
#${workDir} is the working directory
#${monoFile} is the seraut object of monocyte and macrophage

#Load required packages
library(monocle)
library(tidyverse)
library(Seurat)

#Set working directory
setwd(workDir)

#Load the seurat object of monocyte and macrophage
mono.integrated <- readRDS(file =monoFile)

#transfer the seurat obeject to monocle object
RA_matrix <-as(as.matrix(mono.integrated@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-mono.integrated@meta.data
label_new <- as.character(mono.integrated@meta.data[,"assigned_cell_type_byNC"])
label_new[label_new%in%c("CD14_monocyte","CD16_monocyte")] <- "monocyte"
label_new[label_new%in%c("TAM_anti_inflammatory")] <- "SPP1_macrophage"
label_new[label_new%in%c("cDC")] <- "DC_like_monocyte"
sample_ann<- cbind(sample_ann,label_new)
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
RA.cds<- newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())

#trajectories analysis
RA.cds <- estimateSizeFactors(RA.cds) 
RA.cds=detectGenes(RA.cds,min_expr = 0.1)
mono_marker_mono = FindMarkers(mono.integrated,ident.1 = c("CD14_monocyte","CD16_monocyte"), group.by = "assigned_cell_type_byNC" )
cDC_marker_mono = FindMarkers(mono.integrated,ident.1 = c("cDC"), group.by = "assigned_cell_type_byNC" )
early_marker_mono = FindMarkers(mono.integrated,ident.1 = c("early_stage_macrophages"), group.by = "assigned_cell_type_byNC" )
TAM_marker_mono = FindMarkers(mono.integrated,ident.1 = c("TAM_anti_inflammatory"), group.by = "assigned_cell_type_byNC" )
marker_gene=rbind(mono_marker_mono,rbind(cDC_marker_mono,rbind(early_marker_mono,TAM_marker_mono)))
test_ordering_genes=unique(rownames(marker_gene[(marker_gene[,2]> 1 )&(marker_gene[,5]<0.05),]))
RA.cds=setOrderingFilter(RA.cds,ordering_genes = test_ordering_genes) 
RA.cds=reduceDimension(RA.cds,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed")
RA.cds=orderCells(RA.cds,root_state=1)
plot_cell_trajectory(RA.cds,color_by = "label_new")
plot_cell_trajectory(RA.cds,color_by = "Pseudotime")
BEAM_res <- BEAM(RA.cds, branch_point = 1, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

####################################################
#monocle analysis of CD8+ T cells (R version 3.6.1)#
####################################################
#${workDir} is the working directory
#${CD8Tfile} is the seraut object of CD8+ T cells
#${CD8Tmarkers} is the marker gene of CD8+ T cells

#Load required packages
library(monocle)
library(tidyverse)
library(Seurat)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the seurat object of CD8+ T cells
CD8T.integrated <- readRDS(file = CD8Tfile)

#transfer the seurat obeject to monocle object
RA_matrix <-as(as.matrix(CD8T.integrated@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-CD8T.integrated@meta.data
sample_ann<- cbind(sample_ann,label)
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
RA.cds<- newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())

#trajectories analysis
RA.cds <- estimateSizeFactors(RA.cds) 
RA.cds <- estimateDispersions(RA.cds)
RA.cds=detectGenes(RA.cds,min_expr = 1)
CD8T.markers <- readRDS(file = CD8Tmarkers)
test_ordering_genes=unique(CD8T.markers[(CD8T.markers[,"p_val_adj"]<0.05)&(CD8T.markers[,"avg_logFC"]>log(1.8)&(CD8T.markers[,"cluster"]%in%c("0","1","3"))),"gene"])
RA.cds=setOrderingFilter(RA.cds,ordering_genes = test_ordering_genes) 
RA.cds=reduceDimension(RA.cds,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 
RA.cds=orderCells(RA.cds)
RA.cds=orderCells(RA.cds,root_state=4)
plot_cell_trajectory(RA.cds,color_by = "assigned_cell_type_byNC",cell_size = 0.001)
plot_cell_trajectory(RA.cds,color_by = "Pseudotime",cell_size = 0.001)
BEAM_res <- BEAM(RA.cds, branch_point = 1, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

####################################################
#monocle analysis of CD4+ T cells (R version 3.6.1)#
####################################################
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
plot_cell_trajectory(RA.cds,color_by = "Pseudotime",)
BEAM_res <- BEAM(RA.cds, branch_point = 2, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

#############################
#cytoTRACE (R version 4.0.3)#
#############################
#${pythonDir} is the directory of python
#${workDir} is the working directory
#${cancerFile} is the seurat object of cancer cell
#${clinicalFile} is the file of clinical information
#${smp} is the sample ID

#all sample from GGO
#Load required packages
library("reticulate")
Sys.setenv(RETICULATE_PYTHON=pythonDir)
library(CytoTRACE)
library("Seurat")
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the seurat object of cancer cell
cancer.integrated <- readRDS(file =cancerFile)

#Load the clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)

#make input file for cytoTRACE
GGO <- clinical[clinical[,"class"]%in%c("pGGO","s25GGO","s50GGO","s75GGO","s100GGO","HiDenGGO"),"rawDataID"]
#s25GGO,s50GGO,s75GGO,s100GGO and HiDenGGO are called 
#GGO25, GGO50, GGO75, GGO100 and dGGO in the manuscript
cancer.integrated <- subset(cancer.integrated,cells=rownames(cancer.integrated@meta.data)[(as.character(cancer.integrated@meta.data[,"orig.ident"])%in%GGO)])
marrow_10x_expr <- as.matrix(cancer.integrated@assays$RNA@counts)
marrow_10x_pheno <- as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])
names(marrow_10x_pheno) <- rownames(cancer.integrated@meta.data)

#run cytoTRACE
results <- CytoTRACE(marrow_10x_expr,ncores = 32)
cols<-brewer.pal(9,"Blues")[9:2]
pal<-colorRampPalette(cols)
mycolors<-pal(100)
plotCytoTRACE(results, phenotype = marrow_10x_pheno,colors=mycolors)

#cytoTRACE for each sample
#Load required packages
Sys.setenv(RETICULATE_PYTHON=pythonDir)
library(CytoTRACE)
library("Seurat")

#Set working directory
setwd(workDir)

#Load the seurat object of cancer cell
cancer.integrated <- readRDS(file =cancerFile)

#only clusters with more than 5% and 20 cells are used for cytoTRACE analysis
cluster_n <- table(cancer.integrated@meta.data[as.character(cancer.integrated@meta.data[,"orig.ident"])==smp,
"integrated_snn_res.0.25"])
cluster_freq <- cluster_n/sum(cluster_n)
cluster_sel <- names(cluster_n)[(cluster_n>=20)&(cluster_freq>=0.05)]
cancer.integrated <- subset(cancer.integrated,cells=rownames(cancer.integrated@meta.data)[(as.character(cancer.integrated@meta.data[,"orig.ident"])==smp)&(as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])%in%cluster_sel)])

#make the input for cytoTRACE
marrow_10x_expr <- as.matrix(cancer.integrated@assays$RNA@counts)
marrow_10x_pheno <- as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])
names(marrow_10x_pheno) <- rownames(cancer.integrated@meta.data)

##run cytoTRACE
results <- CytoTRACE(marrow_10x_expr,ncores = 32)
plotCytoTRACE(results, phenotype = marrow_10x_pheno)

###############
#RNA velocity #
###############
#${cellRangerDir} is the directory of cellranger result
#${gftFile} is the gtf file of gene
#${pythonDir} is the direction of python
#${smp} is sample ID
#${workDir} is the working directory
#${cancerFile} is the seurat object of cancer cell
#${loomFile} is the loom file of cellranger result

#generation of loom file
velocyto run10x -@ 32 -m repeat_msk.gtf ${cellRangerDir} ${gftFile}

#RNA velocity for each sample
#save h5 file
#Load required packages
library("reticulate")
library("ggplot2")
library("SingleCellExperiment")
library("Seurat")
library(sceasy)
use_python(pythonDir)

#Set working directory
setwd(workDir)

#Load the seurat object of cancer cell
cancer.integrated <- readRDS(file =cancerFile)

#only clusters with more than 5% and 20 cells are used for RNA velocity analysis
cluster_n <- table(cancer.integrated@meta.data[as.character(cancer.integrated@meta.data[,"orig.ident"])==smp,"integrated_snn_res.0.25"])
cluster_freq <- cluster_n/sum(cluster_n)
cluster_sel <- names(cluster_n)[(cluster_n>=20)&(cluster_freq>=0.05)]
cancer.integrated <- subset(cancer.integrated,cells=rownames(cancer.integrated@meta.data)[(as.character(cancer.integrated@meta.data[,"orig.ident"])==smp)&(as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])%in%cluster_sel)])

#save h5 file
sceasy::convertFormat(cancer.integrated, from="seurat", to="anndata",
                       outFile='cancer.all.anndata.h5ad')

#make input file for RNA velocity
import scvelo as scv
scv.settings.presenter_view = True  
import matplotlib.pyplot as pl
import scanpy as sc
import loompy 
import numpy as np
import pandas as pd

adatal = scv.read(loomFile, sparse=True, cache=True) 
adatal.var_names_make_unique() 
adatas = scv.read('cancer.integrated.h5ad', sparse=True, cache=True)
adatas.var_names_make_unique() 
gene1=adatal.var_names.to_list()
gene2=adatas.var_names.to_list()
intergene = list(set(gene1).intersection(set(gene2)))
intersample=adatas.obs_names.to_list()
intersampletmp=[intersample[0].split('_')[0]]
for i in intersample[1:]:     
   intersampletmp.append(i.split('_')[0])
intersample1=[]
for i in intersampletmp:     
   intersample1.append("run_count_"+smp+":"+i+"x")  
intersample_cluster = adatas.obs["integrated_snn_res.0.25"].to_list()
adatal_sub = adatal[intersample1,intergene]
adatal_sub.obsm['X_pca'] = adatas.obsm['X_pca']
adatal_sub.obsm['X_tsne'] = adatas.obsm['X_tsne']
adatal_sub.obsm['X_umap'] = adatas.obsm['X_umap']
num_abc = pd.Series(intersample_cluster, index=intersample1)
adatal_sub.obs['clusters'] = num_abc
adatal_sub.write('adatal_sub.h5ad', compression='gzip')

#run RNA velocity
import scvelo as scv 
scv.settings.presenter_view = True  
import matplotlib.pyplot as pl
import os

os.chdir(workDir)
adata = scv.read('adatal_sub.h5ad', sparse=True, cache=True)
adata.var_names_make_unique()
adata
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=40, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata) 
adata.write('adata_dynamic.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(adata,basis="tsne", color="clusters",save="cancer"+ smp +".dynamic.velocity.tsne") 
scv.pl.velocity_embedding_stream(adata,basis="umap", color="clusters",save="cancer"+ smp +".dynamic.velocity.umap") 
