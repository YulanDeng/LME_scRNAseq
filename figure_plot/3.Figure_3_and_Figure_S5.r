## Script for plot Figure 3 and Figure S5
##Yulan Deng, last updated 2022-6-27

#############################
#Figure 3C (R version 4.0.4)#
#############################
#${workDir} is the working directory
#${survivalFileTCGA} is the TCGA data used for KM plot
#${survivalFileGSE31210} is the GSE31210 data used for KM plot

#Load required packages
library("survminer")

#Set working directory
setwd(workDir)

#Load TCGA data used for KM plot
load(survivalFileTCGA)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="c2signature_TCGA_DFS")

#Load SE31210 data used for KM plot
load(survivalFileGSE31210)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="c2signature_GSE31210_DFS")

##############################
#Figure S5C (R version 4.0.4)#
##############################
#${workDir} is the working directory
#${survivalFileGSE140343} is the GSE140343 data used for KM plot

#Load required packages
library("survminer")

#Set working directory
setwd(workDir)

#Load GSE140343 data used for KM plot
load(survivalFileGSE140343)

#KM plot
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="c2signature_GSE140343_DFS")

#############################################
#Figure 3B,Figure S5A-B(python vesion 3.8.5)#
#############################################
#${smp} is sample ID
#${workDir} is the working directory
import scvelo as scv 
scv.settings.presenter_view = True  
import matplotlib.pyplot as pl
import os
os.chdir(workDir + smp)
adata = scv.read('adata_dynamic.h5ad', sparse=True, cache=True)
scv.pl.velocity_embedding_stream(adata,basis="umap", color="clusters",save="cancer"+ smp +".dynamic.velocity.umap.png",palette=['#66C2A5','#FC8D62','#8DA0CB','#E78AC3'],legend_loc='none') 