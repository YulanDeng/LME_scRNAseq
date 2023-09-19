## Script for plot Figure 6 and Figure S10
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

####################################
#Figure 6b-d,S10b (R version 4.0.3)#
####################################
#${workDir} is the working directory
#${TCGAclinicalDFS} is the DFS of TCGA data
#${TCGAclinicalOS} is the OS of TCGA data
#${TCGAclinical} is the clinical data of TCGA
#${EcotyeFile} is file for Ecotyper

#Load required packages
library(factoextra)
library(cluster)
library(survival)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library("survminer")

#Set working directory
setwd(workDir)

#Load the clinical information
load(TCGAclinicalDFS)
load(TCGAclinicalOS)
load(TCGAclinical)
TCGAclinicallabel <- substr(TCGAclinical[,1],14,15)
TCGAclinicalsub <- TCGAclinical[TCGAclinicallabel=="01",]
TCGAclinicalsub[,1] <- substr(TCGAclinicalsub[,1],9,12)

#Load the result of Ecotyper
abundance <- read.table(file=EcotyeFile,sep="\t",stringsAsFactors=F,header=T)

#Figure S10b
abundance_df <- as.data.frame(t(abundance))
abundance_scale_df <- scale(abundance_df)
dis_obj <- factoextra::get_dist(abundance_scale_df,method="pearson")
f2 <- fviz_nbclust(abundance_scale_df, kmeans, method = "silhouette", nstart = 25,diss=dis_obj)
f2

#Figure 6b
label <- km$"cluster"
labelF <- factor(as.character(label))
annotation_col = data.frame(
	cluster = labelF)
	rownames(annotation_col) = names(label)
	
cell_state_colo <- c("#F8766D", "#00BA38","#619CFF")
names(cell_state_colo) <- as.character(1:3)

ann_colors = list(
         cluster = cell_state_colo)
	
pheatmap(abundance[,order(labelF)],border_color=NA, 
	cluster_rows=F,cluster_cols=F,
	scale ="none",
	color=colorRampPalette(viridis(8))(50),
	show_rownames = T, show_colnames = F,
	annotation_col = annotation_col,
	annotation_colors = ann_colors)

#Figure 6c
labelN <- label
label_abu_T <- substr(names(labelN),14,15)
sampleAssign_abu_T <- labelN[label_abu_T=="01"]
names(sampleAssign_abu_T) <- substr(names(sampleAssign_abu_T),9,12)

intersmp_abu_OS <- intersect(names(sampleAssign_abu_T),TCGAclinicalOS[,1])
TCGAclinicalOSsubsmp <- TCGAclinicalOS[match(intersmp_abu_OS,TCGAclinicalOS[,1]),]
abu_OS <- sampleAssign_abu_T[intersmp_abu_OS]
abu_OS <- as.character(abu_OS)
yOS <-Surv(TCGAclinicalOSsubsmp[,3],TCGAclinicalOSsubsmp[,2])
resLME <- coxph(yOS~ abu_OS) 
Clusters <- as.character(abu_OS)
cliniT_abu_OS <- data.frame(recurrenceTime=TCGAclinicalOSsubsmp[,3],
label=TCGAclinicalOSsubsmp[,2],Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT_abu_OS)
Diff_Survival<-survdiff(yOS~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)
ggsurvplot(Sur_Curve,data=cliniT_abu_OS,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="Survival Probability(OS)")

#Figure 6d
intersmp_abu_DFS <- intersect(names(sampleAssign_abu_T),TCGAclinicalDFS[,1])
TCGAclinicalDFSsubsmp <- TCGAclinicalDFS[match(intersmp_abu_DFS,TCGAclinicalDFS[,1]),]
abu_DFS <- sampleAssign_abu_T[intersmp_abu_DFS]
abu_DFS <- as.character(abu_DFS)
yDFS <-Surv(TCGAclinicalDFSsubsmp[,3],TCGAclinicalDFSsubsmp[,2])
resLME <- coxph(yDFS~ abu_DFS) 
Clusters <- as.character(abu_DFS)
cliniT_abu_DFS <- data.frame(recurrenceTime=TCGAclinicalDFSsubsmp[,3],
label=TCGAclinicalDFSsubsmp[,2],Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT_abu_DFS)
Diff_Survival<-survdiff(yDFS~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)
ggsurvplot(Sur_Curve,data=cliniT_abu_DFS,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="Survival Probability(DFS)")

#################################
#Figure S10l-m (R version 4.0.4)#
#################################
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