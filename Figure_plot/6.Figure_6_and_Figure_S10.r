## Script for plot Figure 6 and Figure S10
##Yulan Deng, last updated 2024-1-5
##my e-mail:kndeajs@163.com

############################################
#Figure 6c, Figure S10b-c (R version 4.0.4)#
############################################
#${workDir} is the working directory
#${TCGAclinicalDFS} is the DFS of TCGA data
#${TCGAclinicalOS} is the OS of TCGA data
#${abundance_df.TCGA} is file abundance of Ecotyper for TCGA datasets

#Load required packages
library(glmnet)
library(survival)
library("survminer")

#Set working directory
setwd(workDir)

#Load the clinical information
load(TCGAclinicalDFS)
load(TCGAclinicalOS)
load(abundance_df.TCGA)

#get the abundance of ecotype for tumor samples from TCGA
label_abu_T <- substr(rownames(abundance_df),14,15)
sampleAssign_abu_T <- abundance_df[label_abu_T=="01",]
rownames(sampleAssign_abu_T) <- substr(rownames(sampleAssign_abu_T),9,12)

##using DFS data to train elastic net regularized Cox regression model
intersmp_abu_DFS <- intersect(rownames(sampleAssign_abu_T),TCGAclinicalDFS[,1])
y <- TCGAclinicalDFS[match(intersmp_abu_DFS,TCGAclinicalDFS[,1]),3:2]
x <- sampleAssign_abu_T[intersmp_abu_DFS,]
colnames(y) <- c("time","status")
y <- as.matrix(y)
rownames(y) <- NULL
x <- as.matrix(x)
rownames(x) <- NULL
yDFS=Surv(y[,"time"],y[,"status"])
fit = glmnet(x, yDFS, family = "cox",alpha=0.5,penalty.factor=c(1,1,10,0,1,1))
set.seed(9)
cvfit <- cv.glmnet(x, y, family = "cox",alpha=0.5, penalty.factor=c(1,1,10,0,1,1),type.measure = "C",nfolds = 10)

#Figure S10b
plot(cvfit)

#Figure S10c
eff <- coef(fit, s = 0.05)
scoreDFS <- apply(x,1,function(zz) sum(zz*eff[,1]))
TCGAclinicalDFSsubsmp <- TCGAclinicalDFS[match(intersmp_abu_DFS,TCGAclinicalDFS[,1]),]
yDFS <-Surv(TCGAclinicalDFSsubsmp[,3],TCGAclinicalDFSsubsmp[,2])
Clusters <- rep(1,length(scoreDFS))
Clusters[scoreDFS>quantile(scoreDFS,probs=2/3)] <- 3
Clusters[(scoreDFS<=quantile(scoreDFS,probs=2/3))&(scoreDFS>quantile(scoreDFS,probs=1/3))] <- 2
cliniT_abu_DFS <- data.frame(recurrenceTime=TCGAclinicalDFSsubsmp[,3],
label=TCGAclinicalDFSsubsmp[,2],Clusters=Clusters,stringsAsFactors=F)
##5-year survival
cliniT_abu_DFS[cliniT_abu_DFS[,"recurrenceTime"]>60,"recurrenceTime"] <- 60
cliniT_abu_DFS[cliniT_abu_DFS[,"recurrenceTime"]>60,"label"] <- 0
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT_abu_DFS)
Diff_Survival<-survdiff(yDFS~Clusters)
ggsurvplot(Sur_Curve,data=cliniT_abu_DFS,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="Survival Probability", pval = TRUE)

#Figure 6c
intersmp_abu_OS <- intersect(rownames(sampleAssign_abu_T),TCGAclinicalOS[,1])
TCGAclinicalOSsubsmp <- TCGAclinicalOS[match(intersmp_abu_OS,TCGAclinicalOS[,1]),]
xOS <- sampleAssign_abu_T[intersmp_abu_OS,]
OS_score <- apply(xOS,1,function(zz) sum(zz*eff[,1]))
yOS <-Surv(TCGAclinicalOSsubsmp[,3],TCGAclinicalOSsubsmp[,2])
Clusters <- rep(1,length(OS_score))
Clusters[OS_score>quantile(OS_score,probs=2/3)] <- 3
Clusters[(OS_score<=quantile(OS_score,probs=2/3))&(OS_score>quantile(OS_score,probs=1/3))] <- 2
cliniT_abu_OS <- data.frame(recurrenceTime=TCGAclinicalOSsubsmp[,3],
label=TCGAclinicalOSsubsmp[,2],Clusters=Clusters,stringsAsFactors=F)
#5-year survival
cliniT_abu_OS[cliniT_abu_OS[,"recurrenceTime"]>60,"recurrenceTime"] <- 60
cliniT_abu_OS[cliniT_abu_OS[,"recurrenceTime"]>60,"label"] <- 0
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT_abu_OS)
Diff_Survival<-survdiff(yOS~Clusters)
ggsurvplot(Sur_Curve,data=cliniT_abu_OS,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="Survival Probability", pval = TRUE)

############################################
#Figure 6b, Figure S10d-e (R version 4.0.3)#
############################################
#${workDir} is the working directory
#${abundance_df.TCGA} is file abundance of Ecotyper for TCGA datasets
#${eff} is the coefficients of elastic net regularized Cox regression model
#${rank_data.file} is the rank file of Ecotyper
#${ecotype.file} describes the composition of cell states for each ecotype
#${pathD} is the result directory for Ecotyper

#Load required packages
library(RColorBrewer)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load the Ecotype abundance
abundance <- read.table(file=abundance_df.TCGA,sep="\t",stringsAsFactors=F,header=T)
abundance_df <- as.data.frame(t(abundance))

#Load coefficients of elastic net regularized Cox regression model
load(eff)

#calculation of LIME score and group patients into top,median and lower tertiles
label_abu_T <- substr(rownames(abundance_df),14,15)
sampleAssign_abu_T <- abundance_df[label_abu_T=="01",]
rownames(sampleAssign_abu_T) <- substr(rownames(sampleAssign_abu_T),9,12)
score <- apply(sampleAssign_abu_T,1,function(zz) sum(zz*eff[,1]))
Clusters <- rep(1,length(score))
Clusters[score>quantile(score,probs=2/3)] <- 3
Clusters[(score<=quantile(score,probs=2/3))&(score>quantile(score,probs=1/3))] <- 2
names(Clusters) <- rownames(sampleAssign_abu_T)

#Figure 6b top
labelF <- factor(as.character(Clusters),levels=c("1","2","3"))
annotation_col = data.frame(
	cluster = labelF,
	score=score)
rownames(annotation_col) = names(Clusters)
cell_state_colo <- brewer.pal(8, "Set3")[1:3]
names(cell_state_colo) <- as.character(1:3)
ann_colors = list(
         cluster = cell_state_colo)
pheatmap(t(sampleAssign_abu_T)[,order(score,decreasing=F)],border_color=NA, 
	cluster_rows=F,cluster_cols=F,
	scale ="none",
	color=colorRampPalette(viridis(8))(50),
	show_rownames = T, show_colnames = F,
	annotation_col = annotation_col,
	annotation_colors = ann_colors)

##Figure 6b bottom
rank_data <- read.table(file=rank_data.file,sep="\t",stringsAsFactors=F,header=T)
ecotype <- read.table(file=ecotype.file,sep="\t",stringsAsFactors=F,header=T)
abundance_list <- lapply(seq(nrow(rank_data)),function(x){
	cell_type=rank_data[x,1]
	abundance_tmp_df=read.table(file =paste0(pathD,
	cell_type,"/state_abundances.txt"),sep="\t",stringsAsFactors=F, header=T)
	return(abundance_tmp_df)
})
names(abundance_list) <- rank_data[,1]
cell_state_abundance <- sapply(seq(nrow(ecotype)),function(x){
	abundance_tmp <- abundance_list[ecotype[x,"CellType"]][[1]][ecotype[x,"State"],]
	return(unlist(abundance_tmp))
})
colnames(cell_state_abundance) <- ecotype[,"ID"]
label_abu_T2 <- substr(rownames(cell_state_abundance),14,15)
sampleAssign_abu_T2 <- cell_state_abundance[label_abu_T2=="01",]
rownames(sampleAssign_abu_T2) <- substr(rownames(sampleAssign_abu_T2),9,12)
annotation_row = data.frame(
	feature_cluster = factor(ecotype[,"Ecotype"]))
	rownames(annotation_row) = ecotype[,"ID"]
annotation_col = data.frame(
	cluster = labelF,
	score=score)
	rownames(annotation_col) = names(Clusters)
pheatmap(t(sampleAssign_abu_T2)[,order(score,decreasing=F)],border_color=NA, 
	cluster_rows=F,cluster_cols=F,
	scale ="row",
	annotation_col = annotation_col,
	annotation_row = annotation_row,
	color=colorRampPalette(c(brewer.pal(9, "Blues")[8],"white",brewer.pal(9, "Reds")[8]))(50),
	show_rownames = T, show_colnames = F)

##Figure S10d
colo <- brewer.pal(8, "Set1")[1:3]
plot(sampleAssign_abu_T[,3],sampleAssign_abu_T[,6],
			xlab="LIME03",
			ylab="LIME06",sub=
			paste0("cor=",
			round(cor(sampleAssign_abu_T[,3],sampleAssign_abu_T[,6],
			method = "spearman"),3),
			";p=",round(cor.test(sampleAssign_abu_T[,3],sampleAssign_abu_T[,6],method = "spearman")$"p.value",4)),pch=19,
			col=colo[Clusters])
			abline(lm(sampleAssign_abu_T[,6]~sampleAssign_abu_T[,3]), lwd=2, col="black")

##Figure S10e
colo <- brewer.pal(8, "Set1")[1:3]
plot(sampleAssign_abu_T[,5],sampleAssign_abu_T[,6],
			xlab="LIME05",
			ylab="LIME06",sub=
			paste0("cor=",
			round(cor(sampleAssign_abu_T[,5],sampleAssign_abu_T[,6],
			method = "spearman"),3),
			";p=",round(cor.test(sampleAssign_abu_T[,5],sampleAssign_abu_T[,6],method = "spearman")$"p.value",4)),pch=19,
			col=colo[Clusters])
			abline(lm(sampleAssign_abu_T[,6]~sampleAssign_abu_T[,5]), lwd=2, col="black")


##################################
#Figure S10f, l (R version 4.0.3)#
##################################
#${workDir} is the working directory
#${LUAD.tmb_table.cBio.maftools} is the mutation file from TCGA
#${TCGAclinical} is the clinical data of TCGA
#${TCGAclinicalOS} is the OS of TCGA data
#${abundance_df.TCGA} is file abundance of Ecotyper for TCGA datasets
#${eff} is the coefficients of elastic net regularized Cox regression model

#Load required packages
library(survival)
library(RColorBrewer)
library(ggpubr)
library("survminer")

#Set working directory
setwd(workDir)

#Load the clinical information and  mutation from TCGA, as well as ecotype abundance
load(LUAD.tmb_table.cBio.maftools)
load(TCGAclinicalOS)
load(eff)
load(TCGAclinical)

#Figure S10f
TCGAclinicallabel <- substr(TCGAclinical[,1],14,15)
TCGAclinicalsub <- TCGAclinical[TCGAclinicallabel=="01",]
TCGAclinicalsub[,1] <- substr(TCGAclinicalsub[,1],9,12)
abundance <- read.table(file=abundance_df.TCGA,sep="\t",stringsAsFactors=F,header=T)
abundance_df <- as.data.frame(t(abundance))
abundance_scale_df <- scale(abundance_df)
labelN <- apply(abundance_df,1,function(zz) sum(zz*eff[,1]))
names(labelN) <- rownames(abundance_scale_df)
label_abu_T <- substr(names(labelN),14,15)
sampleAssign_abu_T <- labelN[label_abu_T=="01"]
names(sampleAssign_abu_T) <- substr(names(sampleAssign_abu_T),9,12)
##TMB
tmb_v <- unlist(tmb_table[,"total_perMB"])
names(tmb_v) <- unlist(tmb_table[,"Tumor_Sample_Barcode"])
labelTMB_T <- substr(names(tmb_v),14,15)
names(tmb_v) <- substr(names(tmb_v),9,12)
interTMB <- intersect(names(sampleAssign_abu_T),names(tmb_v))
sampleAssign_abu_T_tmb <- sampleAssign_abu_T[interTMB]
tmb_v_tmb <- tmb_v[interTMB]
Clusters <- rep(1,length(sampleAssign_abu_T_tmb))
Clusters[sampleAssign_abu_T_tmb>quantile(sampleAssign_abu_T_tmb,probs=2/3)] <- 3
Clusters[(sampleAssign_abu_T_tmb<=quantile(sampleAssign_abu_T_tmb,probs=2/3))&(sampleAssign_abu_T_tmb>quantile(sampleAssign_abu_T_tmb,probs=1/3))] <- 2
TMB_df <- data.frame(
LME_class=factor(as.character(Clusters),levels=c("1","2","3")),
LMB=tmb_v_tmb,stringsAsFactors=F)
my_comparisons <- list()
my_comparisons[[1]] <- c("1","3")
p<-ggboxplot(TMB_df, x = "LME_class", y = "LMB", color="LME_class",outlier.shape = NA,
palette =  brewer.pal(8, "Set2")[1:3],legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("TMB/Mb")+
		  stat_compare_means(comparisons = my_comparisons,size = 8)
p

#Figure S10l
EGFR_df <-all_mut@"data"[(unlist(all_mut@"data"[,"Hugo_Symbol"])=="EGFR")&(unlist(all_mut@"data"[,"COSMIC_overlapping_mutations"])!=""),]
EGFR_smp <- unique(unlist(EGFR_df[,"Tumor_Sample_Barcode"]))
EGFRsample <- substr(as.character(EGFR_smp),9,12)
intersmp_abu_OS <- intersect(intersect(names(sampleAssign_abu_T),TCGAclinicalOS[,1]),EGFRsample)
TCGAclinicalOSsubsmp <- TCGAclinicalOS[match(intersmp_abu_OS,TCGAclinicalOS[,1]),]
abu_OS <- sampleAssign_abu_T[intersmp_abu_OS]
Clusters <- rep(1,length(abu_OS))
Clusters[abu_OS>quantile(abu_OS,probs=2/3)] <- 3
Clusters[(abu_OS<=quantile(abu_OS,probs=2/3))&(abu_OS>quantile(abu_OS,probs=1/3))] <- 2
yOS <-Surv(TCGAclinicalOSsubsmp[,3],TCGAclinicalOSsubsmp[,2])
Clusters1 <- as.numeric(as.character(Clusters)=="1")
cliniT_abu_OS <- data.frame(recurrenceTime=TCGAclinicalOSsubsmp[,3],
label=TCGAclinicalOSsubsmp[,2],Clusters=Clusters1,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT_abu_OS)
Diff_Survival<-survdiff(yOS~Clusters1)
ggsurvplot(Sur_Curve,data=cliniT_abu_OS,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=20,ylab="Survival Probability", pval = TRUE)
