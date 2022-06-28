## Script for plot Figure 6 and Figure S8
##Yulan Deng, last updated 2022-6-28

#############################
#Figure 6A (R version 4.0.3)#
#############################
#${workDir} is the working directory
#${targetRegionFile} is the target region of WES
#${clinicalFile} is the clinical information
#${smp} is sample ID
#${annovarDir} is the directory of ANNOVAR results
#${netMHCpanDir} is the directory of netMHCpan results

#Load required packages
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#Set working directory
setwd(workDir)

##TMB
hg38 <- read.table(file=targetRegionFile,sep="\t",quote = "",stringsAsFactors=F,header=F)
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")
clinical <- cbind(clinical,clinical_label)
clinical[clinical[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical[clinical[,"class"]=="s25GGO","class"] <- "GGO25"
clinical[clinical[,"class"]=="s50GGO","class"] <- "GGO50"
clinical[clinical[,"class"]=="s75GGO","class"] <- "GGO75"
clinical[clinical[,"class"]=="s100GGO","class"] <- "GGO100"
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NoAden"
colo <- c(brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")
anno_list <- lapply(smp,function(x){
	flnm <- paste0(annovarDir,x,".hg38_multianno.txt")
	fl <- read.table(file=flnm,sep="\t",quote = "",stringsAsFactors=F,header=T)
	fl <- cbind(fl,rep(x,nrow(fl)))
	fl <- fl[fl[,6]%in%c(  "exonic","exonic;splicing"),]
	fl <- fl[!(fl[,9]%in%c("synonymous SNV",".","unknown")),]
	colnames(fl)[12] <- "Tumor_Sample"
	return(fl)
})
anno_df <- do.call(rbind,anno_list)
no_mut <- sapply(anno_list,nrow)
names(no_mut) <- c(smp,smp2,smp3)
TMB <- no_mut*10^6/sum(hg38[,3]-hg38[,2])
TMB_df <- data.frame(stage=factor(clinical[match(names(no_mut),clinical[,"clinical_label"]),"class"],
levels=c("pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")),
TMB=TMB,stringsAsFactors=F)
no_mut_list <- split(TMB,
f=factor(clinical[match(names(no_mut),clinical[,"clinical_label"]),"class"],
levels=c( "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")))
p<-ggboxplot(TMB_df, x = "stage", y = "TMB", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("Mutations per megabase")
p

#neoantigen
anno_list <- lapply(smp,function(y){
	flnm <- paste0(netMHCpanDir,y,"_A/",y,"_NetChopcut2netMHCpan.xls")
	flnm2 <- paste0(netMHCpanDir,y,"_A/",y,"NetChopOut_pep.txt")
	netMHCpanRes <- read.table(file =flnm,sep="\t",quote ="",stringsAsFactors=F,fill=T,header=F)
	netChop_anno <- read.table(file =flnm2,sep=" ",quote ="",stringsAsFactors=F,fill=T,header=T)
	hla <- setdiff(unlist(netMHCpanRes[1,]),"")
	netMHCshort <- netMHCpanRes[-c(1:2),]
	binding_list <- lapply(seq(length(hla)),function(x){
		i=x
		idRank=3+(i-1)*6+4
		indexbind <- which(as.numeric(netMHCshort[,idRank]) < 2)
		if(length(indexbind)>0)
		{
			res_df <- data.frame(Peptide=netMHCshort[indexbind,2],
			core=netMHCshort[indexbind,3+(i-1)*6+1],icore=netMHCshort[indexbind,3+(i-1)*6+2],
			EL_score=netMHCshort[indexbind,3+(i-1)*6+3], EL_Rank=netMHCshort[indexbind,3+(i-1)*6+4],
			HLA=hla[x],stringsAsFactors=F)
			return(res_df)
		}
		
	})
	binding_df <- do.call(rbind,binding_list)
	binseq <- match(binding_df[,"Peptide"],netChop_anno[,"peptide"])
	binding_seq_df <- cbind(binding_df,netChop_anno[binseq,],stringsAsFactors=F)
	return(binding_seq_df)
})
names(anno_list) <- c(clinical[match(smp,clinical[,"clinical_label"]),"rawDataID"],
anno_df <- do.call(rbind,anno_list)
no_pep <- sapply(anno_list,function(x) length(unique(x[,"header"])))
names(no_pep) <- c(smp,smp2,smp3)
pep_df <- data.frame(stage=factor(clinical[match(names(no_pep),clinical[,"clinical_label"]),"class"],
levels=c( "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")),
pep=no_pep,stringsAsFactors=F)
no_pep_list <- split(no_pep,
f=factor(clinical[match(names(no_pep),clinical[,"clinical_label"]),"class"],
levels=c( "pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")))
p<-ggboxplot(pep_df, x = "stage", y = "pep", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("Number of predicted binding epitopes")
p

###############################
#Figure 6B-C (R version 4.0.3)#
###############################
#${workDir} is the working directory
#${targetRegionFile} is the target region of WES
#${clinicalFile} is the clinical information
#${smp} is sample ID
#${annovarDir} is the directory of ANNOVAR results
#${cgc} is list of CGC genes

#Load required packages
library(gplots)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")
colo <- c(brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("pGGO","HiDenGGO","s25GGO","s50GGO","s75GGO","s100GGO",
"Solid1","Solid3","SolidN")
colo2 <- brewer.pal(6, "Set3")

#load mutation reuslts
anno_list <- lapply(smp,function(x){
	flnm <- paste0(annovarDir,x,".hg38_multianno.txt")
	fl <- read.table(file=flnm,sep="\t",quote = "",stringsAsFactors=F,header=T)
	fl <- cbind(fl,rep(x,nrow(fl)))
	fl <- fl[fl[,6]%in%c(  "exonic","exonic;splicing"),]
	fl <- fl[!(fl[,9]%in%c("synonymous SNV",".","unknown")),]
	colnames(fl)[12] <- "Tumor_Sample"
	return(fl)
})
anno_df <- do.call(rbind,anno_list)
anno_df <- unique(anno_df[,c("Gene.refGene","ExonicFunc.refGene","Tumor_Sample")])
anno_cgc_df <- anno_df[anno_df[,"Gene.refGene"]%in%cgc,]
cgc_gene <- names(table(anno_cgc_df[,"Gene.refGene"]))[table(anno_cgc_df[,"Gene.refGene"])>0]
cgc_gene <- cgc_gene[order(table(anno_cgc_df[anno_cgc_df[,"Gene.refGene"]%in%cgc_gene,"Gene.refGene"]),decreasing=T)]

#make the matrix  for heatmap
mt <- matrix(0,nrow=length(cgc_gene),ncol=length(anno_list))
rownames(mt) <- cgc_gene
colnames(mt) <- names(anno_list)
for(i in seq(nrow(anno_cgc_df)))
{
	if(anno_cgc_df[i,"Gene.refGene"]%in%cgc_gene)
	{
		if(anno_cgc_df[i,"ExonicFunc.refGene"]=="nonsynonymous SNV")
		{
			mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]] <- mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]]+1
		}else if(anno_cgc_df[i,"ExonicFunc.refGene"]=="nonframeshift substitution"){
			mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]] <- mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]]+2
		}else if(anno_cgc_df[i,"ExonicFunc.refGene"]=="stoploss"){
			mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]] <- mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]]+4
		}else if(anno_cgc_df[i,"ExonicFunc.refGene"]=="frameshift substitution"){
			mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]] <- mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]]+8
		}else if(anno_cgc_df[i,"ExonicFunc.refGene"]=="stopgain"){
			mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]] <- mt[anno_cgc_df[i,"Gene.refGene"],anno_cgc_df[i,"Tumor_Sample"]]+16
		}
	}
}
mt <- mt[,
order(as.numeric(factor(clinical[match(names(anno_list),clinical[,"clinical_label"]),"class"],
levels=c( "pGGO","HiDenGGO","s25GGO","s50GGO","s75GGO","s100GGO","Solid1","Solid3","SolidN"))),decreasing=F)]
mt <- mt[apply(mt,1,sum)>1,]
mt <- mt[order(apply(mt,1,function(x) sum(x>0)),decreasing=T),]

#plot heatmap
heatmap.2 (mt,Rowv = F,Colv=F,dendrogram ="none",scale = "none",
breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,7.5,8.5,15.5,16.5,32),
col=c("white",colo2[1],colo2[4],"black",colo2[5],"black",colo2[3],"black",colo2[6],"black"),
trace="none",
ColSideColors=colo[clinical[match(colnames(mt),clinical[,"clinical_label"]),"class"]],
key = TRUE,keysize = 1.5,density.info="none")

##Figure 6C
anno_df_EGFR <- anno_df[anno_df[,"Gene.refGene"]=="EGFR",]
fractionBelow75 <- sum(clinical[match(unique(anno_df_EGFR[,"Tumor_Sample"]),clinical[,"clinical_label"]),"class"]%in%c( "pGGO","HiDenGGO","s25GGO","s50GGO"))/sum(clinical[match(names(anno_list),clinical[,"clinical_label"]),"class"]%in%c( "pGGO","HiDenGGO","s25GGO","s50GGO"))
#<75
statBelow75 <- c(fractionBelow75,1-fractionBelow75)
names(statBelow75) <- c("EGFRmut","EGFRwt")
pie(statBelow75,col=brewer.pal(6, "Set3")[5:4],main="solid portion < 75%")
##>75
fractionAbove75 <- sum(clinical[match(unique(anno_df_EGFR[,"Tumor_Sample"]),clinical[,"clinical_label"]),"class"]%in%c( "s75GGO","s100GGO",
"Solid1","Solid3","SolidN"))/sum(clinical[match(names(anno_list),clinical[,"clinical_label"]),"class"]%in%c( "s75GGO","s100GGO",
"Solid1","Solid3","SolidN"))
statAbove75 <- c(fractionAbove75,1-fractionAbove75)
names(statAbove75) <- c("EGFRmut","EGFRwt")
pie(statAbove75,col=brewer.pal(6, "Set3")[5:4],main="solid portion > 75%")

#############################
#Figure 6D (R version 4.0.3)#
#############################
#${workDir} is the working directory
#${trainLabel} is the label of training datasets
#${trainPredict} is the prediction of training datasets

#Load required packages
library(pROC)

#Set working directory
setwd(workDir)

#Load the label and prediction of training datasets
train_label <- read.table(file = rainLabel,sep=",",quote = "",stringsAsFactors=F,header=T)
train_predict <- read.table(file = trainPredict,sep=",",quote = "",stringsAsFactors=F,header=F)

#calculate ROC and best threshold
rocobj <- roc(train_label[,1], train_predict[,1])
roc_result <- coords(rocobj, "best",best.weights=c(0.1, 0.3))

#plot
plot(rocobj,
     legacy.axes = TRUE,
     main="best_threshold_for_AUC",
     thresholds=0.81595, 
     print.thres=0.81595) 

#############################
#Figure 6F (R version 4.0.3)#
#############################
#${inputTCRAI} is the input of model on our datasets
#${resTCRAI} is the result of model on our datasets
#${PredictFile} is the prediction of our datasets
#${clinicalFile} is the clinical information
#${netMHCpanDir} is the directory of netMHCpan result
#${smp} is sample ID
#${immuneFile} is the seurat object of immune cell
#${TRUSTdir} is the directory of TRUST result

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(gplots)
library(ggseqlogo)

#Set working directory
setwd(workDir)

#Load the result of model on our datasets
TCRAIinput <- read.table(file = inputTCRAI,sep=",",quote = "",stringsAsFactors=F,header=T)
TCRAIres <- read.table(file = resTCRAI,sep=",",quote = "",stringsAsFactors=F,header=T)
TCRAIpredict <- read.table(file = PredictFile,sep=",",quote = "",stringsAsFactors=F,header=F)
TCRAIres_predict <- cbind(TCRAIres,TCRAIpredict)
inputlabel <- paste(TCRAIinput[,"TRB_cdr3"],TCRAIinput[,"TRB_v_gene"],TCRAIinput[,"TRB_j_gene"],sep="_")
prelabel <- paste(TCRAIres_predict[,"TRB_cdr3"],TCRAIres_predict[,"TRB_v_gene"],
TCRAIres_predict[,"TRB_j_gene"],sep="_")
inputPredict <- TCRAIres_predict[match(inputlabel,prelabel),"V1"]
TCRAIinput_predict <- cbind(TCRAIinput,inputPredict)

#Load the clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")

#Load the result of neoantigen
anno_list <- lapply(smp,function(y){
	flnm <- paste0(netMHCpanDir,y,"_A/",y,"_NetChopcut2netMHCpan.xls")
	flnm2 <- paste0(netMHCpanDir,y,"_A/",y,"NetChopOut_pep.txt")
	netMHCpanRes <- read.table(file =flnm,sep="\t",quote ="",stringsAsFactors=F,fill=T,header=F)
	netChop_anno <- read.table(file =flnm2,sep=" ",quote ="",stringsAsFactors=F,fill=T,header=T)
	hla <- setdiff(unlist(netMHCpanRes[1,]),"")
	netMHCshort <- netMHCpanRes[-c(1:2),]
	binding_list <- lapply(seq(length(hla)),function(x){
		i=x
		idRank=3+(i-1)*6+4
		indexbind <- which(as.numeric(netMHCshort[,idRank]) < 2)
		if(length(indexbind)>0)
		{
			res_df <- data.frame(Peptide=netMHCshort[indexbind,2],
			core=netMHCshort[indexbind,3+(i-1)*6+1],icore=netMHCshort[indexbind,3+(i-1)*6+2],
			EL_score=netMHCshort[indexbind,3+(i-1)*6+3], EL_Rank=netMHCshort[indexbind,3+(i-1)*6+4],
			HLA=hla[x],stringsAsFactors=F)
			return(res_df)
		}
		
	})
	binding_df <- do.call(rbind,binding_list)
	binseq <- match(binding_df[,"Peptide"],netChop_anno[,"peptide"])
	binding_seq_df <- cbind(binding_df,netChop_anno[binseq,],stringsAsFactors=F)
	return(binding_seq_df)
})
names(anno_list) <- clinical[match(smp,clinical[,"clinical_label"]),"rawDataID"]
EGFR_anno_list <- lapply(anno_list,function(x) x[x[,"peptide"]=="KITDFGRAK",])
EGFR_anno_list <- EGFR_anno_list[sapply(EGFR_anno_list,nrow)>0]

#calculate the clone frrequencies for neoantigen-recognizing T cells
smp_list=names(EGFR_anno_list)
egfr_clone <- lapply(smp_list,function(smp){
	TCRAIres_predict_bind <- TCRAIinput_predict[TCRAIinput_predict[,"Sample"]==smp,]
	if(nrow(TCRAIres_predict_bind)>3)
	{
		res =sum(TCRAIres_predict_bind[TCRAIres_predict_bind[,"inputPredict"]>0.81595,"Clonality"])
		return(res)
	}
})
names(egfr_clone) <- smp_list
egfr_clone <- unlist(egfr_clone)
no_egfr_clone <- lapply(setdiff(names(anno_list),smp_list),function(smp){
	TCRAIres_predict_bind <- TCRAIinput_predict[TCRAIinput_predict[,"Sample"]==smp,]
	res =sum(TCRAIres_predict_bind[TCRAIres_predict_bind[,"inputPredict"]>0.81595,"Clonality"])
	if(nrow(TCRAIres_predict_bind)>3)
	{
		res =sum(TCRAIres_predict_bind[TCRAIres_predict_bind[,"inputPredict"]>0.81595,"Clonality"])
		return(res)
	}
})
names(no_egfr_clone) <- setdiff(names(anno_list),smp_list)
no_egfr_clone <- unlist(no_egfr_clone)

#plot the result
clone_df <- data.frame(clone_freq=c(egfr_clone,no_egfr_clone),
neoantigen=c(rep("EGFR L858R neoantigen",length(egfr_clone)),
rep("others",length(no_egfr_clone))),
stage=clinical[match(c(names(egfr_clone),names(no_egfr_clone)),clinical[,"rawDataID"]),"class"],stringsAsFactors=F)
p<-ggplot(clone_df)+geom_boxplot(aes(neoantigen,clone_freq),
col=brewer.pal(6, "Greys")[c(4,4)])+
geom_point(data=clone_df,aes(neoantigen,clone_freq),
col=colo[as.character(clone_df[,"stage"])],
position="jitter",pch=16,cex=1)+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
panel.grid.major=element_line(size=0.1,linetype =1,color="white"),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=90)) + ylab("frequency of predicted TCRs that bind to neoantigen from EGFR L858R")
p

#plot motif
immu.integrated <- readRDS(file = immuneFile)
setwd(TRUSTdir)
fl <- dir()
fl <- grep("_barcode_report.tsv",fl,value=T)
fl_list <- sapply(fl,function(x){
	fl_tmp <- read.table(file =x,sep="\t",stringsAsFactors=F, header=F)
	fl_tmp <- fl_tmp[(fl_tmp[,2]=="abT")&(fl_tmp[,3]!="*"),]
	fl_tmp_df <- data.frame(V1=fl_tmp[,1],
	V2=sapply(strsplit(fl_tmp[,3],","),function(y) y[1]),
	V3=sapply(strsplit(fl_tmp[,3],","),function(y) y[3]),
	V4=sapply(strsplit(fl_tmp[,3],","),function(y) y[6]),stringsAsFactors=F)
	fl_tmp_df[,2] <- sapply(strsplit(fl_tmp_df[,2],"*",fixed=T),function(y) y[1])
	fl_tmp_df[,3] <- sapply(strsplit(fl_tmp_df[,3],"*",fixed=T),function(y) y[1])
	smp <- sub("_barcode_report.tsv", "", x)
	print(smp)
	anno <- subset(immu.integrated,
	cells=rownames(immu.integrated@meta.data)[(as.character(immu.integrated@meta.data[,"orig.ident"])==smp)&(as.character(immu.integrated@meta.data[,"majorLabel"])=="CD8T")])
	smp_index <- strsplit(rownames(anno@meta.data)[1],"_")[[1]][2]
	fl_tmp_df[,1] <- paste(fl_tmp_df[,1],smp_index,sep="_")
	fl_tmp_df <- fl_tmp_df[fl_tmp_df[,1]%in%rownames(anno@meta.data),]

	if(nrow(fl_tmp_df)>0)
	{
		label <- paste(fl_tmp_df[,4],fl_tmp_df[,2],fl_tmp_df[,3],sep="_")
		fl_tmp_label_df <- cbind(fl_tmp_df,label)
		cell_state <- anno@meta.data[fl_tmp_df[,1],"minorLabel"]
		fl_tmp_label_df <- cbind(fl_tmp_label_df,cell_state)
		res_df <- data.frame(TRB_cdr3=fl_tmp_label_df[,"V4"],
		TRB_v_gene=fl_tmp_label_df[,"V2"],TRB_j_gene=fl_tmp_label_df[,"V3"],
		pmhc_code="EGFRpL858R",id="KITDFGRAK",Sample=smp,cell_state=cell_state,
		barcode=fl_tmp_label_df[,"V1"],
		stringsAsFactors=F)
		rownames(res_df) <- NULL
		return(res_df)
	}
})
fl_df <- do.call(rbind,fl_list)
rownames(fl_df) <- NULL
fl_df_label <- paste(fl_df[,1],fl_df[,2],fl_df[,3],sep="_")
predict_v1 <- TCRAIinput_predict[match(fl_df_label,inputlabel),"inputPredict"]
fl_df <- cbind(fl_df,predict_v1)
fl_df_pos <- fl_df[fl_df[,"predict_v1"]>0.81595,]
for(i in 13:15)
{
	t1 <- fl_df_pos_bind[nchar(fl_df_pos_bind[,1])==i,1]
	p1 = ggseqlogo(t1 )
	print(p1)
}

#############################
#Figure 6G (R version 4.0.3)#
#############################
#${clinicalFile} is the clinical information
#${resTCRAI} is the result of model on our datasets
#${immuneFile} is the seurat object of immune cell
#${inputTCRAI} is the input of model on our datasets
#${TRUSTdir} is the directory of TRUST result
#${smp} is sample ID

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(gplots)
library(pheatmap)

#Set working directory
setwd(workDir)

#Load clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical[clinical[,"class"]=="s25GGO","class"] <- "GGO25"
clinical[clinical[,"class"]=="s50GGO","class"] <- "GGO50"
clinical[clinical[,"class"]=="s75GGO","class"] <- "GGO75"
clinical[clinical[,"class"]=="s100GGO","class"] <- "GGO100"
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")

#Load the result of model
load(resTCRAI)

#Load the TCR information
immu.integrated <- readRDS(file = immuneFile)
EGFR_smp_list <- c("WZH","HKJ","ZQZ","XBC","ZL","LTJ","WHZ","YFQ","WLB","CXM","BYQ",
"LYQ","WXY","CSM")
fl_df_pos_bind[,"Sample"] <- factor(fl_df_pos_bind[,"Sample"],levels=EGFR_smp_list)
smp <- table(fl_df_pos_bind[,"Sample"])
stat <- immu.integrated@meta.data[,c("orig.ident","majorLabel")]
stat_CD8T <- stat[(as.character(stat[,"majorLabel"])=="CD8T")&(as.character(stat[,"orig.ident"])%in%names(smp)),]
fraction_list <- split(smp/table(stat_CD8T)[names(smp),],
f=factor(clinical[match(names(smp),clinical[,"rawDataID"]),"class"]))
TCRAIinput <- read.table(file = inputTCRAI,sep=",",quote = "",stringsAsFactors=F,header=T)
fraction_list2 <- split(smp/table(TCRAIinput[,"Sample"])[names(smp)],
f=factor(clinical[match(names(smp),clinical[,"rawDataID"]),"class"]))
fl <- paste0(TRUSTdir,names(smp),"_barcode_report.tsv")
fl_list <- lapply(fl,function(x){
	fl_tmp <- read.table(file =x,sep="\t",stringsAsFactors=F, header=F)
	fl_tmp <- fl_tmp[(fl_tmp[,2]=="abT")&(fl_tmp[,3]!="*"),]
	fl_tmp_df <- data.frame(V1=fl_tmp[,1],
	V2=sapply(strsplit(fl_tmp[,3],","),function(y) y[1]),
	V3=sapply(strsplit(fl_tmp[,3],","),function(y) y[3]),
	V4=sapply(strsplit(fl_tmp[,3],","),function(y) y[6]),stringsAsFactors=F)
	fl_tmp_df[,2] <- sapply(strsplit(fl_tmp_df[,2],"*",fixed=T),function(y) y[1])
	fl_tmp_df[,3] <- sapply(strsplit(fl_tmp_df[,3],"*",fixed=T),function(y) y[1])
	return(fl_tmp_df)
})
names(fl_list) <- names(smp)

#plot the fraction divided by TRUST
TCR_index <- sapply(seq(nrow(fl_df_pos_bind)),function(x){
	smpTmp <- fl_df_pos_bind[x,"Sample"]
	fl_tmp_df <- fl_list[smpTmp][[1]]
	labelVJ <- table(paste(fl_tmp_df[,2],fl_tmp_df[,3],sep="_"))
	NCDR3 <- table(fl_tmp_df[,4])
	res <- 1
	if(labelVJ[paste(fl_df_pos_bind[x,"TRB_v_gene"],fl_df_pos_bind[x,"TRB_j_gene"],sep="_")]>1){
		res <- 2
	}
	
	if(NCDR3[fl_df_pos_bind[x,"TRB_cdr3"]]>1){
		res <- 3
	}
	return(res)
})
fl_df_pos_bind <- cbind(fl_df_pos_bind,TCR_index)
smp <- table(fl_df_pos_bind[,c("Sample","TCR_index")])
stat <- immu.integrated@meta.data[,c("orig.ident","majorLabel")]
stat_CD8T <- stat[(as.character(stat[,"majorLabel"])=="CD8T")&(as.character(stat[,"orig.ident"])%in%rownames(smp)),]
fraction <- apply(smp,1,sum)/as.numeric(table(TCRAIinput[,"Sample"])[rownames(smp)])
clinical_label <- clinical[match(rownames(smp),clinical[,"rawDataID"]),"class"]
colo <- c(brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("pGGO","dGGO","GGO25","GGO50",
	"GGO75","GGO100","Solid1","Solid3","SolidN")
xindex <- rep(0,length(clinical_label))
xindex[clinical_label=="pGGO"] <- 1
xindex[clinical_label=="dGGO"] <- 2
xindex[clinical_label=="GGO25"] <- 3
xindex[clinical_label=="GGO50"] <- 4
xindex[clinical_label=="GGO75"] <- 5
xindex[clinical_label=="Solid1"] <- 6
xindex[clinical_label=="SolidN"] <- 7
plot(jitter(xindex),fraction,col=colo[clinical_label],pch=19,axes = F)
axis(2)
axis(1,labels=c("pGGO","dGGO","GGO25","GGO50",
	"GGO75","Solid1","SolidN"),at=1:7)
box()

#plot the fraction divided by all the CD8+ T cells
fraction2 <- apply(smp,1,sum)/as.numeric(table(stat_CD8T[,"orig.ident"])[rownames(smp)])
clinical_label2 <- clinical[match(rownames(smp),clinical[,"rawDataID"]),"class"]
xindex2 <- rep(0,length(clinical_label2))
xindex2[clinical_label2=="pGGO"] <- 1
xindex2[clinical_label2=="dGGO"] <- 2
xindex2[clinical_label2=="GGO25"] <- 3
xindex2[clinical_label2=="GGO50"] <- 4
xindex2[clinical_label2=="GGO75"] <- 5
xindex2[clinical_label2=="Solid1"] <- 6
xindex2[clinical_label2=="SolidN"] <- 7
plot(jitter(xindex2),fraction2,col=colo[clinical_label2],pch=19,axes = F)
axis(2)
axis(1,labels=c("pGGO","dGGO","GGO25","GGO50",
	"GGO75","Solid1","SolidN"),at=1:7)
box()

################################
#Figure S8B-C (R version 4.0.3)#
################################
#${clinicalFile} is the clinical information
#${ANNOVARdir} is the directory of ANNOVAR result
#${ggoFile} is the seurat object of all the cells
#${EcoFile} is the results of Ecotyper
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${cellStateDir} is the result directory of cell state

#${resTCRAI} is the result of model on our datasets

#${inputTCRAI} is the input of model on our datasets

#${smp} is sample ID

#${clinicalFile} is the file of clinical information of bulk RNaseq


#Load required packages
options(stringsAsFactors=F)
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(Seurat)

#Set working directory
setwd(workDir)

#Load the clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")
colo <- c(brewer.pal(8, "Blues")[c(2,3,4,5,7,8)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("pGGO","dGGO","GGO25","GGO50","GGO75","GGO100","Solid1","Solid3","SolidN")

#Load the result of mutation
anno_list <- lapply(smp,function(x){
	flnm <- paste0(ANNOVARdir,x,".hg38_multianno.txt")
	fl <- read.table(file=flnm,sep="\t",quote = "",stringsAsFactors=F,header=T)
	fl <- cbind(fl,rep(x,nrow(fl)))
	fl <- fl[fl[,6]%in%c(  "exonic","exonic;splicing"),]
	fl <- fl[!(fl[,9]%in%c("synonymous SNV",".","unknown")),]
	colnames(fl)[12] <- "Tumor_Sample"
	return(fl)
})
anno_df <- do.call(rbind,anno_list)
mutGene_sample <- unique(anno_df[,c("Gene.refGene","Tumor_Sample")])
mutGene_sample_stage <- clinical[match(mutGene_sample[,"Tumor_Sample"],clinical[,"clinical_label"]),c("class","rawDataID")]
mutGene_sample_stage_df <- cbind(mutGene_sample,mutGene_sample_stage,stringsAsFactors=F)
EGFR_sample <- unique(mutGene_sample_stage_df[(mutGene_sample_stage_df[,"Gene.refGene"]=="EGFR"),"rawDataID"])

#Load the ecotypr information
ggo.integrated <- readRDS(file=ggoFile)
ecotye <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)
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
clin_meta <- ecotype.integrate@meta.data
clini_class <- clinical[match(clin_meta[,"orig.ident"],clinical[,"rawDataID"]),"class"]
eco_class <- ecotye[match(clin_meta[,"Ecotype"],ecotye[,"ID"]),"Ecotype"]
clin_meta <- cbind(cbind(clin_meta,clini_class,stringsAsFactors=F),eco_class,stringsAsFactors=F)

#Caluculation the assocation between EGFR mutation and LME fraction
wes_sample <- unique(mutGene_sample_stage_df[,"rawDataID"])
E_fraction <- sapply(paste0("E",1:6),function(x) {
	res <- sapply(wes_sample,function(y){
		sel <- clin_meta[clin_meta[,"orig.ident"]==y,"eco_class"]
		fraction_tmp <- sum(sel==x)/length(sel)
		return(fraction_tmp)
	})
	return(res)
})
EGFR_label <- as.numeric(rownames(E_fraction)%in%EGFR_sample)
EGFR_label[EGFR_label==1] <- "EGFR mutation"
EGFR_label[EGFR_label==0] <- "EGFR WT"
E_fraction_list <- lapply(1:6,function(x)
{
	res <- split(E_fraction[,x],f=factor(EGFR_label))
	return(res)
})

#plot association between EGFR mutation and LME03
EcoEGFR_df <- data.frame(EGFR_mutation=factor(EGFR_label),
Eco=E_fraction[,3],stringsAsFactors=F)
colo=brewer.pal(3, "Set1")[c(2,1)]
names(colo) <- c("EGFR WT","EGFR mutation")
my_comparisons <- list( c("EGFR WT","EGFR mutation"))
p<-ggboxplot(EcoEGFR_df, x = "EGFR_mutation", y = "Eco", color="EGFR_mutation",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab("LME3 fraction")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#Caluculation the assocation between EGFR mutation and cell state
celltype_label <- sapply(strsplit(as.character(clin_meta[,"Ecotype"]),"_",fixed=T),function(x) x[1])
z=unique(celltype_label)[9]
clin_meta_sub <- clin_meta[celltype_label==z,]
E_fraction_cell <- sapply(unique(clin_meta_sub[,"Ecotype"]),function(x) {
	res <- sapply(wes_sample,function(y){
		sel <- clin_meta_sub[clin_meta_sub[,"orig.ident"]==y,"Ecotype"]
		fraction_tmp <- sum(sel==x)/length(sel)
		return(fraction_tmp)
	})
	return(res)
})
dim(E_fraction_cell)
E_fraction_cell_list <- lapply(unique(clin_meta_sub[,"Ecotype"]),function(x)
{
	res <- split(E_fraction_cell[,x],f=factor(EGFR_label))
	return(res)
})
i=3;wilcox.test(E_fraction_cell_list[[i]]["EGFR WT"][[1]],E_fraction_cell_list[[i]]["EGFR mutation"][[1]])

#plot association between EGFR mutation and cell state with significance
my_comparisons <- list( c("EGFR WT","EGFR mutation"))
EcoEGFR_df <- data.frame(EGFR_mutation=factor(EGFR_label),
Eco=E_fraction_cell[,i],stringsAsFactors=F)
EcoEGFR_df <- na.omit(EcoEGFR_df)
colo=brewer.pal(3, "Set1")[c(2,1)]
names(colo) <- c("EGFR WT","EGFR mutation")
p<-ggboxplot(EcoEGFR_df, x = "EGFR_mutation", y = "Eco", color="EGFR_mutation",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(color ="black",size=30,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=30),
		  axis.text.y=element_text(color ="black",size=30)) + ylab(unique(clin_meta_sub[,"Ecotype"])[i])+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p


