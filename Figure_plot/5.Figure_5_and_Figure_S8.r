## Script for plot Figure 5 and Figure S8
##Yulan Deng, last updated 2023-9-18

#################################
#Figure 5a,S8d (R version 4.0.3)#
#################################
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

#Figure S8d
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

##############################
#Figure S8c (R version 4.0.3)#
##############################
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