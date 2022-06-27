##Script for survival analysis of cancer cluster and LME 
##Yulan Deng, last updated 2022-6-27

################################################
#survival analysis of cancer cluster c2 in TCGA#
################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data of TCGA
#${expresion} is the RNA expression matrix of TCGA, scaled to log2(TPM)
#${cancer_marker} is the marker genes of cancer clusters 

#Load required packages
library(survival)
library(forestplot)

#Set working directory
setwd(workDir)

#Load expresion matrix and clinical data of TCGA
load(clinicalFile)
load(expresion)

#Load marker genes of cancer clusters
cancer.markers <- readRDS(file = cancer_marker)

#selection of signature genes for cancer cluster c2
gene2 <- cancer.markers[(cancer.markers[,"p_val"]<0.05)&(cancer.markers[,"cluster"]=="2")&(cancer.markers[,"avg_logFC"]>log(1.5)),] 
gene2i <- gene2[gene2[,"gene"]%in%rownames(logexpmtTPMsubsmp),]
gene2i <- gene2i[order(gene2i[,"avg_logFC"],decreasing=T),]

#basic processing of TCGA clinical information
TCGAclinicallabel <- substr(TCGAclinical[,1],14,15)
TCGAclinicalsub <- TCGAclinical[TCGAclinicallabel=="01",]
TCGAclinicalDFS <- TCGAclinicalsub[,c("Sample.ID","Disease.Free.Status","Disease.Free..Months.",
"Diagnosis.Age","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Sex")]
TCGAclinicalDFS[,1] <- substr(TCGAclinicalDFS[,1],9,12)
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,3]),]
TCGAclinicalDFS[,2] <- as.numeric(sapply(strsplit(TCGAclinicalDFS[,2],":",fixed=T),function(x) x[1]))
TCGAclinicalDFS <- TCGAclinicalDFS[TCGAclinicalDFS[,3]!=0,]
colnames(TCGAclinicalDFS)[c(4,5)] <- c("Age","Stage")
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,"Stage"]),]
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IB","Stage IA"),"Stage"] <- "Stage I"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IIB","Stage IIA"),"Stage"] <- "Stage II"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IIIA","Stage IIIB"),"Stage"] <- "Stage III"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage I","Stage II"),"Stage"] <- "Stage I+II"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage III","Stage IV"),"Stage"] <- "Stage III+IV"
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,"Age"]),]
intersmpDFS <- intersect(colnames(logexpmtTPMsubsmp),TCGAclinicalDFS[,1])
TCGAclinicalDFSsub <- TCGAclinicalDFS[match(intersmpDFS,TCGAclinicalDFS[,1]),]
logexpmtTPMsubsmpDFS <- logexpmtTPMsubsmp[,intersmpDFS]
yDFS <-Surv(TCGAclinicalDFSsub[,3],TCGAclinicalDFSsub[,2])

#calculation of c2 signature score
mtZ <- t(apply(logexpmtTPMsubsmpDFS[gene2i[1:40,"gene"],],1,function(y) {
		val <- unlist(y)
		res1 <- (val-mean(val))/sd(val)
		return(res1)
	}))
	score <- apply(mtZ,2,mean)
cluster <- rep(0,length(score))
cluster[score > median(score)] <- 1

#multi-variate cox analysis
resM <-coxph(yDFS~score+TCGAclinicalDFSsub[,"Stage"]+TCGAclinicalDFSsub[,"Age"]+TCGAclinicalDFSsub[,"Sex"])

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("c2_signature","stage","age","sex")
forestplot(
c("c2_signature","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  forest_data[,"high"],
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "TCCA DFS"
)

#save data for KM plot
cliniT <- data.frame(recurrenceTime=TCGAclinicalDFSsub[,3],
label=TCGAclinicalDFSsub[,2],Clusters=cluster,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(yDFS~cluster)
P_Value<-pchisq(Diff_Survival$chisq,length(table(cluster))-1,lower.tail=F)
save(P_Value,Sur_Curve,cliniT,file="c2signature_TCGA_DFS.40.RData")

####################################################
#survival analysis of cancer cluster c2 in GSE31210#
####################################################
#${workDir} is the working directory
#${expresionTCGA} is the RNA expression matrix of TCGA, scaled to log2(TPM)
#${cancer_marker} is the marker genes of cancer clusters 
#${GSE31210File} is the matrix file downloaded from GEO
#${GSE31210ProbeFile} is the probe information downloaded from GEO
#${GSE31210FileHeader} is the header of matrix file downloaded from GEO

#Load required packages
library(survival)
library(forestplot)

#Set working directory
setwd(workDir)

#selection of signature genes for cancer cluster c2
load(expresionTCGA)
cancer.markers <- readRDS(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/20210831/cancer.markers.0.25.rds")
gene2 <- cancer.markers[(cancer.markers[,"p_val"]<0.05)&(cancer.markers[,"cluster"]=="2")&(cancer.markers[,"avg_logFC"]>log(1.5)),] #146
gene2i <- gene2[gene2[,"gene"]%in%rownames(logexpmtTPMsubsmp),]
gene2i <- gene2i[order(gene2i[,"avg_logFC"],decreasing=T),]
gene2signature <- gene2i[1:40,"gene"]

#the basic processing for expression matrix and clinical information of GSE31210
GSE31210_mt <- read.table(file=GSE31210File,sep="\t",comment.char = "!",stringsAsFactors=F,header=T,row.names=1)
GSE31210_mtn <- apply(GSE31210_mt,2,as.numeric)
rownames(GSE31210_mtn) <- rownames(GSE31210_mt)
GSE31210_probe <- read.table(file=GSE31210ProbeFile,sep="\t",comment.char = "#",stringsAsFactors=F,header=T,fill=T)
GSE31210_clinical <- t(read.table(file=GSE31210FileHeader,sep="\t",comment.char = "#",quote = "",stringsAsFactors=F,header=F))
colnames(GSE31210_clinical) <- c("smpID","tissue","age","gender","smoking","bi",
"pathological_stage","pstage","gene_alteration","myc","myc_copy",
"relapse","days before relapse/censor","months before relapse/censor",
"death","days before death/censor","months before death/censor")
intersmp <- intersect(GSE31210_clinical[,"smpID"],colnames(GSE31210_mtn))#226
GSE31210_mtnsub <- GSE31210_mtn[,intersmp]
GSE31210_clinicalsub <- GSE31210_clinical[match(intersmp,GSE31210_clinical[,1]),] 
inter <- intersect(GSE31210_probe[,1],rownames(GSE31210_mtnsub))#54675
gene <- unique(GSE31210_probe[,2]) 
gene <- setdiff(gene,"")
gene_probe <- sapply(gene,function(x){
	probe_id <- GSE31210_probe[GSE31210_probe[,2]==x,1]
	print(x)
	if(length(probe_id)==1)
	{
		return(probe_id)
	}else{
		Mean <- apply(GSE31210_mtnsub[probe_id,],1,mean)
		return(probe_id[which(Mean==max(Mean))[1]])
	}
})
gene_probe_df <- data.frame(probe=gene_probe,gene=gene,stringsAsFactors=F)
inter_probe <- intersect(gene_probe,rownames(GSE31210_mtnsub))
GSE31210_mtnsub <- GSE31210_mtnsub[inter_probe,]
rownames(GSE31210_mtnsub) <- gene_probe_df[match(inter_probe,gene_probe_df[,1]),2]
GSE31210_mtnsub <- GSE31210_mtnsub[apply(GSE31210_mtnsub,1,function(x) sum(is.na(x)))==0,]
GSE31210_clinicalsub <- as.data.frame(GSE31210_clinicalsub)
GSE31210_clinicalsub[,"death"] <- as.numeric(as.character(GSE31210_clinicalsub[,"death"])=="death: dead")
GSE31210_clinicalsub[,"relapse"] <- as.numeric(as.character(GSE31210_clinicalsub[,"relapse"])=="relapse: relapsed")
GSE31210_clinicalsub[,"months before death/censor"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsub[,"months before death/censor"]),": "),function(x) x[2]))
GSE31210_clinicalsub[,"months before relapse/censor"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsub[,"months before relapse/censor"]),": "),function(x) x[2]))
GSE31210_clinicalsub[,"gender"] <- sapply(strsplit(as.character(GSE31210_clinicalsub[,"gender"]),": "),function(x) x[2])
GSE31210_clinicalsub[,"age"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsub[,"age"]),": "),function(x) x[2]))
GSE31210_clinicalsub[,"pathological_stage"] <- sapply(strsplit(as.character(GSE31210_clinicalsub[,"pathological_stage"]),": "),function(x) x[2])
GSE31210_clinicalsub[GSE31210_clinicalsub[,"pathological_stage"]%in%c("IA","IB"),"pathological_stage"] <- "I"
GSE31210_clinicalsubDFS <- GSE31210_clinicalsub[!is.na(GSE31210_clinicalsub[,"months before relapse/censor"]),]
yDFS <-Surv(GSE31210_clinicalsubDFS[,"months before relapse/censor"],GSE31210_clinicalsubDFS[,"relapse"])
GSE31210_mtnsubDFS <- GSE31210_mtnsub[,GSE31210_clinicalsubDFS[,"smpID"]]

#Calculation of c2 signature
mtZ <- t(apply(GSE31210_mtnsubDFS[gene2i[1:40,"gene"],],1,function(y) {
		val <- unlist(y)
		res1 <- (val-mean(val))/sd(val)
		return(res1)
	}))
	score <- apply(mtZ,2,mean)
cluster <- rep(0,length(score))
cluster[score > median(score)] <- 1

#multi-variate cox analysis
resM <-coxph(yDFS~score+GSE31210_clinicalsubDFS[,"pathological_stage"]+GSE31210_clinicalsubDFS[,"age"]+GSE31210_clinicalsubDFS[,"gender"])

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("c2_signature","stage","age","sex")
forestplot(
c("c2_signature","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  forest_data[,"high"],
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "GSE31210 DFS"
)

#save data for KM plot
cliniT <- data.frame(recurrenceTime=GSE31210_clinicalsubDFS[,"months before relapse/censor"],
label=GSE31210_clinicalsubDFS[,"relapse"],Clusters=cluster,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(yDFS~cluster)
P_Value<-pchisq(Diff_Survival$chisq,length(table(cluster))-1,lower.tail=F)
save(P_Value,Sur_Curve,cliniT,file="c2signature_GSE31210_DFS.40.RData")

#####################################################
#survival analysis of cancer cluster c2 in GSE140343#
#####################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data of GSE140343
#${expresionTCGA} is the RNA expression matrix of TCGA, scaled to log2(TPM)
#${cancer_marker} is the marker genes of cancer clusters 
#${GSE140343_Expression_matrix} is the expression matrix file in GSE140343

#Load required packages
library(survival)
library(forestplot)

#Set working directory
setwd(workDir)

#selection of c2 signature
load(clinicalFile)
load(expresionTCGA)
cancer.markers <- readRDS(file = cancer_marker)
gene2 <- cancer.markers[(cancer.markers[,"p_val"]<0.05)&(cancer.markers[,"cluster"]=="2")&(cancer.markers[,"avg_logFC"]>log(1.5)),]
gene2i <- gene2[gene2[,"gene"]%in%rownames(logexpmtTPMsubsmp),]
gene2i <- gene2i[order(gene2i[,"avg_logFC"],decreasing=T),]
gene2signature <- gene2i[1:40,"gene"]

#the basic processing of clinical information and expression matrix of GSE140343
GSE140343_mt <- read.table(file=GSE140343,
sep="\t",header=T,row.names=1)
intersmp <- intersect(clinical_GSE140343[,1],colnames(GSE140343_mt))
GSE140343_mtnsub <- GSE140343_mt[,intersmp]
GSE140343_clinicalsub <- clinical_GSE140343[match(intersmp,clinical_GSE140343[,1]),] 
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IA1","IA2",
"IA3","IB"),"Tumor.stage.A.B."] <- "I"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IIA","IIB"),"Tumor.stage.A.B."] <- "II"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IIIA","IIIB"),"Tumor.stage.A.B."] <- "III"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IVA"),"Tumor.stage.A.B."] <- "IV"
GSE140343_clinicalsub <- GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("I","II"),]
GSE140343_clinicalsubDFS <- GSE140343_clinicalsub[!is.na(GSE140343_clinicalsub[,"DFS.month"]),]
yDFS <-Surv(GSE140343_clinicalsubDFS[,"DFS.month"],GSE140343_clinicalsubDFS[,"DFS"])
GSE140343_mtnsubDFS <- GSE140343_mtnsub[,GSE140343_clinicalsubDFS[,1]]

#Calculation of c2 signature
mtZ <- t(apply(GSE140343_mtnsubDFS[gene2i[1:40,"gene"],],1,function(y) {
		val <- unlist(y)
		res1 <- (val-mean(val))/sd(val)
		return(res1)
	}))
	score <- apply(mtZ,2,mean)

cluster <- rep(0,length(score))
cluster[score >= median(score)] <- 1

#multi-variate cox analysis
resM <-coxph(yDFS~score+GSE140343_clinicalsubDFS[,"Tumor.stage.A.B."]+GSE140343_clinicalsubDFS[,"Age"]+as.character(GSE140343_clinicalsubDFS[,"Gender..1.Male.2.Female."]))

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("c2_signature","stage","age","sex")
forestplot(
c("c2_signature","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  sapply(forest_data[,"high"],function(x) min(x,20)),
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "GSE140343 DFS"
)

#save data for KM plot
cliniT <- data.frame(recurrenceTime=GSE140343_clinicalsubDFS[,"DFS.month"],
label=GSE140343_clinicalsubDFS[,"DFS"],Clusters=cluster,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(yDFS~cluster)
P_Value<-pchisq(Diff_Survival$chisq,length(table(cluster))-1,lower.tail=F)
save(P_Value,Sur_Curve,cliniT,file="c2signature_GSE140343_DFS.40.RData")

###################################################
#survival analysis of cancer cluster LME01 in TCGA#
###################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data for TCGA
#${ecotypeFile} is the fraction of LME

#Load required packages
library(survival)
library(forestplot)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the clinical information for TCGA
load(clinicalFile)

#Load the fraction of LME
fraction <- read.table(file=ecotypeFile,sep="\t",stringsAsFactors=F,header=T)

#Basic processing of TCGA clinical information and LME fraction
TCGAclinicallabel <- substr(TCGAclinical[,1],14,15)
TCGAclinicalsub <- TCGAclinical[TCGAclinicallabel=="01",]
TCGAclinicalsub[,1] <- substr(TCGAclinicalsub[,1],9,12)
labelT <- substr(colnames(fraction),14,15)
fractionsubsmp <- fraction[labelT=="01",]
colnames(fractionsubsmp) <- substr(colnames(fractionsubsmp),9,12)
TCGAclinicalDFS <- TCGAclinicalsub[,c("Sample.ID","Disease.Free.Status","Disease.Free..Months.",
"Diagnosis.Age","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Sex")]
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,3]),]
TCGAclinicalDFS[,2] <- as.numeric(sapply(strsplit(TCGAclinicalDFS[,2],":",fixed=T),function(x) x[1]))
TCGAclinicalDFS <- TCGAclinicalDFS[TCGAclinicalDFS[,3]!=0,]
colnames(TCGAclinicalDFS)[c(4,5)] <- c("Age","Stage")
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,"Stage"]),]
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IB","Stage IA"),"Stage"] <- "Stage I"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IIB","Stage IIA"),"Stage"] <- "Stage II"
TCGAclinicalDFS[TCGAclinicalDFS[,"Stage"]%in%c("Stage IIIA","Stage IIIB"),"Stage"] <- "Stage III"
TCGAclinicalDFS <- TCGAclinicalDFS[!is.na(TCGAclinicalDFS[,"Age"]),]
intersmp <- intersect(intersect(colnames(fractionsubsmp),TCGAclinicalDFS[,1]),
TCGAclinicalsub[TCGAclinicalsub[,"Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"]%in%c("Stage I",
"Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"),1])
TCGAclinicalDFSsubsmp <- TCGAclinicalDFS[match(intersmp,TCGAclinicalDFS[,1]),]
fractionsubsmpDFS <- fractionsubsmp[,intersmp]
yDFS <-Surv(TCGAclinicalDFSsubsmp[,3],TCGAclinicalDFSsubsmp[,2])

#multi-variate cox analysis
resM <-coxph(yDFS~(as.numeric(unlist(fractionsubsmpDFS[1,])> median(as.numeric(unlist(fractionsubsmpDFS[1,])))))+TCGAclinicalDFSsubsmp[,"Stage"]+as.numeric(TCGAclinicalDFSsubsmp[,"Age"]>60 )+TCGAclinicalDFSsubsmp[,"Sex"])

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("E1","stage","age","sex")
forestplot(
c("E1","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  forest_data[,"high"],
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "TCCA DFS stage I/II"
)

#######################################################
#survival analysis of cancer cluster LME01 in GSE31210#
#######################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data for TCGA
#${ecotypeFile} is the fraction of LME
#${GSE31210File} is the matrix file downloaded from GEO
#${GSE31210ProbeFile} is the probe information downloaded from GEO
#${GSE31210FileHeader} is the header of matrix file downloaded from GEO

#Load required packages
library(survival)
library(forestplot)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the fraction of LME
fraction <- read.table(file=ecotypeFile,sep="\t",stringsAsFactors=F,header=T)

#basic processing for clinical information of GSE31210
GSE31210_clinical <- t(read.table(file=GSE31210FileHeader,sep="\t",comment.char = "#",quote = "",stringsAsFactors=F,header=F))
colnames(GSE31210_clinical) <- c("smpID","tissue","age","gender","smoking","bi",
"pathological_stage","pstage","gene_alteration","myc","myc_copy",
"relapse","days before relapse/censor","months before relapse/censor",
"death","days before death/censor","months before death/censor")
intersmp <- intersect(GSE31210_clinical[,"smpID"],colnames(fraction))#226
GSE31210_mtnsub <- fraction[,intersmp]
GSE31210_clinicalsub <- GSE31210_clinical[match(intersmp,GSE31210_clinical[,1]),] 
GSE31210_clinicalsub <- as.data.frame(GSE31210_clinicalsub)
GSE31210_clinicalsub[,"death"] <- as.numeric(as.character(GSE31210_clinicalsub[,"death"])=="death: dead")
GSE31210_clinicalsub[,"relapse"] <- as.numeric(as.character(GSE31210_clinicalsub[,"relapse"])=="relapse: relapsed")
GSE31210_clinicalsub[,"months before death/censor"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsub[,"months before death/censor"]),": "),function(x) x[2]))
GSE31210_clinicalsub[,"months before relapse/censor"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsub[,"months before relapse/censor"]),": "),function(x) x[2]))
GSE31210_clinicalsubDFS <- GSE31210_clinicalsub[(!is.na(GSE31210_clinicalsub[,"months before relapse/censor"])),]
yDFS <-Surv(GSE31210_clinicalsubDFS[,"months before relapse/censor"],GSE31210_clinicalsubDFS[,"relapse"])
GSE31210_mtnsubDFS <- GSE31210_mtnsub[,GSE31210_clinicalsubDFS[,"smpID"]]
GSE31210_clinicalsubDFS[,"pathological_stage"] <- as.character(GSE31210_clinicalsubDFS[,"pathological_stage"] )
GSE31210_clinicalsubDFS[GSE31210_clinicalsubDFS[,"pathological_stage"] %in%c("pathological stage: IA","pathological stage: IB"),
"pathological_stage"] = "pathological stage: I"
GSE31210_clinicalsubDFS[,"age"] <- as.numeric(sapply(strsplit(as.character(GSE31210_clinicalsubDFS[,"age"]),": "),function(x) x[2]))

#multi-variate cox analysis
resM <-coxph(yDFS~(as.numeric(unlist(GSE31210_mtnsubDFS[1,])> median(as.numeric(unlist(GSE31210_mtnsubDFS[1,])))))+GSE31210_clinicalsubDFS[,"pathological_stage"]+as.numeric(GSE31210_clinicalsubDFS[,"age"]>60 )+GSE31210_clinicalsubDFS[,"gender"])

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("E1","stage","age","sex")
forestplot(
c("E1","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  forest_data[,"high"],
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "GSE31210 DFS stage I/II"
)

########################################################
#survival analysis of cancer cluster LME01 in GSE140343#
########################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data of GSE140343
#${ecotypeFile} is the fraction of LME

#Load required packages
library(survival)
library(forestplot)

#Set working directory
setwd(workDir)

#Load the fraction of LME
fraction <- read.table(file=ecotypeFile,sep="\t",stringsAsFactors=F,header=T)

#Load the clinical information of GSE140343
load(clinicalFile)

#Basci processing for clinical information of GSE140343
intersmp <- intersect(clinical_GSE140343[,1],colnames(fraction))
fraction <- fraction[,intersmp]
GSE140343_clinicalsub <- clinical_GSE140343[match(intersmp,clinical_GSE140343[,1]),] 
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IA1","IA2",
"IA3","IB"),"Tumor.stage.A.B."] <- "I"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IIA","IIB"),"Tumor.stage.A.B."] <- "II"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IIIA","IIIB"),"Tumor.stage.A.B."] <- "III"
GSE140343_clinicalsub[GSE140343_clinicalsub[,"Tumor.stage.A.B."]%in%c("IVA"),"Tumor.stage.A.B."] <- "IV"
GSE140343_clinicalsubDFS <- GSE140343_clinicalsub[!is.na(GSE140343_clinicalsub[,"DFS.month"]),]
yDFS <-Surv(GSE140343_clinicalsubDFS[,"DFS.month"],GSE140343_clinicalsubDFS[,"DFS"])
Clusters <- rep(1,nrow(GSE140343_clinicalsubDFS))
Clusters[unlist(fraction[1,GSE140343_clinicalsubDFS[,1]])>=quantile(unlist(fraction[1,GSE140343_clinicalsubDFS[,1]]),probs=0.8)] <- 2

#multi-variate cox analysis
resM <-coxph(yDFS~Clusters+GSE140343_clinicalsubDFS[,"Tumor.stage.A.B."]+GSE140343_clinicalsubDFS[,"Age"]+as.character(GSE140343_clinicalsubDFS[,"Gender..1.Male.2.Female."]))

#forest plot
forest_data <- summary(resM)$"conf.int"[,c(1,3,4)]
colnames(forest_data) <- c("coef", "low" ,"high")
rownames(forest_data) <- c("LME1","stage","age","sex")
forestplot(
c("LME1","stage","age","sex"),
  forest_data[,"coef"],
  forest_data[,"low"],
  forest_data[,"high"],
  zero = 1,
  cex  = 2,
  lineheight = "auto",
  xlab = "GSE140343 DFS"
)

#save data for KM plot
cliniT <- data.frame(recurrenceTime=GSE140343_clinicalsubDFS[,"DFS.month"],
label=GSE140343_clinicalsubDFS[,"DFS"],Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(yDFS~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)
save(P_Value,Sur_Curve,cliniT,file="E1_signature_GSE140343_DFS.RData")

########################################################
#survival analysis of cancer cluster LME04 in GSE135222#
########################################################
#${workDir} is the working directory
#${clinicalFile} is the clinical data of GSE135222
#${ecotypeFile} is the fraction of LME
#${sampleID} is the file of matching sample ID
#${SraRunTable} is the table downloaded from SRA, with SRA sample ID.
#${GEOtable} is the table downloaded from GEO, with GEO sample ID.

#Load required packages
library(survival)
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the clinical information of GSE135222
load(clinicalFile)

#Load the ecotype fraction
fraction <- read.table(file=ecotypeFile,sep="\t",stringsAsFactors=F,header=T)

#match sample ID
sample_match <- read.table(file=sampleID,sep="\t",stringsAsFactors=F,header=F)
smpId1 <- read.table(file=SraRunTable,sep=",",quote =,stringsAsFactors=F,header=T)
clini1 <- read.table(file=GEOtable,stringsAsFactors=F,header=T) 
smpId1_clin <- cbind(smpId1[,c("Run","Sample.Name")],
clini1[match(smpId1[,"GEO_Accession..exp."],clini1[,"Sample_id"]),"Response"])
colnames(smpId1_clin) = c("SRR","GSM","response")

#cox analysis
fraction <- fraction[,smpId1_clin[match(sample_match[,1],smpId1_clin[,2]),1]]
ypfs <-Surv(Jung_clinic[sample_match[,2],"PFS"],Jung_clinic[sample_match[,2],"PD_Event.1_Censoring.0"])
coxph(ypfs~as.numeric(unlist(fraction[4,])>=quantile(unlist(fraction[4,]),probs=0.75)))

#save data for KM plot
Clusters <- rep(2,ncol(fraction))
Clusters[unlist(fraction[4,])>=quantile(unlist(fraction[4,]),probs=0.75)] <- 1
cliniT <- data.frame(recurrenceTime=Jung_clinic[sample_match[,2],"PFS"],
label=Jung_clinic[sample_match[,2],"PD_Event.1_Censoring.0"],Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(ypfs~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)
save(P_Value,Sur_Curve,cliniT,file="E4_signature_GSE135222.RData")