##Script for survival analysis of cancer cluster and LME 
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

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

