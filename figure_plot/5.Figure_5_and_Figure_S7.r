## Script for plot Figure 5 and Figure S7
##Yulan Deng, last updated 2023-4-20

#############################
#Figure 5A (R version 4.0.3)#
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

################################
#Figure S7B-C (R version 4.0.3)#
################################
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
