## Script for plot Figure 4 and Figure S5-7
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

###############################
#Figure 4a-b (R version 4.0.3)#
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

##Figure 4b
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
