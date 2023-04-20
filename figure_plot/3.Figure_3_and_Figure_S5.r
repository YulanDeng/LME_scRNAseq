## Script for plot Figure 3 and Figure S5
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

#############################
#Figure 3A (R version 4.0.3)#
#############################
#${pythonFile} is the executable file for python
#${workDir} is the working directory
#${clinicalFile} is the file of clinical information
#${rank_dataFile} is the parameter for the appropriate number of cell states
#${ggoFile} is the seurat object of all the cell
#${cellStateDir} is the result directory of cell state
#${sliceFile} is the slice file for circos
#${ecoFile} is the result of Ecotyper
#${circos.conf.file} is the configuration file for circos
#${figureName} is the name for circos output figure

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(cluster)
library("reticulate")
use_python(pythonFile)

#Set working directory
setwd(workDir)

#Load input files
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NoAden"
clinical[clinical[,"class"]%in%c("pGGO","HiDenGGO","s25GGO","s50GGO","s75GGO","s100GGO"),"class"] <- "GGO"
clinical[clinical[,"class"]%in%c("Solid1","Solid3","SolidN"),"class"] <- "Solid"
rank_data <- read.table(file=rank_dataFile,sep="\t",stringsAsFactors=F,header=T)

#Exclude epithelium
rank_data <- rank_data[-6,]

ggo.integrated <- readRDS(file=ggoFile)
ggo.integrated$"stage" <- factor(clinical[match(as.character(ggo.integrated@meta.data[,"orig.ident"]),as.character(clinical[,"rawDataID"])),"class"],levels=c("Normal","GGO","Solid","NoAden"))

#pearson residue
#significance KW test
res_list <-  lapply(seq(nrow(rank_data)),function(x){

	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	ecotyper_df=read.table(file =paste0(cellStateDir,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	label <- rep("unassigned",nrow(ggo.integrated@meta.data))
	names(label) = rownames(ggo.integrated@meta.data)
	NcellState <- unique(ecotyper_df[,"State"])
	for(kk in NcellState)
	{
		label[intersect(ecotyper_df[ecotyper_df[,2]==kk,1],names(label))] = kk
	}
	ggo.integrated$ecoRecovery = factor(label,levels=c(NcellState,"unassigned"))

	stage_list <- split(ggo.integrated@meta.data[,c("ecoRecovery","orig.ident")],
	f=ggo.integrated@meta.data[,"stage"]) 

	stage_fraction_list <- lapply(stage_list,function(xx) {
		mt_tmp <- table(xx)
		mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"unassigned"),]
		fraction_tmp <- apply(mt_tmp,2,function(y) y/sum(y))
		fraction_tmp[is.na(fraction_tmp)] <- 0
		return(fraction_tmp)
	})
	names(stage_fraction_list) <- names(stage_list)
	res <- sapply(stage_fraction_list,function(yy) apply(yy,1,mean))
	res <- res[,c("Normal","GGO","Solid")]
	res <- t(apply(res,1,function(xx) (xx-min(xx))/(max(xx)-min(xx))))
	
	stage_fraction_list_p <- sapply(seq(length(NcellState)),function(y){
		resy <- lapply(stage_fraction_list[c("Normal","GGO","Solid")],function(z) z[y,])
		res_df <- data.frame(fraction=unlist(resy),
		stage=rep(c("Normal","GGO","Solid"),sapply(stage_fraction_list[c("Normal","GGO","Solid")],ncol)),stringsAsFactors=F)
		stat <- kruskal.test(fraction~stage, data = res_df)
		return(stat$"p.value")
	})
	
	res <- sapply(stage_fraction_list,function(yy) apply(yy,1,mean))
	res <- res[,c("Normal","GGO","Solid")]
	res <- res - res[,"Normal"]
	res <- t(apply(res,1,function(xx) xx/max(abs(xx))))
	
	##fraction
	n_cell_state <- sapply(stage_list[c("Normal","GGO","Solid")],function(xx) {
		mt_tmp <- table(xx)
		mt_tmp <- mt_tmp[setdiff(rownames(mt_tmp),"unassigned"),]
		res <- apply(mt_tmp,1,sum)
		return(res)
	})
	n_cell_state_fraction <- n_cell_state/sum(n_cell_state)
	res_list_tmp <- list()
	res_list_tmp[[1]] <- stage_fraction_list_p
	res_list_tmp[[2]] <- res
	res_list_tmp[[3]] <- n_cell_state_fraction
	return(res_list_tmp)	
})

slices_res <- read.table(file=sliceFile,sep=" ",stringsAsFactors=F,header=F)

#genes.txt,genes3-5.txt
gene_colo_fun <- function(x)
{
	fraction_v_cut <- cut(x[,4],breaks=seq(from=-1.1,to=1,length.out=10),
	labels=c(paste0("blue",4:1),"white",paste0("red",1:4)))
	return(as.character(fraction_v_cut))
}

#normal
normal_v <- unlist(sapply(sapply(res_list,function(y) y[2]),function(x) x[,1]))
normal_hm <- cbind(slices_res[,1:3],normal_v)
normal_hm <- cbind(normal_hm,paste0("color=",gene_colo_fun(x=normal_hm)))
write.table(normal_hm, file="normalHm.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

#G
G_v <- unlist(sapply(sapply(res_list,function(y) y[2]),function(x) x[,2]))
G_hm <- cbind(slices_res[,1:3],eG_v)
G_hm <- cbind(G_hm,paste0("color=",gene_colo_fun(x=G_hm)))
write.table(G_hm, file="GHm.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

#S
S_v <- unlist(sapply(sapply(res_list,function(y) y[2]),function(x) x[,3]))
S_hm <- cbind(slices_res[,1:3],S_v)
S_hm <- cbind(S_hm,paste0("color=",gene_colo_fun(x=S_hm)))
write.table(S_hm, file="SHm.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

gene_size_fun <- function(yy)
{
	size_v_cut <- cut(yy,breaks=c(0,0.01,0.05,0.1,0.15,0.2,0.25,1),
	labels=c(1:7)*10)
	return(as.character(size_v_cut))
}

#solid
gene1_v <- paste0("color=",gene_colo_fun(x=S_hm),",glyph_size=",
gene_size_fun(yy=unlist(lapply(sapply(res_list,function(y) y[3]),function(x) x[,"Solid"]))))
gene1 <- cbind(slices_res[,1:3],cbind(rep(0,nrow(slices_res)),gene1_v))
write.table(gene1, file="genes.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

#GGO
gene3_v <- paste0("color=",gene_colo_fun(x=G_hm),",glyph_size=",
gene_size_fun(yy=unlist(lapply(sapply(res_list,function(y) y[3]),function(x) x[,"GGO"]))))
gene3 <- cbind(slices_res[,1:3],cbind(rep(0,nrow(slices_res)),gene3_v))
write.table(gene3, file="genes3.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

#normal
gene4_v <- paste0("color=",gene_colo_fun(x=normal_hm),",glyph_size=",
gene_size_fun(yy=unlist(lapply(sapply(res_list,function(y) y[3]),function(x) x[,"Normal"]))))
gene4 <- cbind(slices_res[,1:3],cbind(rep(0,nrow(slices_res)),gene4_v))
write.table(gene4, file="genes4.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )

#genes2.txt
sig <- unlist(lapply(res_list,function(y) y[1]))
fdr <- p.adjust(sig,method="fdr")
gene2 <- cbind(slices_res[,1:3],fdr)
gene2 <- gene2[gene2[,4]<0.1,]
fdr_str <- rep(".",length(gene2[,4]))
fdr_str[gene2[,4]<0.05] <- "*"
gene2_s <- cbind(gene2[,1:3],fdr_str)
write.table(gene2_s, file="genes2.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )
}

rownames(slices_res) <- paste0(slices_res[,1],"_S0",slices_res[,3]/100)

ecotype <- read.table(file=ecoFile,sep="\t",stringsAsFactors=F,header=T)

ecotype <- ecotype[ecotype[,1]!="Epithelial.cells",]

colo=c("55,126,184","77,175,74","228,26,28","255,127,0","152,78,163","166,86,40")
eco <- unique(ecotype[,"Ecotype"])
ribbon_list <- lapply(seq(length(eco)),function(x){
	print(x)
	cell_type=ecotype[ecotype[,"Ecotype"]==eco[x],"ID"]
	cellCombine <- t(combn(cell_type,2))
	res_df <- cbind(slices_res[cellCombine[,1],1:3],slices_res[cellCombine[,2],1:3])
	coloD <- rep(paste0("color=(",colo[x],",0.5),twist=0"),nrow(res_df))
	res_df <- cbind(res_df,coloD)
	res_df[,2] <- res_df[,2]+25
	res_df[,5] <- res_df[,5]+25
	res_df[,3] <- res_df[,3]-25
	res_df[,6] <- res_df[,6]-25
	return(res_df)
})
ribbon_df <- do.call(rbind,ribbon_list)
write.table(ribbon_df, file="ribbon.twist.txt", sep=" " , quote=FALSE, row.names=F, col.names=F )
###
/usr/bin/perl circos -conf ${circos.conf.file} \
-outputdir ${workDir} -outputfile ${figureName}.svg -svg


##############################
#Figure S5D (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${EcoMoFile} is the Ecotyper result of monocyte and macrophage
#${EcoBFile} is the Ecotyper result of B cells
#${EcoCD4Tfile} is the Ecotyper result of CD4+ T cells
#${EcoCD8Tfile} is the Ecotyper result of CD8+ T cells
#${clinicalFile} is the file of clinical information 

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

#Set working directory
setwd(workDir)

fraction_CD8T <- read.table(file=EcoCD8Tfile,sep="\t",stringsAsFactors=F,header=T)
fraction_CD8T <- fraction_CD8T[,!(colnames(fraction_CD8T)%in%c("X10C","X10F","X16C","X16F"))]
fraction_B <- read.table(file=EcoBFile,sep="\t",stringsAsFactors=F,header=T)
fraction_B <- fraction_B[,!(colnames(fraction_B)%in%c("X10C","X10F","X16C","X16F"))]
fraction_Macro <- read.table(file=EcoMoFile,sep="\t",stringsAsFactors=F,header=T)
fraction_Macro <- fraction_Macro[,!(colnames(fraction_Macro)%in%c("X10C","X10F","X16C","X16F"))]
fraction_CD4T <- read.table(file=EcoCD4Tfile,sep="\t",stringsAsFactors=F,header=T)
fraction_CD4T <- fraction_CD4T[,!(colnames(fraction_CD4T)%in%c("X10C","X10F","X16C","X16F"))]
clinical_info <- read.table(file=clinicalFile,
sep="\t",stringsAsFactors=F,header=T)
clinical_info[,1] <- paste0("X",clinical_info[,1],"C")
clinical_info[clinical_info[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical_info[clinical_info[,"class"]=="s25GGO","class"] <- "GGO25"
clinical_info[clinical_info[,"class"]=="s50GGO","class"] <- "GGO50"
clinical_info[clinical_info[,"class"]=="s75GGO","class"] <- "GGO75"
clinical_label <- rep("Normal",ncol(fraction_CD8T))
clinical_label[clinical_label=="HiDenGGO"] <- "dGGO"
clinical_label[clinical_label=="s25GGO"] <- "GGO25"
clinical_label[clinical_label=="s50GGO"] <- "GGO50"
clinical_label[clinical_label=="s75GGO"] <- "GGO75"
names(clinical_label) <- colnames(fraction_CD8T)
clinical_label[clinical_info[,1]] <- clinical_info[,2]
n_clinical_label <- table(clinical_label)
clinical_labelN <- paste0(clinical_label,"(n=",n_clinical_label[clinical_label],")")
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(3,4,5,6,7)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")

#monocyte and macrophage 
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_Macro[3,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of Macro_S03")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#B cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_B[1,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of B_S01")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

##CD4+ T cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_CD4T[1,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of CD4T_S01")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#CD8+ T cells
mono_S3_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
fraction=unlist(fraction_CD8T[2,]),stringsAsFactors=F)
my_comparisons <- list( c("Normal(n=38)", "SolidN(n=5)"))
p<-ggboxplot(mono_S3_df, x = "stage", y = "fraction", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of CD8T S02")+
		  stat_compare_means(comparisons = my_comparisons,
		  symnum.args =list(cutpoints = c(0,  0.001, 0.01, 0.05,0.1, 1), symbols = c( "***", "**", "*", ".","ns")),size = 8)
p

#############################
#Figure 2A (R version 4.1.1)#
#############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${clinicalFile} is file of clinical information
#${EpiFile} is the seurat object of all the Epithelium

#Load required packages
library(reticulate)
use_python(pythonDir,required=T)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library('ggrastr')
library(RColorBrewer)

#Set working directory
setwd(workDir)

#Load the clinical data
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
#Load the seurat object of epithelium
epi.integrated <- readRDS(file = EpiFile)

#label epithlial cells with copy number
label <- as.character(epi.integrated@meta.data[,"orig.ident"])
names(label) <- rownames(epi.integrated@meta.data)
label2 <- clinical[match(label,clinical[,"rawDataID"]),"class"]
names(label2) <- rownames(epi.integrated@meta.data)
label2[label2!="Normal"] <- "precancerous"
label2[grep("c",as.character(epi.integrated@meta.data[,"combAnno"]))] <- "cancer"
epi.integrated$clinical <- label2

#plot for benign patients
nor <- subset(epi.integrated,cells=rownames(epi.integrated@meta.data)[as.character(epi.integrated@meta.data[,"orig.ident"])%in%clinical[clinical[,"class"]=="Normal","rawDataID"]])
b_embed1 <- nor@"reductions"[["tsne"]]@"cell.embeddings"
b_df1 <- data.frame(tSNE_1=b_embed1[,"tSNE_1"],tSNE_2=b_embed1[,"tSNE_2"],stringsAsFactors=F)
p1 <- ggplot(b_df1, aes(tSNE_1, tSNE_2))+
geom_point( size = 0.01,col = brewer.pal(6, "Set3")[5]) +
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p1,scale = 0.25,dpi=300)

#plot for cancer patients
ca <- subset(epi.integrated,cells=rownames(epi.integrated@meta.data)[as.character(epi.integrated@meta.data[,"orig.ident"])%in%clinical[clinical[,"class"]!="Normal","rawDataID"]])
b_embed2 <- ca@"reductions"[["tsne"]]@"cell.embeddings"
b_df2 <- data.frame(tSNE_1=b_embed2[,"tSNE_1"],tSNE_2=b_embed2[,"tSNE_2"],colo=ca@meta.data[,"clinical"],stringsAsFactors=F)
p2 <- ggplot(b_df2, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point(size = 0.01) +scale_colour_manual(values=brewer.pal(7, "Set3")[c(4,7)]) +
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p2,scale = 0.25,dpi=300) 