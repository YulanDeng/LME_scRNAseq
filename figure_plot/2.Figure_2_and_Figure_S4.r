## Script for plot Figure 2 and Figure S4
##Yulan Deng, last updated 2022-6-27

#############################
#Figure 1E (R version 4.1.1)#
#############################
#${workDir} is the working directory
#${pythonDir} is directory of python
#${cancerFile} is the seurat object of cancer cell

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

#Load the seurat object of cancer cell
cancer.integrated <- readRDS(file =cancerFile)

b_embed <- cancer.integrated@"reductions"[["tsne"]]@"cell.embeddings"
colo <- as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])
b_df <- data.frame(tSNE_1=b_embed[,"tSNE_1"],tSNE_2=b_embed[,"tSNE_2"],colo=colo,stringsAsFactors=F)
p <- ggplot(b_df, aes(tSNE_1, tSNE_2, colour = colo))+
geom_point() +scale_colour_manual(values=brewer.pal(8, "Set2")[1:8])+
theme(panel.grid.major=element_line(colour=NA),
panel.background=element_rect(fill="transparent",colour=NA),
plot.background=element_rect(fill="transparent",colour=NA),
panel.grid.minor=element_blank())
rasterize(p,scale = 0.25,dpi=300)

#############################
#Figure 1F (R version 4.0.3)#
#############################
#${workDir} is the working directory
#${cancerFile} is the seurat object of cancer cell
#${ecoFile} is the cancer cell state information from Ecotyper

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(riverplot)
library(RColorBrewer)
library(gplots)
library(fgsea)

#Set working directory
setwd(workDir)

#Load the seurat object of cancer cell
cancer.integrated <- readRDS(file=cancerFile)
#Load the cancer cell state information from Ecotyper
ecotyper_anno <- read.table(file =ecoFile,sep="\t",stringsAsFactors=F, header=T)

#relabel cluster names
label <- rep("unassigned",nrow(cancer.integrated@meta.data))
names(label) = rownames(cancer.integrated@meta.data)
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S01",1],names(label))] = "S01"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S02",1],names(label))] = "S02"
label[intersect(ecotyper_anno[ecotyper_anno[,2]=="S03",1],names(label))] = "S03"
cancer.integrated$ecoRecovery = label
label_new <- as.character(cancer.integrated@meta.data[,"integrated_snn_res.0.25"])
label_new[label_new=="0"] <- "cancer_c0_SCGB3A2"
label_new[label_new=="1"] <- "cancer_c1_NEAT1"
label_new[label_new=="2"] <- "cancer_c2_MT2A"
label_new[label_new=="3"] <- "cancer_c3_FTH1"
label_new[label_new=="4"] <- "cancer_c4_AGER"
label_new[label_new=="5"] <- "cancer_c5_KRT5"
label_new[label_new=="6"] <- "cancer_c6_LAMP3"
label_new[label_new=="7"] <- "cancer_c7_DNAH12"
cancer.integrated$label_new = label_new
stat <- table(cancer.integrated@meta.data[,c("label_new","ecoRecovery")])

#Set the parameter for river plot
edge_cn <- matrix(unlist(strsplit(outer(c("cancer_c0_SCGB3A2","cancer_c1_NEAT1","cancer_c2_MT2A",
"cancer_c3_FTH1","cancer_c4_AGER","cancer_c5_KRT5","cancer_c6_LAMP3","cancer_c7_DNAH12"),paste0("S0",1:3),paste)," ")),
ncol=2,byrow=T)
edges <- data.frame(N1=edge_cn[,1],N2=edge_cn[,2],Value=apply(edge_cn,1,function(x) stat[x[1],x[2]]),stringsAsFactors=FALSE)
nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),
x=rep(c(1,2), c(8,3)), y=c(1:8,1:3))
cols <- c(brewer.pal(9, "Pastel1")[c(2,2,3,9,9,9,9,1)],brewer.pal(3, "Pastel1"))
names(cols) <- unique(c(edges$N1, edges$N2))
style <- sapply(nodes$ID, function(id) list(col=cols[ id ]), simplify=FALSE)
r <- makeRiver(nodes=nodes, edges=edges, styles=style)
d <- list(srt=0, textcex=1)
#river plot
plot(r, plot_area=1, nodewidth=10, default_style=d)

##############################
#Figure S4A (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${smp} is the sample ID
#${heatmapCode} is the code to plot heatmap provided by copyKAT
#To avoid ambiguity of function name, change the name into 'heatmap.3'
#${copyKATfile} is the result of copyKAT
#${meta.dataFile} is the meta data of all the cells

#${cancerFile} is the seurat object of cancer cell
#${ecoFile} is the cancer cell state information from Ecotyper

#Load required packages
library(RColorBrewer)
library(cluster)

#Set working directory
setwd(workDir)

#source the code to plot heatmap provided by copyKAT
source(heatmapCode)
#Load the result of copyKAT
copykat.test <- readRDS(file = copyKATfile)
#Load the meta.data of copyKAT
meta.data <- read.table(file =meta.dataFile,sep="\t",quote = "",stringsAsFactors=F, header=T,row.names=1)

#get the copy number prediction from copykat
pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

#make the matrix for heatmap
colnames(CNA.test) <- gsub(".","-",colnames(CNA.test),fixed=T)
epi <- rownames(meta.data)[(as.character(meta.data[,"orig.ident"])==smp)&(as.character(meta.data[,"assigned_cell_type"])=="Epithelial")]
tumor.cells <- intersect(as.character(pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]),
epi) 
pred.test <- pred.test[tumor.cells,]

#Set the parmeter for heatmap
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]
cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

#plot heatmap
heatmap.3(t(CNA.test[,tumor.cells]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =32, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")

##############################
#Figure S4A (R version 4.1.1)#
##############################
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

##############################
#Figure S4C (R version 4.0.3)#
##############################
#${workDir} is the working directory
#${clinicalFile} is file of clinical information of bulk RNAseq
#${EcoFile} is the recovery of cancer cell states in bulk RNAseq

#Load required packages
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

#Set working directory
setwd(workDir)

#Load the recovery of cancer cell states in bulk RNAseq
fraction_epi <- read.table(file=EcoFile,sep="\t",stringsAsFactors=F,header=T)

#Load the clinical infromation of bulk RNAseq
clinical_info <- read.table(file=clinicalFile,
sep="\t",stringsAsFactors=F,header=T)

#match the sample between cancer cell state and clinical information
fraction_epi <- fraction_epi[,!(colnames(fraction_epi)%in%c("X10C","X10F","X16C","X16F"))]
clinical_info[,1] <- paste0("X",clinical_info[,1],"C")
clinical_info[clinical_info[,"class"]=="HiDenGGO","class"] <- "dGGO"
clinical_info[clinical_info[,"class"]=="s25GGO","class"] <- "GGO25"
clinical_info[clinical_info[,"class"]=="s50GGO","class"] <- "GGO50"
clinical_info[clinical_info[,"class"]=="s75GGO","class"] <- "GGO75"
clinical_label <- rep("Normal",ncol(fraction_epi))
names(clinical_label) <- colnames(fraction_epi)
clinical_label[clinical_info[,1]] <- clinical_info[,2]
n_clinical_label <- table(clinical_label)
clinical_labelN <- paste0(clinical_label,"(n=",n_clinical_label[clinical_label],")")

#set the parameter for box plot
colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[c(3,4,5,6,7)],
brewer.pal(6, "OrRd")[c(2,4,6)])
names(colo) <- c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")
Epi_S03_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","pGGO(n=7)","dGGO(n=9)","GGO25(n=5)","GGO50(n=5)",
	"GGO75(n=1)","Solid1(n=3)","Solid3(n=3)","SolidN(n=5)")),
Epi_S03=unlist(fraction_epi[3,]),stringsAsFactors=F)

#boxplot
p<-ggboxplot(Epi_S03_df, x = "stage", y = "Epi_S03", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("the fraction of Epi S03")
p