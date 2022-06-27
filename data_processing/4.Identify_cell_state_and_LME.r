## Script for the annotation of non-epithlium from scRNAseq
##Yulan Deng, last updated 2022-6-23

##########################################
#Input file for ecotype (R version 3.6.1)#
##########################################
#${workDir} is the working directory
#${exampleFile} is the example file of Ecotype
#${immuneFile} is the seurat object for immune annotation
#${allCellFile} is the seurat object for major cell type annotation
#${clinicalFile} is the file for clinical information, 
#where the column named 'class' represents clinical stage

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)

#Set working directory
setwd(workDir)

#Load example file of Ecotyper 
example_mt <- read.table(file =exampleFile,sep="\t",
stringsAsFactors=F, header=T)

#Load seurat object for immune annotation and major cell type annotation
immu.integrated <- readRDS(file = immuneFile)
ggo.integrated <- readRDS(file = allCellFile)

#Load clinical information
clinical <- read.table(file =clinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical_sub <- clinical[!(clinical[,"class"]%in%c("Dermoid","Sarcomatoid")),]
#datasets from patients of non-adenocarcima are exluded for ecotyper analysis
GGO <- clinical_sub[,"rawDataID"]

#select the cell type and genes for ecotyper analysis
clinical_sub <- clinical[!(clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid")),]
GGO <- clinical_sub[,"rawDataID"]
cellType_all <- as.character(ggo.integrated@meta.data[,"assigned_cell_type"])
names(cellType_all) <- rownames(ggo.integrated@meta.data)
cellType_immune <- as.character(immu.integrated@meta.data[,"majorLabel"])
names(cellType_immune) <- rownames(immu.integrated@meta.data)
cellType_all[intersect(names(cellType_all),names(cellType_immune)[grep("B_sub",cellType_immune)])] <- "PCs"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="CD8T"])] <- "CD8.T.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="CD4T"])] <- "CD4.T.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="MAST"])] <- "Mast.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="NK"])] <- "NK.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="B"])] <- "B.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="mono"])] <- "Monocytes.and.Macrophages"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="neu"])] <- "PMNs"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="DC"])] <- "Dendritic.cells"
cellType_all[cellType_all=="Epithelial"] <- "Epithelial.cells"
cellType_all[cellType_all=="Endothelial"] <- "Endothelial.cells"
cellType_all[cellType_all=="Fibroblasts"] <- "Fibroblasts"
ggo.integrated <- subset(ggo.integrated,cells=rownames(ggo.integrated@meta.data)[(as.character(ggo.integrated@meta.data[,"orig.ident"])%in%GGO)&(cellType_all%in%c("PCs","CD8.T.cells","CD4.T.cells","Mast.cells",
"NK.cells","B.cells","Monocytes.and.Macrophages","PMNs","Dendritic.cells","Epithelial.cells", "Endothelial.cells",
"Fibroblasts"))])
cellType_all <- as.character(ggo.integrated@meta.data[,"assigned_cell_type"])
names(cellType_all) <- rownames(ggo.integrated@meta.data)
cellType_immune <- as.character(immu.integrated@meta.data[,"majorLabel"])
names(cellType_immune) <- rownames(immu.integrated@meta.data)
cellType_all[intersect(names(cellType_all),names(cellType_immune)[grep("B_sub",cellType_immune)])] <- "PCs"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="CD8T"])] <- "CD8.T.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="CD4T"])] <- "CD4.T.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="MAST"])] <- "Mast.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="NK"])] <- "NK.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="B"])] <- "B.cells"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="mono"])] <- "Monocytes.and.Macrophages"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="neu"])] <- "PMNs"
cellType_all[intersect(names(cellType_all),names(cellType_immune)[cellType_immune=="DC"])] <- "Dendritic.cells"
cellType_all[cellType_all=="Epithelial"] <- "Epithelial.cells"
cellType_all[cellType_all=="Endothelial"] <- "Endothelial.cells"
cellType_all[cellType_all=="Fibroblasts"] <- "Fibroblasts"
intergene <- intersect(example_mt[,1],rownames(ggo.integrated@assays$RNA@counts))
count_mt <- as.matrix(ggo.integrated@assays$RNA@counts)

#normalize count matrix to TPM matrix
lib_size=colSums(count_mt)
TPM_dataM <- t(t(as.matrix(count_mt))/lib_size[colnames(count_mt)]*10^6)
TPM_dataM <- cbind(rownames(TPM_dataM),TPM_dataM)
colnames(TPM_dataM)[1] <- "Gene"

#save TPM data and annotation data
write.table(TPM_dataM, file="scRNA_data.txt", sep="\t" , quote=FALSE, row.names=F, col.names=T )
Sample <- as.character(ggo.integrated@meta.data[,"orig.ident"])
anno_df <- data.frame(ID=names(cellType_all),CellType=cellType_all,Sample=Sample,stringsAsFactors=F)
write.table(anno_df, file="scRNA_annotation.txt", sep="\t" , quote=FALSE, row.names=F, col.names=T )

################################
#config file for ecotype (.yml)#
################################
#${workDir} is the working directory
default :
  Input :    
    Discovery dataset name : "discovery_scRNA_GGO"
    Expression matrix : "scRNA_data.txt"    
    Annotation file : "cRNA_annotation.txt" 
    Annotation file column to scale by : NULL
    Annotation file column(s) to plot : []    
  Output :
    Output folder : workDir
  Pipeline settings :
    #Pipeline steps:
    #   step 1 (extract cell type specific genes)
    #   step 2 (cell state discovery on correrlation matrices)
    #   step 3 (choosing the number of cell states)
    #   step 4 (extracting cell state information)
    #   step 5 (cell state re-discovery in expression matrices)
    #   step 6 (extracting information for re-discovered cell states)
    #   step 7 (cell state QC filter)
    #   step 8 (ecotype discovery)
    Pipeline steps to skip : [] 
    Filter non cell type specific genes : True
    Number of threads : 340
    Number of NMF restarts : 50
    Maximum number of states per cell type : 10
    Cophenetic coefficient cutoff : 0.975

####################
#run Ecotype (bash)#
####################
Rscript EcoTyper_discovery_scRNA.R -c config_discovery_scRNA.yml

#################################################################
#select the best number of candidate ecotypes  (R version 4.0.3)#
#################################################################
#Initially, k=2 is the best number of candidate ecotypes.
#However, only one ecotype with more than three cell states was identified.
#Therefore,modified ecotypes_scRNA.R file of Ecotype
#${workDir} is the working directory

#Load required packages
suppressPackageStartupMessages({
library(cluster)
library(ggplot2)
library(viridis)
library(reshape2)
source("lib/misc.R")
source("lib/heatmaps.R")
})

#Set working directory
setwd(workDir)

#set directory for file
dataset = "discovery_scRNA_GGO"
fractions = "Cell_type_specific_genes"
key_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
states_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery")
output_dir = file.path("../EcoTyper", dataset, fractions, "Ecotypes", "discovery")
dir.create(output_dir, recursive = T, showWarning = F) 

#load rank data, which records the best number of candidate cell states for each cell type.
key = read.delim(file.path(key_dir, "rank_data.txt"))

#Load the results of cell state
all_mapping = NULL
all_classes = NULL
all_classes_filt = NULL	
for(cell_type in key[,1])
{	
	#print(cell_type) 
	n_states = key[key[,1] == cell_type,2]

	mapping_path = file.path(states_dir, cell_type, n_states, 'mapping_to_initial_states.txt')
	classes_path = file.path(states_dir, cell_type, n_states, 'initial_state_assignment.txt')
	classes_filt_path = file.path(states_dir, cell_type, n_states, 'state_assignment.txt')

	mapping = read.delim(mapping_path)
	mapping = mapping[,c("State", "InitialState")] 
	mapping$CellType = cell_type
	all_mapping = rbind(all_mapping, mapping)

	classes = read.delim(classes_path)
	clinical = read_clinical(classes$ID, dataset = dataset)
	classes$Sample = clinical$Sample
	
	classes = as.data.frame(table(classes$Sample, classes$State))
	colnames(classes) = c("ID", "State", "Freq")
	splits = split(classes, classes$ID)
	classes = do.call(rbind, lapply(splits, function(spl)
	{
		spl$Frac = spl$Freq / sum(spl$Freq)
		spl$Max= ifelse((!is.na(max(spl$Freq))) & (max(spl$Freq) > 0) & (max(spl$Freq) == spl$Freq), 1, 0)
		spl
		}))
	classes$CellType = cell_type
	all_classes = rbind(all_classes, classes)
}
all_mapping$InitialID = paste(all_mapping$CellType, all_mapping$InitialState, sep = "_")
all_mapping$ID = paste(all_mapping$CellType, all_mapping$State, sep = "_")
write.table(all_mapping, file.path(output_dir, "mapping_all_states.txt"), sep = "\t", row.names = F)

#calculate Jaccard coefficient for the co-occurrence of cell states
casted = dcast(all_classes, ID ~ CellType + State, value.var = "Max")
#jaccard = cor(casted[,-1], use = "complete.obs")
clusters = t(casted[,-1])
clusters[is.na(clusters)] = 0
clusters = clusters[match(all_mapping$InitialID, rownames(clusters)),]
rownames(clusters) = all_mapping$ID
colnames(clusters) = casted[,1]
write.table(clusters, file.path(output_dir, "binary_classification_all_states.txt"), sep = "\t")
jaccard = matrix(NA, nrow(clusters), nrow(clusters))
for(i in 1:(nrow(clusters)))
{
	for(j in (i):(nrow(clusters)))
	{
		int <- sum(clusters[i,] & clusters[j,])
		idx <- int / (sum(clusters[i,]) + sum(clusters[j,]) - int)
		p = 1 - phyper(int, sum(clusters[i,]), ncol(all_classes) - sum(clusters[i,]), sum(clusters[j,]))
		if(is.na(p) | (p > 0.1)){
			#idx = 0
		}
		jaccard[i, j] = jaccard[j, i] = idx
	}
}
jaccard[is.na(jaccard)] = 0
rownames(jaccard) = colnames(jaccard) = rownames(clusters)
write.table(jaccard, file.path(output_dir, "jaccard_matrix.txt"), sep = "\t")

#Hierarchical Clustering of Jaccard coefficient and calculation ofsilhouette coefficient 
hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{
	#gap = clusGap(jaccard, hclusCut, 50, B = 100) 
	silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
		sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
		tmp <<- summary(sil)
		tmp$avg.width
	}))

	g2 <- ggplot(silh, aes(x = K, y = Silhouette)) + 
		geom_point() +
		#ylab("Sparseness (coef)") + 
		geom_line() + 
		geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") + 
		theme_bw() +
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 1) 
		#facet_wrap(CellType, ncol = 4) + 
		
	plot(g2)
	silh[which.max(silh$Silhouette), 1]	
}
hc = hclust(as.dist(1-jaccard), method = "average")
n_clust = choose_clusters(jaccard, "jaccard", range = 3:(nrow(jaccard) - 1))
#k=2 is exluded when enumeration for number of candidate ecotypes,
#therefore, the best resolution was select,that is ,k=12.  
clust =  hclusCut(jaccard, n_clust)$cluster
sil = silhouette(clust, as.dist(1-jaccard))
avg_silhouette = summary(sil)
write.table(avg_silhouette$avg.width, file.path(output_dir, "silhouette_initial.txt"), sep = "\t", row.names = F)

#assign cell states for each ecotype at appropriate number of candidate ecotype
top_ann = as.data.frame(t(sapply(rownames(jaccard), function(x) {
	s= strsplit(x, "_")[[1]]
	c(paste0(s[-length(s)],collapse = "_"), s[length(s)])
})))
colnames(top_ann) = c("CellType","State")
top_ann$InitialEcotype = as.factor(sprintf("IE%02d", clust))
write.table(top_ann, file.path(output_dir, "initial_ecotypes.txt"), sep = "\t")
top_ann$ID = rownames(top_ann)
top_ann = top_ann[order(top_ann$InitialEcotype),]
write.table(top_ann, file.path(output_dir, "ecotypes.txt"), sep = "\t", row.names = F)
jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]
top_ann$"Cell type" = top_ann$CellType
diag(jaccard) = 1
tb = table(top_ann$InitialEcotype)
tb = tb[tb > 2]
top_ann = top_ann[top_ann$InitialEcotype %in% names(tb),]
nm = unique(top_ann$InitialEcotype)
mapping = sprintf("E%d", 1:length(nm))
names(mapping) = nm
top_ann$Ecotype = mapping[as.character(top_ann$InitialEcotype)]
top_ann$Ecotype = ecotype_to_factor(top_ann$Ecotype)
top_ann = top_ann[order(top_ann$Ecotype),]
write.table(top_ann, file.path(output_dir, "ecotypes.txt"), sep = "\t", row.names = F)
jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]
sil <- silhouette(as.numeric(as.character(gsub("E", "", as.character(top_ann$Ecotype)))), as.dist(1-jaccard))
avg_silhouette <<- summary(sil)
write.table(avg_silhouette$avg.width, file.path(output_dir, "silhouette.txt"), sep = "\t", row.names = F)

################################
#recovery of LME in bulk RNAseq#
################################
#align by HISAT2
#${refDir} is the directory for hg38 reference
#${fastqFile1} is the fastq file for read 1
#${fastqFile2} is the fastq file for read 2
#${smp} is the sampleID
#${gtfFile} is the hg38 reference gtf file
#${pin} is the directory of gtf from stringtie
#${clinicalFile} is the file of clinical information
#${bulkRNAexpFile} is the annotation file for bulk RNAseq
#${bulkRNAannoFile} is the expression file for bulk RNAseq
hisat2 -p 8 --dta -t --rna-strandness RF -x ${refDir}  -1 ${fastqFile1} \
-2 ${fastqFile2} -S ${smp}.sam;
samtools view -@ 8 -bS ${smp}.sam > ${smp}.bam;
samtools sort -@ 8 -o ${smp}.sort.bam ${smp}.bam;
rm ${smp}.sam;
rm ${smp}.bam;

#assemble by stringtie
#first step
stringtie ${smp}.sort.bam -p 8 -G ${gtfFile} -l ${smp} >${smp}.gtf
#merge gtf
for j in ${File[@]}; do echo "${smp}.gtf" >> mergelist.txt; done
stringtie --merge -p 8 -G ${gtfFile} mergelist.txt >stringtie_merged.gtf
#second step
stringtie ${smp}.sort.bam -e -p 8 -G stringtie_merged.gtf -l ${smp} >${smp}.gtf

#make input file for Ecotype
setwd(pin)
File <- grep("gtf$",dir(),value=T)
File <- setdiff(File,"stringtie_merged.gtf")
#Load the gtf from stringtie
AllCount <- lapply(File,function(x){
	count <- read.table(file=x,sep="\t",stringsAsFactors=F,header=F)
	pos1 <- which(count[,3]=="transcript")
	subx <- count[pos1,]
	locus <- paste0("chr",subx[,1],":",subx[,4],"-",subx[,5])
	v9 <- gsub(";"," ",subx[,9]) 
	v9 <- strsplit(v9,"\\s")
	geneID <- sapply(v9,function(x) x[2])
	geneName <- sapply(v9,function(x) x[8])
	transcriptID <- sapply(v9,function(x) x[5])
	count <- sapply(v9,function(x) round(as.numeric(x[11])))
	TPM <- sapply(v9,function(x) as.numeric(x[17]))
	result <- data.frame(geneID=geneID,geneName=geneName,transcriptID=transcriptID,
	locus=locus,count=count,TPM=TPM,stringsAsFactors=F)
	return(result)	
})
names(AllCount) <- File

#common transcript
intertran <- AllCount[[1]][,3]
for(i in 2:length(AllCount))
{
intertran <- intersect(intertran,AllCount[[i]][,3])
}

#make TPM matrix
TPMMat <- sapply(AllCount,function(x) x[match(intertran,x[,3]),6])
rownames(TPMMat) <- intertran
colnames(TPMMat) <- sapply(strsplit(colnames(TPMMat),".",fixed=T),function(x) x[1])
clinical_match <- read.table(file=clinicalFile,sep="\t",stringsAsFactors=F,header=T)
colnames(TPMMat) <- clinical_match[match(colnames(TPMMat),clinical_match[,"Samplename"]),"Sampleid"]
TPMMatM <- cbind(rownames(TPMMat),TPMMat)
colnames(TPMMatM)[1] <- "Gene"
anno <- data.frame(ID=colnames(TPMMat),Tissue="Tumor",stringsAsFactors=F)
write.table(TPMMatM, file=bulkRNAexpFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
write.table(anno, file=bulkRNAannoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${bulkRNAexpFile} is the annotation file for bulk RNAseq
#${bulkRNAannoFile} is the expression file for bulk RNAseq
##${workDir} is the working directory
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${bulkRNAannoFile} \
-m ${bulkRNAexpFile} -o ${workDir}

##################################################
#recovery of LME in TCGA dataset(R version 3.6.1)#
##################################################
#make input file
#${expressionFile} is the expression file of TCGA LUAD downloaded from GDAC Firehose.
#${TCGAannoFile} is the annotation file for TCGA
#${TCGAexpFile} is the expression file for TCGA
#${workDir} is the working directory

#Set working directory
setwd(workDir)

#Load expresion file of TCGA
geneExpress <- read.table(file=expressionFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
#calculate TPM
expmt <- geneExpress[-1,seq(from=3,to=ncol(geneExpress),by=3)]
expmtTPM <- apply(expmt,2,function(x) as.numeric(x)*10^6)
rownames(expmtTPM) <- geneExpress[-1,1]
gene <- sapply(strsplit(rownames(expmtTPM),"|",fixed=T),function(x) x[1])
expmtTPM <- expmtTPM[!(gene%in%c("?","SLC35E2")),]
rownames(expmtTPM) <- sapply(strsplit(rownames(expmtTPM),"|",fixed=T),function(x) x[1])
expmtTPM_geneName <- cbind(rownames(expmtTPM),expmtTPM)
colnames(expmtTPM_geneName)[1] <- "Gene"
write.table(expmtTPM_geneName, file=TCGAannoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
ann_data = data.frame(ID=colnames(expmtTPM_geneName)[-1],Tissue="Tumor",stringsAsFactors=F)
write.table(ann_data, file=TCGAexpFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${TCGAannoFile} is the annotation file for TCGA
#${TCGAexpFile} is the expression file for TCGA
##${workDir} is the working directory
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${TCGAannoFile} \
-m ${TCGAexpFile} -o ${workDir}

###############################################
#recovery of LME in GSE31210 (R version 3.6.1)#
###############################################
#make input file
#${expressionFile} is the expression file of GSE31210 downloaded from GEO.
#${GSE31210annoFile} is the annotation file for GSE31210
#${GSE31210expFile} is the expression file for GSE31210
#${GSE31210probeFile} is the probe file for GSE31210
#${workDir} is the working directory

#Set working directory
setwd(workDir)

#Load expresion file of GSE31210
GSE31210_mt <- read.table(file=expressionFile,sep="\t",comment.char = "!",stringsAsFactors=F,header=T,row.names=1)
GSE31210_mtn <- apply(GSE31210_mt,2,as.numeric)
rownames(GSE31210_mtn) <- rownames(GSE31210_mt)

#Load probe file of GSE31210
GSE31210_probe <- read.table(file=GSE31210probeFile,sep="\t",comment.char = "#",
stringsAsFactors=F,header=T,fill=T)
GSE31210_mtnsub <- GSE31210_mtn

#Match probe value to gebe
inter <- intersect(GSE31210_probe[,1],rownames(GSE31210_mtnsub))\
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
expmtTPM_geneName <- cbind(rownames(GSE31210_mtnsub),GSE31210_mtnsub)
colnames(expmtTPM_geneName)[1] <- "Gene"
write.table(expmtTPM_geneName, file=GSE31210expFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

ann_data = data.frame(ID=colnames(expmtTPM_geneName)[-1],Tissue="Tumor",stringsAsFactors=F)
write.table(ann_data, file=GSE31210annoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${GSE31210annoFile} is the annotation file for GSE31210
#${GSE31210expFile} is the expression file for GSE31210
##${workDir} is the working directory
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${GSE31210annoFile} \
-m ${GSE31210expFile} -o ${workDir}

################################################
#recovery of LME in GSE140343 (R version 3.6.1)#
################################################
#make input file
#${gtfDir} is the directory of gtf file of GSE140343.
#${GSE140343annoFile} is the annotation file for GSE140343
#${GSE140343expFile} is the expression file for GSE140343
#${GSE140343clinicalFile} is the clinical file for GSE140343 obtain from supplementary table 
#${GSE140343SampleMatchFile} is the file recording matched sample names
#${workDir} is the working directory

#Set directory of gft file
setwd(gtfDir)
fl <- dir()
label <- grep("gtf$",fl,value=T)
label_end <- sapply(label,function(x) substr(x,nchar(x)-4,nchar(x))=="T.gtf")

#Load and process gtf files
gtf_list <- lapply(label,function(x){
	res <- read.table(file = paste0(gtfDir,
	x),sep="\t",stringsAsFactors=F,header=F)
	res_info <- res[res[,3]=="transcript",9]
	res_info_list <- strsplit(res_info," |;")
	ref_gene_index <- sapply(res_info_list,function(x) "ref_gene_name"%in%x)
	res_ref_info_list <- res_info_list[ref_gene_index]
	res_ref_info_list_len <- res_ref_info_list[sapply(res_ref_info_list,length)==17]
	gene_name <- sapply(res_ref_info_list_len,function(x) x[8])
	tran_id <- sapply(res_ref_info_list_len,function(x) x[5])
	TPM <- as.numeric(sapply(res_ref_info_list_len,function(x) x[17]))
	res_df <- data.frame(gene_name=gene_name,tran_id=tran_id,
	TPM=TPM,stringsAsFactors=F)
	return(res_df)
})
gtf_listT <- gtf_list[label_end] 

#choose transcript with highest expression for each gene 
all_tran_id <- unlist(lapply(gtf_listT,function(x) x[,"tran_id"]))
tab_tran_id <- table(all_tran_id)
com_tran_id <- names(tab_tran_id)[tab_tran_id==length(gtf_listT)]
com_tran_info <- gtf_listT[[1]][gtf_listT[[1]][,"tran_id"]%in%com_tran_id,]
tran_mt <- sapply(gtf_listT,function(x) x[match(com_tran_id,x[,"tran_id"]),"TPM"])
rownames(tran_mt) <- com_tran_id
colnames(tran_mt) <- label[label_end]
gene <- unique(com_tran_info[,"gene_name"])
gene_mt <- c()
for(i in gene)
{
	tran_in <- com_tran_info[com_tran_info[,"gene_name"]==i,"tran_id"]
	if(length(tran_in)==1)
	{
		gene_mt <- rbind(gene_mt,tran_mt[tran_in,])
	}else{
		tran_mean <- apply(tran_mt[tran_in,],1,function(y) mean(y[!is.na(y)])) 
		gene_mt <- rbind(gene_mt,tran_mt[tran_in[which(tran_mean==max(tran_mean))[1]],])
	}
	print(i)
} 
rownames(gene_mt) <- gene
colnames(gene_mt) <- sapply(strsplit(colnames(gene_mt),"_"),function(x) x[1])


#Set working directory
setwd(workDir)

#Load clinical information
clinical <- read.table(file = GSE140343clinicalFile,sep="\t",quote = "",stringsAsFactors=F,header=T)
clinical[,1] <- sub("_","",clinical[,1])
clinical_GSE140343 <- clinical

#Load matched information for sample name
sample_match <- read.table(file = "GSE140343_series_matrix.txt",sep="\t",quote = "",stringsAsFactors=F,header=F)
sample_matchT <- t(sample_match)
TPMMat <- gene_mt
colnames(TPMMat) <- sample_matchT[match(colnames(gene_mt),sample_matchT[,2]),1]
colnames(TPMMat) <- sub("RNA-seq_","",colnames(TPMMat))
TPMMatM <- cbind(rownames(TPMMat),TPMMat)
colnames(TPMMatM)[1] <- "Gene"
anno <- data.frame(ID=colnames(TPMMat),Tissue="Tumor",stringsAsFactors=F)
write.table(TPMMatM, file=GSE140343expFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
write.table(anno, file=GSE140343annoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${GSE140343annoFile} is the annotation file for GSE140343
#${GSE140343expFile} is the expression file for GSE140343
#${workDir} is the working directory
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${GSE140343annoFile} \
-m ${GSE140343expFile} -o ${workDir}

################################################
#recovery of LME in GSE166449 (R version 3.6.1)#
################################################
#the expression of GSE166449 is downloaded from SRA
#convert sra file to fastq file
#${path} is the directory for oufputfile
#${sraFile} is the sra file
fasterq-dump --split-files --outdir ${path} ${sraFile}

#align by hisat2
#${refDir} is the directory for hg38 reference
#${fastqFile1} is the fastq file for read 1
#${fastqFile2} is the fastq file for read 2
#${smp} is the sampleID
/NAS/dyl/software/hisat2/hisat2-2.2.0/hisat2 -p 20 --dta -t -x ${refDir} -q -1 ${fastqFile1} \
-2 ${fastqFile2} -S ${smp}.sam;
samtools view -@ 20 -bS  ${smp}.sam | samtools sort -@ 20 -o ${smp}.bam 

#assemble by stringtie
#${gtfFile} is the hg38 reference gtf file
#${smp} is the sampleID
stringtie -A ${smp} .tab -C ${smp}.cov_refs.gtf -e -b ${smp}.ctab \
-p 10 -G ${gtfFile} -o ${smp}.stringtie.gtf ${smp}.bam

#input for Ecotyper
#${workDir} is the working directory
#${pin} is the directory of gtf from stringtie
#${GSE166449annoFile} is the annotation file for GSE166449
#${GSE166449expFile} is the expression file for GSE166449
setwd(pin)
fl1 <- dir()
fl1s <- setdiff(fl1,"GSE166449.txt")
fl1sd <- sapply(strsplit(fl1s,".",fixed=T),function(x) x[1])
smp <- fl1sd
setwd(workDir)
tpm_list <- lapply(smp,function(x){
	fltmp <- read.table(file=paste0(pin,x,".tab"),
	sep="\t",stringsAsFactors=F,header=T)
	return(fltmp)
})
gene_list <- unlist(lapply(tpm_list,function(x) setdiff(unique(x[,2]),c("-","."))))
geneinter <- names(table(gene_list))[table(gene_list)==length(smp)]
tpm_mt <- sapply(tpm_list,function(x){
	res <- x[match(geneinter,x[,2]),"TPM"]
	return(res)
}) 
rownames(tpm_mt) <- geneinter
colnames(tpm_mt) <- smp
expmtTPM_geneName <- cbind(rownames(tpm_mt),tpm_mt)
colnames(expmtTPM_geneName)[1] <- "Gene"
setwd(workDir)
write.table(expmtTPM_geneName, file=GSE126044expFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
ann_data = data.frame(ID=colnames(expmtTPM_geneName)[-1],Tissue="Tumor",stringsAsFactors=F)
write.table(ann_data, file=GSE166449annoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${workDir} is the working directory
#${GSE166449annoFile} is the annotation file for GSE166449
#${GSE166449expFile} is the expression file for GSE166449
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${GSE166449annoFile} \
-m ${GSE166449expFile} -o ${workDir}

################################################
#recovery of LME in GSE135222 (R version 3.6.1)#
################################################
#the expression of GSE135222 is downloaded from SRA
#convert sra file to fastq file
#${path} is the directory for oufputfile
#${sraFile} is the sra file
fasterq-dump --split-files --outdir ${path} ${sraFile}

#align by hisat2
#${refDir} is the directory for hg38 reference
#${fastqFile1} is the fastq file for read 1
#${fastqFile2} is the fastq file for read 2
#${smp} is the sampleID
/NAS/dyl/software/hisat2/hisat2-2.2.0/hisat2 -p 20 --dta -t -x ${refDir} -q -1 ${fastqFile1} \
-2 ${fastqFile2} -S ${smp}.sam;
samtools view -@ 20 -bS  ${smp}.sam | samtools sort -@ 20 -o ${smp}.bam 

#assemble by stringtie
#${gtfFile} is the hg38 reference gtf file
#${smp} is the sampleID
stringtie -A ${smp} .tab -C ${smp}.cov_refs.gtf -e -b ${smp}.ctab \
-p 10 -G ${gtfFile} -o ${smp}.stringtie.gtf ${smp}.bam

#input for Ecotyper
#${workDir} is the working directory
#${pin} is the directory of gtf from stringtie
#${GSE135222annoFile} is the annotation file for GSE135222
#${GSE135222expFile} is the expression file for GSE135222
setwd(pin)
fl1 <- dir()
fl1s <- setdiff(fl1,"GSE135222.txt")
fl1sd <- sapply(strsplit(fl1s,".",fixed=T),function(x) x[1])
smp <- fl1sd
setwd(workDir)
tpm_list <- lapply(smp,function(x){
	fltmp <- read.table(file=paste0(pin,x,".tab"),
	sep="\t",stringsAsFactors=F,header=T)
	return(fltmp)
})
gene_list <- unlist(lapply(tpm_list,function(x) setdiff(unique(x[,2]),c("-","."))))
geneinter <- names(table(gene_list))[table(gene_list)==length(smp)]
tpm_mt <- sapply(tpm_list,function(x){
	res <- x[match(geneinter,x[,2]),"TPM"]
	return(res)
}) 
rownames(tpm_mt) <- geneinter
colnames(tpm_mt) <- smp
expmtTPM_geneName <- cbind(rownames(tpm_mt),tpm_mt)
colnames(expmtTPM_geneName)[1] <- "Gene"
setwd(workDir)
write.table(expmtTPM_geneName, file=GSE135222expFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
ann_data = data.frame(ID=colnames(expmtTPM_geneName)[-1],Tissue="Tumor",stringsAsFactors=F)
write.table(ann_data, file=GSE135222annoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${workDir} is the working directory
#${GSE135222annoFile} is the annotation file for GSE135222
#${GSE135222expFile} is the expression file for GSE135222
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${GSE135222annoFile} \
-m ${GSE135222expFile} -o ${workDir}

###############################################
#recovery of LME in GSE126044(R version 3.6.1)#
###############################################
#the expression and clinical file of GSE126044 are obtained from R package PD1highCD8Tscore
#${workDir} is the working directory
#${GSE126044annoFile} is the annotation file for GSE126044
#${GSE126044expFile} is the expression file for GSE126044
#${GSE126044File} is the R data from R package PD1highCD8Tscore
setwd(workDir)
load(GSE126044File)
expmtTPM_geneName <- cbind(rownames(Cho_RNA),Cho_RNA)
colnames(expmtTPM_geneName)[1] <- "Gene"
write.table(expmtTPM_geneName, file=GSE126044expFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )
ann_data = data.frame(ID=colnames(expmtTPM_geneName)[-1],Tissue="Tumor",stringsAsFactors=F)
write.table(ann_data, file=GSE126044annoFile, sep="\t" , quote=FALSE, row.names=F, col.names=T )

#run Ecotype
#${workDir} is the working directory
#${GSE126044annoFile} is the annotation file for GSE126044
#${GSE126044expFile} is the expression file for GSE126044
Rscript EcoTyper_recovery_bulk.R -d discovery_scRNA_GGO ${GSE126044annoFile} \
-m ${GSE126044expFile} -o ${workDir}