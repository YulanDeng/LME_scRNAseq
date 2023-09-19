## Script for plot Figure 2 and Figure S3
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

#############################
#Figure 2a (R version 4.0.3)#
#############################
#${clinicalFile} is the file of clinical information
#${fractionFile} is the result of Ecotyper
#${workDir} is the working directory

#Load required packages
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)

#Set working directory
setwd(workDir)

#Load the clinical information and rename
clinical_info <- read.table(file=clinicalFile,sep="\t",stringsAsFactors=F,header=T)
clinical_info[,1] <- paste0("X",clinical_info[,1],"C")
clinical_info[clinical_info[,"class"]%in%c("pGGO","HiDenGGO","s25GGO","s50GGO"),"class"] <- "<50GGO"
clinical_info[clinical_info[,"class"]%in%c("s75GGO","s100GGO","Solid1","Solid3"),"class"] <- ">50GGO"
clinical_info[clinical_info[,"class"]%in%c("s75GGO","Solid1","Solid3"),"class"] <- "Solid"

#Load the result of Ecotyper
fraction <- read.table(file=fractionFile,sep="\t",stringsAsFactors=F,header=T)
fraction <- fraction[,!(colnames(fraction)%in%c("X10C","X10F","X16C","X16F"))]


clinical_label <- rep("Normal",ncol(fraction))
names(clinical_label) <- colnames(fraction)
clinical_label[clinical_info[,1]] <- clinical_info[,2]
n_clinical_label <- table(clinical_label)
clinical_labelN <- paste0(clinical_label,"(n=",n_clinical_label[clinical_label],")")

fraction_list <- lapply(seq(nrow(fraction)),function(x){
	res <- split(unlist(fraction[x,]),
	f=factor(clinical_labelN,
	levels=c("Normal(n=38)","<50GGO(n=26)",">50GGO(n=7)","SolidN(n=5)"  )))
	return(res)
})

colo <- c(brewer.pal(3, "Set1")[3],brewer.pal(8, "Blues")[4],
brewer.pal(6, "OrRd")[c(4,6)])

names(colo) <- c("Normal(n=38)","<50GGO(n=26)",">50GGO(n=7)","SolidN(n=5)")

#LIME01
E1_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","<50GGO(n=26)",">50GGO(n=7)","SolidN(n=5)")),
E1=unlist(fraction[1,]),stringsAsFactors=F)

#kruskal.test
kruskal.test(E1~stage, data = E1_df)
#p-value = 2.947e-09

combn_pair <- combn(unique(clinical_labelN),2)
#my_comparisons <- lapply(seq(ncol(combn_pair)),function(x) combn_pair[,x])
my_comparisons <- list()
my_comparisons[[1]] <-c("Normal(n=38)","<50GGO(n=26)")
my_comparisons[[2]] <-c("Normal(n=38)",">50GGO(n=7)")
my_comparisons[[3]] <-c("Normal(n=38)","SolidN(n=5)")
my_comparisons[[4]] <-c("<50GGO(n=26)",">50GGO(n=7)")

p<-ggboxplot(E1_df, x = "stage", y = "E1", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("Fraction of LME01")+
		  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",size = 8)
p

#LIME05
E5_df <- data.frame(stage=factor(clinical_labelN,
levels=c("Normal(n=38)","<50GGO(n=26)",">50GGO(n=7)","SolidN(n=5)")),
E5=unlist(fraction[5,]),stringsAsFactors=F)
#kruskal.test
kruskal.test(E5~stage, data = E5_df)
#p-value = 3.564e-07

combn_pair <- combn(unique(clinical_labelN),2)
#my_comparisons <- lapply(seq(ncol(combn_pair)),function(x) combn_pair[,x])
my_comparisons <- list()
my_comparisons[[1]] <-c("Normal(n=38)","<50GGO(n=26)")
my_comparisons[[2]] <-c("Normal(n=38)",">50GGO(n=7)")
my_comparisons[[3]] <-c("Normal(n=38)","SolidN(n=5)")

p<-ggboxplot(E5_df, x = "stage", y = "E5", color="stage",outlier.shape = NA,
palette = colo,legend = "none")+geom_jitter(width=0.1,pch=16,cex=3.3,col="black")+xlab(NULL)+
theme(panel.background=element_rect(fill='transparent', color='black'),
axis.ticks.length = unit(.25, "cm"),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
		  axis.text.x=element_text(angle=45,color ="black",size=18,hjust=0.5,vjust=0.5),
		  axis.title.y=element_text(size=18),
		  axis.text.y=element_text(color ="black",size=18)) + ylab("Fraction of LME05")+
		  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",size = 8)
p

#############################
#Figure 2e (R version 4.0.3)#
#############################
#${cellChatFile} is the cellchat result
#${workDir} is the working directory
#${moEcotyperFile} is the results of Ecotyper for macrophage
#${markerDir} is the directory for marker genes
#${EcotyperFile} is the results of Ecotyper 
#${EcotyperDir} is the directory for Ecotyper results

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)

#Set working directory
setwd(workDir)


cell1="Mast.cells"
state1=4
cell2="Monocytes.and.Macrophages"
state2=3
cellchat=cellchat
rank_data=rank_data
ggo.integrated=ggo.integrated

load(cellChatFile)
	nrow(nicheNet_cellChat_con_df)
	res_df <- nicheNet_cellChat_con_df[order(nicheNet_cellChat_con_df[,"cellChat"],decreasing=T),]

	LR <- colnames(cellchat@net[[2]][1, ,])
	e1cd42sig <- LR[cellchat@net[[2]][paste0(cell1,"_S0",state1) ,paste0(cell2,"_S0",state2),]<0.05]
	
	ecotyper_mo_index <- rank_data[rank_data[,1]==cell2,2]
	ecotyper_mo <- read.table(file =paste0("/NAS/dyl/software/ecotyper/ecotyper-master/EcoTyper/discovery_scRNA_GGO_50_975_200_TPM/Cell_type_specific_genes/Cell_States/discovery/",
	cell2,"/",ecotyper_mo_index,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)

	ecotyper_epi_index <- rank_data[rank_data[,1]==cell1,2]
	ecotyper_epi <- read.table(file =moEcotyperFile,sep="\t",stringsAsFactors=F, header=T)
	
	marker.gene_list <-lapply(seq(length(unique(ecotyper_mo[,2]))),function(x){
		marker.gene <- readRDS(file=paste0(markerDir,cell2,".S0",x,".markerGenes.rds"))
		return(marker.gene)
	})

	sub.integrated <- subset(ggo.integrated,
	cells=c(intersect(rownames(ggo.integrated@meta.data),ecotyper_mo[ecotyper_mo[,2]==paste0("S0",state2),1]),
	intersect(rownames(ggo.integrated@meta.data),ecotyper_epi[ecotyper_epi[,2]==paste0("S0",state1),1])))
	DefaultAssay(sub.integrated)<-"RNA"
	expression = t(as.matrix(GetAssayData(sub.integrated,assay = "RNA")))

	sample_info = as.data.frame(sub.integrated@meta.data,stringsAsFactors =F)
	##Step 1: Define expressed genes in sender and receiver cell populations
	cell_type <- rep("unknown",nrow(sample_info))
	names(cell_type) <- rownames(sample_info)
	inter1 <- intersect(ecotyper_epi[ecotyper_epi[,2]==paste0("S0",state1),1],names(cell_type))
	cell_type[inter1] <- paste0(cell1,"_S0",state1)
	inter2 <- intersect(ecotyper_mo[ecotyper_mo[,2]==paste0("S0",state2),1],names(cell_type))
	cell_type[inter2] <- paste0(cell2,"_S0",state2)
	sample_info <- cbind(sample_info,cell_type)
	sub.integrated$minorLabel <- cell_type

	CAF_ids = rownames(sample_info)[sample_info[,"cell_type"]==paste0(cell1,"_S0",state1)] 

	malignant_ids = rownames(sample_info)[sample_info[,"cell_type"]==paste0(cell2,"_S0",state2)] 

	expressed_genes_sender = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 2] %>% names()

	expressed_genes_receiver = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()

	##Step 2: Define the gene set of interest and a background of genes
	geneset_oi = unique(unlist(lapply(marker.gene_list[state2],function(x){
		res <- rownames(x)[(x[,2]>log(1.5))&(x[,5]<0.05)&(x[,3]>0.5)]
		return(res)
	})))

	background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

	##Step 3: Define a set of potential ligands
	lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

	ligands = lr_network %>% pull(from) %>% unique()
	expressed_ligands = intersect(ligands,expressed_genes_sender)

	receptors = lr_network %>% pull(to) %>% unique()
	expressed_receptors = intersect(receptors,expressed_genes_receiver)

	lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 

	potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

	##Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
	ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

	t1 <- ligand_activities %>% arrange(-pearson)
	#best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
	best_upstream_ligands = c("ANXA1","CSF1","IL13")

	#Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
	active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()

	active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

	active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)

	active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.01)
	
	order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
	order_targets = active_ligand_target_links_df$target %>% unique()
	vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
	
	p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(paste0("Prioritized ",cell1,"_S0",state1,"-ligands"),paste0("feature genes in ",cell2,"_S0",state2), color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
	p_ligand_target_network

	ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

	order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
	
	order_targets = active_ligand_target_links_df$target %>% unique()
	vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

	vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

	# get the ligand-receptor network of the top-ranked ligands
	lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
	best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

	# get the weights of the ligand-receptor interactions as used in the NicheNet model
	lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

	# convert to a matrix
	lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
	lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

	# perform hierarchical clustering to order the ligands and receptors
	dist_receptors = dist(lr_network_top_matrix, method = "binary")
	hclust_receptors = hclust(dist_receptors, method = "ward.D2")
	order_receptors = hclust_receptors$labels[hclust_receptors$order]

	dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
	hclust_ligands = hclust(dist_ligands, method = "ward.D2")
	order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

	vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
	p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(paste0("Prioritized ",cell1,"_S0",state1,"-ligands"),paste0("Receptors expressed by ",cell2,"_S0",state2), color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
	p_ligand_receptor_network

	ligand_aupr_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
	vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
	
	p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot(paste0("Prioritized ",cell1,"_S0",state1,"-ligands"),"Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson\n(target gene prediction ability)")
	p_ligand_aupr
	
	ecotype <- read.table(file=EcotyperFile,sep="\t",stringsAsFactors=F,header=T)

ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
	print(x)
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	pathD=EcotyperDir
	ecotyper_df=read.table(file =paste0(pathD,
	cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	
	abundance <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/state_abundances.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	
	log2TPMf <- read.table(file=paste0(pathD,cell_type,"/expression_top_genes_scaled.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	sample_inter <- intersect(colnames(log2TPMf),colnames(abundance))
	cell_state <- sapply(sample_inter,function(x) rownames(abundance)[which(abundance[,x]==max(abundance[,x]))[1]])
	
	log2TPMf <- log2TPMf[,sample_inter]

	geneTopFC <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/gene_info.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	matchingInitial <- read.table(file=paste0(pathD,cell_type,"/",choose_rank,"/mapping_to_initial_states.txt"),
	sep="\t",stringsAsFactors=F,header=T)
	gene_top <- lapply(seq(nrow(matchingInitial)),function(y){
		geneFCtmp <- geneTopFC[geneTopFC[,"InitialState"]==matchingInitial[y,2],c("Gene","MaxFC")]
		gene_inter_tmp <- intersect(geneFCtmp[,"Gene"],rownames(log2TPMf))[1:50]
		meanstate <- apply(log2TPMf[gene_inter_tmp,cell_state==matchingInitial[y,1]],1,mean)
		meanothers <- apply(log2TPMf[gene_inter_tmp,cell_state!=matchingInitial[y,1]],1,mean)
		geneFCtmpOrder <- meanstate- meanothers
		names(geneFCtmpOrder) <- gene_inter_tmp
		return(names(geneFCtmpOrder)[geneFCtmpOrder>0.5])
	})
	names(gene_top) <- matchingInitial[,1]
	
	smp_list <- lapply(seq(nrow(matchingInitial)),function(y){
		res1 <- colnames(log2TPMf)[cell_state==matchingInitial[y,1]]
		return(res1)
	})
	names(smp_list) <- matchingInitial[,1]
	
	res <- list()
	res[[1]] <- gene_top
	res[[2]] <- smp_list
	return(res)
})

names(ecotyper_list) <- rank_data[,1]

cell_type_list <- unique(rank_data[,1])

cell_state_df <- ecotype[ecotype[,"Ecotype"]=="E1",]

sample_list <- lapply(seq(nrow(cell_state_df)),function(x){
		return(ecotyper_list[cell_state_df[x,1]][[1]][[2]][cell_state_df[x,2]][[1]])
	})
names(sample_list) <- cell_state_df[,"ID"]

ecotyper.integrated <- subset(ggo.integrated,cells=unlist(sample_list))
cell_state_label <- rep(names(sample_list),sapply(sample_list,length))
names(cell_state_label) <- unlist(sample_list)
ecotyper.integrated$"ID" <- cell_state_label[rownames(ecotyper.integrated@meta.data)]

DefaultAssay(ecotyper.integrated)<-"RNA"

VlnPlot(ecotyper.integrated, 
features = unique(res_df[,1]), 
ncol=3,   pt.size = 0,group.by="ID",cols=brewer.pal(8, "Set1"))

VlnPlot(ecotyper.integrated, 
features = order_receptors, 
ncol=3,   pt.size = 0,group.by="ID",cols=brewer.pal(8, "Set1"))
