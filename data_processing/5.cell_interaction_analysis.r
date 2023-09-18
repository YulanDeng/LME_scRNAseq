## Script for the transcriptional trajectories analysis
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

######################################
#cell chat analysis (R version 4.0.3)#
######################################
#${workDir} is the working directory
#${ggoFile} is the seurat object for all the cells
#${EcotypeFile} is the result of Ecotype assignment
#${Rank_data} is the appropriate parameter for each cell state
#${CellStateDir} is the directory for the result of cell state

#Load required packages
library(CellChat)
library(patchwork)
library(Seurat)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

#Set working directory
setwd(workDir)

#Load the seurat object of all the cells
ggo.integrated <- readRDS(file=ggoFile)
#Load the result of Ecotype assignment
ecotye <- read.table(file=EcotypeFile,sep="\t",stringsAsFactors=F,header=T)
#Load the appropriate parameter for each cell state
rank_data <- read.table(file=Rank_data,sep="\t",stringsAsFactors=F,header=T)
#Load the result of cell state
ecotyper_list <- lapply(seq(nrow(rank_data)),function(x){
	cell_type=rank_data[x,1]
	choose_rank=rank_data[x,2]
	res=read.table(file =paste0(CellStateDir,cell_type,"/",choose_rank,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
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

#run cell chat
cellchat <- createCellChat(object = ecotype.integrate, group.by = "Ecotype")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 16) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat) 
cellchat <- aggregateNet(cellchat)

#plot the result of cell Chat for each LME
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E1","ID"],ecotye[ecotye[,"Ecotype"]=="E1","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E1","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E1",top = 0.6)
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E2","ID"],ecotye[ecotye[,"Ecotype"]=="E2","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E2","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E2",top = 0.6)
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E3","ID"],ecotye[ecotye[,"Ecotype"]=="E3","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E3","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E3",top = 0.6)
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E4","ID"],ecotye[ecotye[,"Ecotype"]=="E4","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E4","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E4",top = 0.6)
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E5","ID"],ecotye[ecotye[,"Ecotype"]=="E5","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E5","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E5",top = 0.6)
netVisual_circle(cellchat@net$count[ecotye[ecotye[,"Ecotype"]=="E6","ID"],ecotye[ecotye[,"Ecotype"]=="E6","ID"]],
vertex.weight = as.numeric(table(cellchat@idents)[ecotye[ecotye[,"Ecotype"]=="E6","ID"]]),weight.scale = T,
label.edge= F, title.name = "Number of interactions in E6",top = 0.6)

##################################################################################
#Validate the results for CellChat by NicheNet for each ecotype (R version 4.0.3)#
##################################################################################
#${cell1} is sender cell
#${state1} is the cell state of sender cell
#${cell2} is the receiver cell
#${state2} is the cell state of receiver cell
#${cellChatFile} is the results of CellChat
#${seuratFile} is the seurat object for all cells
#${rankFile} is the rank data of ecotype
#${ecotypeDir} is the directory for ecotye results
#${cellStateMarkersDir} is the directory of marker genes for cell states

#Load required packages
library(nichenetr)
library(tidyverse)
library("Seurat")
library(RColorBrewer)
library(cowplot)
library(ggpubr)

##load required files
load(cellChatFile)
ggo.integrated <- readRDS(file=seuratFile)
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
ligand_tf_matrix = readRDS(file="ligand_tf_matrix.rds")
sig_network = readRDS(file="signaling_network.rds")
gr_network = readRDS(file="gr_network.rds")
rank_data <- read.table(file=rankFile,sep="\t",stringsAsFactors=F,header=T)

#Validate the results for CellChat by NicheNet
nichenet_cellchat_fun <- function(cell1,state1,cell2,state2,cellchat=cellchat,rank_data=rank_data,ggo.integrated=ggo.integrated){
	LR <- colnames(cellchat@net[[2]][1, ,])
	e1cd42sig <- LR[cellchat@net[[2]][paste0(cell1,"_S0",state1) ,paste0(cell2,"_S0",state2),]<0.05]
	
	ecotyper_mo_index <- rank_data[rank_data[,1]==cell2,2]
	ecotyper_mo <- read.table(file =paste0(ecotypeDir,
	cell2,"/",ecotyper_mo_index,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)

	ecotyper_epi_index <- rank_data[rank_data[,1]==cell1,2]
	ecotyper_epi <- read.table(file =paste0(ecotypeDir,cell1,"/",ecotyper_epi_index,"/state_assignment.txt"),sep="\t",stringsAsFactors=F, header=T)
	
	marker.gene_list <-lapply(seq(length(unique(ecotyper_mo[,2]))),function(x){
		marker.gene <- readRDS(file=paste0(cellStateMarkersDir,cell2,".S0",x,".markerGenes.rds"))
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
	geneset_oi = unique(unlist(lapply(marker.gene_list[setdiff(seq(length(marker.gene_list)),state2)],function(x){
		res <- rownames(x)[(x[,2]>log(1.5))&(x[,5]<0.05)]
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

	ligand_activities %>% arrange(-pearson)
	best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

	#Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
	active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()

	active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

	active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)

	active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.01)

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

	nicheNet_res <- matrix(unlist(strsplit(outer(rownames(vis_ligand_receptor_network),colnames(vis_ligand_receptor_network),paste)," ")),byrow=T,ncol=2)

	e1cd42sig_list <- strsplit(e1cd42sig,"_")

	nicheNet_res_index <- apply(nicheNet_res,1,function(x) {
		cellchat_in <- sapply(e1cd42sig_list,function(y) all(x%in%y))
		return(which(cellchat_in))
	})

	nicheNet_res_sel <- nicheNet_res[sapply(nicheNet_res_index,length)>0,]
	nicheNet_res_sel_value <- apply(nicheNet_res_sel,1,function(x) vis_ligand_receptor_network[x[1],x[2]])
	nicheNet_res_index_sel <- sapply(nicheNet_res_index[sapply(nicheNet_res_index,length)>0],function(x) e1cd42sig[x[1]])
	cellchat_res_sel_value <- cellchat@net[[1]][paste0(cell1,"_S0",state1) ,paste0(cell2,"_S0",state2),nicheNet_res_index_sel]
	nicheNet_cellChat_con_df <- data.frame(ligand=nicheNet_res_sel[,2],receptor=nicheNet_res_sel[,1],
	nicheNet=nicheNet_res_sel_value,cellChat=cellchat_res_sel_value,stringsAsFactors=F)
	save(nicheNet_cellChat_con_df,file=paste0(cell1,"_S0",state1,"@",paste0(cell2,"_S0",state2),
	"@nicheNet_cellChat.RData"))
}


######################################################################
#NicheNat analysis from endothelium to CD8+ T cells (R version 4.0.3)#
######################################################################
#${workDir} is the working directory
#${BEAMfile} is the result from BEAM
#${cellChatFile} is result from cell chat
#${endoEcoFile} is the Ecotyper result for endothelium
#${CD8TecoFile} is the Ecotyper result for CD8+ T cells
#${CD8Tmarker} is the marker genes of CD8+ T cells
#${ggoFile} is the seurat object of all the cells

#Load required packages
library(nichenetr)
library(tidyverse)
library("Seurat")
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#Set working directory
setwd(workDir)

#Load the result from BEAM
BEAM_res<- readRDS(file = BEAMfile)

#Load the result from cell chat
load(cellChatFile)
#only interaction between endothelium cell state S01 and CD8+ T cell state S04 is considered
LR <- colnames(cellchat@net[[2]][1, ,])
e1cd84sig <- LR[cellchat@net[[2]]["Endothelial.cells_S01" ,"CD8.T.cells_S04",]<0.05]

#Load the result from BEAM
ecotyper_endo <- read.table(file =endoEcoFile,sep="\t",stringsAsFactors=F, header=T)

#the Ecotyper result for endothelium
ecotyper_CD8T <- read.table(file =CD8TecoFile,sep="\t",stringsAsFactors=F, header=T)

#Load the marker genes of CD8+ T cells
CD8T.markers <- readRDS(file = CD8Tmarker)

#Load the seurat obeject of all the cells
ggo.integrated <- readRDS(file=ggoFile)

#Load the reference data for NicheNet
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
ligand_tf_matrix = readRDS(file="ligand_tf_matrix.rds")
sig_network = readRDS(file="signaling_network.rds")
gr_network = readRDS(file="gr_network.rds")

#make the input data
sub.integrated <- subset(ggo.integrated,
cells=c(intersect(rownames(ggo.integrated@meta.data),ecotyper_CD8T[ecotyper_CD8T[,2]=="S04",1]),
intersect(rownames(ggo.integrated@meta.data),ecotyper_endo[ecotyper_endo[,2]=="S01",1])))
DefaultAssay(sub.integrated)<-"RNA"
expression = t(as.matrix(GetAssayData(sub.integrated,assay = "RNA")))
sample_info = as.data.frame(sub.integrated@meta.data,stringsAsFactors =F)

##Step 1: Define expressed genes in sender and receiver cell populations
cell_type <- rep("unknown",nrow(sample_info))
names(cell_type) <- rownames(sample_info)
inter1 <- intersect(ecotyper_endo[ecotyper_endo[,2]=="S01",1],names(cell_type))
cell_type[inter1] <- "endo_S1"
inter2 <- intersect(ecotyper_CD8T[ecotyper_CD8T[,2]=="S04",1],names(cell_type))
cell_type[inter2] <- "CD8T_S4"
sample_info <- cbind(sample_info,cell_type)
sub.integrated$minorLabel <- cell_type
CAF_ids = rownames(sample_info)[sample_info[,"cell_type"]=="endo_S1"] 
malignant_ids = rownames(sample_info)[sample_info[,"cell_type"]=="CD8T_S4"] 
expressed_genes_sender = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 2] %>% names()
expressed_genes_receiver = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()

##Step 2: Define the gene set of interest and a background of genes
geneset_oi = intersect(row.names(subset(BEAM_res,qval < 0.05)),
CD8T.markers[(CD8T.markers[,"p_val_adj"]<0.05)&(CD8T.markers[,"avg_logFC"]>log(1.8)&(CD8T.markers[,"cluster"]%in%c("3"))),"gene"])
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
ligand_activities %>% arrange(-pearson)
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

#Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()

#write the result
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.01)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized endothelial-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized endothelial-ligands","Treg genes in CD8 T cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
best_upstream_ligands=c("IL7","CCL21")
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
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized ligands","Receptors expressed by CD8T", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
ligands_all = c("IL7","CCL21") 
targets_all = c("GZMB","CTLA4","NR3C1","PRDM1","RUNX2")
active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network) 
# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
output_path = workDir
write_output = TRUE 
# change to TRUE for writing output
# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}
# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}
# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")
if(write_output){
bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}

######################################################################
#NicheNat analysis from endothelium to CD4+ T cells (R version 4.0.3)#
######################################################################
#${workDir} is the working directory
#${CD4Tmarker} is the marker genes for CD4+ T cells
#${BEAMfile} is the result of BEAM
#${cellChatFile} is the result of cell chat analysis
#${CD4TecoFile} is the Ecotyper result for CD4+ T cell
#${endoEcoFile} is the Ecotyper result for endothelium
#${ggoFile} is the seurat obeject for all the cell

#Load required packages
library(nichenetr)
library(tidyverse)
library("Seurat")
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#Set working directory
setwd(workDir)

#Load the marker gene of CD4+ T cells
load(CD4Tmarker)

#Load the result of BEAM
BEAM_res<- readRDS(file = BEAMfile)

#Load the result of cell chat
load(cellChatFile)
#only the interaction between endothelium S01 and CD4+ T cell S02 is considered
LR <- colnames(cellchat@net[[2]][1, ,])
e1cd42sig <- LR[cellchat@net[[2]]["Endothelial.cells_S01" ,"CD4.T.cells_S02",]<0.05]

#Load the ecotype result for CD4+ T cells
ecotyper_CD4T <- read.table(file =CD4TecoFile,sep="\t",stringsAsFactors=F, header=T)

#Load the ecotype result for endothelium
ecotyper_endo <- read.table(file =endoEcoFile,sep="\t",stringsAsFactors=F, header=T)

#Load the seurat object for all the cells
ggo.integrated <- readRDS(file=ggoFile)

#Load the reference data for NicheNet
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
ligand_tf_matrix = readRDS(file="ligand_tf_matrix.rds")
sig_network = readRDS(file="signaling_network.rds")
gr_network = readRDS(file="gr_network.rds")

#make the input data
sub.integrated <- subset(ggo.integrated,
cells=c(intersect(rownames(ggo.integrated@meta.data),ecotyper_CD4T[ecotyper_CD4T[,2]=="S02",1]),
intersect(rownames(ggo.integrated@meta.data),ecotyper_endo[ecotyper_endo[,2]=="S01",1])))
DefaultAssay(sub.integrated)<-"RNA"
expression = t(as.matrix(GetAssayData(sub.integrated,assay = "RNA")))
sample_info = as.data.frame(sub.integrated@meta.data,stringsAsFactors =F)

##Step 1: Define expressed genes in sender and receiver cell populations
cell_type <- rep("unknown",nrow(sample_info))
names(cell_type) <- rownames(sample_info)
inter1 <- intersect(ecotyper_endo[ecotyper_endo[,2]=="S01",1],names(cell_type))
cell_type[inter1] <- "endo_S1"
inter2 <- intersect(ecotyper_CD4T[ecotyper_CD4T[,2]=="S02",1],names(cell_type))
cell_type[inter2] <- "CD4T_S2"
sample_info <- cbind(sample_info,cell_type)
sub.integrated$minorLabel <- cell_type
CAF_ids = rownames(sample_info)[sample_info[,"cell_type"]=="endo_S1"] 
malignant_ids = rownames(sample_info)[sample_info[,"cell_type"]=="CD4T_S2"] 
expressed_genes_sender = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 2] %>% names()
expressed_genes_receiver = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()

##Step 2: Define the gene set of interest and a background of genes
geneset_oi = intersect(row.names(subset(BEAM_res,qval < 0.05)),
rownames(CD4_marker_S01)[(CD4_marker_S01[,2]>log(2))&(CD4_marker_S01[,5]<0.05)])
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
ligand_activities %>% arrange(-pearson)
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

#Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()

#write the result
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.01)
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized endothelial-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized endothelial-ligands","Treg genes in CD4 T cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
best_upstream_ligands=c("IL7","TGFB1","APOE")
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
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized ligands","Receptors expressed by CD4T", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
ligands_all = "IL7" 
targets_all = c("BATF","BATF","IL2RA","CTLA4")
active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network) 
# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
output_path = workDir
write_output = TRUE # change to TRUE for writing output
# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}
# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}
# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")
if(write_output){
bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}
