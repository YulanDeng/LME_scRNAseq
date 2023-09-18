## Script for the transcriptional trajectories analysis
##Yulan Deng, last updated 2023-9-18
##my e-mail:kndeajs@163.com

################################################################
#monocle3 analysis of monocyte and macrophage (R version 4.1.1)#
################################################################
#${workDir} is the working directory
#${monoFile} is the seraut object of monocyte and macrophage

#Load required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(magrittr)

#Set working directory
setwd(workDir)

#Load the seurat object of monocyte and macrophage
moMacCellState <- readRDS(file =monoFile)

#transfer the seurat obeject to monocle object
cds <- as.cell_data_set(moMacCellState)
cds <- cluster_cells(cds,k=50)

p1 <- plot_cells(cds, show_trajectory_graph = FALSE, label_cell_groups = TRUE)
p2 <- plot_cells(cds, color_cells_by = "cellState", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
wrap_plots(p1, p2)

#tissue-resident Macrophage and doublets
cds_P1 <- subset(as.Seurat(cds, assay = NULL), monocle3_clusters %in% c(1,2,4,7))

#trajectories analysis
cds <- as.cell_data_set(cds_P1)
cds <- learn_graph(cds,close_loop = FALSE)

plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = TRUE, group_label_size = 6, label_leaves = FALSE, label_branch_points = FALSE)

#set monocyte as root
cds <- order_cells(cds, root_cells = rownames(pData(cds))[as.character(pData(cds)[,"monocle3_clusters"])=="7"])

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)

moMac_P1 <- as.Seurat(cds, assay = NULL)

FeaturePlot(moMac_P1, "monocle3_pseudotime")

plot_cells(cds, color_cells_by = "cellState", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)

#####################################################
#monocle3 analysis of CD8+ T cells (R version 4.1.1)#
#####################################################
#${workDir} is the working directory
#${CD8Tfile} is the seraut object of CD8+ T cells

#Load required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(magrittr)

#Set working directory
setwd(workDir)

#Load the seurat object of CD8+ T cells
CD8TCellState <- readRDS(file = CD8Tfile)

#transfer the seurat obeject to monocle object
cds <- as.cell_data_set(CD8TCellState)
cds <- cluster_cells(cds,k=50)

p1 <- plot_cells(cds, show_trajectory_graph = FALSE, label_cell_groups = TRUE)
p2 <- plot_cells(cds, color_cells_by = "cellState", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
wrap_plots(p1, p2)

plot_cells(cds, show_trajectory_graph = FALSE,color_cells_by = "partition", label_cell_groups = TRUE)

#exclude outliers
cds_P1 <- subset(as.Seurat(cds, assay = NULL),  monocle3_clusters %in% c(1,2,3,4,5,6,7))

#trajectories analysis
cds <- as.cell_data_set(cds_P1)
cds <- learn_graph(cds,close_loop = FALSE)

plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

#set naive CD8 T cells as root
cds <- order_cells(cds, root_cells = rownames(pData(cds))[as.character(pData(cds)[,"monocle3_clusters"])=="4"])

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)

CD8T_P1 <- as.Seurat(cds, assay = NULL)

FeaturePlot(CD8T_P1, "monocle3_pseudotime")

plot_cells(cds, color_cells_by = "cellState", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)

#####################################################
#monocle3 analysis of CD4+ T cells (R version 4.1.1)#
#####################################################
#${workDir} is the working directory
#${CD4Tfile} is the seraut object of CD4+ T cells

#Load required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(magrittr)

#Set working directory
setwd(workDir)

#Load the seurat object and cell state of CD4+ T cells
CD4TCellState <- readRDS(file = CD4Tfile)

#transfer the seurat obeject to monocle object
cds <- as.cell_data_set(CD4TCellState)
cds <- cluster_cells(cds)

p1 <- plot_cells(cds, show_trajectory_graph = FALSE, label_cell_groups = TRUE)
p2 <- plot_cells(cds, color_cells_by = "cellState", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
wrap_plots(p1, p2)

plot_cells(cds, show_trajectory_graph = FALSE,color_cells_by = "partition", label_cell_groups = TRUE)

cds_P1 <- subset(as.Seurat(cds, assay = NULL),  monocle3_partitions == 1)

cds <- as.cell_data_set(cds_P1)
cds <- learn_graph(cds,close_loop = FALSE)

plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

CD4T_P1 <- as.Seurat(cds, assay = NULL)


plot_cells(cds, label_groups_by_cluster = TRUE, group_label_size = 6, label_leaves = FALSE, label_branch_points = FALSE)

#exclude doublets
cds <- order_cells(cds, root_cells = rownames(pData(cds))[as.character(pData(cds)[,"monocle3_clusters"])%in%c("23","6","18","11","24","5","14","22")])

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)


CD4T_P1 <- as.Seurat(cds, assay = NULL)

FeaturePlot(CD4T_P1, "monocle3_pseudotime")

plot_cells(cds, color_cells_by = "cellState", label_cell_groups = FALSE, label_leaves = FALSE, 
    label_branch_points = FALSE,label_roots = FALSE)

