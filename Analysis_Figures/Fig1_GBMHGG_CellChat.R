#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2022/08/31
# Description: CellChat analysis on 13384 CAR T trial tumor samples
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(CellChat)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/Analysis_Figures/CellChat_heatmap.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Importing data and creating inputs
#==============================================================================#

#seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

# Unable to include M6 in the comparative analysis for CD3
seurat_object <- subset(seurat_object, subset = celltype == "M6", invert = T)

# M6 and L9 are only in SPP1 surv bottom
#seurat_object <- subset(seurat_object, subset = celltype %in% c("M6", "L9"), invert = T)

table(seurat_object$celltype)

#seurat_object$celltype <- factor(seurat_object$celltype,
#                                 levels = c(paste0("M", seq(1, 9, by = 1)),
#                                            "B1", "N1",
#                                            paste0("L", seq(1, 10, by = 1)),
#                                            paste0("F", seq(1, 3, by = 1))))

setdiff(Idents(seurat_object), unique(seurat_object$celltype))

seurat_object$celltype <- factor(seurat_object$celltype,
                                 levels = sort(unique(seurat_object$celltype)))
Idents(seurat_object) <- seurat_object$celltype

# Adding a variable
seurat_object$SPP1_surv_extremes <- ifelse(seurat_object$UPN %in% c(129, 146, 185, 224, 228, 234, 237, 248), "bottom",
                                           ifelse(seurat_object$UPN %in% c(109, 141, 181, 223, 265, 301), "top", NA))

# Inputs for CellChat:
# Log-normalized counts and cell labels (cluster/celltype)
# Creating CC objects for CD3_high/low, merging for comparative analysis, as
# well as an object for all samples together
all_counts <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
all_CC <- createCellChat(object = all_counts,
                                meta = seurat_object@meta.data,
                                group.by = "celltype")

## CD3_high
CD3high_seurat_object <- subset(seurat_object, subset = CD3_high == "TRUE")
CD3high_counts <- GetAssayData(CD3high_seurat_object, assay = "RNA",
                               slot = "data")
CD3high_CC <- createCellChat(object = CD3high_counts,
                             meta = CD3high_seurat_object@meta.data,
                             group.by = "celltype")

CD3low_seurat_object <- subset(seurat_object, subset = CD3_high == "FALSE")
CD3low_counts <- GetAssayData(CD3low_seurat_object, assay = "RNA",
                              slot = "data")
CD3low_CC <- createCellChat(object = CD3low_counts,
                            meta = CD3low_seurat_object@meta.data,
                            group.by = "celltype")

setdiff(levels(CD3low_seurat_object$celltype), unique(CD3high_seurat_object$celltype))

# Adding metadata
CD3high_CC <- addMeta(CD3high_CC, meta = CD3high_seurat_object@meta.data)
CD3high_CC <- setIdent(CD3high_CC, ident.use = "celltype")
# Dropping an unused ident
CD3high_CC@idents <- factor(CD3high_CC@idents, levels = unique(CD3high_CC@idents))
CD3low_CC <- addMeta(CD3low_CC, meta = CD3low_seurat_object@meta.data)
CD3low_CC <- setIdent(CD3low_CC, ident.use = "celltype")

#all_CC <- addMeta(all_CC, meta = seurat_object@meta.data)
#all_CC <- setIdent(all_CC, ident.use = "celltype")

setdiff(levels(CD3high_CC@idents), unique(CD3high_CC@idents))

# Creating CC objects for response/non-response
# CR_SD
CR_SD_seurat_object <- subset(seurat_object, subset = binary_response == "CR_SD")
CR_SD_counts <- GetAssayData(CR_SD_seurat_object, assay = "RNA",
                             slot = "data")
CR_SD_CC <- createCellChat(object = CR_SD_counts, meta = CR_SD_seurat_object@meta.data, group.by = "celltype")

PD_seurat_object <- subset(seurat_object, subset = binary_response == "PD")
PD_counts <- GetAssayData(PD_seurat_object, assay = "RNA",
                                     slot = "data")
PD_CC <- createCellChat(object = PD_counts,
                        meta = PD_seurat_object@meta.data,
                        group.by = "celltype")

# Adding metadata
CR_SD_CC <- addMeta(CR_SD_CC, meta = CR_SD_seurat_object@meta.data)
CR_SD_CC <- setIdent(CR_SD_CC, ident.use = "celltype")
# Dropping an unused ident
PD_CC <- addMeta(PD_CC, meta = PD_seurat_object@meta.data)
PD_CC <- setIdent(PD_CC, ident.use = "celltype")

setdiff(unique(PD_seurat_object$celltype), unique(CR_SD_seurat_object$celltype))
setdiff(unique(CR_SD_seurat_object$celltype), unique(PD_seurat_object$celltype))


setdiff(levels(all_CC@idents), unique(CR_SD_CC@idents))
setdiff(levels(all_CC@idents), unique(PD_CC@idents))
setdiff(unique(PD_CC@idents), unique(CR_SD_CC@idents))
setdiff(unique(CR_SD_CC@idents), unique(PD_CC@idents))

setdiff(levels(PD_CC@idents), unique(PD_CC@idents))
setdiff(levels(CR_SD_CC@idents), unique(CR_SD_CC@idents))

## For SPP1/survival extremes
SPP1survtop_seurat_object <- subset(seurat_object, subset = SPP1_surv_extremes == "top")
SPP1survtop_counts <- GetAssayData(SPP1survtop_seurat_object, assay = "RNA",
                                   slot = "data")
SPP1survtop_CC <- createCellChat(object = SPP1survtop_counts,
                                 meta = SPP1survtop_seurat_object@meta.data,
                                 group.by = "celltype")

SPP1survbottom_seurat_object <- subset(seurat_object, subset = SPP1_surv_extremes == "bottom")
SPP1survbottom_counts <- GetAssayData(SPP1survbottom_seurat_object, assay = "RNA",
                                      slot = "data")
SPP1survbottom_CC <- createCellChat(object = SPP1survbottom_counts,
                                    meta = SPP1survbottom_seurat_object@meta.data,
                                    group.by = "celltype")

# Adding metadata
SPP1survtop_CC <- addMeta(SPP1survtop_CC, meta = SPP1survtop_seurat_object@meta.data)
SPP1survtop_CC <- setIdent(SPP1survtop_CC, ident.use = "celltype")

SPP1survbottom_CC <- addMeta(SPP1survbottom_CC, meta = SPP1survbottom_seurat_object@meta.data)
SPP1survbottom_CC <- setIdent(SPP1survbottom_CC, ident.use = "celltype")

setdiff(unique(SPP1survbottom_seurat_object$celltype), unique(SPP1survtop_seurat_object$celltype))
setdiff(unique(SPP1survtop_seurat_object$celltype), unique(SPP1survbottom_seurat_object$celltype))

# Dropping an unused ident
SPP1survtop_CC@idents <- factor(SPP1survtop_CC@idents, levels = unique(SPP1survtop_CC@idents))
SPP1survbottom_CC <- addMeta(SPP1survbottom_CC, meta = SPP1survbottom_seurat_object@meta.data)
SPP1survbottom_CC <- setIdent(SPP1survbottom_CC, ident.use = "celltype")

setdiff(levels(SPP1survbottom_CC@idents), unique(SPP1survtop_CC@idents))

all_CC <- addMeta(all_CC, meta = seurat_object@meta.data)
all_CC <- setIdent(all_CC, ident.use = "celltype")

#==============================================================================#
# CellChat analysis
#==============================================================================#

# Setting the ligand-receptor database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

CC_objects <- list("all_CC" = all_CC,
                   "CD3high_CC" = CD3high_CC,
                   "CD3low_CC" = CD3low_CC)

#CC_objects <- list("all_CC" = all_CC,
#                   "CR_SD_CC" = CR_SD_CC,
#                   "PD_CC" = PD_CC)

#CC_objects <- list("all_CC" = all_CC,
#                   "SPP1survtop_CC" = SPP1survtop_CC,
#                   "SPP1survbottom_CC" = SPP1survbottom_CC)

CC_objects <- lapply(CC_objects, function(xx){
  # Set the used database in the object
  xx@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # Subset the expression data of signaling genes for saving computation cost
  xx <- subsetData(xx) # This step is necessary even if using the whole database
  
  xx <- identifyOverExpressedGenes(xx)
  xx <- identifyOverExpressedInteractions(xx)
  
  # Project gene expression data onto PPI (Optional: when running it, USER should
  # set `raw.use = FALSE` in the function `computeCommunProb()` in order to use
  # the projected data)
  xx <- projectData(xx, PPI.human)
  
  ## Inference of cell-cell communication network
  # Compute the communication probability and infer cellular communication network
  xx <- computeCommunProb(xx) # raw.use = FALSE
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  xx <- filterCommunication(xx, min.cells = 10)
  
  # Infer the cell-cell communication at a signaling pathway level
  unique(xx@idents)
  
  xx <- computeCommunProbPathway(xx)
  
  # Calculate the aggregated cell-cell communication network
  xx <- aggregateNet(xx)
  
  xx
})

# Merging objects for comparative analysis
#CC_objects_compare <- CC_objects[c("CR_SD_CC", "PD_CC")]
CC_objects_compare <- CC_objects[c("CD3high_CC", "CD3low_CC")]
#CC_objects_compare <- CC_objects[c("SPP1survtop_CC", "SPP1survbottom_CC")]
CC_merged_object <- mergeCellChat(CC_objects_compare, add.names = names(CC_objects_compare))

#saveRDS(CC_objects_compare, "/scratch/hnatri/CART/immune_stroma_CC_compare_object_CD3_celltypenames.rds")
#saveRDS(CC_merged_object, "/scratch/hnatri/CART/immune_stroma_CC_merged_object_CD3_celltypenames.rds")
saveRDS(CC_objects_compare, "/scratch/hnatri/CART/immune_stroma_CC_compare_object_response_celltypenames.rds")
saveRDS(CC_merged_object, "/scratch/hnatri/CART/immune_stroma_CC_merged_object_response_celltypenames.rds")
#saveRDS(CC_objects_compare, "/scratch/hnatri/CART/tumors_CC_compare_object_response.rds")
#saveRDS(CC_merged_object, "/scratch/hnatri/CART/tumors_CC_merged_object_response.rds")

q(save="no")

#==============================================================================#
# Inspecting the results
#==============================================================================#

#CC_objects_compare <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_compare_object_CD3_celltypenames.rds")
#CC_merged_object <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_merged_object_CD3_celltypenames.rds")#
#CC_objects_compare <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_compare_object_SPP1survtop.rds
#CC_merged_object <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_merged_object_SPP1survtop.rds
CC_objects_compare <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/CellChat/immune_stroma_CC_compare_object_response_celltypenames.rds")
CC_merged_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/CellChat/immune_stroma_CC_merged_object_response_celltypenames.rds")
#CC_objects_compare <- readRDS("/scratch/hnatri/CART/tumors_CC_compare_object_response.rds")
#CC_merged_object <- readRDS("/scratch/hnatri/CART/tumors_CC_merged_object_response.rds")

plot_levels <- c(paste0("M", seq(1, 9, by = 1)),
                 "B1", "N1",
                 paste0("L", seq(1, 10, by = 1)),
                 paste0("F", seq(1, 3, by = 1)))

CC_objects_compare$CR_SD_CC@meta$celltype <- factor(CC_objects_compare$CR_SD_CC@meta$celltype,
                                                    levels = plot_levels)
CC_objects_compare$PD_CC@meta$celltype <- factor(CC_objects_compare$PD_CC@meta$celltype,
                                                    levels = plot_levels)

#CC_objects_compare$CD3high_CC@meta$celltype <- factor(CC_objects_compare$CD3high_CC@meta$celltype,
#                                                    levels = plot_levels)
#CC_objects_compare$CD3low_CC@meta$celltype <- factor(CC_objects_compare$CD3low_CC@meta$celltype,
#                                                 levels = plot_levels)

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(CC_merged_object, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(CC_merged_object, show.legend = F, group = c(1,2),
                           measure = "weight")
gg1 + gg2

# Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(CC_merged_object, weight.scale = T)
netVisual_diffInteraction(CC_merged_object, weight.scale = T,
                          measure = "weight")

# Heatmaps
Heatmap <- ComplexHeatmap::Heatmap

gg1 <- netVisual_heatmap(CC_merged_object,
                         cluster.rows = T,
                         cluster.cols = T,
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         color.use = immune_fibro_celltype_col)
gg2 <- netVisual_heatmap(CC_merged_object,
                         measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T,
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         color.use = immune_fibro_celltype_col)
gg1 + gg2

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_CD3_heatmaps.pdf"
pdf(file = filename,
    width = 8,
    height = 4)

gg1 + gg2

dev.off()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_nInteractions_heatmap.pdf"
pdf(file = filename,
    width = 4.5,
    height = 4)

gg1

dev.off()

# To better control the node size and edge weights of the inferred networks
# across different datasets, we compute the maximum number of cells per cell
# group and the maximum number of interactions (or interaction weights) across all datasets.
weight.max <- getMaxWeight(CC_objects_compare, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_circle(CC_objects_compare[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(CC_objects_compare)[i]),
                   color.use = immune_fibro_celltype_col)
}

# Need to run netAnalysis_computeCentrality before netAnalysis_signalingRole_scatter

# Compute the network centrality scores
CC_objects_compare <- lapply(CC_objects_compare, function(xx){
  netAnalysis_computeCentrality(xx, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
})

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#netAnalysis_signalingRole_network(CC_objects_compare[[1]], signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Compare the major sources and targets in 2D space
num.link <- sapply(CC_objects_compare, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(CC_objects_compare)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(CC_objects_compare[[i]],
                                               title = names(CC_objects_compare)[i],
                                               weight.MinMax = weight.MinMax,
                                               color.use = immune_fibro_celltype_col)
}

patchwork::wrap_plots(plots = gg)

# Identify signaling changes associated with one cell group
# Visualizing differential outgoing and incoming signaling changes

# Identify signaling groups based on their functional similarity
CC_merged_object <- computeNetSimilarityPairwise(CC_merged_object, type = "functional")
CC_merged_object <- netEmbedding(CC_merged_object, type = "functional")
CC_merged_object <- netClustering(CC_merged_object, type = "functional")

# Visualization in 2D-space
netVisual_embeddingPairwise(CC_merged_object, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
CC_merged_object <- computeNetSimilarityPairwise(CC_merged_object, type = "structural")
CC_merged_object <- netEmbedding(CC_merged_object, type = "structural")
CC_merged_object <- netClustering(CC_merged_object, type = "structural")

# Visualization in 2D-space
netVisual_embeddingPairwise(CC_merged_object, type = "structural", label.size = 3.5)

# We can identify the signaling networks with larger (or less) difference based
# on their Euclidean distance in the shared two-dimensions space. Larger distance
# implies larger difference of the communication networks between two datasets
# in terms of either functional or structure similarity.
rankSimilarity(CC_merged_object, type = "functional")
rankSimilarity(CC_merged_object, type = "functional", comparison1 = c(1,2))

# Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1_CRSD_circle.pdf"
pdf(file = filename,
    width = 8,
    height = 8)

netVisual_aggregate(CC_objects_compare[[1]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[1]))

dev.off()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1_PD_circle.pdf"
pdf(file = filename,
    width = 8,
    height = 8)

netVisual_aggregate(CC_objects_compare[[2]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[2]))

dev.off()

# Heatmap based on a single object 
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  ht[[i]] <- netVisual_heatmap(CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ",
                                                  names(CC_objects_compare)[i]))
  #color.use = tumor_celltype_col)
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  ht[[i]] <- netVisual_heatmap(CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ", names(CC_objects_compare)[i]))
  #color.use = tumor_cluster_col_c[1:13])
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  CC_objects_compare[[i]]@meta$celltype <- factor(CC_objects_compare[[i]]@meta$celltype,
                                                  levels = plot_levels)
  CC_objects_compare[[i]]@idents <- factor(CC_objects_compare[[i]]@idents,
                                           levels = plot_levels)
  
  ht[[i]] <- netVisual_heatmap(object = CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = c('#ffffff','#b2182b'),
                               measure = "count",
                               slot.name = "netP",
                               #remove.isolate = T,
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ", names(CC_objects_compare)[i]),
                               color.use = immune_fibro_celltype_col)
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1_CRSD_PD_heatmap.pdf"
pdf(file = filename,
    width = 7,
    height = 3.5)

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

dev.off()

# Delta heatmap
#ht[[1]]@matrix - ht[[2]]@matrix

delta_ht <- ht[[2]]
delta_ht@matrix <- ht[[2]]@matrix - ht[[1]]@matrix

delta_ht@matrix_color_mapping

mx_cols <- c("#4287F5FF", "#ffffff", "#A84E76FF", "#B12A3DFF", "#B2182BFF")
mx_breaks <- c(-0.3, 0, 0.3, 0.6, 0.9)
col_fun = colorRamp2(mx_breaks, mx_cols)

delta_ht@matrix_color_mapping <- ComplexHeatmap::ColorMapping(col_fun = col_fun,
                                                              breaks = mx_breaks,
                                                              na_col = "#e6e7e8")

ComplexHeatmap::draw(delta_ht)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1_PD_CRSD_delta_heatmap.pdf"
pdf(file = filename,
    width = 4,
    height = 3.5)

ComplexHeatmap::draw(delta_ht)

dev.off()

# SPP1 signaling across whole dataset
#CC_merged_object@meta$celltype <- factor(CC_merged_object@meta$celltype,
#                                                levels = plot_levels)
#CC_merged_object@idents <- factor(CC_merged_object@idents,
#                                         levels = plot_levels)
#
#SPP1_combined_heatmap <- netVisual_heatmap(object = CC_merged_object,
#                                           signaling = pathways.show,
#                                           color.heatmap = c('#ffffff','#b2182b'),
#                                           measure = "count",
#                                           slot.name = "netP",
#                                           #comparison = NULL,
#                                           #remove.isolate = T,
#                                           #cluster.rows = T,
#                                           #cluster.cols = T,
#                                           title.name = paste(pathways.show, "signaling"),
#                                           color.use = immune_fibro_celltype_col)
#
#SPP1_combined_heatmap@matrix_color_mapping
#
#filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1_combined_heatmap.pdf"
#pdf(file = filename,
#    width = 4,
#    height = 3.5)
#
#ComplexHeatmap::draw(SPP1_combined_heatmap)
#
#dev.off()

# Chord diagram
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = tumor_celltype_col)
}

pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = tumor_celltype_col)
}

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = immune_fibro_celltype_col)
}

#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(CC_objects_compare[[1]]@idents)
pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  #group = group.cellType,
  netVisual_chord_cell(CC_objects_compare[[i]],
                       signaling = pathways.show,
                       title.name = paste0(pathways.show, " signaling network - ", names(CC_objects_compare)[i]),
                       color.use = tumor_celltype_col)
}

par(mfrow = c(1, 2), xpd=TRUE)


# Compare the signaling gene expression distribution between different datasets
#CC_merged_object@idents <- factor(factor(CC_merged_object@meta$cluster_cc, levels = paste0("C", c(0, seq(1, 13, by=1)))), levels = paste0("C", c(0, seq(1, 13, by=1))))
CC_merged_object@meta$binary_response = factor(CC_merged_object@meta$binary_response, levels = c("CR_SD", "PD")) # set factor level
plotGeneExpression(CC_merged_object,
                   signaling = "SPP1",
                   split.by = "binary_response",
                   colors.ggplot = T)
#color.use = c("red", "blue"))

plot_features <- c("SPP1", "CD44", "ITGAV", "ITGA4", "ITGA5", "ITGB1", "ITGB5")

VlnPlot(seurat_object,
        features = plot_features,
        group.by = "celltype",
        split.by = "binary_response",
        pt.size = 0,
        slot = "data") &
  theme_bw()

myel_fib <- c(paste0("M", seq(1, 10)), c("N1"), paste0("F", seq(1, 3)))

seurat_object_subset <- subset(seurat_object, subset = celltype %in% myel_fib)

VlnPlot(seurat_object_subset,
        features = plot_features,
        group.by = "celltype",
        split.by = "binary_response",
        pt.size = 0,
        slot = "data",
        ncol = 1,
        cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend() +
  manuscript_theme

p1 <- VlnPlot(seurat_object_subset,
              features = c("SPP1"),
              group.by = "celltype",
              split.by = "binary_response",
              pt.size = 0,
              slot = "data",
              ncol = 1,
              cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend()

seurat_object_fibrob <- subset(seurat_object, subset = celltype %in% c("Fibroblast"))
p2 <- VlnPlot(seurat_object_fibrob,
              features = c("CD44", "ITGAV"),
              group.by = "celltype",
              split.by = "binary_response",
              pt.size = 0,
              slot = "data",
              ncol = 2,
              cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend()

p1 / p2

# Violinplots of SPP1 receptors
receptor_vln <- VlnPlot(seurat_object,
                        features = c("CD44", "ITGAV", "ITGA4", "ITGA5", "ITGB1", "ITGB5"),
                        group.by = "celltype",
                        split.by = "binary_response",
                        pt.size = 0,
                        slot = "data",
                        ncol = 1,
                        cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend() +
  manuscript_theme

receptor_vln

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_immune_fibro_response_SPP1receptors_VlnPlot.pdf"
pdf(file = filename,
    width = 6,
    height = 8)

receptor_vln

dev.off()

# Barplot of pathways
# Extracting data from the plot
SPP1surv_pathways <- rankNet(CC_merged_object,
                             mode = "comparison",
                             stacked = F,
                             do.stat = T,
                             cutoff.pvalue = 0.01,
                             #tol = 1,
                             thresh = 0.01)

SPP1surv_pathways_data <- ggplot_build(SPP1surv_pathways)
SPP1surv_pathways_info_flow_data <- SPP1surv_pathways_data$plot$data

# Constructing the plot
SPP1surv_pathways_info_flow_data$log2_contribution_1 <- log2(SPP1surv_pathways_info_flow_data$contribution+1)
SPP1surv_pathways_info_flow_data$group <- ifelse(SPP1surv_pathways_info_flow_data$group == "SPP1survtop_CC", "SPP1survtop", "SPP1survbottom")

barplot_colors <- c("SPP1survtop" = "aquamarine4", "SPP1survbottom" = "deeppink3")

hist(SPP1surv_pathways_info_flow_data$contribution, breaks = 55, main = "", xlab = "Information flow")
abline(v=1, col = "red")

dev.off()

# Plotting delta
barplot_colors <- c("aquamarine3", "deeppink3")

# MK, VEGF, PlAU, GRN, LAIR1, PARs
SPP1surv_delta_plot <- SPP1surv_pathways_info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = SPP1survbottom-SPP1survtop) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0) %>% # summarise(mind = min(delta))
  filter(abs(delta) > 0.8) %>%
  ggplot(aes(x = reorder(name, -delta), y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("SPP1 low, longer survival", "SPP1 high, lower survival")) +
  theme_bw() + 
  xlab("") +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "SPP1/survival") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
  manuscript_theme
#scale_y_break(c(-15, -45))

SPP1surv_delta_plot

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_tumor_response_delta_infoflow.pdf"
pdf(file = filename,
    width = 5,
    height = 2)

SPP1surv_delta_plot

dev.off()


#==============================================================================#
# Saving the results
#==============================================================================#

CC_merged_object

head(CC_merged_object@netP$CR_SD_CC$prob)

# Using the netVisual_bubble function to pull all LR pair results
lr_all <- netVisual_bubble_edit(CC_merged_object,
                                sources.use = NULL, targets.use = NULL,
                                comparison = c(1, 2), angle.x = 45)


lr_all_data <- ggplot_build(lr_all)
lr_all_data <- lr_all_data$plot$data

write.table(lr_all_data, "/home/hnatri/13384_CART/13384_Tumors/CRSD_PD_CellChat_LR_all_res.tsv",
            quote = F, sep = "\t", row.names = F)
