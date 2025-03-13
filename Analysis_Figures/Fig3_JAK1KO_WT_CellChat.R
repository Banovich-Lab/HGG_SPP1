#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2022/08/31
# Description: CellChat analysis of the JAK1KO/WT YUMM CD45+ data
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
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/Github/Analysis_Figures/CellChat_heatmap.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Importing data and creating inputs
#==============================================================================#

#seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/SPP1_release/JAK1KO_WildType_Control_Cells.rds")
# Not including epithelial cells
# Now /tgen_labs/banovich/BCTCSF/Heini/JAK1KO_scRNAseq/JAK1KO_WT_CTRL_immune_stroma.rds
seurat_object <- readRDS("/scratch/hnatri/CART/ko_wt_data.rds")

#seurat_object <- subset(seurat_object, subset = CD45_Group == "CD45_Pos")
seurat_object <- subset(seurat_object, subset = CART_Group == "Control")

seurat_object$celltype <- factor(seurat_object$celltype, levels = sort(as.character(unique(seurat_object$celltype))))
Idents(seurat_object) <- seurat_object$celltype

seurat_object$celltype <- factor(seurat_object$celltype, levels = sort(unique(seurat_object$celltype)))
Idents(seurat_object) <- seurat_object$celltype

# Inputs for CellChat:
# Log-normalized counts and cell labels (cluster/celltype)
# Creating CC objects for each group, merging for comparative analysis, as
# well as an object for all samples together
all_counts <- LayerData(seurat_object, assay = "RNA", layer = "data")
all_CC <- createCellChat(object = all_counts,
                                meta = seurat_object@meta.data,
                                group.by = "celltype")

# Creating CC objects
WT_seurat_object <- subset(seurat_object, subset = JAK1_Group == "WildType")
WT_counts <- LayerData(WT_seurat_object, assay = "RNA",
                             layer = "data")
WT_CC <- createCellChat(object = WT_counts, meta = WT_seurat_object@meta.data,
                        group.by = "celltype")

JAK1KO_seurat_object <- subset(seurat_object, subset = JAK1_Group == "Jak1_KO")
JAK1KO_counts <- LayerData(JAK1KO_seurat_object, assay = "RNA",
                                     layer = "data")
JAK1KO_CC <- createCellChat(object = JAK1KO_counts,
                            meta = JAK1KO_seurat_object@meta.data,
                            group.by = "celltype")

# Adding metadata
WT_CC <- addMeta(WT_CC, meta = WT_seurat_object@meta.data)
WT_CC <- setIdent(WT_CC, ident.use = "celltype")
# Dropping an unused ident
JAK1KO_CC <- addMeta(JAK1KO_CC, meta = JAK1KO_seurat_object@meta.data)
JAK1KO_CC <- setIdent(JAK1KO_CC, ident.use = "celltype")

all_CC <- addMeta(all_CC, meta = seurat_object@meta.data)
all_CC <- setIdent(all_CC, ident.use = "celltype")

#==============================================================================#
# CellChat analysis
#==============================================================================#

# Setting the ligand-receptor database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

CC_objects <- list("all_CC" = all_CC,
                   "WT_CC" = WT_CC,
                   "JAK1KO_CC" = JAK1KO_CC)

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
CC_objects_compare <- CC_objects[c("WT_CC", "JAK1KO_CC")]
CC_merged_object <- mergeCellChat(CC_objects_compare, add.names = names(CC_objects_compare))

saveRDS(CC_objects_compare, "/scratch/hnatri/CART/JAK1KO_WT_CTRL_CC_compare_object.rds")
saveRDS(CC_merged_object, "/scratch/hnatri/CART/JAK1KO_WT_CTRL_CC_merged_object.rds")
#q(save="no")
#
#CC_objects_compare <- readRDS("/scratch/hnatri/CART/JAK1KO_WT_CTRL_CC_compare_object.rds")
#CC_merged_object <- readRDS("/scratch/hnatri/CART/JAK1KO_WT_CTRL_CC_merged_object.rds")

#CC_objects_compare$WT_CC@meta$celltype <- factor(CC_objects_compare$WT_CC@meta$celltype,
#                                                    levels = levels(CC_objects_compare$JAK1KO_CC@meta$celltype))

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

#CC_merged_object@meta$cluster_cc <- factor(CC_merged_object@meta$cluster_cc, levels = paste0("C", c(0, seq(1, 13, by=1))))

# Heatmaps
Heatmap <- ComplexHeatmap::Heatmap

gg1 <- netVisual_heatmap(CC_merged_object,
                         cluster.rows = T,
                         cluster.cols = T,
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         color.use = jak1_celltype_col)
gg2 <- netVisual_heatmap(CC_merged_object,
                         measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T,
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         color.use = jak1_celltype_col)
gg1 + gg2

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_JAK1KO_WT_heatmaps.pdf"
pdf(file = filename,
    width = 8,
    height = 4)

gg1 + gg2

dev.off()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_JAK1KO_WT_nInteractions_heatmap.pdf"
pdf(file = filename,
    width = 4.5,
    height = 4)

gg1

dev.off()

# Identify and visualize the conserved and context-specific signaling pathways:
# signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv)
# increase, by change their information flow at one condition as compared to
# another condition.
# The top signaling pathways colored red are enriched in CRSD, and these colored green were enriched in PD
gg1 <- rankNet(CC_merged_object, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(CC_merged_object, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#mm_CC_merged_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/CellChat/mm_CC_merged_object.rds")

# Extracting data from the plot
mm <- rankNet(CC_merged_object,
              mode = "comparison",
              stacked = F,
              do.stat = T,
              cutoff.pvalue = 0.01,
              #tol = 1,
              thresh = 0.01)
mm_data <- ggplot_build(mm)
info_flow_data <- mm_data$plot$data

# Constructing the plot
info_flow_data$log2_contribution_1 <- log2(info_flow_data$contribution+1)
#info_flow_data$group <- ifelse(info_flow_data$group == "Jak1_KO", "Jak1 KO", "WT")

barplot_colors <- c("WT_CC" = "aquamarine4", "JAK1KO_CC" = "deeppink3")

hist(info_flow_data$contribution, breaks = 55, main = "", xlab = "Information flow")
abline(v=1, col = "red")

dev.off()

info_flow_data %>%
  ggplot(aes(x = reorder(name, -log2_contribution_1), y = log2_contribution_1, fill = group)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors) +
  theme_bw() + 
  xlab("") + 
  ylab(bquote(log[2]("Information flow + 1"))) +
  labs(fill = "Jak1 status") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9))

barplot_colors <- c("aquamarine3", "deeppink3")

barp <- info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = JAK1KO_CC-WT_CC) %>%
  dplyr::select(name, delta) %>%
  #mutate(log2_delta = log2(delta)) %>%
  mutate(sign = delta > 0) %>%
  filter(abs(delta) > 1) %>%
  ggplot(aes(x = reorder(name, -delta), y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("WT", "JAK1KO")) +
  theme_bw() + 
  xlab("") + 
  #ylab(bquote(log[2]("Information flow + 1"))) +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "Jak1 status") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
  manuscript_theme

barp

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CellChat_pathways.pdf"
pdf(file = filename,
    width = 6,
    height = 2)

barp

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
                   color.use = jak1_celltype_col)
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
                                               color.use = jak1_celltype_col)
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

# Heatmap based on a single object 
pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)

ht_JAK1KO <- netVisual_heatmap(CC_objects_compare[["JAK1KO_CC"]],
                             signaling = pathways.show,
                             color.heatmap = "Reds",
                             #cluster.rows = T,
                             #cluster.cols = T,
                             title.name = paste(pathways.show, "signaling in JAK1KO"),
                             color.use = jak1_celltype_col)

ht_WT <- netVisual_heatmap(CC_objects_compare[["WT_CC"]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling in JAK1KO"),
                               color.use = jak1_celltype_col)

ht_JAK1KO

summary(CC_objects_compare[["WT_CC"]]@data.signaling["Spp1",])
hist(CC_objects_compare[["WT_CC"]]@data.signaling["Spp1",])
hist(CC_objects_compare[["JAK1KO_CC"]]@data.signaling["Spp1",])

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_CellChat_SPP1_heatmap.pdf"
pdf(file = filename,
    width = 5,
    height = 4)

ht_JAK1KO

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
                                                  levels = sort(as.character(unique(CC_objects_compare[[i]]@meta$celltype))))
  CC_objects_compare[[i]]@idents <- factor(CC_objects_compare[[i]]@idents,
                                           levels = sort(as.character(unique(CC_objects_compare[[i]]@idents))))
  
  
  ht[[i]] <- netVisual_heatmap(object = CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ", names(CC_objects_compare)[i]),
                               color.use = jak1_celltype_col)
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

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

