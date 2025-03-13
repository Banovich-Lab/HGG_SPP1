#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2023/07/16
# Description: Plotting infercnv result
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

.libPaths("/home/hnatri/R/library")
#.libPaths("/home/hnatri/R/4.1_libs")

library(Seurat)
library(ggplot2)
#library(stringi)
library(tidyverse)
library(ggpubr)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Parsing results
#==============================================================================#

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")

res_list <- lapply(unique(tumors$celltype), function(ct){
  message(ct)
  readRDS(paste0("/tgen_labs/banovich/BCTCSF/Heini/infercnv/tumors_", ct, "_infercnv_HMMi3_samplelevel_seurat.rds"))
})

plot_features_loss <- colnames(res_list[[1]]@meta.data)[grep("has_loss", colnames(res_list[[1]]@meta.data))]
plot_features_dupli <- colnames(res_list[[1]]@meta.data)[grep("has_dupli", colnames(res_list[[1]]@meta.data))]
# Proportion loss and duplication
featureplot_features_loss <- colnames(res_list[[1]]@meta.data)[grep("proportion_loss", colnames(res_list[[1]]@meta.data))]
featureplot_features_dupli <- colnames(res_list[[1]]@meta.data)[grep("proportion_dupli", colnames(res_list[[1]]@meta.data))]
# Proportion of all CNVs
featureplot_features_cnv <- colnames(res_list[[1]]@meta.data)[grep("proportion_cnv_", colnames(res_list[[1]]@meta.data))]
# Only looking at chrs 1-22
featureplot_features_cnv <- featureplot_features_cnv[grep("chr", featureplot_features_cnv)]
featureplot_features_loss <- featureplot_features_loss[grep("chr", featureplot_features_loss)]
featureplot_features_dupli <- featureplot_features_dupli[grep("chr", featureplot_features_dupli)]
feature_plot_chrs <- c(featureplot_features_dupli, featureplot_features_loss)

# Getting the average proportion of CNVs in each cell
res_list <- lapply(res_list, function(xx){
  xx$average_prop_cnv <- rowMeans(xx@meta.data[,featureplot_features_cnv])
  xx$average_prop_dupli <- rowMeans(xx@meta.data[,featureplot_features_dupli])
  xx$average_prop_loss <- rowMeans(xx@meta.data[,featureplot_features_loss])
  
  xx
})

# Adding info to the complete Seurat object
# Adjust cell barcodes to match
#substring(rownames(tumors@meta.data), 24, 24) <- "-"
rownames(tumors@meta.data) <- gsub("-", "_", rownames(tumors@meta.data))
colnames(tumors@assays$RNA@counts) <- rownames(tumors@meta.data)
colnames(tumors@assays$RNA@data) <- rownames(tumors@meta.data)

tumors[["integrated_sct_umap"]] <- RenameCells(
  object = tumors[["integrated_sct_umap"]],
  new.names = colnames(tumors@assays$RNA@counts))

tumors[["integrated_sct_pca"]] <- RenameCells(
  object = tumors[["integrated_sct_pca"]],
  new.names = colnames(tumors@assays$RNA@counts))

tumors[["integrated_sct"]] <- RenameCells(
  object = tumors[["integrated_sct"]],
  new.names = colnames(tumors@assays$RNA@counts))

metadata_list <- lapply(res_list, function(xx){
  xx@meta.data[,c("celltype", "average_prop_cnv", "average_prop_dupli", "average_prop_loss", featureplot_features_cnv)]
})
metadata <- do.call("rbind", metadata_list)

rownames(metadata) <- gsub("-", "_", rownames(metadata))

tumors$average_prop_cnv <- plyr::mapvalues(colnames(tumors),
                                           from = rownames(metadata),
                                           to = metadata$average_prop_cnv)

tumors$average_prop_cnv <- as.numeric(tumors$average_prop_cnv)

for(feature in featureplot_features_cnv){
  tumors@meta.data[,feature] <- plyr::mapvalues(colnames(tumors),
                                             from = rownames(metadata),
                                             to = metadata[,feature])
  tumors@meta.data[-which(colnames(tumors) %in% rownames(metadata)), feature] <- NA
  tumors@meta.data[,feature] <- as.numeric(tumors@meta.data[,feature])
}

saveRDS(tumors, "/labs/banovich/BCTCSF/Heini/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_inferCNV.rds")

#==============================================================================#
# Plotting
#==============================================================================#

Idents(tumors) <- tumors$cluster
selected <- WhichCells(tumors, idents = 10)
tumor_subset <- subset(tumors, cells = selected)

VlnPlot(tumors, features = "average_prop_cnv",
        group.by = "cluster", pt.size = 0, cols = tumor_cluster_col) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

VlnPlot(tumor_subset, features = "average_prop_cnv",
        group.by = "binary_response") &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

DotPlot(tumors, features = "average_prop_cnv",
        group.by = "UPN") &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

DotPlot(tumor_subset, features = featureplot_features_cnv,
        group.by = "UPN") &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

DotPlot(tumor_subset, features = featureplot_features_cnv,
        group.by = "UPN") &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

VlnPlot(subset(tumors, subset = binary_response %in% c("PD", "CR_SD")), features = "average_prop_cnv",
        group.by = "cluster", split.by = "binary_response", pt.size = 0) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

FeaturePlot(tumors, features = "average_prop_cnv") &
  coord_fixed(ratio=1) &
  theme_bw()

FeaturePlot(tumors, features = "average_prop_cnv",
            split.by = "cluster") &
  coord_fixed(ratio=1) &
  theme_bw()

plot_data <- tumors@meta.data
plot_data <- plot_data[which(plot_data$binary_response %in% c("PD", "CR_SD")),]

ggplot(plot_data, aes(x = binary_response, y = average_prop_cnv, fill = binary_response)) +
  facet_wrap(~plot_data$cluster) +
  geom_boxplot(alpha=0.8) +
  ylim(0, 1) +
  theme_bw() +
  geom_point() +
  stat_compare_means(comparisons = list(c("PD", "CR_SD")), label="p.format", method="wilcox.test")

plot_list <- lapply(featureplot_features_cnv, function(xx){
  ggplot(plot_data, aes_string(x = "binary_response", y = xx, fill = "binary_response")) +
    facet_wrap(~plot_data$cluster, ncol = 5) +
    geom_boxplot(alpha=0.8) +
    ylim(0, 1) +
    theme_bw() +
    geom_point() +
    stat_compare_means(comparisons = list(c("PD", "CR_SD")), label="p.format", method="wilcox.test") +
    ggtitle(xx)
})

#patchwork::wrap_plots(plot_list, nrow = 22)

filename <- "/home/hnatri/CART/Tumors/violinplot_propcnv_by_chr_cluster_response.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 8,
    height = 6)

plot_list

dev.off()

