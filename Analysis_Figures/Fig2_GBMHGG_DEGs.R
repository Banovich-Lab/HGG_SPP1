#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/30/2023
# Description: Cell type markers and DEGs from the GBM scRNA-seq
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(tidyverse)
library(googlesheets4)

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
reduction <- "integrated_sct_umap"

#==============================================================================#
# Import data
#==============================================================================#

# The final object
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#==============================================================================#
# Cell type markers
#==============================================================================#

# Top markers for each cluster
markers <- presto::wilcoxauc(tumors,
                             group_by = "celltype",
                             assay = "data",
                             seurat_assay = "RNA")

sig_markers <- markers %>% filter(auc>0.6, padj<0.01)
sig_markers <- sig_markers[order(sig_markers$auc, decreasing = T),]

sig_markers %>% filter(group == "Lymph1") %>% head(n = 20)
sig_markers %>% filter(group == "Lymph2") %>% head(n = 40)

#==============================================================================#
# DEGs
#==============================================================================#

DEG_response <- lapply(unique(immune_fibro$celltype), function(xx){
  #message(xx)
  data_subset <- subset(immune_fibro, subset = celltype == xx)
  Idents(data_subset) <- data_subset$binary_response
  
  table(data_subset$binary_response)
  if (min(table(data_subset$binary_response))<10){
    return(NULL)
  }
  
  if (all((c("CR_SD", "PD") %in% data_subset$binary_response) == c(T, T))){
    markers <- FindMarkers(data_subset,
                           ident.1 = "CR_SD",
                           ident.2 = "PD",
                           assay = "RNA",
                           verbose = F)
    markers$feature <- rownames(markers)
    markers$celltype <- xx
    
    return(markers)
  } else {
    return(NULL)
  }
})
names(DEG_response) <- unique(immune_fibro$celltype)
DEG_response[sapply(DEG_response, is.null)] <- NULL

DEG_response_df <- as.data.frame(do.call(rbind, DEG_response))

# Distribution of log2FC
hist(DEG_response_df$avg_log2FC)

DEG_response_df_sig <- DEG_response_df %>%
  filter(p_val_adj < 0.01,
         abs(avg_log2FC) > 2)

DEG_response_df %>% filter(feature == "SPP1") %>%
  arrange(celltype)

# Saving to a file
write.table(DEG_response_df_sig,
            "/scratch/hnatri/CART/13384_immune_fibro_response_DEGs_Seurat_sig.tsv",
            sep = "\t", quote = F, row.names = F)

# CD3 score
DEG_CD3 <- lapply(unique(immune_fibro$celltype), function(xx){
  message(xx)
  data_subset <- subset(immune_fibro, subset = celltype == xx)
  Idents(data_subset) <- data_subset$CD3_high_low
  
  table(data_subset$CD3_high_low)
  
  if (min(table(data_subset$CD3_high_low))<10){
    return(NULL)
  }
  
  if (min(table(data_subset$CD3_high_low)<10)){
    return(NULL)
  }
  
  if (all((c("High", "Low") %in% data_subset$CD3_high_low) == c(T, T))){
    markers <- FindMarkers(data_subset,
                           ident.1 = "High",
                           ident.2 = "Low",
                           assay = "RNA",
                           verbose = F)
    markers$feature <- rownames(markers)
    markers$celltype <- xx
    
    return(markers)
  } else {
    return(NULL)
  }
})
names(DEG_CD3) <- unique(immune_fibro$celltype)
DEG_CD3[sapply(DEG_CD3, is.null)] <- NULL

DEG_CD3_df <- as.data.frame(do.call(rbind, DEG_CD3))

# Distribution of log2FC
hist(DEG_CD3_df$avg_log2FC)

DEG_CD3_df_sig <- DEG_CD3_df %>%
  filter(p_val_adj < 0.01,
         abs(avg_log2FC) > 2)

DEG_CD3_df %>% filter(feature == "SPP1") %>%
  arrange(celltype)

# Saving to a file
write.table(DEG_CD3_df_sig,
            "/scratch/hnatri/CART/13384_immune_fibro_CD3high_DEGs_Seurat_sig.tsv",
            sep = "\t", quote = F, row.names = F)

