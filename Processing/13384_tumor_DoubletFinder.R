#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 8/23/2023
# Description: DoubletFinder on batches run without multiplexing
#==============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(googlesheets4)
library(dplyr)
library(DoubletFinder)

#==============================================================================
# Environment
#==============================================================================

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/SingleCellBestPractices/scripts/preprocessing_qc_module.R")

#==============================================================================
# Import data
#==============================================================================

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")

#==============================================================================
# Running DoubletFinder
#==============================================================================

seurat_list <- SplitObject(tumors, split.by = "Batch")
seurat_subset <- seurat_list[grepl("F0", names(seurat_list))]

# Skipping samples with <50 cells
F04739 <- seurat_list["F04739"]
seurat_subset["F04739"] <- NULL

# F04739
seurat_subset_db <- run_doubletfinder(seurat_subset)

# Saving cell names
cellinfo <- lapply(seurat_subset_db, function(xx){
  xx@meta.data %>% 
    rownames_to_column("cellname") %>%
    dplyr::select("doublet_finder", "cellname")
})
cellinfo <- do.call("rbind", cellinfo)

write.table(cellinfo, "/scratch/hnatri/13384_tumor_doubletfinder.tsv",
            quote = F, row.names = F, sep = "\t")

#cellinfo <- read.table("/scratch/hnatri/13384_tumor_doubletfinder.tsv")

# Adding DoubletFinder info

# The final objects
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered.rds")

doublets <- cellinfo[which(cellinfo$doublet_finder=="Doublet"),]$cellname
length(intersect(doublets, colnames(tumors)))
length(intersect(doublets, colnames(immune_fibro)))

tumors$DoubletFinder <- ifelse(colnames(tumors) %in% doublets, "Doublet", "Singlet")
immune_fibro$DoubletFinder <- ifelse(colnames(immune_fibro) %in% doublets, "Doublet", "Singlet")

table(tumors$DoubletFinder)
table(tumors$celltype, tumors$DoubletFinder)

table(immune_fibro$DoubletFinder)
table(immune_fibro$celltype, immune_fibro$DoubletFinder)

tumors <- subset(tumors, subset = DoubletFinder == "Singlet")
immune_fibro <- subset(immune_fibro, subset = DoubletFinder == "Singlet")

saveRDS(tumors, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
saveRDS(immune_fibro, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

