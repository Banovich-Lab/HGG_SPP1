#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 8/23/2023
# Description: Only keeping singlets after running DoubletFinder
#==============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(googlesheets4)
library(dplyr)
library(DoubletFinder)

#==============================================================================
# Import data and filter
#==============================================================================

cellinfo <- read.table("/scratch/hnatri/13384_tumor_doubletfinder.tsv", header = T)

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
