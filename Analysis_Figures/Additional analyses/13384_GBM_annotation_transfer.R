#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 07/22/2024
# Description: Transfering immune+fibroblast annotations to the full object
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(googlesheets4)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
reduction <- "integrated_sct_umap"

#==============================================================================#
# Import data
#==============================================================================#

# The final object
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#==============================================================================#
# Transfer annotations and plot
#==============================================================================#

tumors$celltype2 <- mapvalues(x = colnames(tumors),
                              from = colnames(immune_fibro),
                              to = immune_fibro$celltype)
tumors@meta.data[!(tumors$celltype2 %in% immune_fibro$celltype),]$celltype2 <- "NA"
tumors$celltype2 <- factor(tumors$celltype2,
                           levels = c(paste0("M", seq(1, 9, by = 1)),
                                      "B1", "N1",
                                      paste0("L", seq(1, 10, by = 1)),
                                      paste0("F", seq(1, 3, by = 1)),
                                      "NA"))

tumors$celltype <- factor(tumors$celltype, levels = c("Olig1", "Olig2",
                                                      "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6",
                                                      "Fibroblast", "Lymph1", "Lymph2",
                                                      "Myel1", "Myel2", "Myel3", "Myel4", "Myel5",  "Myel6"))

plot_col <- c(immune_fibro_celltype_col, c("NA" = "gray89"))

tumor_celltype <- DimPlot(tumors,
                          group.by = "celltype2",
                          cols = plot_col,
                          reduction = reduction,
                          label = T,
                          label.box = T,
                          label.size = 3,
                          repel = T,
                          raster = T,
                          raster.dpi = c(1024, 1024),
                          pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  manuscript_theme + 
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

tumor_celltype

barplot1 <- create_barplot(tumors,
                           group_var = "celltype",
                           plot_var = "celltype2",
                           plot_levels = levels(tumors$celltype2),
                           group_levels = levels(tumors$celltype),
                           plot_colors = plot_col,
                           var_names =  c("Frequency", "Full object celltype"),
                           legend_title = "Immune+fibroblast celltype")# + manuscript_theme


barplot1

barplot2 <- create_barplot(tumors,
                           group_var = "celltype2",
                           plot_var = "celltype",
                           plot_levels = levels(tumors$celltype),
                           group_levels = levels(tumors$celltype2),
                           plot_colors = tumor_celltype_col,
                           var_names =  c("Frequency", "Immune+fibroblast celltype"),
                           legend_title = "Full object celltype")# + manuscript_theme


barplot2

barplot1 / barplot2

