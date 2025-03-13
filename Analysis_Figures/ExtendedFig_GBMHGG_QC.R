#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/24/2024
# Description: SFig of scRNA-seq quality metrics
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

#==============================================================================#
# Import data
#==============================================================================#

# The final object
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#==============================================================================#
# Plotting
#==============================================================================#\

phase_p1 <- DimPlot(tumors,
                    group.by = "Phase",
                    #cols = immune_fibro_celltype_col,
                    reduction = "integrated_sct_umap",
                    #label = T,
                    #label.box = T,
                    #label.size = 3,
                    #repel = T,
                    raster = T,
                    raster.dpi = c(1024, 1024),
                    pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  manuscript_theme + 
  #NoLegend() +
  NoAxes() +
  coord_fixed(1)

phase_p1

phase_p2 <- DimPlot(immune_fibro,
                    group.by = "Phase",
                    #cols = immune_fibro_celltype_col,
                    reduction = "umap",
                    #label = T,
                    #label.box = T,
                    #label.size = 3,
                    #repel = T,
                    raster = T,
                    raster.dpi = c(1024, 1024),
                    pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  manuscript_theme + 
  #NoLegend() +
  NoAxes() +
  coord_fixed(1)

phase_p2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/13384_tumor_phase.pdf"
pdf(file = filename,
    width = 7,
    height = 3)

phase_p + phase_p2

dev.off()

f1 <- FeaturePlot(tumors,
                  features = c("S.Score", "G2M.Score",
                               "nCount_SoupX_RNA",
                               "nFeature_SoupX_RNA",
                               "percent.mt_SoupX_RNA",
                               "percent.ribo_SoupX_RNA"),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  keep.scale = "all",
                  ncol = 3,
                  raster = T,
                  raster.dpi = c(2048, 1024),
                  pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  #manuscript_theme &
  NoAxes()
  #NoLegend()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

f1

f2 <- FeaturePlot(immune_fibro,
                  features = c("S.Score", "G2M.Score",
                               "nCount_SoupX_RNA",
                               "nFeature_SoupX_RNA",
                               "percent.mt_SoupX_RNA",
                               "percent.ribo_SoupX_RNA"),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  keep.scale = "all",
                  ncol = 3,
                  raster = T,
                  raster.dpi = c(2048, 1024),
                  pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  #manuscript_theme &
  NoAxes()
  #NoLegend()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

f1 / f2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/13384_tumor_QC.pdf"
pdf(file = filename,
    width = 10,
    height = 9)

f1 / f2

dev.off()

# Violinplots
tumors$celltype <- factor(tumors$celltype, levels = sort(as.character(unique(tumors$celltype))))
v1 <- VlnPlot(tumors,
              features = c("nCount_SoupX_RNA",
                           "nFeature_SoupX_RNA"),
              group.by = "celltype",
              pt.size = 0,
              cols = tumor_celltype_col,
              ncol = 2) &
  theme_classic() &
  manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v1

v2 <- VlnPlot(immune_fibro,
              features = c("nCount_SoupX_RNA",
                           "nFeature_SoupX_RNA"),
              group.by = "celltype",
              pt.size = 0,
              cols = immune_fibro_celltype_col,
              ncol = 2) &
  theme_classic() &
  manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/13384_tumor_QC_vln.pdf"
pdf(file = filename,
    width = 6,
    height = 4)

v1 / v2

dev.off()
