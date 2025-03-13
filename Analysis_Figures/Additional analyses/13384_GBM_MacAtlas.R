#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 07/14/2024
# Description: Lifting annotations from the macrophage atlas
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
library(tidyverse)
library(googlesheets4)
library(survminer) # for theme_classic2()

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import data
#==============================================================================#

# The final object
#tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

# Macrophage annotations from the Seurat object
macs <- readRDS("/scratch/hnatri/CART/mac.atlas.zenodo.200524.rds")

#==============================================================================#
# Label transfer
#==============================================================================#

DefaultAssay(macs) <- "RNA"

#macs[["integrated"]] <- NULL

#macs@assays$RNA@layers$data <- NULL
#macs@assays$RNA@layers$scale.data <- NULL
counts <- GetAssayData(macs, assay = "RNA", layer = "counts")

rownames(counts)
colnames(counts)
dim(counts)

# Zero counts?
colSums(counts)

rm_cells <- ""

rownames(macs)
colnames(macs)

head(macs@meta.data)
# Replacing cell names
rownames(macs@meta.data) <- macs$cellid

identical(rownames(macs@meta.data), colnames(counts))
setdiff(rownames(macs@meta.data), colnames(counts))

#macs[["RNA"]] <- CreateAssay5Object(counts = counts)
#macs[["RNA5"]] <- as(object = macs[["RNA"]], Class = "Assay5")

#DefaultAssay(macs) <- "RNA"
#macs <- NormalizeData(macs)
#macs <- ScaleData(macs)
#macs <- SCTransform(macs,
#                    vst.flavor = "v2",
#                    verbose = F)
macs <- RunPCA(macs)

DefaultAssay(immune_fibro) <- "RNA"
immune_fibro <- NormalizeData(immune_fibro)
anchors <- FindTransferAnchors(reference = macs,
                               query = immune_fibro,
                               dims = 1:20,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors,
                            refdata = macs$short.label,
                            dims = 1:20)
immune_fibro <- AddMetaData(immune_fibro, metadata = predictions)

saveRDS(immune_fibro, "/scratch/hnatri/CART/tumor_immune_fibroblast_MacAtlas.rds")
q(save = "no")

immune_fibro <- readRDS("/scratch/hnatri/CART/tumor_immune_fibroblast_MacAtlas.rds")

#==============================================================================#
# Plotting
#==============================================================================#

# UMAPs
dimplot_celltype <- DimPlot(immune_fibro,
                            group.by = "celltype",
                            cols = immune_fibro_celltype_col,
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

dimplot_macatlas <- DimPlot(immune_fibro,
                            group.by = "short.label",
                            #cols = tumor_celltype_col,
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

dimplot_celltype