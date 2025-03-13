#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 12/16/2024
# Description: Plotting JAK1KO/WT CAR product scRNA-seq
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)
library(presto)
library(ggrepel)
library(tidyverse)

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
reduction <- "umap"

#==============================================================================#
# Import data
#==============================================================================#

load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Mouse_CART.RData")

merged_data <- merge(CART_KO.Tcells, CART_WT.Tcells, merge.dr = T)

# CAR positive cells
counts <- LayerData(merged_data, layer = "counts",
                    features = c("CD19T"))

hist(counts, breaks = 25, main = "CD19T")

carpos <- names(counts[which(counts>0)])
merged_data$CARpos <- ifelse(colnames(merged_data) %in% carpos, "CARpos", "CARneg")
merged_data$CD19T <- counts

# Difference in CAR+ proportion between JAK1/KO and WT?
cont_table <- table(merged_data$CARpos, merged_data$JAK1_Group)

chi2 <- chisq.test(cont_table)

# All cells
combined <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/scRNA_seq_Pool.combined.rds")
combined <- scRNA_seq_Pool.combined
preCAR <- subset(combined, subset = CART_Group == "Control")
#preCAR <- subset(preCAR, subset = seurat_clusters %in% c("T1", "T2", "T3", "T4", "Neu3"))

DefaultAssay(preCAR) <- "RNA"
counts <- LayerData(preCAR, layer = "counts",
                    features = c("CD19T"))

#==============================================================================#
# Plotting
#==============================================================================#

barplot <- merged_data@meta.data %>% group_by(JAK1_Group) %>%
  dplyr::summarise(n = n(),
                   n_pos = length(which(CD19T>0)),
                   pct_pos = n_pos/n*100) %>%
  ggplot(aes(x = JAK1_Group, y = pct_pos)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("") +
  ylab("% CAR+")

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CARpos_pct.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

barplot

dev.off()

# Significance
merged_data@meta.data %>% group_by(JAK1_Group) %>%
  dplyr::summarise(n = n(),
                   n_pos = length(which(CD19T>0)),
                   pct_pos = n_pos/n*100)

data <- matrix(c(2864, 951, 96, 18), nrow = 2)
result <- chisq.test(data)
print(result)

umap_celltype <- DimPlot(merged_data,
                         group.by = "seurat_clusters",
                         #cols = jak1_celltype_col,
                         reduction = "umap",
                         label = T,
                         label.box = T,
                         label.size = 3,
                         repel = T,
                         raster = T,
                         raster.dpi = c(1024, 1024),
                         pt.size = 3) +
  ggtitle("KO+WT") +
  theme_classic() +
#  manuscript_theme + 
  NoLegend() +
  #NoAxes() +
  coord_fixed(1)

umap_celltype

DimPlot(merged_data,
        group.by = "seurat_clusters",
        #cols = jak1_celltype_col,
        reduction = "umap",
        split.by = "JAK1_Group",
        #label = T,
        #label.box = T,
        #label.size = 3,
        #repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
#  manuscript_theme + 
  NoLegend() +
  #NoAxes() +
  coord_fixed(1)

f1 <- FeaturePlot(merged_data,
                  features = c("Cd3e",
                               "CD19T",
                               "Cd44",
                               "Spp1"),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  #keep.scale = "all",
                  ncol = 2,
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

v1 <- VlnPlot(merged_data,
              features = c("Cd3e",
                           "CD19T",
                           "Cd44",
                           "Spp1"),
              group.by = "seurat_clusters",
              pt.size = 0,
              #cols = jak1_celltype_col,
              ncol = 1) &
  theme_classic() &
  #manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v1

v2 <- VlnPlot(merged_data,
              features = c("Cd3e",
                           "CD19T",
                           "Cd44",
                           "Spp1"),
              split.by = "JAK1_Group",
              group.by = "seurat_clusters",
              pt.size = 0,
              layer = "data",
              #cols = jak1_celltype_col,
              ncol = 1) &
  theme_classic() &
  #manuscript_theme &
  #NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v2

v3 <- VlnPlot(merged_data,
              features = c("Cd3e",
                           "CD19T",
                           "Cd44",
                           "Spp1"),
              split.by = "JAK1_Group",
              group.by = "seurat_clusters",
              split.plot = T,
              pt.size = 0,
              layer = "counts",
              #cols = jak1_celltype_col,
              ncol = 1) &
  theme_classic() &
  #manuscript_theme &
  #NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v3

v4 <- VlnPlot(merged_data,
              features = c("Cd3e",
                           "CD19T",
                           "Cd44",
                           "Spp1"),
              split.by = "JAK1_Group",
              group.by = "CARpos",
              split.plot = T,
              pt.size = 0,
              layer = "counts",
              #cols = jak1_celltype_col,
              ncol = 1) &
  theme_classic() &
  #manuscript_theme &
  #NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v4


#filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CAR_Spp1_Cd44_Cd19T.pdf"
#pdf(file = filename,
#    width = 4,
#    height = 8)
#
#f1
#
#dev.off()

# Exhaustion markers

vp <- VlnPlot(merged_data,
        features = c("Cd3e",
                     "CD19T",
                     "Cd44",
                     "Pdcd1"),
        split.by = "JAK1_Group",
        group.by = "CARpos",
        split.plot = T,
        pt.size = 0,
        layer = "counts",
        #cols = jak1_celltype_col,
        ncol = 2,
        cols = c("deeppink3", "aquamarine3")) &
  theme_classic() &
  #manuscript_theme &
  #NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

merged_data <- NormalizeData(merged_data)
Idents(merged_data) <- merged_data$CARpos
markers_ko <- FindMarkers(subset(merged_data, subset = JAK1_Group == "Jak1_KO"),
                          ident.1 = "CARpos",
                          ident.2 = "CARneg",
                          assay = "RNA",
                          verbose = F)
markers_wt <- FindMarkers(subset(merged_data, subset = JAK1_Group == "WildType"),
                          ident.1 = "CARpos",
                          ident.2 = "CARneg",
                          assay = "RNA",
                          verbose = F)

markers_ko %>% filter(rownames(markers_ko) %in% c("Cd44", "Pdcd1"))
markers_wt %>% filter(rownames(markers_wt) %in% c("Cd44", "Pdcd1"))

VlnPlot(merged_data,
        features = c("Cd3e",
                     "CD19T",
                     "Cd44",
                     "Pdcd1",
                     "Tigit",
                     "Havcr2"),
        split.by = "CARpos",
        group.by = "JAK1_Group",
        split.plot = T,
        pt.size = 0,
        layer = "counts",
        #cols = jak1_celltype_col,
        ncol = 2) &
  theme_classic() &
  #manuscript_theme &
  #NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Idents(merged_data) <- merged_data$JAK1_Group
markers_cp <- FindMarkers(subset(merged_data, subset = CARpos == "CARpos"),
                          ident.1 = "Jak1_KO",
                          ident.2 = "WildType",
                          assay = "RNA",
                          verbose = F)
markers_cn <- FindMarkers(subset(merged_data, subset = CARpos == "CARneg"),
                          ident.1 = "Jak1_KO",
                          ident.2 = "WildType",
                          assay = "RNA",
                          verbose = F)

markers_cp %>% filter(rownames(markers_cp) %in% c("Havcr2", "Tigit", "Pdcd1"))
markers_cn %>% filter(rownames(markers_cn) %in% c("Havcr2", "Tigit", "Pdcd1"))

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CAR_Spp1_Cd44_Cd19T_PD1_vlnplot.pdf"
pdf(file = filename,
    width = 6,
    height = 4)

vp

dev.off()
