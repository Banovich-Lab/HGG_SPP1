#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/30/2023
# Description: Generating plots for the 13384 GBM tumor manuscript
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

integrated_seurat <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/Immune_Combined_Human_Mouse_RPCA_SCT.rds")
Labels_JAK1_KO <- as.data.frame(Labels_JAK1_KO)
Labels_JAK1_KO$Labels_JAK1_KO <- as.character(Labels_JAK1_KO$Labels_JAK1_KO)

immune.combined.sct$Mouse_Label <- plyr::mapvalues(x = colnames(immune.combined.sct),
                                                   from = rownames(Labels_JAK1_KO),
                                                   to = Labels_JAK1_KO$Labels_JAK1_KO)
immune.combined.sct$Mouse_Label[-which(colnames(immune.combined.sct) %in% rownames(Labels_JAK1_KO))] <- NA
DefaultAssay(immune.combined.sct) <- "RNA"
gbm_jak1ko <- immune.combined.sct
rm(immune.combined.sct)

#==============================================================================#
# Figure 3: GBM and JAK1KO mouse tumors
#==============================================================================#

#========================
# UMAPs
#========================

DimPlot(gbm_jak1ko,
        group.by = "seurat_clusters",
        label = T,
        ncol = 1,
        cols = integrated_tumor_clusters_col) &
  coord_fixed(ratio=1) &
  theme_bw() &
  NoLegend()

DimPlot(gbm_jak1ko,
        group.by = "Species",
        label = F,
        ncol = 1,
        cols = c("azure3", "deeppink"),
        order = T) &
  coord_fixed(ratio=1) &
  theme_bw()

integrated_clusters <- DimPlot(gbm_jak1ko,
                               group.by = "seurat_clusters",
                               cols = integrated_tumor_clusters_col,
                               reduction = "umap",
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

integrated_clusters

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UMAP_cluster_GBM_JAK1KO.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

integrated_clusters

dev.off()

# Cell type annotations
gbm_jak1ko$celltype <- plyr::mapvalues(x = gbm_jak1ko$seurat_clusters,
                                       from = celltype_annot_integrated$cluster,
                                       to = celltype_annot_integrated$annotation)

unique(gbm_jak1ko$celltype)

integrated_celltypes <- DimPlot(gbm_jak1ko,
                                group.by = "celltype",
                                cols = celltype_annot_integrated_col,
                                reduction = "umap",
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

integrated_celltypes

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UMAP_celltype_GBM_JAK1KO.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

integrated_celltypes

dev.off()

#========================
# Cluster similarity and proportions
#========================

unique(gbm_jak1ko$Species)
p3 <- create_clusterpropplot(seurat_object = gbm_jak1ko,
                             group_var = "Species",
                             group1 = "Mm",
                             group2 = "Hs",
                             plot_var = "celltype",
                             plot_colors = celltype_annot_integrated_col,
                             var_names = c("JAK1KO", "GBM"),
                             legend_title = "") + manuscript_theme

p3

p3 + scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0))

filename <- "/home/hnatri/CART/13384_Tumors/Plots/celltype_prop_scatter_GBM_JAK1KO.pdf"
pdf(file = filename,
    width = 2.75,
    height = 2.75)

p3 + scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0))

dev.off()


#========================
# Feature expression
#========================

plot_features <- c("CD3D", "CD4", "CD68",
                   "C1QA", "C1QB", "C1QC",
                   "SPP1", "NLRP3", "APOE",
                   "ICAM1", "S100A8", "S100A9",
                   "S100A10", "COL1A1", "COL11A1",
                   "ACTA2", "PDGFRB", "PDPN")

# Dotplot
gbm_jak1ko@assays$integrated
DotPlot(gbm_jak1ko,
        features = plot_features,
        group.by = "seurat_clusters",
        assay = "integrated") +
  #cols = c("white", "tomato3")) +
  coord_flip()

# Cluster markers
# Removing ribosomal protein

genes <- rownames(LayerData(gbm_jak1ko, assay = "RNA", layer = "counts"))
genes <- genes[-grep("RPS", genes)]
genes <- genes[-grep("RPL", genes)]
genes <- genes[-grep("MT-", genes)]
gbm_jak1ko_woutrp <- subset(gbm_jak1ko, features = genes)
markers <- presto::wilcoxauc(gbm_jak1ko_woutrp,
                             group_by = "seurat_clusters",
                             assay = "data",
                             seurat_assay = "RNA")

top_markers <- markers %>%  group_by(group) %>% slice_max(order_by = auc, n = 50)
#write.table(top_markers, "/home/hnatri/CART/gbm_jak1ko_top50_markers.tsv", quote = F, row.names = F, sep = "\t")

# T cell markers
gs4_deauth()
canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/186PIcU3Jb0Xm3IuZgENEChluykvwsK4JnHF8FbeJASA/edit?usp=sharing")
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "T cells, gene sets")

DotPlot(immune_fibro,
        features = unique(canonical_markers$RNA),
        group.by = "celltype",
        cols = c("azure", "tomato3")) +
  coord_flip()

