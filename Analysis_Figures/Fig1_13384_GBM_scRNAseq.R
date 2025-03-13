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
library(scProportionTest)

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

# Updating cell type names
tumors$celltype <- plyr::mapvalues(x = tumors$cluster,
                                   from = celltype_annot$cluster,
                                   to = celltype_annot$annotation)

immune_fibro$celltype <- plyr::mapvalues(x = immune_fibro$sub.cluster3,
                                         from = celltype_annot_immune_fibro$cluster,
                                         to = celltype_annot_immune_fibro$annotation)

tumors$celltype <- factor(tumors$celltype, levels = c("Olig1", "Olig2",
                                                      "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6",
                                                      "Fibroblast", "Lymph1", "Lymph2",
                                                      "Myel1", "Myel2", "Myel3", "Myel4", "Myel5",  "Myel6"))

immune_fibro$celltype <- factor(immune_fibro$celltype,
                                levels = c(paste0("M", seq(1, 9, by = 1)),
                                           "B1", "N1",
                                           paste0("L", seq(1, 10, by = 1)),
                                           paste0("F", seq(1, 3, by = 1))))

top_surv <- c(109, 141, 181, 223, 265, 301)
bottom_surv <- c(129, 146, 185, 224, 228, 234, 237, 248)
immune_fibro$SPP1_high_low <- ifelse(immune_fibro$UPN %in% top_surv, "SPP1low",
                                       ifelse(immune_fibro$UPN %in% bottom_surv, "SPP1high", "SPP1middle"))

tumors@meta.data %>%
  as.data.frame() %>%
  dplyr::select("UPN", "Diagnosis.Histology") %>%
  distinct() %>%
  dplyr::select("Diagnosis.Histology") %>%
  table()

# Top markers for each celltype
markers <- presto::wilcoxauc(tumors,
                             group_by = "celltype",
                             assay = "data",
                             seurat_assay = "RNA")
top_markers <- markers %>%  group_by(group) %>% slice_max(order_by = auc, n = 50)

write.table(top_markers, "/home/hnatri/13384_CART/13384_Tumors/tumor_celltype_markers_top50.tsv",
            quote = F, sep = "\t", row.names = F)

markers_tme <- presto::wilcoxauc(immune_fibro,
                                 group_by = "celltype",
                                 assay = "data",
                                 seurat_assay = "RNA")
top_markers_tme <- markers_tme %>%  group_by(group) %>% slice_max(order_by = auc, n = 50)

write.table(top_markers_tme, "/home/hnatri/13384_CART/13384_Tumors/immune_fibro_celltype_markers_top50.tsv",
            quote = F, sep = "\t", row.names = F)

#saveRDS(tumors, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
#saveRDS(immune_fibro, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#========================
# Cell cycle
#========================

VlnPlot(subset(immune_fibro, subset = SPP1_high_low %in% c("SPP1low", "SPP1high")),
        features = c("S.Score", "G2M.Score"),
        group.by = "celltype",
        split.by = "SPP1_high_low",
        split.plot = T,
        pt.size = 0,
        ncol = 1) &
  theme_minimal()

FeaturePlot(subset(immune_fibro, subset = SPP1_high_low %in% c("SPP1low", "SPP1high")),
            features = c("S.Score", "G2M.Score"),
            split.by = "SPP1_high_low",
            ncol = 1) &
  theme_minimal()

VlnPlot(subset(immune_fibro, subset = SPP1_high_low %in% c("SPP1low", "SPP1high")),
        group.by = "celltype",
        features = c("GPX4", "SPP1"),
        split.by = "SPP1_high_low",
        split.plot = T,
        pt.size = 0,
        ncol = 1) &
  theme_minimal()

# Checking for significance
DEG_SPP1high <- lapply(unique(immune_fibro$celltype), function(xx){
  message(xx)
  data_subset <- subset(immune_fibro, subset = celltype == xx)
  Idents(data_subset) <- data_subset$SPP1_high_low
  
  table(data_subset$SPP1_high_low)
  if (min(table(data_subset$SPP1_high_low))<2){
    return(NULL)
  }
  
  if (all((c("SPP1low", "SPP1high") %in% data_subset$SPP1_high_low) == c(T, T))){
    markers <- FindMarkers(data_subset,
                           ident.1 = "SPP1low",
                           ident.2 = "SPP1high",
                           assay = "RNA",
                           verbose = F)
    markers$feature <- rownames(markers)
    markers$celltype <- xx
    
    return(markers)
  } else {
    return(NULL)
  }
})
names(DEG_SPP1high) <- unique(immune_fibro$celltype)
DEG_SPP1high[sapply(DEG_SPP1high, is.null)] <- NULL

DEG_SPP1high_df <- as.data.frame(do.call(rbind, DEG_SPP1high))

# Distribution of log2FC
hist(DEG_SPP1high_df$avg_log2FC)

DEG_SPP1high_df %>% filter(feature == "GPX4") %>%
  arrange(celltype)

group_col <- c("SPP1low" = "aquamarine3",
               "SPP1high" = "deeppink3",
               "SPP1middle" = "gray80")
phase_col <- c("G1" = "lightblue2",
               "G2M" = "azure4",
               "S" = "gray85")

create_barplot(immune_fibro,
               group_var = "SPP1_high_low",
               plot_var = "Phase",
               plot_levels = sort(unique(immune_fibro$Phase)),
               group_levels = sort(unique(immune_fibro$SPP1_high_low)),
               plot_colors = phase_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Phase")

create_barplot(immune_fibro,
               group_var = "celltype",
               plot_var = "Phase",
               plot_levels = sort(unique(immune_fibro$Phase)),
               group_levels = sort(unique(immune_fibro$celltype)),
               plot_colors = phase_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Phase")

#========================
# Cell counts
#========================

celln <- table(tumors$UPN, tumors$celltype)
rowsums <- rowSums(celln)
colsums <- colSums(celln)
celln <- rbind(celln, colsums)
celln <- cbind(celln, c(rowsums, NA))

write.table(celln, "/home/hnatri/CART/13384_Tumors/celln_UPN_celltype.tsv",
            quote = F, sep = "\t", row.names = T)

celln <- table(immune_fibro$UPN, immune_fibro$celltype)
rowsums <- rowSums(celln)
colsums <- colSums(celln)
celln <- rbind(celln, colsums)
celln <- cbind(celln, c(rowsums, NA))

write.table(celln, "/home/hnatri/CART/13384_Tumors/celln_immune_fibro_UPN_celltype.tsv",
            quote = F, sep = "\t", row.names = T)

table(tumors$binary_response, tumors$celltype) %>%
  as.data.frame() %>%
  filter(Var2 %in% paste0("Myel", seq(1, 6))) %>%
  group_by(Var1) %>%
  dplyr::summarise(sum = sum(Freq))

# By batch and UPN
batch_upn_celln <- table(tumors$Batch, tumors$UPN) %>%
  as.data.frame() %>% 
  filter(Freq>0)

write.table(batch_upn_celln, "/home/hnatri/CART/13384_Tumors/celln_tumor_batch_UPN.tsv",
            quote = F, sep = "\t", row.names = F)

#==============================================================================#
# Figure 1: GBM tumors
#==============================================================================#

#========================
# Sample numbers
#========================

# Sample numbers by CD3 score and response
metadata <- tumors@meta.data %>%
  dplyr::select("UPN", "CD3_high_low", "binary_response") %>%
  distinct()

plot_data <- table(metadata$CD3_high_low,
                   metadata$binary_response)

small_heatmap <- ComplexHeatmap::pheatmap(plot_data,
                                          legend = TRUE,
                                          cluster_rows = FALSE,
                                          cluster_cols = FALSE,
                                          display_numbers = plot_data,
                                          color = c("white", "tomato2"),
                                          heatmap_legend_param = list(title = "# Donors",
                                                                      at = c(0, 10, 20)))

filename <- "/home/hnatri/CART/13384_Tumors/Plots/samplen_small_heatmap.pdf"
pdf(file = filename,
    width = 2.75,
    height = 2.75)

small_heatmap

dev.off()

#========================
# UMAPs
#========================

# UMAPs
tumor_celltype <- DimPlot(tumors,
                          group.by = "celltype",
                          cols = tumor_celltype_col,
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

# Plot for presentation, highlighting immune+fibroblast clusters
tumor_celltype_col_cp <- tumor_celltype_col
tumor_celltype_col_cp <- ifelse(names(tumor_celltype_col_cp) %in% c("Olig1", "Olig2",
                                                                    "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5", "Tumor6"),
                                "whitesmoke", tumor_celltype_col_cp)
names(tumor_celltype_col_cp) <- names(tumor_celltype_col)

DimPlot(tumors,
        group.by = "celltype",
        cols = tumor_celltype_col_cp,
        reduction = reduction,
        label = T,
        label.box = T,
        repel = T) +
  ggtitle("") +
  theme_classic() +
  manuscript_theme + 
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

DimPlot(immune_fibro,
        group.by = "sub.cluster3",
        #cols = immune_fibro_celltype_col,
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

immune_fibro_celltype <- DimPlot(immune_fibro,
                                 group.by = "celltype",
                                 cols = immune_fibro_celltype_col,
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

immune_fibro_celltype

tumor_celltype / immune_fibro_celltype

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UMAP_celltype_tumor_immunefib.pdf"
pdf(file = filename,
    width = 3,
    height = 6)

tumor_celltype / immune_fibro_celltype

dev.off()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UMAP_celltype_tumor.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

tumor_celltype

dev.off()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UMAP_celltype_immunefib.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

immune_fibro_celltype

dev.off()

DimPlot(immune_fibro,
        group.by = "celltype",
        split.by = "celltype",
        cols = immune_fibro_celltype_col,
        reduction = "umap",
        #label = T,
        #label.box = T,
        #label.size = 3,
        #repel = T,
        raster = T,
        pt.size = 3,
        ncol = 5) +
  #raster.dpi = c(1024, 1024),
  #pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  #manuscript_theme + 
  #NoLegend() +
  #NoAxes() +
  coord_fixed(1)

#========================
# Cell type proportions by UPN
#========================

barplot1 <- create_barplot(tumors,
                          group_var = "UPN",
                          plot_var = "celltype",
                          plot_levels = levels(tumors$celltype),
                          group_levels = sort(unique(tumors$UPN)),
                          plot_colors = tumor_celltype_col,
                          var_names =  c("Frequency", "UPN"),
                          legend_title = "UPN") + manuscript_theme

filename <- "/home/hnatri/CART/13384_Tumors/Plots/UPN_celltype_barplot.pdf"
pdf(file = filename,
    width = 6,
    height = 2.5)

barplot1

dev.off()

barplot2 <- create_barplot(immune_fibro,
                           group_var = "UPN",
                           plot_var = "celltype",
                           plot_levels = c(paste0("M", seq(1, 9, by = 1)),
                                           "B1", "N1",
                                           paste0("L", seq(1, 10, by = 1)),
                                           paste0("F", seq(1, 3, by = 1))),
                           group_levels = sort(unique(tumors$UPN)),
                           plot_colors = immune_fibro_celltype_col,
                           var_names =  c("Frequency", "UPN"),
                           legend_title = "UPN") + manuscript_theme

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumors_immune_fibro_UPN_celltype_barplot.pdf"
pdf(file = filename,
    width = 8,
    height = 6)

barplot1 / barplot2

dev.off()

#========================
# Feature expression
#========================

# CNVs
cnvs <- FeaturePlot(infercnv,
                    features = "average_prop_cnv",
                    cols = c("gray89", "tomato3"),
                    #order = T,
                    #keep.scale = "all",
                    ncol = 1,
                    raster = T,
                    raster.dpi = c(2048, 1024),
                    pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes()

cnvs

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CNV_FeaturePlot.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

cnvs

dev.off()

# Suppl. Fig. features for the whole object
suppl_features <- c("PTPRC", "CD14", "CD3E", "OLIG1", "NLGN1", "COL1A1")

tumor_sf <- FeaturePlot(tumors,
                        features = suppl_features,
                        #group.by = "celltype",
                        cols = c("gray89", "tomato3"),
                        order = T,
                        keep.scale = "feature",
                        ncol = 3,
                        raster = T) &
  #raster.dpi = c(1024, 1024),
  #pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

tumor_sf

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_global_FeaturePlots.pdf"
pdf(file = filename,
    width = 6,
    height = 3)

tumor_sf

dev.off()

v1 <- VlnPlot(tumors,
              features = suppl_features,
              group.by = "celltype",
              pt.size = 0,
              cols = tumor_celltype_col,
              ncol = 3) &
  theme_classic() &
  manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v1

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_global_VlnPlots.pdf"
pdf(file = filename,
    width = 8,
    height = 3)

v1

dev.off()

plot_features <- c("PTMA", "PFN1", "CFL1", "TMSB4X", "TPT1", "TMSB10", "MIF",
                   "PDPN", "NLRP3", "IL1B", "CCL4", "S100A8", "S100A9",
                   "S100A10", "TYROBP", "CD68", "ICAM1", "C1QA", "C1QB", "C1QC",
                   "CD74", "AREG", "CD4", "APOE", "FABP5", "SPP1", "CD274",
                   "CD96", "PTPRC", "CEMIP2", "KLRD1", "CD8A", "NKG7", "IL32",
                   "CD3D", "BTG1", "IFITM2", "ITM2A", "SELL", "GZMB", "CD79A",
                   "ACTA2", "PDGFRB", "COL1A1", "CD163", "MRC1", "ITGAM", "CD14",
                   "CD279", "PDCD1", "TREM2", "TMEM119", "P2RY12", "CX3CR1")

# CPVL, FLT3, CLEC9A, C1orf54, CSF1R, CD1C, FCER1A, LTB, LAG3, BTLA, TIGIT,
# HAVCR2, ICOS, IL10RA, CD80, CD86, FOXP3, TMSB4X, SELL, LY6C2 (?), CD160, GZMB,
# ROCR (?), AREG, KLRG1, CDH1 (E-Cadherin), TACR1 (NK1.1, NK1R1), NKX1-2 (NK1.2),
# NCAM (CD56), KLRB1 (CD161), CD27 (TNFRSF7), TNFSF7 (?), CD74, MIF, IL1-B
additional_features <- c("CPVL", "FLT3", "CLEC9A", "C1orf54", "CSF1R", "CD1C",
                         "FCER1A", "LTB", "LAG3", "BTLA", "TIGIT", "HAVCR2",
                         "ICOS", "IL10RA", "CD80", "CD86", "FOXP3", "TMSB4X",
                         "SELL", "CD160", "GZMB", "AREG", "KLRG1", "CDH1",
                         "TACR1", "KLRB1", "CD27", "CD74", "MIF", "IL1B")

# CD15 (FUT4),CD45 (PTPRC),CD66B (not in our data), CD62L (SELL), Arg1 and PDL-1 (CD274)
# Dotplot
immune_fibro_dotplot <- DotPlot(immune_fibro,
                                features = plot_features,
                                group.by = "celltype",
                                cols = c("azure", "tomato3")) +
  #coord_flip() +
  manuscript_theme
  #RotatedAxis()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot.pdf"
pdf(file = filename,
    width = 4,
    height = 7)

immune_fibro_dotplot

dev.off()

DotPlot(immune_fibro,
  features = additional_features,
  group.by = "celltype",
  cols = c("azure", "tomato3")) +
  coord_flip() +
  manuscript_theme +
  RotatedAxis()

# seurat_object, plot_features, group_var, group_colors, column_title, km=5, row.order = NULL
immune_fibro_dotplot <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro,
                                               plot_features = sort(plot_features),
                                               group_var = "celltype",
                                               group_colors = immune_fibro_celltype_col,
                                               column_title = "",
                                               row_km = 5,
                                               col.order = plot_features)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot_heatmap_clustercelltypes.pdf"
pdf(file = filename,
    width = 10,
    height = 4)

immune_fibro_dotplot

dev.off()

celltype_order <- c(paste0("M", seq(1, 9, by = 1)),
                    "N1", "B1",
                    paste0("L", seq(1, 10, by = 1)),
                    paste0("F", seq(1, 3, by = 1)))

immune_fibro_dotplot_unclustered <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro,
                                                          plot_features = plot_features,
                                                          group_var = "celltype",
                                                          group_colors = immune_fibro_celltype_col,
                                                          column_title = "",
                                                          col.order = plot_features,
                                                          row.order = celltype_order)

immune_fibro_dotplot_clustered <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro,
                                                          plot_features = plot_features,
                                                          group_var = "celltype",
                                                          group_colors = immune_fibro_celltype_col,
                                                          column_title = "",
                                                          row.order = NULL,
                                                          col.order = NULL)

immune_fibro_dotplot_clustergenes <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro,
                                                                      plot_features = plot_features,
                                                                      group_var = "celltype",
                                                                      group_colors = immune_fibro_celltype_col,
                                                                      column_title = "",
                                                                      #row_km=5,
                                                                      #col.order = plot_features,
                                                                      row.order = celltype_order)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot_clusternone.pdf"
pdf(file = filename,
    width = 10,
    height = 4)

immune_fibro_dotplot_unclustered

dev.off()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot_clusterall.pdf"
pdf(file = filename,
    width = 10,
    height = 4)

immune_fibro_dotplot_clustered

dev.off()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot_clustegenes.pdf"
pdf(file = filename,
    width = 10,
    height = 4)

immune_fibro_dotplot_clustergenes

dev.off()

# Without T cells
immune_fibro_noT <- subset(immune_fibro, subset = celltype %in% paste0("L", seq(1, 10)), invert = T)

plot_features_noT <- c("PTMA", "PFN1", "CFL1", "TMSB4X", "TPT1", "TMSB10", "MIF",
                   "PDPN", "NLRP3", "IL1B", "CCL4", "S100A8", "S100A9",
                   "S100A10", "TYROBP", "CD68", "ICAM1", "C1QA", "C1QB", "C1QC",
                   "CD74", "AREG", "CD4", "APOE", "FABP5", "SPP1", "CD274",
                   "PTPRC", "BTG1", "IFITM2", "SELL", "GZMB", "CD79A",
                   "ACTA2", "PDGFRB", "COL1A1", "CD163", "MRC1", "ITGAM", "CD14",
                   "CD279", "TREM2", "TMEM119", "P2RY12", "CX3CR1")

immune_fibro_dotplot_noT <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro_noT,
                                                              plot_features = plot_features_noT,
                                                              group_var = "celltype",
                                                              group_colors = immune_fibro_celltype_col,
                                                              column_title = "",
                                                              row_km = 2,
                                                              col_km = 3,
                                                              row.order = NULL,
                                                              col.order = NULL)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_dotplot_noTcells.pdf"
pdf(file = filename,
    width = 8.5,
    height = 3.5)

immune_fibro_dotplot_noT

dev.off()

tumor_dotplot <- create_dotplot_heatmap_horizontal(seurat_object = tumors,
                                                              plot_features = plot_features,
                                                              group_var = "celltype",
                                                              group_colors = tumor_celltype_col,
                                                              column_title = "",
                                                              row_km = 2,
                                                              col_km = 3,
                                                              row.order = NULL,
                                                              col.order = NULL)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_mainmarkers_dotplot.pdf"
pdf(file = filename,
    width = 8.5,
    height = 3.5)

immune_fibro_dotplot_noT

dev.off()

# T cell markers
gs4_deauth()
canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/186PIcU3Jb0Xm3IuZgENEChluykvwsK4JnHF8FbeJASA/edit?usp=sharing")
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "T cells, gene sets")
canonical_markers <- canonical_markers %>% filter(!is.na(Bigger_gene_sets))



d1 <- DotPlot(subset(immune_fibro, subset = celltype %in% paste0("L", seq(1, 10))),
        features = unique(canonical_markers$RNA),
        group.by = "celltype",
        cols = c("azure", "tomato3")) +
  coord_flip() +
  manuscript_theme

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_immune_fibro_Tcellmarkers_dotplot.pdf"
pdf(file = filename,
    width = 4,
    height = 7)

d1

dev.off()

f1 <- FeaturePlot(immune_fibro,
                  features = sort(plot_features),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  keep.scale = "feature",
                  ncol = 8,
                  raster = T,
                  raster.dpi = c(1024, 1024),
                  pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

f1

filename <- "/scratch/hnatri/tumor_immune_fibro_featureplots.pdf"
pdf(file = filename,
    width = 12,
    height = 16)

f1

dev.off()

f2 <- FeaturePlot(immune_fibro,
                  features = sort(unique(canonical_markers$RNA)),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  keep.scale = "all",
                  ncol = 5,
                  raster = T) &
  #raster.dpi = c(1024, 1024),
  #pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

f2

filename <- "/scratch/hnatri/tumor_immune_fibro_Tcellmarkers_featureplots.pdf"
pdf(file = filename,
    width = 12,
    height = 16)

f2

dev.off()

FeaturePlot(immune_fibro,
            features = additional_features,
            #group.by = "celltype",
            cols = c("gray89", "tomato3"),
            order = T,
            keep.scale = "all",
            ncol = 6,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes() &
  NoLegend()


presentation_features <- c("MIF", "NLRP3", "IL1B", "CCL4",
                           "S100A8", "S100A9", "S100A10", "TYROBP",
                           "C1QA", "C1QB", "C1QC", "CD74", "AREG",
                           "APOE", "SPP1", "PTPRC", "CD8A", "NKG7",
                           "IL32", "CD3D", "GZMB", "CD79A", "ACTA2",
                           "COL1A1")

FeaturePlot(immune_fibro,
            features = presentation_features,
            #group.by = "celltype",
            cols = c("gray89", "tomato3"),
            order = T,
            keep.scale = "all",
            ncol = 6,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes() &
  NoLegend()

#========================
# Cell type proportions
#========================

# Barplots
tumor_response_bp <- create_barplot(seurat_object = tumors,
                                    group_var = "binary_response",
                                    plot_var = "celltype",
                                    plot_levels = levels(tumors$celltype),
                                    group_levels = c("CR_SD", "PD"),
                                    plot_colors = tumor_celltype_col,
                                    var_names =  c("Frequency (%)", ""),
                                    legend_title = "Celltype")

tme_response_bp <- create_barplot(seurat_object = immune_fibro,
                                  group_var = "binary_response",
                                  plot_var = "celltype",
                                  plot_levels = levels(immune_fibro$celltype),
                                  group_levels = c("CR_SD", "PD"),
                                  plot_colors = immune_fibro_celltype_col,
                                  var_names =  c("Frequency (%)", ""),
                                  legend_title = "Celltype")

tme_cd3_bp <- create_barplot(seurat_object = immune_fibro,
                             group_var = "CD3_high_low",
                             plot_var = "celltype",
                             plot_levels = levels(immune_fibro$celltype),
                             group_levels = c("High", "Low"),
                             plot_colors = immune_fibro_celltype_col,
                             var_names =  c("Frequency (%)", ""),
                             legend_title = "Celltype")

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/tumor_tme_celltype_prop_barplots.pdf"
pdf(file = filename,
    width = 8,
    height = 5)

tumor_response_bp | tme_response_bp | tme_cd3_bp

dev.off()


# Cluster proportion scatterplots

# Full object
# Cell type proportions by CD3 score
# Inputs
# seurat_object = Seurat object
# group_var = e.g. CD3_high_low
# group1 = e.g. High
# group2 = e.g. Low
# plot_var = e.g. celltype
# plot_colors = A named vector of colors to use in the plot, corresponding to plot_var (e.g. cell type colors)
# var_names = Used as axis titles, c("group2", "group1")
# legend_title = Legend title, ("" for no title)
p1 <- create_clusterpropplot(seurat_object = tumors,
                             group_var = "CD3_high_low",
                             group2 = "High",
                             group1 = "Low",
                             plot_var = "celltype2",
                             plot_colors = tumor_celltype2_col,
                             var_names = c("CD3 low", "CD3 high"),
                             legend_title = "") + manuscript_theme

p2 <- create_clusterpropplot(seurat_object = tumors,
                             group_var = "binary_response",
                             group2 = "CR_SD",
                             group1 = "PD",
                             plot_var = "celltype2",
                             plot_colors = tumor_celltype2_col,
                             var_names = c("PD", "CR/SD"),
                             legend_title = "") + manuscript_theme

p1 | p2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_celltype_prop_scatterplot.pdf"
pdf(file = filename,
    width = 5.5,
    height = 2.75)

p1 | p2

dev.off()

# Immune + fibroblast object
p1 <- create_clusterpropplot(seurat_object = immune_fibro,
                             group_var = "CD3_high_low",
                             group1 = "High",
                             group2 = "Low",
                             plot_var = "celltype",
                             plot_colors = immune_fibro_celltype_col,
                             var_names = c("CD3 high", "CD3 low"),
                             legend_title = "") + manuscript_theme

p2 <- create_clusterpropplot(seurat_object = immune_fibro,
                             group_var = "binary_response",
                             group1 = "CR_SD",
                             group2 = "PD",
                             plot_var = "celltype",
                             plot_colors = immune_fibro_celltype_col,
                             var_names = c("CR/SD", "PD"),
                             legend_title = "") + manuscript_theme

p1 | p2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_immune_fibro_celltype_prop_scatterplot.pdf"
pdf(file = filename,
    width = 5.5,
    height = 2.75)

p1 | p2

dev.off()

# Using scProportionTest
prop_test <- sc_utils(tumors)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "High", sample_2 = "Low",
  sample_identity = "CD3_high_low")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_celltype_prop_forest_CD3.pdf"
pdf(file = filename,
    width = 2.5,
    height = 4)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# Response
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "CR_SD", sample_2 = "PD",
  sample_identity = "binary_response")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_celltype_prop_forest_response.pdf"
pdf(file = filename,
    width = 2.5,
    height = 4)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# For immune+fibroblasts
prop_test <- sc_utils(immune_fibro)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "High", sample_2 = "Low",
  sample_identity = "CD3_high_low")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_immune_fibro_celltype_prop_forest_CD3.pdf"
pdf(file = filename,
    width = 1.8,
    height = 4)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# Response
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "CR_SD", sample_2 = "PD",
  sample_identity = "binary_response")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_immune_fibro_celltype_prop_forest_response.pdf"
pdf(file = filename,
    width = 1.8,
    height = 4)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# Violinplot of cell type proportions
prop_table <- as.data.frame(table(immune_fibro$celltype, immune_fibro$UPN))
colnames(prop_table) <- c("Celltype", "UPN", "Freq")
prop_table <- spread(prop_table, Celltype, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- pivot_longer(prop_table, cols = names(prop_table)[2:length(names(prop_table))],
                           names_to = "Celltype", values_to = "Freq")
prop_table$CD3 <- mapvalues(prop_table$UPN,
                            from = immune_fibro$UPN,
                            to = immune_fibro$CD3_high_low)

prop_table$Celltype <- factor(prop_table$Celltype,
                              levels = c(paste0("M", seq(1, 11, by = 1)),
                                         paste0("L", seq(1, 10, by = 1)),
                                         paste0("F", seq(1, 3, by = 1))))

celltype_prop_vln <- ggplot(prop_table, aes(x = CD3, y = Freq, fill = Celltype)) +
  geom_violin() + 
  geom_jitter(alpha = 1, size = 0.6) +
  scale_fill_manual(name = "Celltype", values = immune_fibro_celltype_col, guide="none") + 
  theme_bw() +
  facet_wrap(~Celltype) +
  ylab("Cell type proportion (%)")
  xlab("CD3 score") +
  theme(panel.background = element_rect(colour = "gray99")) +
  theme(legend.position = "none")

celltype_prop_vln             

#========================
# Differential expression violin plots
#========================

vlnplot_de_features <- c("APOE", "C1QA", "C1QB", "C1QC", "COL1A1",
                         "S100A8", "S100A9", "SPP1" )

VlnPlot(immune_fibro,
        group.by = "celltype",
        split.by = "CD3_high_low",
        features = vlnplot_de_features,
        pt.size = 0,
        ncol = 2,
        cols = CD3_high_low_col) &
  theme_minimal() &
  manuscript_theme

DotPlot(immune_fibro,
        group.by = "celltype",
        split.by = "CD3_high_low",
        features = vlnplot_de_features,
        #pt.size = 0,
        #ncol = 2,
        cols = c("azure", "tomato3"))
#theme_minimal() &
#manuscript_theme

#========================
# SPP1+ myeloid clusters vs. others
#========================

unique(immune_fibro$celltype)
myeloid <- subset(immune_fibro, subset = celltype %in% paste0("M", seq(1, 9)))

# Top markers for each cluster
markers <- presto::wilcoxauc(myeloid,
                             group_by = "celltype",
                             assay = "data",
                             seurat_assay = "RNA")
#top_markers <- markers %>%  group_by(group) %>% slice_max(order_by = auc, n = 4)

sig_markers <- markers %>% filter(auc>0.6, padj<0.01)


