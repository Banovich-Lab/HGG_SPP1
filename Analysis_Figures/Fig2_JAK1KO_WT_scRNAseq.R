#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 01/29/2024
# Description: Replotting the JAK1KO scRNAseq data
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
library(survival)

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

# Seurat data
mm_ko_cells <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/CellChat/Jak1_KO_Control_Cells.rds")
mm_wt_cells <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/CellChat/WildType_Control_Cells.rds")

combined <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/scRNA_seq_Pool.combined.rds")
combined <- scRNA_seq_Pool.combined

# Object after removing unwanted clusters
# ko_wt_data <- readRDS("/scratch/hnatri/CART/ko_wt_data.rds")

combined_CD45pos <- subset(combined, subset = CD45_Group == "CD45_Pos")

#ko_wt_data <- merge(Jak1_KO_Control_Cells, WildType_Control_Cells, merge.dr = c("umap"))
ko_wt_data <- combined_CD45pos

# Excluding MT genes
assay_data <- LayerData(ko_wt_data, layer = "counts", assay = "RNA")
assay_data <- assay_data[-grep("mt-", rownames(assay_data)),]
ko_wt_data <- subset(ko_wt_data, features = rownames(assay_data))

# Cluster markers
#markers <- presto::wilcoxauc(ko_wt_data,
#                             group_by = "celltype",
#                             assay = "data",
#                             seurat_assay = "RNA")
#
#top_markers <- markers %>% group_by(group) %>% slice_max(order_by = auc, n = 50)
#write.table(top_markers, "/home/hnatri/CART/13384_Tumors/JAK1KO_WT_CART_CTRL_top50_markers.tsv", quote = F, row.names = F, sep = "\t")

# Converting mouse gene names to human
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")

convert_mouse_to_human <- function(gene_list){
  gene_names <- as.data.frame(matrix(nrow = length(gene_list),
                                     ncol = 2))
  colnames(gene_names) <- c("mouse", "human")
  rownames(gene_names) <- gene_list
  gene_names$mouse <- gene_list
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      
      if(length(human_genes)==0){
        gene_names[gene, "human"] <- NA
      } else if (length(human_genes)>1){
        #  human_genes <- paste0(human_genes, collapse = ", ")
        bind_df <- data.frame("mouse" = rep(gene, times = length(human_genes)),
                              "human" = human_genes)
        gene_names <- rbind(gene_names, bind_df)
      } else {
        gene_names[gene, "human"] <- human_genes
      }
    }
  }
  return(gene_names)
}

gene_names <- convert_mouse_to_human(rownames(ko_wt_data@assays$RNA))

# Keeping mouse genes with a single human ortholog
gene_names <- gene_names %>%
  group_by(mouse) %>%
  filter(!is.na(human),
         n() == 1) %>%
  ungroup()

DefaultAssay(ko_wt_data) <- "RNA"
assay_data <- GetAssayData(ko_wt_data, slot = "counts")
assay_data <- assay_data[which(rownames(assay_data) %in% gene_names$mouse),]
new_names <- rownames(assay_data)
new_names <- mapvalues(x = new_names,
                       from = gene_names$mouse,
                       to = gene_names$human)
rownames(assay_data) <- new_names

ko_wt_data[["RNA_human"]] <- CreateAssayObject(assay_data,
                                               min.cells = 0,
                                               min.features = 0)
DefaultAssay(ko_wt_data) <- "RNA_human"

# Updating cell type annotations
gs4_deauth()
cluster_annotations  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
cluster_annotations <- read_sheet(cluster_annotations, sheet = "Cluster annotations, JAK mouse")

ko_wt_data$celltype <- mapvalues(ko_wt_data$seurat_clusters,
                                 from = cluster_annotations$cluster,
                                 to = cluster_annotations$annotation)

# Excluding Neut4
ko_wt_data <- subset(ko_wt_data, subset = celltype == "Neu4", invert = T)

#saveRDS(ko_wt_data, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/SPP1_release/JAK1KO_WT_CART_CTRL.rds")
#ko_wt_data <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/SPP1_release/JAK1KO_WT_CART_CTRL.rds")

# Dropping epithelial cells
ko_wt_data <- subset(ko_wt_data, subset = celltype %in% paste0("Epi", seq(1, 5)), invert = T)

#bp <- create_barplot(seurat_object = combined,
#                     plot_var = "celltype",
#                     group_var = "CART_JAK1_Group",
#                     group_levels = unique(combined$CART_JAK1_Group),
#                     plot_levels = sort(unique(combined$celltype)),
#                     plot_colors = jak1_celltype_col,
#                     var_names = c("ct", "group"),
#                     legend_title = "")

# Reclustering
ko_wt_data <- SCTransform(ko_wt_data,
                          vars.to.regress = c("percent.mt"),
                          vst.flavor = "v2")
ko_wt_data <- RunPCA(ko_wt_data,
                     reduction.name = "pca",
                     verbose = F)
pcs <- get_pcs(ko_wt_data, reduction_name = "pca")
ko_wt_data <- RunUMAP(ko_wt_data,
                      reduction = "pca",
                      reduction.name = "umap",
                      dims = 1:min(pcs),
                      return.model = TRUE)
ko_wt_data <- FindNeighbors(ko_wt_data,
                            reduction = "pca",
                            dims = 1:min(pcs),
                            graph.name = c("nn", "snn"))
# resolution 0.2
ko_wt_data <- FindClusters(ko_wt_data,
                           resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                           graph.name = "nn")

saveRDS(ko_wt_data, "/scratch/hnatri/CART/ko_wt_data.rds")
q(save = "no")

#ko_wt_data <- readRDS("/scratch/hnatri/CART/ko_wt_data.rds")
#ko_wt_data <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/JAK1KO_scRNAseq/JAK1KO_WT_CTRL_immune_stroma.rds")

#==============================================================================#
# Figure 2: JAK1KO vs. WT scRNA-seq
#==============================================================================#

#========================
# Quality metrics
#========================

f1 <- FeaturePlot(ko_wt_data,
                  features = c("nCount_RNA",
                               "nFeature_RNA",
                               "percent.mt"),
                  #group.by = "celltype",
                  cols = c("gray89", "tomato3"),
                  order = T,
                  #keep.scale = "all",
                  ncol = 1,
                  raster = T,
                  raster.dpi = c(2048, 1024),
                  pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  manuscript_theme &
  NoAxes()
#NoLegend()
#theme(legend.key.size = unit(0.3, 'cm')) # change legend key size

f1

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_QC_FeaturePlot.pdf"
pdf(file = filename,
    width = 4,
    height = 6)

f1

dev.off()

# Violinplots
ko_wt_data$celltype <- factor(ko_wt_data$celltype, levels = sort(as.character(unique(ko_wt_data$celltype))))
v1 <- VlnPlot(ko_wt_data,
              features = c("nCount_RNA",
                           "nFeature_RNA",
                           "percent.mt"),
              group.by = "celltype",
              pt.size = 0,
              cols = jak1_celltype_col,
              ncol = 1) &
  theme_classic() &
  manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

v1

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_QC_vln.pdf"
pdf(file = filename,
    width = 4,
    height = 5.5)

v1

dev.off()


#========================
# UMAPs
#========================

umap_celltype <- DimPlot(ko_wt_data,
                         group.by = "celltype",
                         cols = jak1_celltype_col,
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

umap_celltype

filename <- "/home/hnatri/CART/13384_Tumors/Plots/JAK1KO_WT_CART_CTRL_celltype_UMAP.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

umap_celltype

dev.off()

umap_celltype_ko <- DimPlot(subset(ko_wt_data, subset = JAK1_Group == "Jak1_KO"),
                         group.by = "celltype",
                         #split.by = "JAK1_Group",
                         cols = jak1_celltype_col,
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

umap_celltype_wt <- DimPlot(subset(ko_wt_data, subset = JAK1_Group == "WildType"),
                            group.by = "celltype",
                            #split.by = "JAK1_Group",
                            cols = jak1_celltype_col,
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

umap_celltype_ko + umap_celltype_wt

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CART_CTRL_split_celltype_UMAP.pdf"
pdf(file = filename,
    width = 6,
    height = 3)

umap_celltype_ko + umap_celltype_wt

dev.off()

#========================
# Cell type proportions
#========================

ko_wt_CART <- subset(ko_wt_data, subset = CART_Group == "CART")
ko_wt_CTRL <- subset(ko_wt_data, subset = CART_Group == "Control")
ko_wt_JAK1KO <- subset(ko_wt_data, subset = JAK1_Group == "Jak1_KO")
ko_wt_WT <- subset(ko_wt_data, subset = JAK1_Group == "WildType")

# Immune + fibroblast object
p1 <- create_clusterpropplot(seurat_object = ko_wt_data,
                             group_var = "JAK1_Group",
                             group1 = "Jak1_KO",
                             group2 = "WildType",
                             plot_var = "celltype",
                             plot_colors = jak1_celltype_col,
                             var_names = c("JAK1 KO", "WT"),
                             legend_title = "") + manuscript_theme

p1

filename <- "/home/hnatri/CART/13384_Tumors/Plots/JAK1KOWT_celltypeprop_scatterplot.pdf"
pdf(file = filename,
    width = 2.75,
    height = 2.75)

p1

dev.off()

# Using scProportionTest
prop_test <- sc_utils(ko_wt_CART)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WildType", sample_2 = "Jak1_KO",
  sample_identity = "JAK1_Group")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CART_celltypeprop_forestplot.pdf"
pdf(file = filename,
    width = 2.25,
    height = 3.5)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

prop_test <- sc_utils(ko_wt_CTRL)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WildType", sample_2 = "Jak1_KO",
  sample_identity = "JAK1_Group")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CTRL_celltypeprop_forestplot.pdf"
pdf(file = filename,
    width = 2.25,
    height = 3.5)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# JAK1KO
prop_test <- sc_utils(ko_wt_JAK1KO)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "CART", sample_2 = "Control",
  sample_identity = "CART_Group")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_CART_CTRL_celltypeprop_forestplot.pdf"
pdf(file = filename,
    width = 2.25,
    height = 3.5)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# WT
prop_test <- sc_utils(ko_wt_WT)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "CART", sample_2 = "Control",
  sample_identity = "CART_Group")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/WT_CART_CTRL_celltypeprop_forestplot.pdf"
pdf(file = filename,
    width = 2.25,
    height = 3.5)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

# Violinplot of cell type proportions
#prop_table <- as.data.frame(table(immune_fibro$celltype, immune_fibro$UPN))
#colnames(prop_table) <- c("Celltype", "UPN", "Freq")
#prop_table <- spread(prop_table, Celltype, Freq)
## Converting to percetange
#prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
#prop_table <- pivot_longer(prop_table, cols = names(prop_table)[2:length(names(prop_table))],
#                           names_to = "Celltype", values_to = "Freq")
#prop_table$CD3 <- mapvalues(prop_table$UPN,
#                            from = immune_fibro$UPN,
#                            to = immune_fibro$CD3_high_low)
#
#prop_table$Celltype <- factor(prop_table$Celltype,
#                              levels = c(paste0("M", seq(1, 11, by = 1)),
#                                         paste0("L", seq(1, 10, by = 1)),
#                                         paste0("F", seq(1, 3, by = 1))))
#
#celltype_prop_vln <- ggplot(prop_table, aes(x = CD3, y = Freq, fill = Celltype)) +
#  geom_violin() + 
#  geom_jitter(alpha = 1, size = 0.6) +
#  scale_fill_manual(name = "Celltype", values = immune_fibro_celltype_col, guide="none") + 
#  theme_bw() +
#  facet_wrap(~Celltype) +
#  ylab("Cell type proportion (%)")
#xlab("CD3 score") +
#  theme(panel.background = element_rect(colour = "gray99")) +
#  theme(legend.position = "none")
#
#celltype_prop_vln

#========================
# Feature expression
#========================

DefaultAssay(ko_wt_data) <- "RNA_human"

ko_wt_data$nCount_RNA_log <- log(ko_wt_data$nCount_RNA)
ko_wt_data$nFeature_RNA_log <- log(ko_wt_data$nFeature_RNA)

# Quality metrics
p1 <- FeaturePlot(ko_wt_data,
                  features = c("percent.mt",
                               "nCount_RNA_log",
                               "nFeature_RNA_log"),
                  ncol = 3,
                  order = T,
                  cols = c("gray89", "tomato3"),
                  raster = T,
                  raster.dpi = c(2048, 1024),
                  pt.size = 4) &
  coord_fixed(ratio=1) &
  theme_classic() &
  manuscript_theme

filename <- "/home/hnatri/CART/13384_Tumors/Plots/JAK1KO_WT_CTRL_QC_features.pdf"
pdf(file = filename,
    width = 7.5,
    height = 2.5)

p1

dev.off()

# Plot features
plot_features <- c("CD3D", "CD4", "CD8A",
                   "ITM2A", "CEMIP2", "BTG1",
                   "IL32", "PFN1", "NKG7",
                   "TMSB4X", "KLRD1",
                   "TPT1", "CD96", "IFITM2",
                   "HMBG1", "PTMA", "TMSB10",
                   "CFL1", "CD68", "MIF",
                   "NLRP3", "APOE", "FUT4",
                   "FTH", "FTL1", "AREG",
                   "GZMB", "CCL4", "CORO1A",
                   "C1QA", "C1QB", "C1QC",
                   "TYROBP", "CD79A", "SPP1",
                   "PTPRC", "SELL",  "ICAM1", 
                   "ARG1", "S100A8", "S100A9",
                   "S100A10", "CD274", "COL1A1",
                   "ACTA2", "PDGFRB", "FABP5",
                   "PDPN", "IL1B", "CD74")

# Dotplot
ko_wt_dotplot <- DotPlot(ko_wt_data,
                         features = sort(plot_features),
                         group.by = "celltype",
                         split.by = "JAK1_Group",
                         cols = c("gray89", "tomato3")) +
  #coord_flip() +
  manuscript_theme +
  RotatedAxis()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/ko_wt_dotplot.pdf"
pdf(file = filename,
    width = 4,
    height = 7)

ko_wt_dotplot

dev.off()

# seurat_object, plot_features, group_var, group_colors, column_title, km=5, row.order = NULL
DefaultAssay(ko_wt_data) <- "RNA_human"
ko_wt_dotplot <- create_dotplot_heatmap(seurat_object = ko_wt_data,
                                        plot_features = plot_features,
                                        group_var = "celltype",
                                        group_colors = jak1_celltype_col,
                                        column_title = "",
                                        col_km = 5,
                                        row.order = plot_features)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_canonical_dotplot_heatmap.pdf"
pdf(file = filename,
    width = 8,
    height = 6.5)

ko_wt_dotplot

dev.off()

subset_data <- subset(ko_wt_data, subset = JAK1_Group == "Jak1_KO")
DefaultAssay(subset_data) <- "RNA_human"
VariableFeatures(subset_data) <- plot_features
subset_data <- NormalizeData(subset_data)
subset_data <- ScaleData(subset_data)

ko_dotplot <- create_dotplot_heatmap(seurat_object = subset_data,
                                     plot_features = plot_features,
                                     group_var = "celltype",
                                     group_colors = jak1_celltype_col,
                                     column_title = "",
                                     col_km = 5,
                                     row.order = plot_features)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_canonical_dotplot_heatmap.pdf"
pdf(file = filename,
    width = 6,
    height = 6.5)

ko_dotplot

dev.off()

subset_data <- subset(ko_wt_data, subset = JAK1_Group == "WildType")
DefaultAssay(subset_data) <- "RNA_human"
VariableFeatures(subset_data) <- plot_features
subset_data <- NormalizeData(subset_data)
subset_data <- ScaleData(subset_data)

wt_dotplot <- create_dotplot_heatmap(seurat_object = subset_data,
                                     plot_features = plot_features,
                                     group_var = "celltype",
                                     group_colors = jak1_celltype_col,
                                     column_title = "",
                                     col_km = 5,
                                     row.order = plot_features)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/WT_canonical_dotplot_heatmap.pdf"
pdf(file = filename,
    width = 6,
    height = 6.5)

wt_dotplot

dev.off()

# Violin plots by cell type
VlnPlot(ko_wt_data,
        group.by = "celltype",
        features = sort(plot_features),
        pt.size = 0,
        cols = jak1_celltype_col,
        ncol = 6,
        log = T) &
  theme_bw() &
  manuscript_theme &
  NoLegend() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

VlnPlot(ko_wt_data,
        group.by = "celltype",
        split.by = "JAK1_Group",
        features = "SPP1",
        split.plot = T,
        pt.size = 0,
        cols = c("deeppink3", "aquamarine3"),
        log = T) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ko_wt_data$celltype <- factor(ko_wt_data$celltype, levels = sort(as.character(unique(ko_wt_data$celltype))))

#DefaultAssay(ko_wt_data) <- "RNA_human"
#ko_wt_data <- NormalizeData(ko_wt_data)

v11 <- VlnPlot(subset(ko_wt_data, CART_Group == "CART"),
        group.by = "celltype",
        split.by = "JAK1_Group",
        features = c("SPP1", "CD44"),
        split.plot = T,
        pt.size = 0,
        ncol = 1,
        cols = c("deeppink3", "aquamarine3"),
        slot = "data") &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CARtreated_SPP1_CD44_vlnplot.pdf"
pdf(file = filename,
    width = 7,
    height = 4)

v11

dev.off()

# Top markers for each cluster
markers <- presto::wilcoxauc(ko_wt_data,
                             group_by = "celltype",
                             assay = "data",
                             seurat_assay = "RNA_human")
top_markers <- markers %>%  group_by(group) %>% slice_max(order_by = auc, n = 4)

ko_wt_dotplot <- create_dotplot_heatmap(seurat_object = ko_wt_data,
                                        plot_features = unique(top_markers$feature),
                                        group_var = "celltype",
                                        group_colors = jak1_celltype_col,
                                        column_title = "")

filename <- "/home/hnatri/CART/13384_Tumors/Plots/JAK1KO_WT_topmarker_dotplot_heatmap_smaller.pdf"
pdf(file = filename,
    width = 7,
    height = 9)

ko_wt_dotplot

dev.off()

# SPP1
v1 <- VlnPlot(ko_wt_CART,
        group.by = "celltype",
        split.by = "JAK1_Group",
        features = "SPP1",
        assay = "RNA_human",
        slot = "data",
        split.plot = T,
        pt.size = 0,
        cols = c("#931A1D", "#3E58A8"),
        log = T) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) &
  ggtitle("SPP1, CAR T")

v2 <- VlnPlot(ko_wt_CTRL,
              group.by = "celltype",
              split.by = "JAK1_Group",
              features = "SPP1",
              assay = "RNA_human",
              slot = "data",
              split.plot = T,
              pt.size = 0,
              cols = c("#931A1D", "#3E58A8"),
              log = T) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) &
  ggtitle("SPP1, CTRL")

v1 / v2

v3 <- VlnPlot(ko_wt_JAK1KO,
              group.by = "celltype",
              split.by = "CART_Group",
              features = "SPP1",
              assay = "RNA_human",
              slot = "data",
              split.plot = T,
              pt.size = 0,
              cols = c("#931A1D", "#3E58A8"),
              log = T) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) &
  ggtitle("SPP1, JAK1KO")

v4 <- VlnPlot(ko_wt_WT,
              group.by = "celltype",
              split.by = "CART_Group",
              features = "SPP1",
              assay = "RNA_human",
              slot = "data",
              split.plot = T,
              pt.size = 0,
              cols = c("#931A1D", "#3E58A8"),
              log = T) &
  theme_bw() &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) &
  ggtitle("SPP1, WT")

v3 / v4

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CART_CTRL_SPP1_vln.pdf"
pdf(file = filename,
    width = 7,
    height = 8)

v1 / v2 / v3 / v4

dev.off()

# Differential expression
ko_wt_CART <- subset(ko_wt_data, subset = CART_Group == "CART")
unique(ko_wt_CART$JAK1_Group)
DefaultAssay(ko_wt_CART) <- "RNA_human"

DEG_JAK1 <- lapply(unique(ko_wt_CART$celltype), function(xx){
  #message(xx)
  data_subset <- subset(ko_wt_CART, subset = celltype == xx)
  Idents(data_subset) <- data_subset$JAK1_Group
  
  if (min(table(data_subset$JAK1_Group))<10){
    return(NULL)
  }
  
  if (all((c("WildType", "Jak1_KO") %in% data_subset$JAK1_Group) == c(T, T))){
    markers <- FindMarkers(data_subset,
                           ident.1 = "WildType",
                           ident.2 = "Jak1_KO",
                           assay = "RNA_human",
                           verbose = F)
    markers$feature <- rownames(markers)
    markers$celltype <- xx
    
    return(markers)
  } else {
    return(NULL)
  }
})
names(DEG_JAK1) <- unique(ko_wt_CART$celltype)
DEG_JAK1[sapply(DEG_JAK1, is.null)] <- NULL

DEG_JAK1_df <- as.data.frame(do.call(rbind, DEG_JAK1))

# Distribution of log2FC
hist(DEG_JAK1_df$avg_log2FC)

DEG_JAK1_df_sig <- DEG_JAK1_df %>%
  filter(p_val_adj < 0.01,
         abs(avg_log2FC) > 2)

spp1_cd44 <- DEG_JAK1_df %>% filter(feature %in% c("SPP1", "CD44")) %>%
  arrange(celltype)

# Saving to a file
write.table(spp1_cd44,
            "/home/hnatri/13384_CART/13384_Tumors/JAK1KO_WT_SPP1_CD44_DE.tsv",
            sep = "\t", quote = F, row.names = F)
