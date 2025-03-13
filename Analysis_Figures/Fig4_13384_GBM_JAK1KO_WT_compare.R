#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 03/07/2024
# Description: GBM and JAK1KO/WT expression correlation
#==============================================================================

suppressMessages({library(Seurat)
                  library(ggplot2)
                  library(data.table)
                  library(dplyr)
                  library(patchwork)
                  library(tidyr)
                  library(stringr)
                  library(RColorBrewer)
                  library(ComplexHeatmap)
                  library(circlize)
                  library(gridExtra)
                  library(ggpubr)
                  library(qdap)
                  library(plyr)
                  library(viridis)
                  library(colourvalues)
                  library(googlesheets4)
                  library(ggrepel)})

#==============================================================================
# Helper functions
#==============================================================================

# Colors and plot functions
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "umap"

#==============================================================================
# Import Seurat objects
#==============================================================================

integrated_seurat <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/Immune_Combined_Human_Mouse_RPCA_SCT.rds")

DefaultAssay(immune.combined.sct) <- "RNA"
seurat_object <- immune.combined.sct

# Adding cell type annotations
seurat_object$celltype <- plyr::mapvalues(x = seurat_object$seurat_clusters,
                                          from = celltype_annot_integrated$cluster,
                                          to = celltype_annot_integrated$annotation)

# Additional colors
tumor_clusters <- as.factor(c(0, seq(1:max(as.numeric(seurat_object$seurat_clusters)))))
tumor_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(tumor_clusters))
names(tumor_cluster_col) <- levels(tumor_clusters)
tumor_cluster_col_c <- tumor_cluster_col
names(tumor_cluster_col_c) <- paste0("C", names(tumor_cluster_col_c))

# WT and GBM
wt_gbm_seurat <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Integration_MmWT_HsGBM/Immune_Combined_Human_Mouse_RPCA_SCT.rds")
wt_gbm_seurat <- immune.combined.sct

# GBM PD and JAK1KO mouse, pretreatment
ko_gbm_pd_seurat <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Integration_Mm_KO_Pretreatment_Hs_PD/Integration_Mouse_KO_pretreatment_Human_PD.rds")
ko_gbm_pd_seurat <- Integration_Mouse_KO_pretreatment_Human_PD

#==============================================================================
# Marker gene expression heatmaps
#==============================================================================

# Canonical markers
gs4_deauth()
canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
canonical_markers_all <- read_sheet(canonical_markers, sheet = "Heatmap genes")

canonical_markers_all <- canonical_markers_all %>% 
  mutate(markers = strsplit(as.character(markers), ", ")) %>% 
  tidyr::unnest(markers)

# Adding cell type annotations from the GBM object
immune_fibro <- immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

immune_fibro@meta.data$cellname <- paste0("Hs_", rownames(immune_fibro@meta.data))
seurat_object$GBM_celltype <- mapvalues(rownames(seurat_object@meta.data),
                                        from = immune_fibro@meta.data$cellname,
                                        to = immune_fibro@meta.data$celltype)
seurat_object$GBM_celltype[-which(seurat_object$GBM_celltype %in% immune_fibro$celltype)] <- NA

DimPlot(seurat_object,
        group.by = "GBM_celltype",
        reduction = reduction,
        label=T,
        cols = immune_fibro_celltype_col) + NoLegend()

DimPlot(seurat_object,
        group.by = "celltype",
        split.by = "GBM_celltype",
        reduction = reduction,
        #label=T,
        cols = celltype_annot_integrated_col,
        ncol = 5) & NoLegend()

DimPlot(seurat_object,
        group.by = "GBM_celltype",
        split.by = "celltype",
        reduction = reduction,
        #label=T,
        cols = immune_fibro_celltype_col,
        ncol = 5) & NoLegend()

# Downsampling for plotting (Seurat DoHeatmap uses ggplot, doesn't work with
# >30k cells. downsample = # of cells per identity class. Same problem with
# ComplexHeatmap).
DimPlot(seurat_object, group.by = "seurat_clusters", reduction = reduction, label=T) + NoLegend()
VlnPlot(seurat_object, group.by = "seurat_clusters", features = c("nCount_RNA", "nFeature_RNA"), pt.size=0)

# 30000
seurat_object_downsample <- subset(seurat_object, downsample = 30000/length(unique(seurat_object@meta.data$seurat_clusters)))

DimPlot(seurat_object_downsample, group.by = "seurat_clusters", reduction = "umap", label=T) + NoLegend()

# Rescaling all features
DefaultAssay(seurat_object_downsample) <- "RNA"
VariableFeatures(seurat_object_downsample) <- rownames(seurat_object_downsample)
seurat_object_downsample <- NormalizeData(seurat_object_downsample)
seurat_object_downsample <- ScaleData(seurat_object_downsample)
#vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

# Top markers for each cluster, using the full object
cluster_markers_rna <- presto::wilcoxauc(seurat_object,
                                         group_by = "seurat_clusters",
                                         #groups_use = "",
                                         assay = "data",
                                         seurat_assay = "RNA")

markers <- presto::top_markers(cluster_markers_rna, n = 10, auc_min = 0.5, pct_in_min = 50, pct_out_max = 100)

top_cluster_markers_rna <- markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

# Top markers and canonical markers
#all_plot_markers <- unique(c(top_cluster_markers_rna, setdiff(canonical_markers_all$markers, c(NA))))
all_plot_markers <- canonical_markers_all$markers

# Gene expression data from the RNA assay
cellpop_rna_all_markers <- LayerData(seurat_object_downsample, layer = "scale.data", assay = "RNA")
cellpop_rna_all_markers <- cellpop_rna_all_markers[which(rownames(cellpop_rna_all_markers) %in% all_plot_markers),]

# Order by cluster, then by response
cellpop_cells_metadata <- data.frame("cell" = rownames(seurat_object_downsample@meta.data),
                                     #"cluster" = seurat_object_downsample@meta.data$seurat_clusters,
                                     "celltype" = seurat_object_downsample@meta.data$celltype,
                                     "dataset" = seurat_object_downsample@meta.data$Species)

cellpop_cells_metadata <- cellpop_cells_metadata[order(cellpop_cells_metadata$celltype, cellpop_cells_metadata$dataset),]
cellpop_rna_all_markers <- as.data.frame(cellpop_rna_all_markers)[cellpop_cells_metadata$cell]

# RNA heatmap
col <- list()
col$celltype <- celltype_annot_integrated_col
col_ha <- HeatmapAnnotation(
  df = data.frame(celltype = cellpop_cells_metadata$celltype,
                  dataset = cellpop_cells_metadata$dataset),
  annotation_height = unit(4, "mm"),
  col = col
  #col = c(col, cont_col)
)

# Only annotating selected genes
labels <- unique(canonical_markers_all$markers)
label_positions <- grep(paste(labels, collapse = "|"), rownames(cellpop_rna_all_markers))
row_ha_labels <- rowAnnotation(link = anno_mark(at = label_positions, 
                                                labels = labels, 
                                                labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

# Quantile-normalizing rows and columns for plotting
cellpop_rna_all_markers_qqnorm <- t(apply(cellpop_rna_all_markers, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))

cellpop_rna_all_markers_qqnorm <- (apply(cellpop_rna_all_markers_qqnorm, 2, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(cellpop_rna_all_markers_qqnorm) <- colnames(cellpop_rna_all_markers)
rownames(cellpop_rna_all_markers_qqnorm) <- rownames(cellpop_rna_all_markers)

# Heatmap colors
col_fun2 <- colorRamp2(quantile(as.matrix(cellpop_rna_all_markers_qqnorm), seq(0, 1, by = 0.25)), viridis(5))

plot_func <- function(){
  hm_rna  <- Heatmap(as.matrix(cellpop_rna_all_markers_qqnorm), 
                     name = "Gene exp", # Title of legend
                     col = viridis(100),
                     use_raster = T,
                     column_title = NULL, row_title = NULL,
                     show_column_names = FALSE,
                     show_row_names = T,
                     cluster_rows = TRUE,
                     row_km = 10,
                     cluster_columns = T,
                     top_annotation = col_ha,
                     height = nrow(cellpop_rna_all_markers)*unit(4, "mm"),
                     row_names_gp = gpar(fontsize = 7),  # Text size for row names
                     row_names_side = "right")
  
  heatmap <- draw(hm_rna)
}

p <- plot_func()

# Saving to a file
filename <- "/home/hnatri/CART/Tumors/GBM_JAK1KO_heatmap_canonical_clustered_cells.pdf"
pdf(file = filename,
    width = 10,
    height = 10)
p
dev.off()

png(file="/home/hnatri/CART/Tumors/GBM_JAK1KO_heatmap_canonical_clustered_cells.png",
    width=1000, height=1000)
p
dev.off()

#==============================================================================
# Pseudobulk
#==============================================================================

DefaultAssay(seurat_object) <- "RNA"
VariableFeatures(seurat_object) <- canonical_markers_all$markers
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object)

seurat_object$celltype_dataset <- paste0(seurat_object$celltype, "-", seurat_object$Species)
#pseudobulk <- AverageExpression(seurat_object, group.by = "celltype_dataset", assay = "RNA", layer = "data")
#pseudobulk <- pseudobulk$RNA
#pseudobulk <- pseudobulk[which(rownames(pseudobulk) %in% canonical_markers_all$markers),]
#
#pseudobulk <- pseudobulk[!is.infinite(rowSums(pseudobulk)),]
#colnames(pseudobulk) <- as.character(colnames(pseudobulk))
lognorm_counts <- LayerData(seurat_object,
                            assay = "RNA",
                            layer = "data",
                            features = intersect(rownames(seurat_object), canonical_markers_all$markers))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$celltype_dataset <- seurat_object$celltype_dataset

lognorm_counts <- lognorm_counts %>% group_by(celltype_dataset) %>%
  summarise_at(vars(-group_cols()), mean) %>%
  as.data.frame()

lognorm_counts <- t(lognorm_counts) %>%
  as.data.frame()

colnames(lognorm_counts) <- lognorm_counts[1,]
lognorm_counts <- lognorm_counts[2:nrow(lognorm_counts),]
pseudobulk <- lognorm_counts
pseudobulk <- sapply(pseudobulk, as.numeric)
rownames(pseudobulk) <- rownames(lognorm_counts)

# Annotations
pseudobulk_metadata <- data.frame("celltype" = sapply(strsplit(colnames(pseudobulk),"-"), `[`, 1),
                                  "dataset" = sapply(strsplit(colnames(pseudobulk),"-"), `[`, 2))

pseudobulk_metadata$dataset <- gsub("Hs", "GBM", pseudobulk_metadata$dataset)
pseudobulk_metadata$dataset <- gsub("Mm", "JAK1KO", pseudobulk_metadata$dataset)

col <- list()
col$celltype <- celltype_annot_integrated_col
col$dataset <- dataset_col
row_ha <- rowAnnotation(
  df = data.frame(celltype = pseudobulk_metadata$celltype,
                  dataset = pseudobulk_metadata$dataset),
  annotation_height = unit(3, "mm"),
  annotation_width = unit(3, "mm"),
  col = col,
  annotation_name_gp = gpar(fontsize = 7)
)

# Quantile-normalizing rows and columns for plotting
pseudobulk_qqnorm <- t(apply(pseudobulk, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
pseudobulk_qqnorm <- (apply(pseudobulk_qqnorm, 2, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(pseudobulk_qqnorm) <- colnames(pseudobulk)
rownames(pseudobulk_qqnorm) <- rownames(pseudobulk)

colnames(pseudobulk_qqnorm) <- gsub("Mm", "JAK1KO", colnames(pseudobulk_qqnorm))
colnames(pseudobulk_qqnorm) <- gsub("Hs", "GBM", colnames(pseudobulk_qqnorm))

pseudobulk_metadata$celltype_dataset <- paste0(pseudobulk_metadata$celltype, "_",
                                               pseudobulk_metadata$dataset)

plot_func <- function(){
  hm_rna  <- ComplexHeatmap::Heatmap(as.matrix(t(pseudobulk_qqnorm)), 
                     name = "Gene_exp", # Name of the heatmap and legend title
                     col = rev(c("#d53e4f", "#f46d43", "#fdae61", "#f7fcca", "#cffcd3", "#abdda4", "#849482")), # "#ddf1da"
                     use_raster = T,
                     column_title = NULL, row_title = NULL,
                     row_split = pseudobulk_metadata$celltype,
                     show_column_names = T,
                     show_row_names = T,
                     cluster_columns = T,
                     column_km = 5,
                     left_annotation = row_ha,
                     #right_annotation = row_ha_labels,
                     height = nrow(t(pseudobulk_qqnorm))*unit(3, "mm"),
                     width = ncol(t(pseudobulk_qqnorm))*unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 7),  # Text size for row names
                     column_names_gp = gpar(fontsize = 7),  # Text size for col names
                     row_names_side = "right")
  
  heatmap <- draw(hm_rna)
}

p <- plot_func()

# Saving to a file
filename <- "/home/hnatri/CART/13384_Tumors/GBM_JAK1KO_heatmap_pseudobulk_canonical_qqnorm_horizontal_2.pdf"
pdf(file = filename,
    width = 12,
    height = 10)
p
dev.off()

#png(file="/home/hnatri/CART/13384_Tumors/GBM_JAK1KO_heatmap_pseudobulk_canonical_qqnorm.png",
#    width=1000, height=1000)
#p
#dev.off()

# Extracting genes in each cluster
# loop to extract genes for each cluster.
for (i in 1:length(column_order(p))){
  if (i == 1) {
      clu <- t(t(row.names(pseudobulk_qqnorm[column_order(p)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("GeneID", "Cluster")
      } else {
          clu <- t(t(row.names(pseudobulk_qqnorm[column_order(p)[[i]],])))
          clu <- cbind(clu, paste("cluster", i, sep=""))
          out <- rbind(out, clu)
          }
}

head(out)
cluster_info <- as.data.frame(out)

#==============================================================================
# Scatterplots of GBM and JAK1KO pseudobulk data
#==============================================================================

# Scatterplot of the pseudobulk data
label_subset <- pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  filter(Hs> 2.5 & Mm > 2.5) %>%
  mutate(label = gene)

# Regression
r2_data <- pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene)

r2 <- round(summary(lm(Hs~Mm, r2_data))$adj.r.squared, digits = 2)

pseudobulk_scatter <- pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene) %>%
  ggplot(aes(y = log2(Hs+1), x = log2(Mm+1), color = celltype, label = label)) +
  geom_point() +
  scale_color_manual(values = celltype_annot_integrated_col) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 3), expand = c(0, 0)) +
  #ylab(expression(log[2]*"(GBM expression + 1)")) + 
  #xlab(expression(log[2]*"(JAK1KO expression + 1)")) + 
  ylab("HGG/GBM average expression") +
  xlab("JAK1/KO average expression") +
  annotate("text", x = 0.25, y = 2.75, label = r2, size = 2) +
  geom_smooth(method='lm', fill = "lightgray", linetype = "dashed",
              alpha = 0.1) +
  labs(color = "Cell type") +
  geom_text_repel(data=label_subset, aes(label=label),
                  size = 2,
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  #arrow = arrow(length = unit(0.010, "npc")),
                  nudge_x = .15,
                  nudge_y = .15) +
  NoLegend() +
  coord_fixed(ratio=1) +
  manuscript_theme

pseudobulk_scatter

# Genes with log2(exp) above threshold in both
select_genes <- pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene) %>%
  filter(Hs > 2.5 & Mm > 2.5) %>%
  dplyr::select(gene) %>%
  unlist() %>% unique()

#==============================================================================
# Comparing GBM and JAK1 WT
#==============================================================================
  
# Heatmap of WT+GBM
DefaultAssay(wt_gbm_seurat) <- "RNA"
VariableFeatures(wt_gbm_seurat) <- canonical_markers_all$markers
wt_gbm_seurat <- NormalizeData(wt_gbm_seurat)
wt_gbm_seurat <- ScaleData(wt_gbm_seurat)

wt_gbm_seurat$celltype_dataset <- paste0(wt_gbm_seurat$seurat_clusters, "-", wt_gbm_seurat$orig.ident)

lognorm_counts <- LayerData(wt_gbm_seurat,
                            assay = "RNA",
                            layer = "data",
                            features = intersect(rownames(wt_gbm_seurat), canonical_markers_all$markers))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$celltype_dataset <- wt_gbm_seurat$celltype_dataset

lognorm_counts <- lognorm_counts %>% group_by(celltype_dataset) %>%
  summarise_at(vars(-group_cols()), mean) %>%
  as.data.frame()

lognorm_counts <- t(lognorm_counts) %>%
  as.data.frame()

colnames(lognorm_counts) <- lognorm_counts[1,]
lognorm_counts <- lognorm_counts[2:nrow(lognorm_counts),]
wt_pseudobulk <- lognorm_counts
wt_pseudobulk <- sapply(wt_pseudobulk, as.numeric)
rownames(wt_pseudobulk) <- rownames(lognorm_counts)

wt_pseudobulk <- wt_pseudobulk[!is.infinite(rowSums(wt_pseudobulk)),]
colnames(wt_pseudobulk) <- as.character(colnames(wt_pseudobulk))

# Annotations
wt_pseudobulk_metadata <- data.frame("celltype" = sapply(strsplit(colnames(wt_pseudobulk),"-"), `[`, 1),
                                     "dataset" = sapply(strsplit(colnames(wt_pseudobulk),"-"), `[`, 2))

wt_pseudobulk_metadata$dataset <- gsub("Hs", "GBM", wt_pseudobulk_metadata$dataset)
wt_pseudobulk_metadata$dataset <- gsub("Mm", "WT", wt_pseudobulk_metadata$dataset)

# Colors
wt_pseudobulk_clusters <- as.factor(c(0, seq(1:length(unique(wt_pseudobulk_metadata$celltype)))))
wt_pseudobulk_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(wt_pseudobulk_clusters))
#names(wt_pseudobulk_cluster_col) <- paste0("g", unique(wt_pseudobulk_clusters))

col <- list()
col$celltype <- wt_pseudobulk_cluster_col
col$dataset <- dataset_col
row_ha <- rowAnnotation(
  df = data.frame(celltype = wt_pseudobulk_metadata$celltype,
                  dataset = wt_pseudobulk_metadata$dataset),
  annotation_height = unit(3, "mm"),
  annotation_width = unit(3, "mm"),
  col = col,
  annotation_name_gp = gpar(fontsize = 7)
  #col = c(col, cont_col)
)

# Quantile-normalizing rows and columns for plotting
wt_pseudobulk_qqnorm <- t(apply(wt_pseudobulk, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
wt_pseudobulk_qqnorm <- (apply(wt_pseudobulk_qqnorm, 2, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(wt_pseudobulk_qqnorm) <- colnames(wt_pseudobulk)
rownames(wt_pseudobulk_qqnorm) <- rownames(wt_pseudobulk)

colnames(wt_pseudobulk_qqnorm) <- gsub("Mm", "WT", colnames(wt_pseudobulk_qqnorm))
colnames(wt_pseudobulk_qqnorm) <- gsub("Hs", "GBM", colnames(wt_pseudobulk_qqnorm))

wt_pseudobulk_metadata$celltype_dataset <- paste0(wt_pseudobulk_metadata$celltype, "_",
                                                  wt_pseudobulk_metadata$dataset)

plot_func <- function(){
  hm_rna  <- Heatmap(as.matrix(t(wt_pseudobulk_qqnorm)), 
                     name = "Gene_exp",
                     col = rev(c("#d53e4f", "#f46d43", "#fdae61", "#f7fcca", "#cffcd3", "#abdda4", "#849482")), # "#ddf1da"
                     use_raster = T,
                     column_title = NULL, row_title = NULL,
                     row_split = wt_pseudobulk_metadata$celltype,
                     show_column_names = T,
                     show_row_names = T,
                     cluster_columns = T,
                     column_km = 10,
                     left_annotation = row_ha,
                     height = nrow(t(wt_pseudobulk_qqnorm))*unit(3, "mm"),
                     width = ncol(t(wt_pseudobulk_qqnorm))*unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 7),
                     column_names_gp = gpar(fontsize = 7),
                     row_names_side = "right")
  
  heatmap <- draw(hm_rna)
}

p <- plot_func()


#==============================================================================
# Scatterplots of GBM and JAK1 WT pseudobulk data
#==============================================================================

# Scatterplot of the pseudobulk data
wt_label_subset <- wt_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(wt_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  filter(Hs > 2.5 & Mm > 2.5) %>%
  mutate(label = gene)

# Regression
r2_data <- wt_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(wt_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene)

r2 <- round(summary(lm(Hs~Mm, r2_data))$adj.r.squared, digits = 2)

wt_pseudobulk_scatter <- wt_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(wt_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene) %>%
  ggplot(aes(y = log2(Hs+1), x = log2(Mm+1), color = celltype, label = label)) +
    geom_point() +
    scale_color_manual(values = wt_pseudobulk_cluster_col) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 3), expand = c(0, 0)) +
  ylab("HGG/GBM average expression") +
  xlab("WT average expression") +
    annotate("text", x = 0.25, y = 2.75, label = r2, size = 2) +
    geom_smooth(method='lm', fill = "lightgray", linetype = "dashed",
                alpha = 0.1) +
    labs(color = "Cell type") +
    geom_text_repel(data=wt_label_subset, aes(label=label),
                    size = 2,
                    min.segment.length = 0, 
                    seed = 42, 
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    #arrow = arrow(length = unit(0.010, "npc")),
                    nudge_x = .15,
                    nudge_y = .15) +
    NoLegend() +
    coord_fixed(ratio=1) +
    manuscript_theme

wt_pseudobulk_scatter

pseudobulk_scatter + wt_pseudobulk_scatter

# Saving to a file
filename <- "/home/hnatri/CART/13384_Tumors/GBM_JAK1KO_WT_pseudobulk_scatter.pdf"
pdf(file = filename,
    width = 5,
    height = 2.75)

pseudobulk_scatter + wt_pseudobulk_scatter

dev.off()

#==============================================================================
# Comparing GBM PD and JAK1KO
#==============================================================================

# Heatmap of GBM PD only and JAK1KO pre-treatment
DefaultAssay(ko_gbm_pd_seurat) <- "RNA"
VariableFeatures(ko_gbm_pd_seurat) <- canonical_markers_all$markers
ko_gbm_pd_seurat <- NormalizeData(ko_gbm_pd_seurat)
ko_gbm_pd_seurat <- ScaleData(ko_gbm_pd_seurat)

ko_gbm_pd_seurat$celltype_dataset <- paste0(ko_gbm_pd_seurat$seurat_clusters, "-", ko_gbm_pd_seurat$orig.ident)

lognorm_counts <- LayerData(ko_gbm_pd_seurat,
                            assay = "RNA",
                            layer = "data",
                            features = intersect(rownames(ko_gbm_pd_seurat), canonical_markers_all$markers))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$celltype_dataset <- ko_gbm_pd_seurat$celltype_dataset

lognorm_counts <- lognorm_counts %>% group_by(celltype_dataset) %>%
  summarise_at(vars(-group_cols()), mean) %>%
  as.data.frame()

lognorm_counts <- t(lognorm_counts) %>%
  as.data.frame()

colnames(lognorm_counts) <- lognorm_counts[1,]
lognorm_counts <- lognorm_counts[2:nrow(lognorm_counts),]
ko_gbm_pd_pseudobulk <- lognorm_counts
ko_gbm_pd_pseudobulk <- sapply(ko_gbm_pd_pseudobulk, as.numeric)
rownames(ko_gbm_pd_pseudobulk) <- rownames(lognorm_counts)

# Annotations
ko_gbm_pd_pseudobulk_metadata <- data.frame("celltype" = sapply(strsplit(colnames(ko_gbm_pd_pseudobulk),"-"), `[`, 1),
                                     "dataset" = sapply(strsplit(colnames(ko_gbm_pd_pseudobulk),"-"), `[`, 2))

ko_gbm_pd_pseudobulk_metadata$dataset <- gsub("Hs", "GBM", ko_gbm_pd_pseudobulk_metadata$dataset)
ko_gbm_pd_pseudobulk_metadata$dataset <- gsub("Mm", "WT", ko_gbm_pd_pseudobulk_metadata$dataset)

# Colors
ko_gbm_pd_pseudobulk_clusters <- as.factor(c(0, seq(1:length(unique(ko_gbm_pd_pseudobulk_metadata$celltype)))))
ko_gbm_pd_pseudobulk_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(ko_gbm_pd_pseudobulk_clusters))
#names(ko_gbm_pd_pseudobulk_cluster_col) <- paste0("g", unique(ko_gbm_pd_pseudobulk_clusters))

col <- list()
col$celltype <- ko_gbm_pd_pseudobulk_cluster_col
col$dataset <- dataset_col
row_ha <- rowAnnotation(
  df = data.frame(celltype = ko_gbm_pd_pseudobulk_metadata$celltype,
                  dataset = ko_gbm_pd_pseudobulk_metadata$dataset),
  annotation_height = unit(3, "mm"),
  annotation_width = unit(3, "mm"),
  col = col,
  annotation_name_gp = gpar(fontsize = 7)
  #col = c(col, cont_col)
)

# Quantile-normalizing rows and columns for plotting
ko_gbm_pd_pseudobulk_qqnorm <- t(apply(ko_gbm_pd_pseudobulk, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
ko_gbm_pd_pseudobulk_qqnorm <- (apply(ko_gbm_pd_pseudobulk_qqnorm, 2, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(ko_gbm_pd_pseudobulk_qqnorm) <- colnames(ko_gbm_pd_pseudobulk)
rownames(ko_gbm_pd_pseudobulk_qqnorm) <- rownames(ko_gbm_pd_pseudobulk)

colnames(ko_gbm_pd_pseudobulk_qqnorm) <- gsub("Mm", "WT", colnames(ko_gbm_pd_pseudobulk_qqnorm))
colnames(ko_gbm_pd_pseudobulk_qqnorm) <- gsub("Hs", "GBM", colnames(ko_gbm_pd_pseudobulk_qqnorm))

ko_gbm_pd_pseudobulk_metadata$celltype_dataset <- paste0(ko_gbm_pd_pseudobulk_metadata$celltype, "_",
                                                         ko_gbm_pd_pseudobulk_metadata$dataset)

plot_func <- function(){
  hm_rna  <- Heatmap(as.matrix(t(ko_gbm_pd_pseudobulk_qqnorm)), 
                     name = "Gene_exp",
                     col = rev(c("#d53e4f", "#f46d43", "#fdae61", "#f7fcca", "#cffcd3", "#abdda4", "#849482")), # "#ddf1da"
                     use_raster = T,
                     column_title = NULL, row_title = NULL,
                     row_split = ko_gbm_pd_pseudobulk_metadata$celltype,
                     show_column_names = T,
                     show_row_names = T,
                     cluster_columns = T,
                     column_km = 10,
                     left_annotation = row_ha,
                     height = nrow(t(ko_gbm_pd_pseudobulk_qqnorm))*unit(3, "mm"),
                     width = ncol(t(ko_gbm_pd_pseudobulk_qqnorm))*unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 7),
                     column_names_gp = gpar(fontsize = 7),
                     row_names_side = "right")
  
  heatmap <- draw(hm_rna)
}

p <- plot_func()

#==============================================================================
# Scatterplots of GBM PD and JAK1KO pseudobulk data
#==============================================================================

# Scatterplot of the pseudobulk data
ko_gbm_pd_label_subset <- ko_gbm_pd_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(ko_gbm_pd_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  filter(Hs > 2.5 & Mm > 2.5) %>%
  mutate(label = gene)

# Regression
r2_data <- ko_gbm_pd_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(ko_gbm_pd_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene)

r2 <- round(summary(lm(Hs~Mm, r2_data))$adj.r.squared, digits = 2)
#   annotate("text", x = 2, y = 9, label = r2, size = 2) +

ko_gbm_pd_pseudobulk_scatter <- ko_gbm_pd_pseudobulk %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = setdiff(colnames(ko_gbm_pd_pseudobulk), "gene"), names_to = "celltype_dataset", values_to ="exp") %>%
  separate(celltype_dataset, into = c("celltype", "dataset"), sep = "-") %>% 
  pivot_wider(id_cols = c("celltype", "gene"), names_from = c("dataset"), values_from = "exp") %>%
  mutate(label = gene) %>%
  ggplot(aes(y = log2(Hs+1), x = log2(Mm+1), color = celltype, label = label)) +
  geom_point() +
  scale_color_manual(values = ko_gbm_pd_pseudobulk_cluster_col) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 3), expand = c(0, 0)) +
  ylab("HGG/GBM PD average expression") +
  xlab("JAK1/KO average expression") +
  annotate("text", x = 0.25, y = 2.75, label = r2, size = 2) +
  geom_smooth(method='lm', fill = "lightgray", linetype = "dashed",
              alpha = 0.1) +
  labs(color = "Cell type") +
  geom_text_repel(data=ko_gbm_pd_label_subset, aes(label=label),
                  size = 2,
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  #arrow = arrow(length = unit(0.010, "npc")),
                  nudge_x = .15,
                  nudge_y = .15) +
  NoLegend() +
  coord_fixed(ratio=1) +
  manuscript_theme

pseudobulk_scatter + wt_pseudobulk_scatter + ko_gbm_pd_pseudobulk_scatter

# Saving to a file
filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/3plots_pseudobulk_scatter.pdf"
pdf(file = filename,
    width = 8.25,
    height = 2.75)

pseudobulk_scatter + wt_pseudobulk_scatter + ko_gbm_pd_pseudobulk_scatter

dev.off()

