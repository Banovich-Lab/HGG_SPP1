#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 10/03/2023
# Description: 13384 tumor data visualization colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Colors and themes
#==============================================================================

# ggplot theme
manuscript_theme <- theme(text = element_text(size = 6),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),  
                          axis.title.x = element_text(size = 6),
                          axis.title.y = element_text(size = 6))

# Colors for plotting
# Define colors for each level of categorical variables

# Clusters
tumor_clusters <- as.factor(c(0, seq(1:20)))
tumor_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(tumor_clusters))
names(tumor_cluster_col) <- levels(tumor_clusters)

# Integrated GBM + JAK1KO mouse data

# Cluster colors
integrated_tumor_clusters <- as.factor(c(0, seq(1:23)))
integrated_tumor_clusters_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(integrated_tumor_clusters))
names(integrated_tumor_clusters_col) <- levels(integrated_tumor_clusters)

# Cell types
# https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing
gs4_deauth()
tumor_tables  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
sheet_names(tumor_tables)
celltype_annot <- read_sheet(tumor_tables, sheet = "Cluster annotations")

tumor_celltype_col <- celltype_annot$color_fig1
names(tumor_celltype_col) <- celltype_annot$annotation
tumor_celltype2_col <- celltype_annot$color_fig1
names(tumor_celltype2_col) <- celltype_annot$annotation2

celltype_annot_immune_fibro <- read_sheet(tumor_tables, sheet = "Cluster annotations, immune+fibroblast")
celltype_annot_immune_fibro$cluster <- as.character(celltype_annot_immune_fibro$cluster)
celltype_annot_immune_fibro$orig_cluster <- as.character(celltype_annot_immune_fibro$orig_cluster)
immune_fibro_celltype_col <- celltype_annot_immune_fibro$color_fig1
names(immune_fibro_celltype_col) <- celltype_annot_immune_fibro$annotation

celltype_annot_integrated <- read_sheet(tumor_tables, sheet = "Cluster annotations, GBM+mouse")
celltype_annot_integrated$cluster <- as.character(celltype_annot_integrated$cluster)
celltype_annot_integrated_col <- celltype_annot_integrated$color_fig1
names(celltype_annot_integrated_col) <- celltype_annot_integrated$annotation

jak1_celltype_annot <- read_sheet(tumor_tables, sheet = "Cluster annotations, JAK mouse")
jak1_celltype_annot$annotation <- jak1_celltype_annot$cluster
jak1_celltype_col <- jak1_celltype_annot$color_fig1
names(jak1_celltype_col) <- jak1_celltype_annot$annotation

# CD3_high_low
CD3_high_low_col <- c("High" = brewer.pal(3, "Set1")[1],
                      "Low" = brewer.pal(3, "Set1")[2],
                      "NA" = "azure3")

# Dataset
dataset_col <- c("GBM" = "ghostwhite",
                 "JAK1KO" = "ivory4",
                 "WT" = "ivory4")

## CD4/CD8
#cd4_cd8_col <- c("CD8" = brewer.pal(3, "Set1")[1],
#                 "CD4" = brewer.pal(3, "Set1")[2])
#
## Response
#response_col <- c("NA" = nord("silver_mine", 5)[3],
#                  "Progression Disease (PD)" = nord("aurora", 5)[1],
#                  "Partial Response (PR)" = nord("aurora", 5)[3],
#                  "Stable Disease (SD)" = nord("aurora", 5)[5],
#                  "Complete Response (CR)" = nord("aurora", 5)[4])

# Batch
batch_col <- colorRampPalette(brewer.pal(10, "Spectral"))(19)
names(batch_col) <- paste0("Batch", seq(1,19))
colScale_batch <- scale_fill_manual(name = "Batch", values = batch_col)

simple_batch_col <- colorRampPalette(brewer.pal(10, "Spectral"))(19)
names(simple_batch_col) <- seq(1,19)

# Diagnosis
diagnosis_col <- c(colorRampPalette(brewer.pal(10, "Spectral"))(10), "#666666")
diagnosis <- c("Glioblastoma, NOS", "Ependymoma, NOS", "Astrocytoma, anaplastic", 
               "Oligodendroglioma, anaplastic", "Primitive neuroectodermal tumor",                            
               "Astrocytoma, NOS", "Fibrillary astrocytoma", "Pleomorphic xanthoastrocytoma",
               "Oligodendroglioma, NOS", "Ependymoma, anaplastic", "NA")

names(diagnosis_col) <- diagnosis
colScale_diagnosis <- scale_fill_manual(name = "Diagnosis", values = diagnosis_col)

# Tumor UPN color
tumor_UPN_col <- colorRampPalette(brewer.pal(11, "Spectral"))(44)
names(tumor_UPN_col) <- sort(c("109", "243", "185", "277", "241", "315", "303", "282", "275", "289", "265", "125", "181", "234", "117", "232", "191", "129", "146", "157", "208", "145", "214", "237", "240", "218", "350", "228", "248", "266", "131", "239", "301", "215", "223", "224", "230", "260", "213", "193", "149", "141", "124", "122"))

# Cell cycle
#cell_cycle_col <- c("G1" = nord("lumina", 5)[1],
#                    "G2M" = nord("lumina", 5)[3],
#                    "S" = nord("lumina", 5)[5])

# TCR clone size
tcr_clonetype_col <- c("Single (0 < X <= 1)" = "lightyellow",
                       "Small (1 < X <= 5)" = "yellow2",
                       "Medium (5 < X <= 20)" = "orange",
                       "Large (20 < X <= 100)" = "red",
                       "Hyperexpanded (100 < X <= 500)" = "darkred")

