#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/15/2024
# Description: Gene Set Enrichment Analysis for JAK1KO/WT tumor CD45+ data
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
library(escape)
library(dittoSeq)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/CART/CART_plot_functions.R")
source("/home/hnatri/CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Import data
#==============================================================================#

# The final object
seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/SPP1_release/JAK1KO_WildType_Control_Cells.rds")

head(seurat_object@meta.data)
unique(seurat_object$celltype)

table(seurat_object$celltype, seurat_object$JAK1_Group)

celltypes_keep <- c("M1", "M2", "M3", "Neut1", "Neut2")

seurat_object <- subset(seurat_object, subset = celltype %in% celltypes_keep)

#seurat_object <- immune_fibro

#==============================================================================#
# Run GSEA
#==============================================================================#

# C2 = curated gene sets
GS <- getGeneSets(species = "Homo sapiens", library ="C2")
GS_CANONICAL <- GS[grep("KEGG|REACTOME|BIOCARTA", names(GS),
                        ignore.case = TRUE)]

res_seurat_object <- enrichIt(obj = seurat_object,
                       gene.sets = GS_CANONICAL,
                       groups = 50, cores = 4)
seurat_object <- AddMetaData(seurat_object, res_seurat_object)

saveRDS(seurat_object, "/scratch/hnatri/CART/mouse_myeloid_GSEA.rds")
saveRDS(res_seurat_object, "/scratch/hnatri/CART/mouse_myeloid_GSEA_res.rds")

q(save="no")

seurat_object <- readRDS("/scratch/hnatri/CART/mouse_myeloid_GSEA.rds")
res_seurat_object <- readRDS("/scratch/hnatri/CART/mouse_myeloid_GSEA_res.rds")

# Only CD45+
seurat_object <- subset(seurat_object, subset = orig.ident %in% c("KO_CD45_Pos", "WT_CD45_Pos"))

# Plotting
seurat_object_res2 <- seurat_object@meta.data %>%
  dplyr::select(c("orig.ident", colnames(res_seurat_object)))
#colnames(myeloid_res2)[ncol(myeloid_res2)] <- "celltype"

#myeloid_res2$binary_response <- factor(myeloid_res2$binary_response, levels = c("CR_SD", "PD"))
output <- data.frame(getSignificance(seurat_object_res2, group = "orig.ident", fit = "Wilcoxon"))
output$pathways <- rownames(output)
output <- output %>% filter(FDR < 0.01) %>% arrange(FDR)

write.table(output, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/JAK1KO_WT_GSEA_sig.tsv",
            quote = F, sep = "\t", row.names = F)

output <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/JAK1KO_WT_GSEA_sig.tsv",
                     header = T)

DefaultAssay(seurat_object) <- "RNA"

dittoHeatmap(seurat_object,
             genes = NULL,
             metas = output_seurat_object$pathways, 
             heatmap.colors = colorRampPalette(c("dodgerblue3", "white", "tomato"))(100),
             annot.colors = c("deeppink3", "aquamarine3"),
             annot.by = "orig.ident",
             order.by = "orig.ident",
             cluster_cols = F,
             fontsize = 7)

png(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_JAK1KOWT_sig.png",
    width=1200, height=300)

dittoHeatmap(seurat_object,
             genes = NULL,
             metas = output_seurat_object$pathways[1:15], 
             heatmap.colors = colorRampPalette(c("dodgerblue3", "white", "tomato"))(100),
             annot.colors = c("deeppink3", "aquamarine3"),
             annot.by = "orig.ident",
             order.by = "orig.ident",
             cluster_cols = F,
             fontsize = 10)

dev.off()

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_JAK1KOWT_sig.pdf",
    width=9, height=3)

dittoHeatmap(seurat_object,
             genes = NULL,
             metas = output_seurat_object$pathways, 
             heatmap.colors = colorRampPalette(c("dodgerblue3", "white", "tomato"))(100),
             annot.colors = c("deeppink3", "aquamarine3"),
             annot.by = "orig.ident",
             order.by = "orig.ident",
             cluster_cols = F,
             fontsize = 5)

dev.off()
