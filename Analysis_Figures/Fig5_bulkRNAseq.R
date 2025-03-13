#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/15/2024
# Description: Plotting the SPP1+CAR bulk-RNAseq
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(tidyverse)
library(DESeq2)
library(googlesheets4)
library(pheatmap)
library(edgeR)

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

# Mouse immune markers
gs4_deauth()
canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
canonical_markers_all <- read_sheet(canonical_markers, sheet = "Mm immune markers")
deconv_markers <- read_sheet(canonical_markers, sheet = "Deconvolution markers")

#count_data_yz <- load("/scratch/hnatri/DESeq2_Expression_ISG15_SPP1.rds")
count_data <- read.csv("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Results_IGC-CM-21417_10202023/DESeq2/gene_count_matrix.csv")

ensemble <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 1)
gene_names <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 2)

rownames(count_data) <- count_data$gene_name
count_data <- count_data %>% dplyr::select(c("CAR_T", "Control", "SPP1", "SPP1_CAR_T"))

cpm_data <- cpm(count_data, log = T)
col_data <- data.frame("group" = colnames(count_data))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~1)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts$gene_name <- gene_names

norm_counts <- cpm(count_data, log = T)
norm_counts <- as.data.frame(norm_counts)
norm_counts$gene_name <- gene_names

plot_count_data <- norm_counts %>% filter(gene_name %in% canonical_markers_all$gene_name)

rownames(plot_count_data) <- plot_count_data$gene_name
plot_count_data <- plot_count_data %>% dplyr::select(c("CAR_T", "Control", "SPP1", "SPP1_CAR_T"))

canonical_markers_all <- canonical_markers_all %>% filter(gene_name %in% rownames(plot_count_data))
annot <- as.data.frame(canonical_markers_all$annotation)
rownames(annot) <- canonical_markers_all$gene_name

plot_count_data <- plot_count_data[match(rownames(annot), rownames(plot_count_data)),]
plot_count_data <- plot_count_data[complete.cases(plot_count_data),]

pheatmap(plot_count_data,
         annotation_row = annot,
         fontsize = 7,
         cluster_rows = F,
         scale = "column")

# Deconvolution markers
ComICS::DCQ_mar

plot_count_data <- norm_counts %>% filter(gene_name %in% ComICS::DCQ_mar$DCQ)

rownames(plot_count_data) <- plot_count_data$gene_name
plot_count_data <- plot_count_data %>% dplyr::select(c("CAR_T", "Control", "SPP1", "SPP1_CAR_T"))

#canonical_markers_all <- canonical_markers_all %>% filter(gene_name %in% rownames(plot_count_data))
#annot <- as.data.frame(canonical_markers_all$annotation)
#rownames(annot) <- canonical_markers_all$gene_name

#plot_count_data <- plot_count_data[match(rownames(annot), rownames(plot_count_data)),]
plot_count_data <- plot_count_data[complete.cases(plot_count_data),]

pheatmap(t(plot_count_data),
         #annotation_row = annot,
         fontsize = 7,
         cluster_rows = T,
         cluster_cols = T)

# Deconvolution markers
plot_count_data <- norm_counts %>% filter(gene_name %in% deconv_markers$mm_gene_name)

deconv_markers$gene_name
norm_counts$gene_name[grep("Dcn", norm_counts$gene_name)]

rownames(plot_count_data) <- plot_count_data$gene_name
plot_count_data <- plot_count_data %>% dplyr::select(c("CAR_T", "Control", "SPP1", "SPP1_CAR_T"))

pheatmap(t(plot_count_data),
         #annotation_row = annot,
         fontsize = 7,
         cluster_rows = T,
         cluster_cols = T)

# Metabolism genes
meta_gene_table  <- gs4_get("https://docs.google.com/spreadsheets/d/16sI_lATPquB5yHQYXtXVs9UORQ6TNL1cqh2Yh3jy_Ho/edit?usp=sharing")
meta_gene_table <- read_sheet(meta_gene_table, sheet = "Sheet1")

# Converting to human gene names
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

norm_counts <- norm_counts %>% group_by(gene_name) %>% 
  filter(n() == 1) %>% 
  ungroup()
gene_names <- convert_mouse_to_human(norm_counts$gene_name)

# Keeping mouse genes with a single human ortholog
gene_names <- gene_names %>%
  group_by(mouse) %>%
  filter(!is.na(human),
         n() == 1) %>%
  ungroup()

norm_counts <- norm_counts[which(norm_counts$gene_name %in% gene_names$mouse),]
new_names <- mapvalues(x = norm_counts$gene_name,
                       from = gene_names$mouse,
                       to = gene_names$human)
norm_counts$human_gene <- new_names

norm_counts <- norm_counts %>% group_by(human_gene) %>% 
  filter(n() == 1) %>% 
  ungroup()

norm_counts <- as.data.frame(norm_counts)
rownames(norm_counts) <- norm_counts$human_gene
plot_count_data <- norm_counts %>% dplyr::select(-c("gene_name", "human_gene"))

hlas <- rownames(plot_count_data)[grep("^HLA", rownames(plot_count_data))]

plot_count_data_select <- plot_count_data[rownames(plot_count_data) %in% sort(c(meta_gene_table$Gene, hlas)),]

plot_count_data_select <- plot_count_data[rownames(plot_count_data) %in% hlas,]

metab_hm <- pheatmap(plot_count_data_select,
                     #annotation_row = annot,
                     fontsize = 8,
                     cluster_rows = T,
                     cluster_cols = T)

metab_hm


filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/JAK1KO_WT_CAR_bulkRNAseq_metabolism_heatmap.pdf"
pdf(file = filename,
    width = 4,
    height = 9)

metab_hm

dev.off()

