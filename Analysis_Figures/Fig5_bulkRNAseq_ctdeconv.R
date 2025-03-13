#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 08/29/2024
# Description: Celltype deconvolution of SPP1+CAR bulk-RNAseq
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(tidyverse)
library(DESeq2)
library(googlesheets4)
library(pheatmap)
library(edgeR)
library(BayesPrism)
library(org.Hs.eg.db)

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

# scRNAseq for reference
#tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
#immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")
#
## Adding immune annotations
#tumors$celltype_combined <- mapvalues(x = rownames(tumors@meta.data),
#                                      from = rownames(immune_fibro@meta.data),
#                                      to = immune_fibro$celltype)
#tumors$celltype_combined <- ifelse(rownames(tumors@meta.data) %in% rownames(immune_fibro@meta.data), tumors$celltype_combined, as.character(tumors$celltype))

combined <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/scRNA_seq_Pool.combined.rds")
seurat_object <- scRNA_seq_Pool.combined

rm(scRNA_seq_Pool.combined)

# Updating cell type annotations
gs4_deauth()
cluster_annotations  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
cluster_annotations <- read_sheet(cluster_annotations, sheet = "Cluster annotations, JAK mouse")

seurat_object$celltype <- mapvalues(seurat_object$seurat_clusters,
                                    from = cluster_annotations$cluster,
                                    to = cluster_annotations$annotation)

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

gene_names <- convert_mouse_to_human(rownames(seurat_object@assays$RNA))

# Keeping mouse genes with a single human ortholog
gene_names <- gene_names %>%
  group_by(mouse) %>%
  filter(!is.na(human),
         n() == 1) %>%
  ungroup()

DefaultAssay(seurat_object) <- "RNA"
assay_data <- LayerData(seurat_object, layer = "counts")
assay_data <- assay_data[which(rownames(assay_data) %in% gene_names$mouse),]
new_names <- rownames(assay_data)
new_names <- mapvalues(x = new_names,
                       from = gene_names$mouse,
                       to = gene_names$human)
rownames(assay_data) <- new_names

seurat_object[["RNA_human"]] <- CreateAssayObject(assay_data,
                                                  min.cells = 0,
                                                  min.features = 0)
DefaultAssay(seurat_object) <- "RNA_human"

# Bulk-RNAseq
count_data <- read.csv("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Results_IGC-CM-21417_10202023/DESeq2/gene_count_matrix.csv")

ensemble <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 1)
gene_names <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 2)

count_data$gene_name <- gene_names

# Removing duplicate gene names
count_data <- count_data %>%
  group_by(gene_name) %>%
  filter(n() == 1) %>%
  ungroup()

gene_names <- count_data$gene_name
count_data <- as.data.frame(count_data)
rownames(count_data) <- gene_names

head(count_data)

gene_names <- convert_mouse_to_human(gene_list = rownames(count_data))
head(gene_names)

# Keeping mouse genes with a single human ortholog
gene_names <- gene_names %>%
  group_by(mouse) %>%
  filter(!is.na(human),
         n() == 1) %>%
  ungroup()

count_data <- count_data[which(rownames(count_data) %in% gene_names$mouse),]
new_names <- mapvalues(x = rownames(count_data),
                       from = gene_names$mouse,
                       to = gene_names$human)
count_data$human_gene <- new_names

count_data <- count_data %>% group_by(human_gene) %>% 
  filter(n() == 1) %>% 
  ungroup()

count_data <- as.data.frame(count_data)
rownames(count_data) <- count_data$human_gene
count_data <- count_data %>% dplyr::select(-c("gene_name", "human_gene", "gene_id"))

#==============================================================================#
# Running BayesPrism
#==============================================================================#

load("/home/hnatri/Tools/BayesPrism/tutorial.dat/tutorial.gbm.rdata")

# Bulk
bk.dat <- t(count_data)

# scRNAseq
#head(sc.dat)
table(cell.type.labels)
table(seurat_object$celltype)

# Combining tumor labels
seurat_object$celltype <- ifelse(seurat_object$celltype %in% paste0("Epi", seq(1,6)), "Tumor", as.character(seurat_object$celltype))

sc.dat <- LayerData(seurat_object, assay = "RNA_human", layer = "counts")
sc.dat <- t(sc.dat)
sc.dat <- as.matrix(sc.dat)

cell.type.labels <- as.character(seurat_object$celltype)

# sc
plot.cor.phi(input=sc.dat,
             input.labels=cell.type.labels, #cell.state.labels,
             title="cell state correlation",
             cexRow=0.2, cexCol=0.2,
             margins=c(2,2))

# Outliers in sc
sc.stat <- plot.scRNA.outlier(input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
                              cell.type.labels=cell.type.labels,
                              species="hs", #currently only human(hs) and mouse(mm) annotations are supported
                              return.raw=TRUE) #return the data used for plotting.

# Outliers in bulk
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE)

# Filtering outlier genes from sc
sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                 input.type="count.matrix",
                                 species="hs", 
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","chrX","chrY"),
                                 exp.cells=5)

length(intersect(colnames(sc.dat.filtered), colnames(bk.dat)))

saveRDS(sc.dat.filtered, "/scratch/hnatri/CART/sc.dat.filtered.rds")
saveRDS(bk.dat, "/scratch/hnatri/CART/bk.dat.rds")

sc.dat.filtered <- readRDS("/scratch/hnatri/CART/sc.dat.filtered.rds")
bk.dat <- readRDS("/scratch/hnatri/CART/bk.dat.rds")

length(intersect(colnames(sc.dat.filtered), colnames(bk.dat)))

# Gene expression concordance between bulk and sc
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                       bulk.input = bk.dat,
                       pdf.prefix = NULL)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/bulk_sc_JAK1KOWT_corr.pdf"
pdf(file = filename,
    width = 12,
    height = 4)

plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bk.dat,
                pdf.prefix = NULL)


dev.off()

# Subset protein coding genes
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                       gene.type = "protein_coding")

# Constructing the BayesPrism object
myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=c("Tumor"),
  outlier.cut=0.01,
  outlier.fraction=0.1)

# Running BayesPrism
bp.res <- run.prism(prism = myPrism, n.cores=4)

theta <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")

#write.table(theta, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bayesprism_res_KOWTscref_Hhgenename.tsv", quote = F, row.names = T, sep = "\t")
#theta <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bayesprism_res_KOWTscref_Hhgenename.tsv", header = T)

# Plotting
plot_cols <- cluster_annotations$color_fig1
names(plot_cols) <- cluster_annotations$annotation
plot_cols <- plot_cols[!(names(plot_cols) %in% paste0("Epi", c(1, 2, 4, 5, 6)))]
names(plot_cols) <- gsub("Epi3", "Tumor", names(plot_cols))

#ct_cols <- as.factor(colnames(theta))
#plot_cols <- colorRampPalette(brewer.pal(11, "Paired"))(length(ct_cols))
#names(plot_cols) <- levels(ct_cols)
#plot_cols["Tumor"] <- "#cc2a56"

pbp <- theta %>% rownames_to_column(var = "group") %>%
  pivot_longer(cols = setdiff(colnames(theta), "group"),
               values_to = "prop",
               names_to = "celltype") %>%
  mutate(group = factor(group, levels = c("SPP1_CAR_T", "CAR_T", "SPP1", "Control")),
         pct = prop*100,) %>%
  ggplot(aes(x = pct, y = group, fill = celltype)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = plot_cols) +
  ylab("") +
  xlab("Proportion (%)") #manuscript_theme

pbp

theta %>% rownames_to_column(var = "group") %>%
  pivot_longer(cols = setdiff(colnames(theta), "group"),
               values_to = "prop",
               names_to = "celltype") %>%
  filter(group == "SPP1_CAR_T") %>%
  print(n=30)

# c("#FFFF99")"

theta %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "celltype") %>%
  ggplot(aes(x = 100*SPP1_CAR_T, y = 100*Control, color = celltype, label = celltype)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = plot_cols) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = "dashed") +
  #theme(panel.background = element_rect(colour = "gray88")) +
  theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm")) +
  geom_text_repel(size = 2.5) +
  #ylim(0, 40) + 
  #xlim(0, 40) +
  theme(legend.position = "none") + 
  ylab("Control") +
  xlab("SPP1 + CAR T")
  #manuscript_theme

theta %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "celltype") %>%
  filter(celltype == "M6")
