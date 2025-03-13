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

source("/home/hnatri/13384_CART/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Import data
#==============================================================================#

# scRNAseq for reference
combined <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/scRNA_seq_Seurat/scRNA_seq_Pool.combined.rds")
seurat_object <- scRNA_seq_Pool.combined

rm(scRNA_seq_Pool.combined)

head(seurat_object@meta.data)

# Updating cell type annotations
gs4_deauth()
cluster_annotations  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
sheet_names(cluster_annotations)
cluster_annotations <- read_sheet(cluster_annotations, sheet = "Cluster annotations, JAK mouse")

seurat_object$celltype <- mapvalues(seurat_object$seurat_clusters,
                                    from = cluster_annotations$cluster,
                                    to = cluster_annotations$annotation)

table(seurat_object$celltype)

table(seurat_object$CART_Group, seurat_object$JAK1_Group)

#seurat_object <- subset(seurat_object, subset = CD45_Group == "CD45_Pos")
#seurat_object <- subset(seurat_object, subset = CART_Group == "Control")
seurat_object <- subset(seurat_object, subset = JAK1_Group == "Jak1_KO")

# Bulk -NAseq
count_data <- read.csv("/tgen_labs/banovich/BCTCSF/Mouse_tumor/Results_IGC-CM-21417_10202023/DESeq2/gene_count_matrix.csv")

ensemble <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 1)
gene_names <- sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 2)

count_data$gene_name <- ensemble

# Removing duplicate gene names
count_data <- count_data %>%
  group_by(gene_name) %>%
  filter(n() == 1) %>%
  ungroup()

gene_names <- count_data$gene_name

count_data <- as.data.frame(count_data)
rownames(count_data) <- gene_names

#plot_count_data <- norm_counts %>% dplyr::select(c("CAR_T", "Control", "SPP1", "SPP1_CAR_T"))

head(count_data)

count_data <- as.data.frame(count_data)
rownames(count_data) <- count_data$gene_name
count_data <- count_data %>% dplyr::select(-c("gene_name", "gene_id"))

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
seurat_object$celltype <- ifelse(seurat_object$celltype %in% paste0("T", seq(1,6)), "Tumor", as.character(seurat_object$celltype))

sc.dat <- LayerData(seurat_object, assay = "RNA", layer = "counts")
sc.dat <- t(sc.dat)
sc.dat <- as.matrix(sc.dat)

cell.type.labels <- as.character(seurat_object$celltype)

sc.dat[1:10, 1:10]

#colnames(sc.dat) <-  sapply(strsplit(colnames(sc.dat), "\\."), `[`, 1)
library(org.Mm.eg.db)
sc_gene_names <- mapIds(org.Mm.eg.db, keys = colnames(sc.dat), keytype = "SYMBOL", column = "ENSEMBL")

colnames(sc.dat) <- sc_gene_names
sc.dat <- sc.dat[,!is.na(colnames(sc.dat))]
sc.dat <- sc.dat[, !duplicated(colnames(sc.dat))]

#sapply(strsplit(count_data$gene_id, split='|', fixed=TRUE), `[`, 1)
colnames(bk.dat) <- sapply(strsplit(colnames(bk.dat), split='.', fixed=TRUE), `[`, 1)
#bk.dat <- bk.dat[,!is.na(colnames(bk.dat))]
bk.dat <- bk.dat[, !duplicated(colnames(bk.dat))]


# sc
plot.cor.phi(input=sc.dat,
             input.labels=cell.type.labels, #cell.state.labels,
             title="cell state correlation",
             cexRow=0.2, cexCol=0.2,
             margins=c(2,2))

# Outliers in sc
sc.stat <- plot.scRNA.outlier(input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
                              cell.type.labels=cell.type.labels,
                              species="mm", #currently only human(hs) and mouse(mm) annotations are supported
                              return.raw=TRUE) #return the data used for plotting.

length(intersect(colnames(bk.dat), colnames(sc.dat)))

# Outliers in bulk
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE)

# Filtering outlier genes from sc
sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                 input.type="count.matrix",
                                 species="mm", 
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","chrX","chrY"),
                                 exp.cells=5)

length(intersect(colnames(sc.dat.filtered), colnames(bk.dat)))

saveRDS(sc.dat.filtered, "/scratch/hnatri/CART/sc.dat.filtered_mmEnsembl.rds")
saveRDS(bk.dat, "/scratch/hnatri/CART/bk.dat_mmEnsembl.rds")

sc.dat.filtered <- readRDS("/scratch/hnatri/CART/sc.dat.filtered_mmEnsembl.rds")
bk.dat <- readRDS("/scratch/hnatri/CART/bk.dat_mmEnsembl.rds")

length(intersect(colnames(sc.dat.filtered), colnames(bk.dat)))

sc.dat.filtered[1:10, 1:10]
bk.dat[1:4, 1:10]

# Gene expression concordance between bulk and sc
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bk.dat,
                pdf.prefix = NULL)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/bulk_sc_JAK1KOWT_corr.pdf"
pdf(file = filename,
    width = 8,
    height = 8)

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

head(theta)


#write.table(theta, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bayesprism_res_KOWTscref_Hhgenename.tsv", quote = F, row.names = T, sep = "\t")
theta <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bayesprism_res_KOWTscref_Hhgenename.tsv", header = T)

# Plotting
plot_cols <- c(tumor_celltype_col, immune_fibro_celltype_col)
plot_cols <- plot_cols[!(names(plot_cols) %in% paste0("Tumor", seq(2,7)))]
names(plot_cols) <- gsub("Tumor1", "Tumor", names(plot_cols))

ct_cols <- as.factor(colnames(theta))
plot_cols <- colorRampPalette(brewer.pal(11, "Paired"))(length(ct_cols))
names(plot_cols) <- levels(ct_cols)

sum(theta["CAR_T",])

summary(theta$Tumor)
summary(theta$Stromal1)

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
  ylim(0, 40) + 
  xlim(0, 40) +
  theme(legend.position = "none") + 
  ylab("Control") +
  xlab("SPP1 + CAR T")
#manuscript_theme

theta %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "celltype") %>%
  filter(celltype == "M6")
