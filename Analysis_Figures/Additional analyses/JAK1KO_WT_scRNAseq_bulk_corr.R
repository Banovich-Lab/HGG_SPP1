#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/02/2024
# Description: JAK1KO/WT scRNAseq and bulk gene expression correlation
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/Github/Analysis_Figures/13384_tumor_ms_themes.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/Github/Analysis_Figures/CellChat_heatmap.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Importing data
#==============================================================================#

sc.dat.filtered <- readRDS("/scratch/hnatri/CART/sc.dat.filtered.rds")
bk.dat <- readRDS("/scratch/hnatri/CART/bk.dat.rds")

#==============================================================================#
# Plotting
#==============================================================================#

gene.shared <- intersect(colnames(sc.input), colnames(bulk.input))

# align reference and mixture
bulk.input <- bulk.input[, gene.shared]
sc.input <- sc.input[, gene.shared]

# gencode.v22.broad.category.txt
# genelist.mm.new.txt
gene.tab.path <- system.file("extdata","gencode.v22.broad.category.txt", package="BayesPrism")	
gene.list <- read.table(gene.tab.path, sep="\t",header=F,stringsAsFactors=F)

intersect(gene.shared, gene.list[,2])

gene.df <- gene.list[match(gene.shared, gene.list[,5]),c(5,9)]
colnames(gene.df) <- c("gene.name", "category")

plot.df <- data.frame(log2.bulk = log2(colSums(bulk.input)),
                      log2.sc = log2(colSums(sc.input)),
                      gene.df)

head(plot.df)

write.table(plot.df, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bulk_sc_corr.tsv",
            sep = "\t", quote = F, row.names = F)

plot.df <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/bulk_sc_corr.tsv",
                      sep = "\t", header = T)

scatterp <- plot.df %>%
  ggplot(aes(x = log2.sc, y = log2.bulk)) +
  geom_point(color = "dodgerblue4") +
  theme_bw() +
  xlim(0, 25) +
  ylim(0, 25) +
  xlab("log2(sc-RNAseq expression)") +
  ylab("log2(bulk-RNAseq expression)") +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = "dashed") +
  coord_fixed()

scatterp

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/YUMM_scRNAseq_bulk_corr_scatter.pdf"
pdf(file = filename,
    width = 4,
    height = 4)

scatterp

dev.off()


