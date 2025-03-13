#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/07/2024
# Description: SPP1 by UPN, survival analysis
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
library(survival)
library(survminer)
library(dittoSeq)
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

#==============================================================================#
# Import data
#==============================================================================#

# The final object
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")
#immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#==============================================================================#
# Pseudobulk gene expression in whole tumor
#==============================================================================#

# Genes to test
gs4_deauth()

gene_table  <- gs4_get("https://docs.google.com/spreadsheets/d/1Io8SYwIcMlUbkg1yLvcX8uwRxv7udK-BcC3_jiJovPk/edit?usp=sharing")
sheet_names(gene_table)
gene_table <- read_sheet(gene_table, sheet = "Sheet1")
head(gene_table)
length(unique(gene_table$feature))

plot_features <- c("PTMA", "PFN1", "CFL1", "TMSB4X", "TPT1", "TMSB10", "MIF",
                   "ARG1", "PDPN", "NLRP3", "IL1B", "CCL4", "S100A8", "S100A9",
                   "S100A10", "TYROBP", "CD68", "ICAM1", "C1QA", "C1QB", "C1QC",
                   "CD74", "AREG", "CD4", "APOE", "FABP5", "SPP1", "CD274",
                   "CD96", "PTPRC", "CEMIP2", "KLRD1", "CD8A", "NKG7", "IL32",
                   "CD3D", "BTG1", "IFITM2", "ITM2A", "SELL", "GZMB", "CD79A",
                   "ACTA2", "PDGFRB", "COL1A1")

test_genes <- c(unique(gene_table$feature), c("SPP1"), plot_features)

lognorm_counts <- LayerData(tumors,
                            assay = "RNA",
                            layer = "data",
                            features = unique(test_genes))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$UPN <- tumors$UPN

lognorm_counts <- lognorm_counts %>% group_by(UPN) %>%
  summarise_at(vars(-group_cols()), mean) %>%
  as.data.frame()

#lognorm_counts_spp1 <- lognorm_counts %>% group_by(UPN) %>%
#  dplyr::summarize(ave_SPP1 = mean(SPP1)) %>%
#  as.data.frame()

lognorm_counts <- t(lognorm_counts) %>%
  as.data.frame()

colnames(lognorm_counts) <- paste0("g", lognorm_counts[1,])
lognorm_counts <- lognorm_counts[2:nrow(lognorm_counts),]

pseudobulk_exp <- lognorm_counts

pseudobulk_exp <- sapply(pseudobulk_exp, as.numeric )
rownames(pseudobulk_exp) <- rownames(lognorm_counts)


metadata <- tumors@meta.data %>%
  dplyr::select("UPN",
                "Death.Status",
                "Survival.time.in.Months.from.surgery",
                "Diagnosis.Histology",
                "Grade",
                "binary_response",
                "EGFR.Amplification",
                "EGFRvIII..Exon.2.7.deletion",
                "EGFR.Missense",
                "TP53.Frameshift",
                "TP53.Missense",
                "TP53.Nonsense",
                "IDH1.Missense",
                "IDH2.Missense",
                "TERT.Promoter.Mutation",
                "PTEN.Codon.Deletion",
                "PTEN.Frameshift",
                "PTEN.Interference.of.splice.acceptor.site.in.Intron.6") %>%
  distinct()

# Looking at SPP1 and different metadata
spp1_pseudobulk_gene_exp <- pseudobulk_exp %>%
  as.data.frame() %>%
  dplyr::filter(row.names(.) %in% "SPP1") %>%
  pivot_longer(cols = colnames(pseudobulk_exp),
               names_to = "UPN",
               values_to = "SPP1") %>%
  mutate(UPN = gsub("g", "", UPN)) %>%
  mutate(SPP1 = as.numeric(SPP1))

# Adding survival info
spp1_exp_survival <- merge(spp1_pseudobulk_gene_exp, metadata, by = "UPN")

spp1_exp_survival$GBM <- ifelse(spp1_exp_survival$Diagnosis.Histology == "Glioblastoma, NOS",
                                "GBM", "HGG")

#spp1_exp_survival$log_SPP1 <- log(spp1_exp_survival$SPP1+1)
spp1_scatter <- spp1_exp_survival %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response",
                    shape = "GBM",
                    label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  #geom_text(data = subset(spp1_exp_survival,
  #                        (Survival.time.in.Months.from.surgery > 15 & SPP1 < 3) |
  #                          (Survival.time.in.Months.from.surgery < 15 & SPP1 > 3)),
  #          aes(label = UPN),
  #          size = 2.5,
  #          nudge_x = 0.1, nudge_y = 1) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  #geom_hline(yintercept = 15, color = "gray89") +
  #geom_vline(xintercept = 50, color = "gray89") +
  #geom_vline(xintercept = 3, color = "gray89") +
  theme_bw() +
  NoLegend() +
  theme_classic()  +
  manuscript_theme

spp1_scatter

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/wholetumor_SPP1_survival_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

spp1_scatter + NoLegend()

dev.off()

# R2 and p for correlation
cor.test(spp1_exp_survival$SPP1, spp1_exp_survival$Survival.time.in.Months.from.surgery,
         method = "pearson")
