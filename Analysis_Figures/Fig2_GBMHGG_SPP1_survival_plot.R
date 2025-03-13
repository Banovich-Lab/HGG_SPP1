#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/07/2024
# Description: SPP1 by UPN, survival analysis, alternative plots
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
#tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

# Adding IDH status
gs4_deauth()
idh_info  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
idh_info <- read_sheet(idh_info, sheet = "Clinical info")

immune_fibro@meta.data$IDH <- mapvalues(immune_fibro$UPN,
                                        from = idh_info$UPN,
                                        to = idh_info$idh_mutated) %>%
  unlist() %>%
  as.numeric()

immune_fibro@meta.data$IDH[-which(immune_fibro@meta.data$UPN %in% idh_info$UPN)] <- NA
immune_fibro$IDH <- gsub(0, "N", immune_fibro$IDH)
immune_fibro$IDH <- gsub(1, "Y", immune_fibro$IDH)

#myeloid$SPP1high <- ifelse(myeloid$celltype %in% c("M1", "M2", "M3", "M5", "M6"), "High", "Low")
myeloid <- subset(immune_fibro, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))
myeloid <- NormalizeData(myeloid)
myeloid <- ScaleData(myeloid)

#==============================================================================#
# Pseudobulk gene expression in the myeloid fraction
#==============================================================================#

# Genes to test
gs4_deauth()

gene_table  <- gs4_get("https://docs.google.com/spreadsheets/d/1Io8SYwIcMlUbkg1yLvcX8uwRxv7udK-BcC3_jiJovPk/edit?usp=sharing")
gene_table <- read_sheet(gene_table, sheet = "Sheet1")

fgf_genes <- row.names(myeloid)[grep("FGF", row.names(myeloid))]

plot_features <- c("PTMA", "PFN1", "CFL1", "TMSB4X", "TPT1", "TMSB10", "MIF",
                   "ARG1", "PDPN", "NLRP3", "IL1B", "CCL4", "S100A8", "S100A9",
                   "S100A10", "TYROBP", "CD68", "ICAM1", "C1QA", "C1QB", "C1QC",
                   "CD74", "AREG", "CD4", "APOE", "FABP5", "SPP1", "CD274",
                   "CD96", "PTPRC", "CEMIP2", "KLRD1", "CD8A", "NKG7", "IL32",
                   "CD3D", "BTG1", "IFITM2", "ITM2A", "SELL", "GZMB", "CD79A",
                   "ACTA2", "PDGFRB", "COL1A1")

test_genes <- c(unique(gene_table$feature), c("SPP1"), plot_features)

# Pseudobulk expression across all myeloid celltypes
pseudobulk_exp <- AverageExpression(myeloid,
                                    group.by = "UPN",
                                    features = test_genes,
                                    assay = "RNA",
                                    layer = "data") # AverageExpression will exponentiate "data"

pseudobulk_exp <- PseudobulkExpression(myeloid,
                                       assays = "RNA",
                                       features = test_genes,
                                       return.seurat = FALSE,
                                       group.by = "UPN",
                                       add.ident = NULL,
                                       layer = "counts",
                                       method = "average",
                                       normalization.method = "LogNormalize",
                                       scale.factor = 10000,
                                       margin = 1,
                                       verbose = TRUE)

myeloid <- NormalizeData(myeloid)
lognorm_counts <- LayerData(myeloid,
                        assay = "RNA",
                        layer = "data",
                        features = unique(test_genes))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$UPN <- myeloid$UPN

lognorm_counts <- lognorm_counts %>% group_by(UPN) %>%
  dplyr::summarize(ave_SPP1 = mean(SPP1)) %>%
  as.data.frame()

# Sum
#agg_counts <- AggregateExpression(myeloid,
#                                  group.by = "UPN",
#                                  assay = "RNA",
#                                  slot = "counts",
#                                  features = unique(test_genes))
#agg_counts <- agg_counts$RNA
#agg_counts <- edgeR::cpm(agg_counts, normalized.lib.sizes = TRUE, log = TRUE)
#agg_counts <- t(agg_counts) %>% as.data.frame()
#agg_counts$UPN <- gsub("g", "", rownames(agg_counts))

metadata <- myeloid@meta.data %>%
  dplyr::select("UPN",
                "Death.Status",
                "Survival.time.in.Months.from.surgery",
                "Diagnosis.Histology",
                "Grade",
                "IDH",
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
#spp1_pseudobulk_gene_exp <- pseudobulk_exp$RNA %>% as.data.frame() %>%
#  dplyr::filter(row.names(.) %in% "SPP1") %>%
#  pivot_longer(cols = colnames(pseudobulk_exp$RNA),
#               names_to = "UPN",
#               values_to = "SPP1") %>%
#  mutate(UPN = gsub("g", "", UPN))

# Adding survival info
spp1_exp_survival <- merge(lognorm_counts, metadata, by = "UPN")

spp1_exp_survival$GBM <- ifelse(spp1_exp_survival$Diagnosis.Histology == "Glioblastoma, NOS",
                                "GBM", "HGG")

#==============================================================================#
# Plotting SPP1 by UPN and group
#==============================================================================#

spp1_scatter <- spp1_exp_survival %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response",
                    shape = "GBM",
                    label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  geom_text(data = subset(spp1_exp_survival,
                          (Survival.time.in.Months.from.surgery > 15 & ave_SPP1 < 3) |
                            (Survival.time.in.Months.from.surgery < 15 & ave_SPP1 > 3)),
            aes(label = UPN),
            size = 2.5,
            nudge_x = 0.1, nudge_y = 1) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  geom_hline(yintercept = 15, color = "gray89") +
  #geom_vline(xintercept = 2.75, color = "gray89") +
  geom_vline(xintercept = 3, color = "gray89") +
  theme_bw() +
  NoLegend() +
  theme_classic()  +
  manuscript_theme

spp1_scatter

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/logSPP1_survival_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

spp1_scatter + NoLegend()

dev.off()

# R2 and p for correlation
cor.test(spp1_exp_survival$ave_SPP1, spp1_exp_survival$Survival.time.in.Months.from.surgery,
         method = "pearson")

# UPNs at the extremes
top_surv <- spp1_exp_survival %>%
  dplyr::filter(Survival.time.in.Months.from.surgery > 15, ave_SPP1 < 3) %>%
  dplyr::select("UPN")

bottom_surv <- spp1_exp_survival %>%
  dplyr::filter(Survival.time.in.Months.from.surgery < 15, ave_SPP1 > 3) %>%
  dplyr::select("UPN")

myeloid$SPP1_surv_extremes <- ifelse(myeloid$UPN %in% top_surv$UPN, "top",
                                     ifelse(myeloid$UPN %in% bottom_surv$UPN, "bottom", NA))

#==============================================================================#
# Subsets
#==============================================================================#

spp1_exp_survival$IDH <- ifelse(is.na(spp1_exp_survival$IDH), "NA", spp1_exp_survival$IDH)

# Shape by IDH status
spp1_scatter_idh <- spp1_exp_survival %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response",
                    shape = "IDH",
                    label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  geom_text(data = subset(spp1_exp_survival,
                          (Survival.time.in.Months.from.surgery > 15 & ave_SPP1 < 3) |
                            (Survival.time.in.Months.from.surgery < 15 & ave_SPP1 > 3)),
            aes(label = UPN),
            size = 2.5,
            nudge_x = 0.1, nudge_y = 1) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  geom_hline(yintercept = 15, color = "gray89") +
  #geom_vline(xintercept = 50, color = "gray89") +
  geom_vline(xintercept = 3, color = "gray89") +
  theme_bw() +
  NoLegend() +
  theme_classic()  +
  manuscript_theme

spp1_scatter_idh

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/logSPP1_IDH_survival_scatter.pdf"
pdf(file = filename,
    width = 4,
    height = 3)

spp1_scatter_idh

dev.off()

# Correlation with Grade IV only
spp1_exp_survival_gradeIV <- spp1_exp_survival %>%
  filter(Grade %in% c("IV-IDH", "IV-ND", "IV"))

spp1_exp_survival_gradeIV_plot <- spp1_exp_survival_gradeIV %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  theme_bw() +
  NoLegend() +
  theme_classic() +
  ggtitle("Grade IV tumors only") +
  manuscript_theme

# R2 and p for correlation
cor.test(spp1_exp_survival_gradeIV$ave_SPP1, spp1_exp_survival_gradeIV$Survival.time.in.Months.from.surgery,
         method = "pearson")

# Correlation with GBM only
spp1_exp_survival_GBM <- spp1_exp_survival %>%
  filter(Diagnosis.Histology %in% c("Glioblastoma, NOS"))

spp1_exp_survival_GBM_plot <- spp1_exp_survival_GBM %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  theme_bw() +
  NoLegend() +
  theme_classic() +
  ggtitle("GBM only") +
  manuscript_theme

# R2 and p for correlation
cor.test(spp1_exp_survival_GBM$ave_SPP1, spp1_exp_survival_GBM$Survival.time.in.Months.from.surgery,
         method = "pearson")

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/logSPP1_survival_scatter_Grade4_GBM.pdf"
pdf(file = filename,
    width = 3.5,
    height = 5)

spp1_exp_survival_gradeIV_plot / spp1_exp_survival_GBM_plot

dev.off()

# Correlation with IDH-mutated only
spp1_exp_survival_IDH <- spp1_exp_survival %>%
  filter(IDH %in% c("Y"))

spp1_exp_survival_IDH_plot <- spp1_exp_survival_IDH %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  theme_bw() +
  NoLegend() +
  theme_classic() +
  ggtitle("IDH-mutated only") +
  manuscript_theme

# R2 and p for correlation
cor.test(spp1_exp_survival_IDH$ave_SPP1, spp1_exp_survival_IDH$Survival.time.in.Months.from.surgery,
         method = "pearson")

# Correlation with IDH WT only
spp1_exp_survival_IDHWT <- spp1_exp_survival %>%
  filter(IDH %in% c("N"))

spp1_exp_survival_IDHWT_plot <- spp1_exp_survival_IDHWT %>%
  ggplot(aes_string(x = "ave_SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("SPP1 expression") +
  theme_bw() +
  NoLegend() +
  theme_classic() +
  ggtitle("IDH WT only") +
  manuscript_theme

# R2 and p for correlation
cor.test(spp1_exp_survival_IDHWT$ave_SPP1, spp1_exp_survival_IDHWT$Survival.time.in.Months.from.surgery,
         method = "pearson")

(spp1_exp_survival_IDH_plot + NoLegend()) / (spp1_exp_survival_IDHWT_plot + NoLegend())

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_survival_scatter_IDH.pdf"
pdf(file = filename,
    width = 3,
    height = 5)

(spp1_exp_survival_IDH_plot + NoLegend()) / (spp1_exp_survival_IDHWT_plot + NoLegend())

dev.off()
