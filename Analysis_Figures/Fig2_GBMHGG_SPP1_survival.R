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

myeloid <- subset(immune_fibro, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))
myeloid$SPP1high <- ifelse(myeloid$celltype %in% c("M6"), "High",
                           ifelse(myeloid$celltype %in% c("M1", "M2", "M3", "M5"), "Int", "Low"))

# M1, M2, M3, M5, and M6
myeloid$SPP1high_low <- ifelse(myeloid$celltype %in% c("M1", "M2", "M3", "M5", "M6"), "High_Int", "Low")

myeloid <- NormalizeData(myeloid)
myeloid <- ScaleData(myeloid)

#==============================================================================#
# Plotting SPP1 by UPN and group
#==============================================================================#

VlnPlot(myeloid,
        features = c("SPP1"),
        split.by = "UPN",
        group.by = "celltype",
        pt.size = 0) +
  coord_flip()

vln_idh <- VlnPlot(subset(myeloid, subset = IDH %in% c("Y", "N")),
        features = c("SPP1"),
        group.by = "IDH",
        #group.by = "celltype",
        pt.size = 0,
        cols = c("aquamarine3", "deeppink3"),
        split.plot = F) &
  theme_classic() +
  NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_vln_IDH_response.pdf"
pdf(file = filename,
    width = 3,
    height = 2.5)

vln_idh

dev.off()


vln_idh <- VlnPlot(subset(myeloid, subset = IDH %in% c("N", "Y")),
        features = c("SPP1"),
        group.by = "IDH",
        split.by = "binary_response",
        slot = "data",
        pt.size = 0,
        cols = c("aquamarine3", "deeppink3"),
        split.plot = F) &
  theme_classic()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_vln_IDH_response.pdf"
pdf(file = filename,
    width = 6,
    height = 2.5)

vln_idh

dev.off()

p1 <- VlnPlot(myeloid,
        features = c("SPP1"),
        split.by = "binary_response",
        group.by = "celltype",
        pt.size = 0,
        cols = c("aquamarine3", "deeppink3"),
        split.plot = T) &
  theme_classic() &
  manuscript_theme

p2 <- VlnPlot(myeloid,
        features = c("SPP1"),
        split.by = "CD3_high_low",
        group.by = "celltype",
        pt.size = 0,
        cols = c("aquamarine3", "deeppink3"),
        split.plot = T) &
  theme_classic() &
  manuscript_theme

p1 / p2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/SPP1_violinplots.pdf"
pdf(file = filename,
    width = 5,
    height = 3)

p1 / p2

dev.off()

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
#pseudobulk_exp <- AverageExpression(myeloid,
#                                    group.by = "UPN",
#                                    features = test_genes,
#                                    assay = "RNA",
#                                    layer = "data") # AverageExpression will exponentiate "data"

lognorm_counts <- LayerData(myeloid,
                            assay = "RNA",
                            layer = "data",
                            features = unique(test_genes))
lognorm_counts <- lognorm_counts %>% t() %>% as.data.frame()
lognorm_counts$UPN <- myeloid$UPN

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

# Looping through each gene and plotting
plot_list <- lapply(test_genes, function(gene){
  message(gene)

  pseudobulk_gene_exp <- pseudobulk_exp %>%
    filter(row.names(.) %in% gene) %>%
    pivot_longer(cols = colnames(pseudobulk_exp),
                 names_to = "UPN",
                 values_to = gene) %>%
    mutate(UPN = gsub("g", "", UPN))
  
  # Adding survival info
  exp_survival <- merge(pseudobulk_gene_exp, metadata, by = "UPN")
  
  # R2 and p for correlation
  cor_res <- cor.test(exp_survival[,gene], exp_survival$Survival.time.in.Months.from.surgery,
           method = "pearson")
  
  cor <- as.numeric(cor_res$estimate)
  pval <- as.numeric(cor_res$p.value)
  xpos <- max(exp_survival[,gene])
  if(xpos < 1){
    xpos <- 1
  } else {
    xpos <- xpos
  }
  
  exp_survival_plot <- exp_survival %>%
    ggplot(aes_q(x = as.name(gene), y = as.name("Survival.time.in.Months.from.surgery"), color = as.name("binary_response"))) + # , label = "UPN"
      geom_point() +
      scale_color_manual(values = c("aquamarine3", "deeppink3")) +
      ylab("Survival (months from surgery)") +
      xlab(paste0("Mean log-normalized ", gene, " expression in myeloid cells")) +
      theme_bw() +
      NoLegend() +
      theme_classic()  +
      manuscript_theme
      

  return(exp_survival_plot)
})

names(plot_list) <- test_genes

# Saving as a pdf
filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/all_genes_survival_scatter.pdf"

pdf(file=filename,
    width = 4.5,
    height = 3)

# One plot per page
for (plot in plot_list) {   
  print(plot)
}

dev.off()

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

write.table(spp1_exp_survival, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1exp_metadata_UPN.tsv",
            quote = F, sep = "\t", row.names = F)

spp1_scatter <- spp1_exp_survival %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response",
                    shape = "GBM",
                    label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  geom_text(data = subset(spp1_exp_survival,
                          (Survival.time.in.Months.from.surgery > 15 & SPP1 < 3) |
                            (Survival.time.in.Months.from.surgery < 15 & SPP1 > 3)),
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

spp1_scatter

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_survival_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

spp1_scatter + NoLegend()

dev.off()

# R2 and p for correlation
cor.test(spp1_exp_survival$SPP1, spp1_exp_survival$Survival.time.in.Months.from.surgery,
         method = "pearson")

spp1_exp_survival$IDH <- ifelse(is.na(spp1_exp_survival$IDH), "NA", spp1_exp_survival$IDH)

# Shape by IDH status
spp1_scatter_idh <- spp1_exp_survival %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response",
                    shape = "IDH",
                    label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  geom_text(data = subset(spp1_exp_survival,
                          (Survival.time.in.Months.from.surgery > 15 & SPP1 < 3) |
                            (Survival.time.in.Months.from.surgery < 15 & SPP1 > 3)),
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

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_IDH_survival_scatter.pdf"
pdf(file = filename,
    width = 4,
    height = 3)

spp1_scatter_idh

dev.off()

# DE of SPP1 between IDH WT or mutated for Extended figure
Idents(myeloid) <- myeloid$IDH
markers <- FindMarkers(myeloid,
                       ident.1 = "Y",
                       ident.2 = "N",
                       assay = "RNA",
                       verbose = F)
markers$feature <- rownames(markers)

markers_sig <- markers %>%
  dplyr::filter(p_val_adj < 0.01, abs(avg_log2FC) > 2)

markers %>% dplyr::filter(feature == "SPP1")

# DE of SPP1 between outcomes in IDH WT and mutated
Idents(myeloid) <- myeloid$binary_response
markers <- FindMarkers(subset(myeloid, subset = IDH == "Y"),
                       ident.1 = "CR_SD",
                       ident.2 = "PD",
                       assay = "RNA",
                       verbose = F)
markers$feature <- rownames(markers)

markers_sig <- markers %>%
  dplyr::filter(p_val_adj < 0.01, abs(avg_log2FC) > 2)

markers %>% dplyr::filter(feature == "SPP1")

markers <- FindMarkers(subset(myeloid, subset = IDH == "N"),
                       ident.1 = "CR_SD",
                       ident.2 = "PD",
                       assay = "RNA",
                       verbose = F)
markers$feature <- rownames(markers)

markers_sig <- markers %>%
  dplyr::filter(p_val_adj < 0.01, abs(avg_log2FC) > 2)

markers %>% dplyr::filter(feature == "SPP1")

# Plotting SPP1 and different metadata
plot_list <- lapply(c("EGFR.Amplification",
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
                      "PTEN.Interference.of.splice.acceptor.site.in.Intron.6"), function(xx){
                      spp1_exp_survival %>%
                        ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery", color = xx)) +
                        geom_point() +
                        scale_color_manual(values = c("deeppink3", "aquamarine3")) +
                        ylab("Survival (months from surgery)") +
                        xlab("SPP1 expression") +
                        theme_bw() +
                        NoLegend() +
                        theme_classic()  +
                        manuscript_theme
                    })

patchwork::wrap_plots(plot_list)

# UPNs at the extremes for SPP1 high vs. low tumor comparison
top_surv <- spp1_exp_survival %>%
  dplyr::filter(Survival.time.in.Months.from.surgery > 15, SPP1 < 3) %>%
  dplyr::select("UPN")

bottom_surv <- spp1_exp_survival %>%
  dplyr::filter(Survival.time.in.Months.from.surgery < 15, SPP1 > 3) %>%
  dplyr::select("UPN")

myeloid$SPP1_surv_extremes <- ifelse(myeloid$UPN %in% top_surv$UPN, "top",
                                     ifelse(myeloid$UPN %in% bottom_surv$UPN, "bottom", NA))

# DEGs between SPP1 top vs. bottom
Idents(myeloid) <- myeloid$SPP1_surv_extremes
markers <- FindMarkers(myeloid,
                       ident.1 = "bottom",
                       ident.2 = "top",
                       assay = "RNA",
                       verbose = F)
markers$feature <- rownames(markers)

markers_sig <- markers %>%
  dplyr::filter(p_val_adj < 0.01, abs(avg_log2FC) > 2)

# Saving to a file
write.table(markers_sig,
            "/scratch/hnatri/CART/13384_myeloid_SPP1surv_top_bottom_DEGs_Seurat_sig.tsv",
            sep = "\t", quote = F, row.names = F)

# Correlation with Grade IV only
spp1_exp_survival_gradeIV <- spp1_exp_survival %>%
  filter(Grade %in% c("IV-IDH", "IV-ND", "IV"))
  
spp1_exp_survival_gradeIV_plot <- spp1_exp_survival_gradeIV %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
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
cor.test(spp1_exp_survival_gradeIV$SPP1, spp1_exp_survival_gradeIV$Survival.time.in.Months.from.surgery,
         method = "pearson")

# Correlation with GBM only
spp1_exp_survival_GBM <- spp1_exp_survival %>%
  filter(Diagnosis.Histology %in% c("Glioblastoma, NOS"))

spp1_exp_survival_GBM_plot <- spp1_exp_survival_GBM %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
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
cor.test(spp1_exp_survival_GBM$SPP1, spp1_exp_survival_GBM$Survival.time.in.Months.from.surgery,
         method = "pearson")

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_survival_scatter_Grade4_GBM.pdf"
pdf(file = filename,
    width = 4.5,
    height = 6)

spp1_exp_survival_gradeIV_plot / spp1_exp_survival_GBM_plot

dev.off()

# Correlation with IDH-mutated only
spp1_exp_survival_IDH <- spp1_exp_survival %>%
  filter(IDH %in% c("Y"))

spp1_exp_survival_IDH_plot <- spp1_exp_survival_IDH %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
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
cor.test(spp1_exp_survival_IDH$SPP1, spp1_exp_survival_IDH$Survival.time.in.Months.from.surgery,
         method = "pearson")

# Correlation with IDH WT only
spp1_exp_survival_IDHWT <- spp1_exp_survival %>%
  filter(IDH %in% c("N"))

spp1_exp_survival_IDHWT_plot <- spp1_exp_survival_IDHWT %>%
  ggplot(aes_string(x = "SPP1", y = "Survival.time.in.Months.from.surgery", color = "binary_response", label = "UPN")) +
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
cor.test(spp1_exp_survival_IDHWT$SPP1, spp1_exp_survival_IDHWT$Survival.time.in.Months.from.surgery,
         method = "pearson")

spp1_exp_survival_IDH_plot + NoLegend() / spp1_exp_survival_IDHWT_plot + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_survival_scatter_IDH.pdf"
pdf(file = filename,
    width = 3,
    height = 5)

(spp1_exp_survival_IDH_plot + NoLegend()) / (spp1_exp_survival_IDHWT_plot + NoLegend())

dev.off()

#==============================================================================#
# Violinplots
#==============================================================================#

# Plotting test features by top/bottom
VlnPlot(myeloid,
        features = test_genes,
        group.by = "SPP1_surv_extremes",
        ncol = 4,
        pt.size = 0) &
  theme_bw()

main_genes <- c("HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "HLA-DOA",
                "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", 
                "HLA-DRA", "HLA-DRB5", "HLA-DRB1", "CCL2", "CCL3", "CCL4",
                "CD86", "CD14", "TLR2", "PDL1", "CD68", "IL1A", "IL6ST")

myeloid$SPP1_surv_extremes <- gsub("top", "SPP1low", myeloid$SPP1_surv_extremes)
myeloid$SPP1_surv_extremes <- gsub("bottom", "SPP1high", myeloid$SPP1_surv_extremes)
myeloid$x <- "x"

vp1 <- VlnPlot(subset(myeloid, subset = SPP1_surv_extremes %in% c("SPP1low", "SPP1high")),
        features = main_genes,
        group.by = "x",
        split.by = "SPP1_surv_extremes",
        ncol = 5,
        split.plot = T,
        pt.size = 0,
        cols = c("orangered1", "royalblue3")) &
  theme_classic2() &
  NoLegend() &
  manuscript_theme &
  xlab("") &
  ylab("")

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_surv_extremes_HLA_vlnplot_recolor.pdf"
pdf(file = filename,
    width = 5,
    height = 4.5)

vp1

dev.off()

extra_genes <- c("USP53", "STAB1", "ST6GAL1", "SLC40A1", "SELENOP", "SDC3",
                 "PLD4", "OLFML3", "OLDML2B", "MERTK", "IL1RA", "GPR34",
                 "FOLR2", "FCGBP", "F13A1", "EPB41L2", "ENPP2", "C3", "ADORA3",
                 "ADAM28", "FTH1", "FTL", "LDLR", "SLC40A")

VlnPlot(subset(myeloid, subset = SPP1_surv_extremes %in% c("SPP1low", "SPP1high")),
        features = extra_genes,
        #group.by = "x",
        group.by = "SPP1_surv_extremes",
        ncol = 5,
        #split.plot = T,
        pt.size = 0,
        cols = c("deeppink3", "aquamarine3")) &
  theme_classic2() &
  NoLegend() &
  manuscript_theme &
  xlab("") &
  ylab("")

# Metabolism genes
meta_gene_table  <- gs4_get("https://docs.google.com/spreadsheets/d/16sI_lATPquB5yHQYXtXVs9UORQ6TNL1cqh2Yh3jy_Ho/edit?usp=sharing")
meta_gene_table <- read_sheet(meta_gene_table, sheet = "Sheet1")
meta_gene_table <- meta_gene_table %>% filter(Plot == "y")

DefaultAssay(myeloid) <- "RNA"
VlnPlot(subset(myeloid, subset = SPP1_surv_extremes %in% c("SPP1low", "SPP1high")),
        features = unique(meta_gene_table$Gene),
        group.by = "x",
        split.by = "SPP1_surv_extremes",
        ncol = 7,
        split.plot = T,
        pt.size = 0,
        cols = c("deeppink3", "aquamarine3")) &
  theme_classic2() &
  NoLegend() &
  manuscript_theme &
  xlab("") &
  ylab("")

# All genes
all_genes <- c("HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "HLA-DOA",
               "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", 
               "HLA-DRA", "HLA-DRB5", "HLA-DRB1", "CCL2", "CCL3", "CCL4",
               "CD86", "CD14", "TLR2", "PDL1", "CD68", "IL1A", "IL6ST",
               "USP53", "STAB1", "ST6GAL1", "SLC40A1", "SELENOP", "SDC3",
               "PLD4", "OLFML3", "OLDML2B", "MERTK", "IL1RA", "GPR34",
               "FOLR2", "FCGBP", "F13A1", "EPB41L2", "ENPP2", "C3", "ADORA3",
               "ADAM28", "FTH1", "FTL", "LDLR", "SLC40A",
               "MCTP1", "PKM", "GAPDH", "LDHA", "MGLL", "AHR", "ECHS1",
               "PTGES2", "APOE", "APOC1", "CD63", "CD36", "ATF1", "C1QA",
               "C1QB", "C1QC", "CD163", "CHI3L1", "IRF3", "FABP4", "FABP5",
               "CTSB", "CTSD", "CTSL", "F13A1", "FOLR2", "GPNMB", "LGALS3",
               "LPL", "LIPA", "MACRO", "MERTK", "MMP7", "MMP9", "MMP12", "MRC1",
               "NR1H3", "NRF1", "NUPR1", "PLA2G7", "RNASE1", "SPARC", "TFDP2",
               "TREM2", "ZEB1", "EPAS1", "CXCR4", "SEMA3A", "NRP1", "BMAL1",
               "ACLY", "SLC7A11", "PTGES", "PTGS2", "ABCA1", "ABCG1", "ADORA2A",
               "NT5E", "ENTPD1", "HIF1A", "FABP5", "ACP5")

all_genes_stats <- markers %>% filter(feature %in% all_genes)
write.table(all_genes_stats, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1highvslow_DEGs_mainfig.tsv",
            quote = F, sep = "\t", row.names = F)

DefaultAssay(myeloid) <- "RNA"
v1 <- VlnPlot(subset(myeloid, subset = SPP1_surv_extremes %in% c("SPP1low", "SPP1high")),
        features = all_genes,
        group.by = "x",
        split.by = "SPP1_surv_extremes",
        ncol = 5,
        split.plot = T,
        pt.size = 0,
        cols = c("orangered1", "royalblue3")) &
  theme_classic2() &
  NoLegend() &
  manuscript_theme &
  xlab("") &
  ylab("") &
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1_surv_extremes_HLA_vlnplot_recolor_moregenes.pdf"
pdf(file = filename,
    width = 5,
    height = 15)

v1

dev.off()

#==============================================================================#
# Testing for significance using Cox's proportional hazards model
#==============================================================================#

cox_res <- coxph(Surv(Survival.time.in.Months.from.surgery, Death.Status) ~ SPP1, data = spp1_exp_survival)
cox_res

summary(cox_res)

fit <- survfit(cox_res)

p1 <- ggforest(cox_res, data = spp1_exp_survival)

dependent <- "Surv(Survival.time.in.Months.from.surgery, Death.Status)"
explanatory <- c("SPP1")

spp1_exp_survival %>% 
  finalfit(dependent, explanatory, add_dependent_label = FALSE)

spp1_exp_survival %>% 
  coxphmulti(dependent, explanatory) %>% 
  cox.zph() %>% 
  {zph_result <<- .} %>% 
  plot(var=1)

#==============================================================================#
# Cell type proportion differences and DEGs between the extremes
#==============================================================================#

immune_fibro$SPP1_surv_extremes <- ifelse(immune_fibro$UPN %in% top_surv$UPN, "top",
                                          ifelse(immune_fibro$UPN %in% bottom_surv$UPN, "bottom", NA))
immune_fibro$SPP1_surv_extremes <- gsub("top", "SPP1low", immune_fibro$SPP1_surv_extremes)
immune_fibro$SPP1_surv_extremes <- gsub("bottom", "SPP1high", immune_fibro$SPP1_surv_extremes)

DimPlot(immune_fibro,
        group.by = "celltype",
        split.by = "SPP1_surv_extremes",
        label = T,
        cols = immune_fibro_celltype_col) +
  theme_classic2() +
  NoLegend() +
  NoAxes() +
  ggtitle("")

prop_test <- sc_utils(immune_fibro)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_2 = "SPP1high", sample_1 = "SPP1low",
  sample_identity = "SPP1_surv_extremes")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/tumor_immune_fibro_celltype_prop_forest_SPP1low_high.pdf"
pdf(file = filename,
    width = 1.8,
    height = 4)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()

dev.off()

