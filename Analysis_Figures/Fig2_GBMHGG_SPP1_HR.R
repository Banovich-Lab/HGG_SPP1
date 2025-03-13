#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/24/2024
# Description: SPP1/survival hazard ratio
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
#library(dittoSeq)
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

# Subsetting out myeloid cell types only
myeloid <- subset(immune_fibro, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))
#lymphoid <- subset(immune_fibro, subset = celltype %in% paste0("L", seq(1, 10)))
#fibroblast <- subset(immune_fibro, subset = celltype %in% paste0("F", seq(1, 3)))
#bcells <- subset(immune_fibro, subset = celltype %in% c("B1"))

myeloid$SPP1high <- ifelse(myeloid$celltype %in% c("M6"), "High",
                           ifelse(myeloid$celltype %in% c("M1", "M2", "M3", "M5"), "Int", "Low"))

# M1, M2, M3, M5, and M6
myeloid$SPP1high_low <- ifelse(myeloid$celltype %in% c("M1", "M2", "M3", "M5", "M6"), "High_Int", "Low")

# Genes to test
gs4_deauth()

gene_table  <- gs4_get("https://docs.google.com/spreadsheets/d/1Io8SYwIcMlUbkg1yLvcX8uwRxv7udK-BcC3_jiJovPk/edit?usp=sharing")
sheet_names(gene_table)
gene_table <- read_sheet(gene_table, sheet = "Sheet1")

test_genes <- c(unique(gene_table$feature), c("SPP1"))

#==============================================================================#
# Calculating hazard ratios for gene expression
#==============================================================================#

hr_test_genes <- c("SPP1", "C1QA", "C1QB", "C1QC", "CD68", "CD14",
                   "IL1B", "MIF", "CD74",  "AREG", "TYROBP")

#myeloid$x <- "x"
#VlnPlot(myeloid,
#        group.by = "x",
#        features = hr_test_genes,
#        split.by = "binary_response",
#        split.plot = T,
#        pt.size = 0,
#        ncol = 4,
#        cols = c("aquamarine3", "deeppink3")) &
#  manuscript_theme
#
## SPP1 HR in each compartment
#object_list <- list("myeloid" = myeloid,
#                    "lymphoid" = lymphoid,
#                    "bcells" = bcells,
#                    "fibroblast" = fibroblast)
#pseudobulk_list <- lapply(names(object_list), function(xx){
#  obj <- object_list[[xx]]
#  pseudobulk_exp <- AverageExpression(obj,
#                                      group.by = "UPN",
#                                      features = "SPP1",
#                                      assay = "RNA",
#                                      slot = "data")
#  pseudobulk_exp <- pseudobulk_exp$RNA %>%
#    #t() %>%
#    as.data.frame() %>%
#    #rownames_to_column(var = "UPN") %>%
#    pivot_longer(cols = colnames(pseudobulk_exp$RNA),
#                 names_to = "UPN",
#                 values_to = "SPP1") %>%
#    mutate(UPN = gsub("g", "", UPN))
#  
#  exp_survival <- merge(pseudobulk_exp, metadata, by = "UPN")
#  exp_survival$compartment <- xx
#  
#  exp_survival
#})
#
#head(pseudobulk_list[[1]])
#names(pseudobulk_list) <- names(object_list)
#
#df <- do.call("rbind", pseudobulk_list) %>%
#  dplyr::select(c("Survival.time.in.Months.from.surgery", 
#                  "Death.Status",
#                  "UPN",
#                  "SPP1",
#                  "compartment")) %>%
#  pivot_wider(#cols = c("Survival.time.in.Months.from.surgery", 
#              #         "Death.Status",
#                       #"SPP1",
#              #         "UPN"),
#    values_from = "SPP1",
#    names_from = "compartment")
#
#head(df)
#
## Fit survival data using the Kaplan-Meier method
#surv_object <- Surv(time = df$Survival.time.in.Months.from.surgery,
#                    event = df$Death.Status)
#surv_object 
#
## Fit a Cox proportional hazards model for a single gene
#fit.coxph <- coxph(surv_object ~ myeloid + lymphoid + bcells + fibroblast,
#                   data = df, na.action = na.omit)
#p2 <- ggforest(fit.coxph, data = df)
#
#
#p1 + p2
#
#filename <- "/home/hnatri/CART/13384_Tumors/Plots/GBM_SPP1_HR_all_compartments.pdf"
#pdf(file = filename,
#    width = 5,
#    height = 10)
#p1 / p2
#dev.off()

#==============================================================================#
# Calculating hazard ratios for cell type proportions
#==============================================================================#

prop_table <- as.data.frame(table(immune_fibro@meta.data$celltype, as.character(immune_fibro@meta.data$UPN)))
colnames(prop_table) <- c("celltype", "UPN", "Freq")
prop_table <- spread(prop_table, celltype, Freq)
# Converting to percentage
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, celltype, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

prop_table <- prop_table %>% pivot_wider(names_from = celltype,
                                         values_from = Freq)
prop_table %>% filter(UPN == 109) %>%
  dplyr::select(!(UPN)) %>%
  as.numeric() %>% sum()

# Adding survival data
metadata <- immune_fibro@meta.data %>%
  dplyr::select("UPN",
                "Death.Status",
                "Survival.time.in.Months.from.surgery",
                "binary_response",
                "Grade",
                "Diagnosis.Histology",
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

prop_survival <- merge(prop_table, metadata, by = "UPN")

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = prop_survival$Survival.time.in.Months.from.surgery,
                    event = prop_survival$Death.Status)
surv_object 

# Fit a Cox proportional hazards model for a single gene
fit.coxph <- coxph(surv_object ~ M1 + M2 + M3 + M4 + M5 + M6 + M7 + M8 + M9 + N1 +
                     L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L8 + L9 + L10 + B1 +
                     F1 + F2 + F3,
                   data = prop_survival)

p1 <- ggforest(fit.coxph, data = prop_survival)

p1

# Proportion of myeloid cells out of all immune cells
immune <- subset(immune_fibro, subset = celltype %in% c(c("N1", paste0("M", seq(1, 9))), paste0("L", seq(1, 10))))
immune$lineage <- ifelse(immune$celltype %in% c("N1", paste0("M", seq(1, 9))), "myeloid", "lymphoid")

prop_table <- as.data.frame(table(immune@meta.data$lineage, as.character(immune@meta.data$UPN)))
colnames(prop_table) <- c("lineage", "UPN", "Freq")
prop_table <- spread(prop_table, lineage, Freq)
# Converting to percentage
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, lineage, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

prop_table <- prop_table %>% pivot_wider(names_from = lineage,
                                         values_from = Freq)

# Adding survival data
prop_survival <- merge(prop_table, metadata, by = "UPN")

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = prop_survival$Survival.time.in.Months.from.surgery,
                    event = prop_survival$Death.Status)
surv_object 

# Fit a Cox proportional hazards model for a single gene
fit.coxph <- coxph(surv_object ~ myeloid + lymphoid,
                   data = prop_survival)

p2 <- ggforest(fit.coxph, data = prop_survival)

p2

# SPP1 high vs. low myeloid clusters
prop_table <- as.data.frame(table(myeloid@meta.data$SPP1high, as.character(myeloid@meta.data$UPN)))
colnames(prop_table) <- c("SPP1exp", "UPN", "Freq")
prop_table <- spread(prop_table, SPP1exp, Freq)
# Converting to percentage
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, SPP1exp, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

prop_table <- prop_table %>% pivot_wider(names_from = SPP1exp,
                                         values_from = Freq)

# Adding survival data
prop_survival <- merge(prop_table, metadata, by = "UPN")

write.table(prop_survival, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1prop_metadata.tsv",
            quote = F, sep = "\t", row.names = F)

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = prop_survival$Survival.time.in.Months.from.surgery,
                    event = prop_survival$Death.Status)
surv_object 

# Fit a Cox proportional hazards model for a single gene
fit.coxph <- coxph(surv_object ~ High + Int + Low,
                   data = prop_survival)

p3 <- ggforest(fit.coxph, data = prop_survival)

p3

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highintlow_HR.pdf"
pdf(file = filename,
    width = 5,
    height = 3)

p3

dev.off()

#==============================================================================#
# SPP1 high and intermediate vs. low myeloid clusters
#==============================================================================#

#myeloid$SPP1high_low <- ifelse(myeloid$SPP1high %in% c("High", "Int"), "High_Int", "Low")

prop_table <- as.data.frame(table(myeloid@meta.data$SPP1high_low, as.character(myeloid@meta.data$UPN)))
colnames(prop_table) <- c("SPP1high_low", "UPN", "Freq")
prop_table <- spread(prop_table, SPP1high_low, Freq)
# Converting to percentage
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, SPP1high_low, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

prop_table %>% filter(UPN == 146) %>%
  dplyr::select(Freq) %>%
  unlist() %>% sum()

prop_table <- prop_table %>% pivot_wider(names_from = SPP1high_low,
                                         values_from = Freq)

# Adding survival data
prop_survival <- merge(prop_table, metadata, by = "UPN")

prop_survival %>%
  ggplot(aes(x = High_Int, y = Survival.time.in.Months.from.surgery, label = UPN)) +
  geom_point() +
  geom_text(hjust=0, vjust=0) +
  theme_bw() +
  xlab("SPP1 high/intermediate proportion")

prop_survival$GBM <- ifelse(prop_survival$Diagnosis.Histology == "Glioblastoma, NOS",
                            "GBM", "HGG")

spp1highint_scatter <- prop_survival %>%
  ggplot(aes_string(x = "High_Int", y = "Survival.time.in.Months.from.surgery",
                    color = "binary_response", shape = "GBM", label = "UPN")) +
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  #geom_text(data = subset(prop_survival,
  #                        (Survival.time.in.Months.from.surgery > 15 & High_Int < 0.2) |
  #                          (Survival.time.in.Months.from.surgery < 15 & High_Int > 0.2)),
  #          aes(label = UPN),
  #          size = 2.5,
  #          nudge_x = 0.02, nudge_y = 2) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Survival (months from surgery)") +
  xlab("Proportion of SPP1 high/int. myeloid cells") +
  geom_hline(yintercept = 15, color = "gray89") +
  #geom_vline(xintercept = 0.2, color = "gray89") +
  theme_bw() +
  NoLegend() +
  theme_classic()  +
  manuscript_theme

spp1highint_scatter

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highint_survival_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

spp1highint_scatter + NoLegend()

dev.off()

# Significance, cor = -0.2968635, p-value = 0.07438
cor.test(prop_survival$High_Int, prop_survival$Survival.time.in.Months.from.surgery)

#prop_survival <- prop_survival %>%
#  filter(!(UPN == 265))

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = prop_survival$Survival.time.in.Months.from.surgery,
                    event = prop_survival$Death.Status)
surv_object 

# Fit a Cox proportional hazards model for a single gene
fit.coxph <- coxph(surv_object ~ High_Int + Low,
                   data = prop_survival)

p5 <- ggforest(fit.coxph, data = prop_survival)

p5

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highint_survival_HR_pct.pdf"
pdf(file = filename,
    width = 7,
    height = 2)

p5

dev.off()

cor.test(spp1_exp_survival$SPP1, spp1_exp_survival$Survival.time.in.Months.from.surgery,
         method = "pearson")

# GBM only
prop_survival_gbm <- prop_survival %>% filter(Diagnosis.Histology == "Glioblastoma, NOS")

surv_object <- Surv(time = prop_survival_gbm$Survival.time.in.Months.from.surgery,
                    event = prop_survival_gbm$Death.Status)

fit.coxph <- coxph(surv_object ~ High_Int + Low,
                   data = prop_survival_gbm)


p6 <- ggforest(fit.coxph, data = prop_survival)

p6

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highint_survival_HR_GBMonly.pdf"
pdf(file = filename,
    width = 4,
    height = 2)

p6

dev.off()

# Grade IV only
prop_survival_gradev <- prop_survival %>% filter(Grade == "IV")

surv_object <- Surv(time = prop_survival_gradev$Survival.time.in.Months.from.surgery,
                    event = prop_survival_gradev$Death.Status)

fit.coxph <- coxph(surv_object ~ High_Int + Low,
                   data = prop_survival_gradev)


p7 <- ggforest(fit.coxph, data = prop_survival)

p7

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highint_survival_HR_Grade4only.pdf"
pdf(file = filename,
    width = 4,
    height = 2)

p7

dev.off()

