#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/17/2024
# Description: Lymphoid/myeloid proportions
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
library(survminer) # for theme_classic2()

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
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered.rds")

# Keeping immune cells only
immune_fibro <- subset(immune_fibro, subset = celltype %in% c("F1", "F2", "F3", "B1"), invert = T)
immune_fibro$lm <- ifelse(immune_fibro$celltype %in% paste0("L", seq(1, 10)), "L", "M")

top <- c("109", "265", "181", "301", "223", "141")
bottom <- c("275", "228", "237", "234", "224", "129", "248", "185", "146")

immune_fibro$SPP1_surv_extremes <- ifelse(immune_fibro$UPN %in% top, "top",
                                          ifelse(immune_fibro$UPN %in% bottom, "bottom", NA))

prop_table <- as.data.frame(table(immune_fibro@meta.data[,"lm"], as.character(immune_fibro@meta.data[,"SPP1_surv_extremes"])))
colnames(prop_table) <- c("lm", "SPP1_surv_extremes", "Freq")
prop_table <- spread(prop_table, lm, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, lm, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

plot_colors <- c("#4CA3B6", "#9ACA3C")

barplot <- ggplot(prop_table, aes(x=SPP1_surv_extremes, y = Freq, fill = lm)) +
  geom_bar(stat="identity", position='stack', width = 0.8) +
  scale_fill_manual(name = "", values = plot_colors) +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #NoLegend() +
  manuscript_theme

barplot

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1hvsl_myel_lymph_prop_barplot.pdf"
pdf(file = filename,
    width = 3,
    height = 4)

barplot

dev.off()

# X2 test
head(prop_table)

data <- t(matrix(prop_table$Freq, nrow = 2))
chisq.test(data)
