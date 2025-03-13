#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/16/2025
# Description: Plotting SPP1 high/med myeloid proportion and SPP1 expression
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(ggplot2)
library(patchwork)
library(tidyr)
library(ggrepel)
library(tidyverse)

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

spp1_exp <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1exp_metadata_UPN.tsv",
                       sep = "\t", header = T)

spp1high_prop <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1prop_metadata.tsv",
                            sep = "\t", header = T)
spp1high_prop$SPP1high_int <- spp1high_prop$High + spp1high_prop$Int

spp1_prop_exp <- merge(spp1_exp, spp1high_prop, by = "UPN")

spp1_prop_exp_scatter <- spp1_prop_exp %>%
  #filter(SPP1high_int < 25) %>%
  ggplot(aes(y = SPP1, x = SPP1high_int,
             color = binary_response.x, shape = GBM, label = UPN)) + # label = UPN
  geom_point() +
  #geom_text(nudge_x = 2, nudge_y = 2) +
  #geom_text(aes(label = UPN),
  #          size = 2.5,
  #          nudge_x = 0.02, nudge_y = 0.02) +
  scale_color_manual(values = c("aquamarine3", "deeppink3")) +
  ylab("Ave. SPP1 expression in myeloid cells") +
  xlab("Proportion of SPP1 high/int. myeloid cells") +
  #geom_smooth(method='lm', formula= y~x) +
  #geom_hline(yintercept = 15, color = "gray89") +
  #geom_vline(xintercept = 0.2, color = "gray89") +
  theme_bw() +
  NoLegend() +
  theme_classic()
  #manuscript_theme

spp1_prop_exp_scatter

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/SPP1highint_survival_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

spp1_prop_exp_scatter + NoLegend()

dev.off()

# R2 and p for correlation
#test_data <- spp1_prop_exp %>%
#  filter(SPP1high_int < 25)

cor.test(spp1_prop_exp$SPP1, spp1_prop_exp$SPP1high_int,
         method = "pearson")
