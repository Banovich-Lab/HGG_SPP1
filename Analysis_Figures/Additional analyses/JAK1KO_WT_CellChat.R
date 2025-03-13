#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2022/08/31
# Description: CellChat analysis on 13384 CAR T trial tumor samples
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(CellChat)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/CART/13384_Tumors/SPP1_ms/13384_tumor_ms_themes.R")
source("/home/hnatri/CART/13384_Tumors/SPP1_ms/cellchat_heatmap.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import data
#==============================================================================#

# Cell chat objects
cc_mm_ko_ctrl <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/CellChat/cellchat_KO_Control.rds")
cc_mm_ko_ctrl
cellchat_KO_Control_CD45_Neg

cellchat_KO_Control_CD45_Neg@options$datatype <- "RNA"
cellchat_KO_Control_CD45_Neg@images <- list()

cc_mm_wt_ctrl <- load("/tgen_labs/banovich/BCTCSF/Mouse_tumor/CellChat/cellchat_WildType_Control.rds")
cc_mm_wt_ctrl
cellchat_WildType_Control_CD45_Neg

cellchat_WildType_Control_CD45_Neg@options$datatype <- "RNA"
cellchat_WildType_Control_CD45_Neg@images <- list()

# Seurat objects
#mm_ko_ctrl_cells <- load("/labs/banovich/BCTCSF/Mouse_tumor/CellChat/Jak1_KO_Control_Cells.rds")
#Jak1_KO_Control_Cells
#head(Jak1_KO_Control_Cells@meta.data)
#unique(Jak1_KO_Control_Cells$orig.ident)

#==============================================================================#
# CellChat analysis
#==============================================================================#

# Creating the merged object for comparative analysis
head(cellchat_KO_Control_CD45_Neg@meta)
unique(cellchat_KO_Control_CD45_Neg@meta$JAK1_Group)
head(cellchat_WildType_Control_CD45_Neg@meta)
unique(cellchat_WildType_Control_CD45_Neg@meta$JAK1_Group)

CC_objects_compare <- list("Jak1_KO" = cellchat_KO_Control_CD45_Neg,
                           "WT" = cellchat_WildType_Control_CD45_Neg)
CC_merged_object <- mergeCellChat(CC_objects_compare, add.names = names(CC_objects_compare))

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(CC_merged_object, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(CC_merged_object, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Compute the network centrality scores
CC_objects_compare <- lapply(CC_objects_compare, function(xx){
  netAnalysis_computeCentrality(xx, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
})

saveRDS(CC_merged_object, "/scratch/hnatri/CART/mm_CC_merged_object.rds")

q(save = "no")

# Compare the major sources and targets in 2D space
num.link <- sapply(CC_objects_compare, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(CC_objects_compare)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(CC_objects_compare[[i]],
                                               title = names(CC_objects_compare)[i],
                                               weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

# Identify and visualize the conserved and context-specific signaling pathways:
# signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv)
# increase, by change their information flow at one condition as compared to
# another condition.
# The top signaling pathways colored red are enriched in CRSD, and these colored green were enriched in PD
gg1 <- rankNet(CC_merged_object, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(CC_merged_object, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# Extracting data from the plot
gg2_data <- ggplot_build(gg2)
info_flow_data <- gg2_data$plot$data

# Constructing the plot
info_flow_data$log2_contribution_1 <- log2(info_flow_data$contribution+1)
#info_flow_data$group <- ifelse(info_flow_data$group == "Jak1_KO", "Jak1 KO", "WT")

barplot_colors <- c("WT" = "aquamarine4", "Jak1_KO" = "deeppink3")

hist(info_flow_data$contribution, breaks = 55, main = "", xlab = "Information flow")
abline(v=1, col = "red")

dev.off()

info_flow_data %>%
  ggplot(aes(x = reorder(name, -log2_contribution_1), y = log2_contribution_1, fill = group)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors) +
  theme_bw() + 
  xlab("") + 
  ylab(bquote(log[2]("Information flow + 1"))) +
  labs(fill = "Jak1 status") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9))

barplot_colors <- c("aquamarine3", "deeppink3")

barp <- info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = Jak1_KO-WT) %>%
  dplyr::select(name, delta) %>%
  #mutate(log2_delta = log2(delta)) %>%
  mutate(sign = delta > 0) %>%
  ggplot(aes(x = reorder(name, -delta), y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("CR/SD", "PD")) +
  theme_bw() + 
  xlab("") + 
  #ylab(bquote(log[2]("Information flow + 1"))) +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "Jak1 status") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9))

filename <- "/home/hnatri/CART/13384_Tumors/Plots/JAK1KO_WT_CellChat_pathways.pdf"
pdf(file = filename,
    width = 4,
    height = 2.5)

barp

dev.off()

# Heatmaps
Heatmap <- ComplexHeatmap::Heatmap

gg1 <- netVisual_heatmap(CC_merged_object,
                         cluster.rows = T,
                         cluster.cols = T)
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         #color.use = jak1_celltype_col)
gg2 <- netVisual_heatmap(CC_merged_object,
                         measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T)
                         #color.use = immune_fibro_celltype_col[-which(names(immune_fibro_celltype_col) %in% c("M8"))])
                         #color.use = jak1_celltype_col)
gg1 + gg2

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_response_heatmaps.pdf"
pdf(file = filename,
    width = 8,
    height = 4)

gg1 + gg2

dev.off()