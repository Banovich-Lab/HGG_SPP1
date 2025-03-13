#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/10/2024
# Description: Plotting information flow
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
library(ggbreak)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/Analysis_Figures/CellChat_heatmap.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Importing CellChat objects
#==============================================================================#

CC_CD3_objects_compare <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_compare_object_CD3_celltypenames.rds")
CC_CD3_merged_object <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_merged_object_CD3_celltypenames.rds")

CC_response_objects_compare <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_compare_object_response_celltypenames.rds")
CC_response_merged_object <- readRDS("/scratch/hnatri/CART/immune_stroma_CC_merged_object_response_celltypenames.rds")

CC_tumor_response_objects_compare <- readRDS("/scratch/hnatri/CART/tumor_CC_compare_object_response.rds")
CC_tumor_response_merged_object <- readRDS("/scratch/hnatri/CART/tumor_CC_merged_object_response.rds")

mm_CC_merged_object <- readRDS("/scratch/hnatri/CART/mm_CC_merged_object.rds")

CC_response_objects_compare <- CC_tumor_response_objects_compare
CC_response_merged_object <- CC_tumor_response_merged_object

# Identify and visualize the conserved and context-specific signaling pathways:
# signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv)
# increase, by change their information flow at one condition as compared to
# another condition.
# The top signaling pathways colored red are enriched in CRSD, and these colored green were enriched in PD

# CR/SD vs. PD

# Extracting data from the plot
response <- rankNet(CC_response_merged_object,
                    mode = "comparison",
                    stacked = F,
                    do.stat = T,
                    cutoff.pvalue = 0.01,
                    #tol = 1,
                    thresh = 0.01)

response_data <- ggplot_build(response)
response_info_flow_data <- response_data$plot$data

# Constructing the plot
response_info_flow_data$log2_contribution_1 <- log2(response_info_flow_data$contribution+1)
response_info_flow_data$group <- ifelse(response_info_flow_data$group == "CR_SD_CC", "CR_SD", "PD")

barplot_colors <- c("CR_SD" = "aquamarine4", "PD" = "deeppink3")

hist(response_info_flow_data$contribution, breaks = 55, main = "", xlab = "Information flow")
abline(v=1, col = "red")

dev.off()

# Plotting delta
barplot_colors <- c("aquamarine3", "deeppink3")

# MK, VEGF, PlAU, GRN, LAIR1, PARs
response_delta_plot <- response_info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = PD-CR_SD) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0) %>% # summarise(mind = min(delta))
  filter(abs(delta) > 0.8) %>%
  ggplot(aes(x = reorder(name, -delta), y = delta, fill = sign)) +
    geom_bar(position = "dodge", stat="identity") +
    scale_fill_manual(values = barplot_colors, labels=c("CR_SD", "PD")) +
    theme_bw() + 
    xlab("") +
    ylab(expression(Delta ~ "information flow")) +
    labs(fill = "Response") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
    manuscript_theme
    #scale_y_break(c(-15, -45))

response_delta_plot

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_tumor_response_delta_infoflow.pdf"
pdf(file = filename,
    width = 5,
    height = 2)

response_delta_plot

dev.off()

# Adding CD3 high vs. low

# Extracting data from the plot
cd3 <- rankNet(CC_CD3_merged_object,
               mode = "comparison",
               stacked = F,
               do.stat = T,
               cutoff.pvalue = 0.01,
               #tol = 1,
               thresh = 0.01)

cd3_data <- ggplot_build(cd3)
cd3_info_flow_data <- cd3_data$plot$data

# Constructing the plot
cd3_info_flow_data$log2_contribution_1 <- log2(cd3_info_flow_data$contribution+1)
cd3_info_flow_data$group <- ifelse(cd3_info_flow_data$group == "CD3low_CC", "CD3_low", "CD3_high")

barplot_colors <- c("CD3_high" = "aquamarine4", "CD3_low" = "deeppink3")

hist(cd3_info_flow_data$contribution, breaks = 55, main = "", xlab = "Information flow")
abline(v=1, col = "red")

dev.off()

# Extracting and sorting data to match the GBM response result
plot_pathways <- cd3_info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = PD-CR_SD) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0) %>%
  filter(abs(delta) > 0.8) %>%
  arrange(-delta)

cd3_info_flowplot <- cd3_info_flow_data %>%
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = CD3_low-CD3_high) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0)

# Lock in factor level order
cd3_info_flowplot <- cd3_info_flowplot[match(plot_pathways$name, cd3_info_flowplot$name),]
cd3_info_flowplot$name <- factor(cd3_info_flowplot$name, levels = cd3_info_flowplot$name)

unique(cd3_info_flowplot$sign)
cd3_info_flowplot %>% dplyr::filter(is.na(sign))

# Plotting delta
barplot_colors <- c("aquamarine3", "deeppink3")

cd3_delta_plot <- cd3_info_flowplot %>% 
  filter(!is.na(name)) %>%
  ggplot(aes(x = name, y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("CD3 high", "CD3 low")) + # labels=c("CR/SD", "PD")
  theme_bw() + 
  xlab("") + 
  #ylab(bquote(log[2]("Information flow + 1"))) +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "CD3 score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
  manuscript_theme

cd3_delta_plot

# Adding mouse data

# Extracting data from the plot
mm <- rankNet(mm_CC_merged_object,
              mode = "comparison",
              stacked = F,
              do.stat = T,
              cutoff.pvalue = 0.01,
              #tol = 1,
              thresh = 0.01)
mm_data <- ggplot_build(mm)
mm_data_info_flow <- mm_data$plot$data

mm_data_info_flow$log2_contribution_1 <- log2(mm_data_info_flow$contribution+1)

mm_data_info_flowplot <- mm_data_info_flow %>%
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = Jak1_KO-WT) %>%
  dplyr::select(name, delta) %>%
  #mutate(log2_delta = log2(delta)) %>%
  mutate(sign = delta > 0)

# Lock in factor level order
mm_data_info_flowplot <- mm_data_info_flowplot[match(plot_pathways$name, mm_data_info_flowplot$name),]
mm_data_info_flowplot$name <- factor(mm_data_info_flowplot$name, levels = mm_data_info_flowplot$name)

unique(mm_data_info_flowplot$sign)
mm_data_info_flowplot %>% dplyr::filter(is.na(sign))

mm_delta_plot <- mm_data_info_flowplot %>% 
  filter(!is.na(name)) %>%
  ggplot(aes(x = name, y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("WT", "Jak1 KO")) + # labels=c("CR/SD", "PD")
  theme_bw() + 
  xlab("") + 
  #ylab(bquote(log[2]("Information flow + 1"))) +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "Jak1 status") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
  manuscript_theme

mm_delta_plot

plot_list <- list(response_delta_plot, cd3_delta_plot, mm_delta_plot)

p1 <- patchwork::wrap_plots(plot_list, ncol = 1)
p2 <- response_delta_plot / cd3_delta_plot / mm_delta_plot


filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_response_CD3_JAK1KOWT_delta_infoflow.pdf"
pdf(file = filename,
    width = 4,
    height = 5)

patchwork::wrap_plots(plot_list, ncol = 1)

dev.off()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_response_delta_infoflow.pdf"
pdf(file = filename,
    width = 4,
    height = 2)

response_delta_plot

dev.off()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_CD3_delta_infoflow.pdf"
pdf(file = filename,
    width = 4,
    height = 2)

cd3_delta_plot

dev.off()

filename <- "/home/hnatri/CART/13384_Tumors/Plots/CellChat_JAK1KOWT_delta_infoflow.pdf"
pdf(file = filename,
    width = 4,
    height = 2)

mm_delta_plot

dev.off()

# All top pathways for CD3h vs l
cd3_info_flowplot <- cd3_info_flow_data %>%
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = CD3_low-CD3_high) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0)

# Lock in factor level order
#cd3_info_flowplot <- cd3_info_flowplot[match(plot_pathways$name, cd3_info_flowplot$name),]
cd3_info_flowplot$name <- factor(cd3_info_flowplot$name, levels = cd3_info_flowplot$name)

unique(cd3_info_flowplot$sign)
cd3_info_flowplot %>% dplyr::filter(is.na(sign))

# Plotting delta
barplot_colors <- c("aquamarine3", "deeppink3")

# Plotting delta
barplot_colors <- c("aquamarine3", "deeppink3")

cd3_delta_plot <- cd3_info_flow_data %>%
  #mutate_at(c("contribution"), function(x) ifelse(is.infinite(x), 0, x)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(any(contribution>1)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(any(pvalues<0.01)) %>%
  ungroup() %>% 
  dplyr::select(name, contribution, group) %>%
  tidyr::pivot_wider(names_from = "group", values_from = "contribution") %>%
  mutate(delta = CD3_low-CD3_high) %>%
  dplyr::select(name, delta) %>%
  mutate(sign = delta > 0) %>% # summarise(mind = min(delta))
  filter(abs(delta) > 0.8) %>%
  ggplot(aes(x = reorder(name, -delta), y = delta, fill = sign)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = barplot_colors, labels=c("CD3_high", "CD3_low")) +
  theme_bw() + 
  xlab("") +
  ylab(expression(Delta ~ "information flow")) +
  labs(fill = "Response") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.9)) +
  manuscript_theme
#scale_y_break(c(-15, -45))

cd3_delta_plot

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/CellChat_tumor_CD3_allpathways_delta_infoflow.pdf"
pdf(file = filename,
    width = 5,
    height = 2)

cd3_delta_plot

dev.off()

#####

## Identify the upregulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabilities

plot_stroma <- paste0("F", seq(1, 3))
plot_myeloid <- c(paste0("M", seq(1, 9)), "N1")
plot_lymphoid <- paste0("L", seq(1, 9))

# Filtering based on pathway name
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
all_interactions <- CellChatDB$interaction
interactions_keep <- filter(all_interactions, grepl(paste(c("SPP1",
                                                            "COLLAGEN",
                                                            "LAMININ",
                                                            "VISFANTIN"), collapse="|"), pathway_name))

all_interactions

# Filtering based on a ligand/receptor
# ligand.symbol
interactions_keep <- filter(all_interactions, grepl(paste(c("SPP1", "APOE",
                                                            "C1QA", "C1QB",
                                                            "C1QC", "COL1A1",
                                                            "S100A8", "S100A9"), collapse="|"), ligand.symbol))

pairlr_use <- data.frame("interaction_name" = unique(interactions_keep$interaction_name))

CC_response_merged_object@LR$CR_SD_CC$LRsig

# Upregulated/downregulated in PD
gg1 <- netVisual_bubble(CC_response_merged_object,
                        sources.use = plot_myeloid,
                        targets.use = plot_stroma,
                        comparison = c(1, 2),
                        max.dataset = 2,
                        thresh = 0.01,
                        pairLR.use = pairlr_use, # a df with one column of either "interaction_name" or "pathway_name"
                        title.name = paste0("Up-regulated signaling in ", names(CC_response_objects_compare)[2]),
                        angle.x = 45,
                        remove.isolate = T)
gg2 <- netVisual_bubble(CC_response_merged_object,
                        sources.use = plot_myeloid,
                        targets.use = plot_stroma,
                        comparison = c(1, 2),
                        max.dataset = 1,
                        thresh = 0.01,
                        pairLR.use = pairlr_use,
                        title.name = paste0("Down-regulated signaling in ", names(CC_response_objects_compare)[2]),
                        angle.x = 45,
                        remove.isolate = T)
gg1 + gg2

filename <- "/home/hnatri/CART/13384_Tumors/SPP1_ms/CellChat_bubble_SPP1.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 10,
    height = 2.5)

gg1

dev.off()

gg1 <- netVisual_bubble(CC_response_merged_object,
                        sources.use = plot_lymphoid,
                        targets.use = plot_stroma,
                        comparison = c(1, 2),
                        max.dataset = 2,
                        title.name = paste0("Up-regulated signaling in ", names(CC_response_objects_compare)[2]),
                        angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(CC_response_merged_object,
                        sources.use = plot_lymphoid,
                        targets.use = plot_stroma,
                        comparison = c(1, 2),
                        max.dataset = 1,
                        title.name = paste0("Down-regulated signaling in ", names(CC_response_objects_compare)[2]),
                        angle.x = 45, remove.isolate = T)
gg1 / gg2

gg1 + manuscript_theme + ylab("")

filename <- "/home/hnatri/CART/13384_Tumors/SPP1_ms/CellChat_bubble_SPP1.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 8.5,
    height = 2.2)

gg1 + manuscript_theme + ylab("")

dev.off()

# Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(CC_objects_compare, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC_objects_compare)[i]))
}

# Heatmap based on a single object 
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  ht[[i]] <- netVisual_heatmap(CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ",
                                                  names(CC_objects_compare)[i]))
  #color.use = tumor_celltype_col)
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  ht[[i]] <- netVisual_heatmap(CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               #cluster.rows = T,
                               #cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ", names(CC_objects_compare)[i]))
  #color.use = tumor_cluster_col_c[1:13])
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

tumor_celltype_col[-which(names(tumor_celltype_col) %in% seurat_object$celltype)]
unique(seurat_object$celltype)[-which(unique(seurat_object$celltype) %in% names(tumor_celltype_col))]

tumor_celltype_col <- tumor_celltype_col[!names(tumor_celltype_col) %in% c('Myel7','Tumor6')]
tumor_celltype_col <- tumor_celltype_col[unique(names(tumor_celltype_col))]

pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CC_objects_compare)) {
  ht[[i]] <- netVisual_heatmap(object = CC_objects_compare[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               cluster.rows = T,
                               cluster.cols = T,
                               title.name = paste(pathways.show, "signaling ", names(CC_objects_compare)[i]),
                               color.use = tumor_celltype_col)
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = tumor_celltype_col)
}

pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = tumor_celltype_col)
}

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_aggregate(CC_objects_compare[[i]],
                      signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show, names(CC_objects_compare)[i]),
                      color.use = tumor_celltype_col)
}

#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(CC_objects_compare[[1]]@idents)
pathways.show <- c("LAMININ") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  #group = group.cellType,
  netVisual_chord_cell(CC_objects_compare[[i]],
                       signaling = pathways.show,
                       title.name = paste0(pathways.show, " signaling network - ", names(CC_objects_compare)[i]),
                       color.use = tumor_celltype_col)
}

par(mfrow = c(1, 2), xpd=TRUE)

# compare all the interactions sending from one cell type to another
for (i in 1:length(CC_objects_compare)) {
  netVisual_chord_gene(CC_objects_compare[[i]],
                       sources.use = c("C9"),
                       targets.use = c("C3", "C6", "C7", "C10"),
                       lab.cex = 0.5,
                       #small.gap = 0.2,
                       #big.gap = 2,
                       title.name = paste0("Signaling from fibroblasts - ", names(CC_objects_compare)[i]))
}

# Compare all the interactions sending from fibroblast to inflammatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_chord_gene(CC_objects_compare[[i]],
                       sources.use = c("C9"),
                       targets.use = c("C3", "C6", "C7", "C10"),
                       title.name = paste0("Signaling received by myeloid cells - ", names(CC_objects_compare)[i]),
                       legend.pos.x = 10)
}

# Show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(CC_objects_compare)) {
  netVisual_chord_gene(CC_objects_compare[[i]],
                       sources.use = c("Fibroblasts", "Olig1", "Olig2", "Olig3"),
                       targets.use = c("Glia", "MDSC1", "MDSC2", "TAM1", "Neut", "Myel"),
                       slot.name = "netP",
                       #small.gap = 8,
                       big.gap = 30,
                       title.name = paste0("Signaling pathways from stroma to myeloid cells - ", names(CC_objects_compare)[i]),
                       legend.pos.x = 10)
}

# Compare the signaling gene expression distribution between different datasets
#CC_merged_object@idents <- factor(factor(CC_merged_object@meta$cluster_cc, levels = paste0("C", c(0, seq(1, 13, by=1)))), levels = paste0("C", c(0, seq(1, 13, by=1))))
CC_merged_object@meta$binary_response = factor(CC_merged_object@meta$binary_response, levels = c("CR_SD", "PD")) # set factor level
plotGeneExpression(CC_merged_object,
                   signaling = "SPP1",
                   split.by = "binary_response",
                   colors.ggplot = T)
#color.use = c("red", "blue"))

plot_features <- c("SPP1", "CD44", "ITGAV", "ITGA4", "ITGA5", "ITGB1", "ITGB5")

seurat_object$celltype <- factor(seurat_object$celltype,
                                 levels = c("Olig1", "Olig2", "Olig3", "Fibroblast",
                                            "Tumor1", "Tumor2", "Tumor3", "Tumor4", "Tumor5",
                                            "Teff", "NK", "Glia", "MDSC1", "MDSC2", "TAM1", "Neut", "Myel"))

VlnPlot(seurat_object,
        features = plot_features,
        group.by = "celltype",
        split.by = "binary_response",
        pt.size = 0,
        slot = "data") &
  theme_bw()

seurat_object_subset <- subset(seurat_object, subset = celltype %in% c("MDSC1", "MDSC2", "Microglia"))
p1 <- VlnPlot(seurat_object_subset,
              features = c("SPP1"),
              group.by = "celltype",
              split.by = "binary_response",
              pt.size = 0,
              slot = "data",
              ncol = 1,
              cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend()

seurat_object_fibrob <- subset(seurat_object, subset = celltype %in% c("Fibroblast"))
p2 <- VlnPlot(seurat_object_fibrob,
              features = c("CD44", "ITGAV"),
              group.by = "celltype",
              split.by = "binary_response",
              pt.size = 0,
              slot = "data",
              ncol = 2,
              cols = c("aquamarine3", "deeppink3")) &
  theme_bw() &
  xlab("") &
  ylab("Log-normalized expression") &
  NoLegend()

p1 / p2


