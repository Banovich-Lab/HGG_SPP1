#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/15/2024
# Description: Gene Set Enrichment Analysis for 13384 tumor myeloid cells
# with high/low SPP1/survival
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
library(UpSetR)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables and helper functions
#==============================================================================#

set.seed(1234)

getSignificance <- function(enriched, group = NULL,
                            gene.sets = NULL,
                            fit = NULL) {
  fit <- match.arg(fit,  choices = c("T.test", "ANOVA", "Wilcoxon", "LR", "KW"))
  group2 <- enriched[,group]
  gr_names <- unique(group2)
  if (!is.null(gene.sets)) {
    input <- enriched[,colnames(enriched) %in% gene.sets]
  } else {
    input <- select_if(enriched, is.numeric)
  }
  medians <- get.medians(input, group2)
  output <- NULL
  if (fit == "T.test" || fit == "Wilcoxon" || fit == "LR") {
    if (length(unique(group2)) != 2) {
      message("Ensure the group selection has only two levels for T.test 
                fit") 
    } else {
      if (fit == "T.test") {
        out <- lapply(input, function(x) t.test(x ~ group2))
        stat <- "T"
      } else if (fit == "Wilcoxon") {
        out <- lapply(input, function(x) wilcox.test(x ~ group2))
        stat <- "W"
      }  else if (fit == "LR") {
        group2 <- ifelse(group2 == gr_names[1], 0,1)
        out <- lapply(input, function(x) glm(group3 ~ x, family = "binomial"))
        out <- lapply(out, function(x) tidy(x)[2,])
        stat <- "L"
      }
      for (i in seq_along(out)) {
        df <- out[[i]]
        mat <- c(df$statistic, df$p.value)
        output <- rbind(output,mat)
      }
      output <- as.data.frame(output)
      colnames(output) <- c(paste0(stat, ".statistic"), "p.value")
    }
  } else if (fit == "ANOVA") {
    if (length(unique(group2)) <= 2) {
      message("Ensure the group selection has more than two levels 
                for ANOVA fit") }
    out <- lapply(input, function(x) aov(x ~ group2))
    for (i in seq_along(out)) {
      df <- out[[i]]
      tukey <- TukeyHSD(df)
      ind.p.values <- melt(tukey$group2[,4])
      names.ind.p.values <- gsub("-", "v", rownames(ind.p.values))
      names.ind.p.values <- paste0(names.ind.p.values,".p.value")
      fval <- summary(df)[[1]]$'F value'[[1]]
      pval <- summary(df)[[1]]$'Pr(>F)'[[1]]
      output <- rbind(output, c(fval, pval, t(ind.p.values)))
    }
    output <- as.data.frame(output)
    colnames(output) <- c("f.value", "p.value", names.ind.p.values)
  } else if (fit == "KW") {
    if (length(unique(group2)) <= 2) {
      message("Ensure the group selection has more than two levels 
                for Kruskal-Wallis test")}
    out <- lapply(input, function(x) kruskal.test(x ~ group2))
    out.ind <- lapply(input, function(x) pairwise.wilcox.test(x, group2, p.adjust.method = "BH"))
    for (i in seq_along(out)) {
      ind.p.values <- na.omit(melt(out.ind[[i]]$p.value))
      names.ind.p.values <- paste0(ind.p.values$Var1, "v", ind.p.values$Var2)
      names.ind.p.values <- paste0(names.ind.p.values,".p.value")
      ind.p.values <- ind.p.values[,3]
      Chi.squared <- out[[i]]$statistic 
      pval <- out[[i]]$p.value
      output <- rbind(output, c(Chi.squared, pval, t(ind.p.values)))
    }
    output <- as.data.frame(output)
    colnames(output) <- c("Chi.square", "p.value", names.ind.p.values)
  }
  rownames(output) <- colnames(input)
  output$FDR <- p.adjust(output$p.value) 
  output <- cbind.data.frame(output, medians)
  return(output)
}

get.medians<- function(input, group2) {
  input <- cbind.data.frame(group2, input)
  num <- ncol(input)-1
  med <- input %>%
    group_by(group2) %>%
    dplyr::summarise(across(seq_len(all_of(num)), median))
  med <- as.data.frame(med[,seq_len(num) + 1])
  rownames(med) <- paste0("median.", unique(group2))
  med <- as.data.frame(t(med))
  return(med)
}

#==============================================================================#
# Import data
#==============================================================================#

# The final object
tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

#myeloid <- subset(immune_fibro, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))
#myeloid$SPP1high <- ifelse(myeloid$celltype %in% c("M6"), "High",
#                           ifelse(myeloid$celltype %in% c("N1", "M7"), "Low", "Int"))

top <- c("109", "265", "181", "301", "223", "141")
bottom <- c("228", "237", "234", "224", "129", "248", "185", "146") #"275", 

tumors$SPP1_surv_extremes <- ifelse(tumors$UPN %in% top, "top",
                                    ifelse(tumors$UPN %in% bottom, "bottom", NA))

immune_fibro$SPP1_surv_extremes <- ifelse(immune_fibro$UPN %in% top, "top",
                                          ifelse(immune_fibro$UPN %in% bottom, "bottom", NA))

seurat_object <- tumors

#==============================================================================#
# Run GSEA
#==============================================================================#

# C2 = curated gene sets,
# H = Hallmark
GS <- getGeneSets(species = "Homo sapiens", library = c("C2", "H"))
GS_CANONICAL <- GS[grep("KEGG|REACTOME|BIOCARTA|HALLMARK", names(GS),
                        ignore.case = TRUE)]

res <- enrichIt(obj = seurat_object,
                gene.sets = GS_CANONICAL,
                groups = 50, cores = 4)
seurat_object <- AddMetaData(seurat_object, res)

#saveRDS(seurat_object, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_C2_H.rds")
#saveRDS(res, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_res_C2_H.rds")

#q(save="no")

#==============================================================================#
# Plot results
#==============================================================================#

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_immune_fibro_GSEA_C2_H.rds")
#res_tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_immune_fibro_GSEA_res_GSEA_C2_H.rds")

tumors$SPP1_surv_extremes <- ifelse(tumors$UPN %in% top, "top",
                                    ifelse(tumors$UPN %in% bottom, "bottom", NA))

#res_tumors$SPP1_surv_extremes <- ifelse(res_tumors$UPN %in% top, "top",
#                                        ifelse(res_tumors$UPN %in% bottom, "bottom", NA))

colnames(tumors@meta.data)

# Plotting
tumors_res2 <- tumors@meta.data %>%
  dplyr::select(c("orig.ident", "SPP1_surv_extremes", colnames(tumors@meta.data)[95:length(colnames(tumors@meta.data))]))
#colnames(myeloid_res2)[ncol(myeloid_res2)] <- "celltype"

#myeloid_res2$binary_response <- factor(myeloid_res2$binary_response, levels = c("CR_SD", "PD"))
tumors_res2 <- tumors_res2 %>% filter(SPP1_surv_extremes %in% c("top", "bottom"))
output_tumors <- data.frame(getSignificance(tumors_res2, group = "SPP1_surv_extremes", fit = "Wilcoxon"))
output_tumors$pathways <- rownames(output_tumors)
output_tumors <- output_tumors %>% filter(FDR < 0.01) %>% arrange(FDR)

head(output_tumors)

write.table(output_tumors, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1top_bottom_GSEA_sig_v2.tsv",
            quote = F, row.names = F, sep = "\t")

output1 <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1top_bottom_GSEA_sig.tsv",
                      sep = "\t", header = T)

setdiff(output1$pathways, output_tumors$pathways)

#output_tumors <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1top_bottom_GSEA_sig.tsv",
#                            header = T)
#
DefaultAssay(tumors) <- "RNA"

tumors_subset <- subset(tumors, subset = SPP1_surv_extremes %in% c("top", "bottom"))
myeloid <- subset(tumors, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))

# CR/SD vs. PD
tumors_res3 <- tumors@meta.data %>%
  dplyr::select(c("orig.ident", "binary_response", colnames(res_tumors)))

tumors_res3 <- tumors_res3 %>% filter(binary_response %in% c("CR_SD", "PD"))
output_CRSD_PD <- data.frame(getSignificance(tumors_res3, group = "binary_response", fit = "Wilcoxon"))
output_CRSD_PD$pathways <- rownames(output_CRSD_PD)
output_CRSD_PD <- output_CRSD_PD %>% filter(FDR < 0.01) %>% arrange(FDR)

write.table(output_CRSD_PD, "/home/hnatri/CART/13384_Tumors/SPP1_ms/CRSD_PD_GSEA_sig.tsv",
            quote = F, row.names = F, sep = "\t")

# CD3 high vs. low
tumors_res4 <- tumors@meta.data %>%
  dplyr::select(c("orig.ident", "CD3_high_low", colnames(res_tumors)))

tumors_res4 <- tumors_res4 %>% filter(CD3_high_low %in% c("High", "Low"))
output_CD3high_low <- data.frame(getSignificance(tumors_res4, group = "CD3_high_low", fit = "Wilcoxon"))
output_CD3high_low$pathways <- rownames(output_CD3high_low)
output_CD3high_low <- output_CD3high_low %>% filter(FDR < 0.01) %>% arrange(FDR)

write.table(output_CD3high_low, "/home/hnatri/CART/13384_Tumors/SPP1_ms/CD3high_low_GSEA_sig.tsv",
            quote = F, row.names = F, sep = "\t")

# Significant in all comparisons
input_list <- list("CD3high vs. low" = output_CD3high_low$pathways,
                   "CR/SD vs. PD" = output_CRSD_PD$pathways,
                   "SPP1high vs. SPP1low" = output_tumors$pathways)

upset(fromList(input_list), order.by = "freq")

pdf(file = "/home/hnatri/CART/13384_Tumors/SPP1_ms/GSEA_upset.pdf",
    width=5, height=2.5)

upset(fromList(input_list), order.by = "freq")

dev.off()

# Barplot of selected pathways for SPP1 high vs. low
head(output_tumors)

gs4_deauth()
selected_pathways  <- gs4_get("https://docs.google.com/spreadsheets/d/1jagg8T5KwYnx68PiruLfa4aOz67BG-dxFeA-pp5tZwQ/edit?usp=sharing")
sheet_names(selected_pathways)
selected_pathways <- read_sheet(selected_pathways, sheet = "Final selected paths")
head(selected_pathways)

setdiff(selected_pathways$pathways, output_tumors$pathways)
output_tumors$type <- mapvalues(x = output_tumors$pathways,
                                from = selected_pathways$pathways,
                                to = selected_pathways$Type)

delta_plot <- output_tumors %>%
  filter(pathways %in% selected_pathways$pathways) %>% # other_pathways
  mutate(delta = median.bottom - median.top,
         sign = sign(delta),
         signstr = if_else(sign == 1, "SPP1high", "SPP1low")) %>%
  ggplot(aes(x = delta, y = reorder(pathways, delta), fill = signstr)) +
    geom_bar(stat = "identity") +
    #geom_col(width = 0.85) +
    scale_fill_manual(values = c("orangered1", "royalblue3")) +
    theme_classic() +
    manuscript_theme +
    ylab("") +
    xlab(expression(Delta ~ "median enrichment score")) +
    #coord_flip() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
      theme(
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
        #plot.title = element_text(size = 15, face = "bold"),
        strip.text.x = element_text(angle = 270, face = "bold"),
        strip.placement = "outside",
        axis.title.y = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
        axis.title.x = element_blank(),
        #axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
      )

delta_plot

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_SPP1highlow_selected_paths_faceted.pdf",
    width=7, height=7)

delta_plot

dev.off()

# Antigen processing and ECM
unique(selected_pathways$Type)
selected_pathways_subset <- selected_pathways %>% filter(Type %in% c("ECM/Antigen processing"))

delta_plot_1 <- output_tumors %>%
  filter(pathways %in% selected_pathways_subset$pathways) %>% # other_pathways
  mutate(delta = median.bottom - median.top,
         sign = sign(delta),
         signstr = if_else(sign == 1, "SPP1high", "SPP1low")) %>%
  ggplot(aes(x = delta, y = reorder(pathways, delta), fill = signstr)) +
  geom_bar(stat = "identity") +
  #geom_col(width = 0.85) +
  scale_fill_manual(values = c("orangered1", "royalblue3")) +
  theme_classic() +
  manuscript_theme +
  ylab("") +
  xlab(expression(Delta ~ "median enrichment score")) +
  #coord_flip() +
  facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    #plot.title = element_text(size = 15, face = "bold"),
    strip.text.x = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.y = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
    axis.title.x = element_blank(),
    #axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
  )

delta_plot_1

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_SPP1highlow_antigenproc_ECM.pdf",
    width=7, height=4)

delta_plot_1

dev.off()

# Metabolism pathways
metab_pathways <- selected_pathways %>%
  filter(Type == "Metabolism")

metab_delta_plot <- output_tumors %>%
  filter(pathways %in% metab_pathways$pathways) %>%
  mutate(delta = median.bottom - median.top,
         sign = sign(delta),
         signstr = if_else(sign == 1, "SPP1high", "SPP1low")) %>%
  ggplot(aes(x = delta, y = reorder(pathways, delta), fill = signstr)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("orangered1", "royalblue3")) +
  theme_classic() +
  manuscript_theme +
  ylab("") +
  xlab(expression(Delta ~ "median enrichment score"))

metab_delta_plot

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_SPP1highlow_metabolism_paths.pdf",
    width=5.6, height=3)

metab_delta_plot

dev.off()

selected_pathways_subset2 <- selected_pathways %>% filter(Type %in% c("Interleukin pathways upregulated in SPP-1 high",
                                                                      "Interleukin",
                                                                      "Ligand-receptor"))

delta_plot_2 <- output_tumors %>%
  filter(pathways %in% selected_pathways_subset2$pathways) %>% # other_pathways
  mutate(delta = median.bottom - median.top,
         sign = sign(delta),
         signstr = if_else(sign == 1, "SPP1high", "SPP1low")) %>%
  ggplot(aes(x = delta, y = reorder(pathways, delta), fill = signstr)) +
  geom_bar(stat = "identity") +
  #geom_col(width = 0.85) +
  scale_fill_manual(values = c("orangered1", "royalblue3")) +
  theme_classic() +
  manuscript_theme +
  ylab("") +
  xlab(expression(Delta ~ "median enrichment score")) +
  #coord_flip() +
  facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    #plot.title = element_text(size = 15, face = "bold"),
    strip.text.x = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.y = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
    axis.title.x = element_blank(),
    #axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
  )

delta_plot_2

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_SPP1highlow_interleukins_ligandreceptor.pdf",
    width=5.2, height=5)

delta_plot_2

dev.off()
