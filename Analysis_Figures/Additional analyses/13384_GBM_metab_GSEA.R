#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 10/31/2024
# Description: Parsing and plotting HGG tumor escape/GSEA results
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
library(survminer) # for theme_classic2()

#==============================================================================#
# Helper functions
#==============================================================================#

# getSignificance function from the previous version of escape
getSignificance <- function(enriched, 
                            group = NULL,
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
  new.grouping <- unlist(med[,1])
  med <- as.data.frame(med[,seq_len(num) + 1])
  rownames(med) <- paste0("median.", new.grouping)
  med <- as.data.frame(t(med))
  return(med)
}

#==============================================================================#
# Import data
#==============================================================================#

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_C2_H.rds")
res_tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_res_C2_H.rds")
tumors <- subset(tumors, subset = celltype %in% paste0("Tumor", seq(1,7)))
tumors_res <- tumors@meta.data %>%
  dplyr::select(c("orig.ident", "SPP1_surv_extremes", colnames(res_tumors)))
tumors_res <- tumors_res %>% filter(SPP1_surv_extremes %in% c("top", "bottom"))

#output_tumors <- read.table("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1top_bottom_GSEA_sig.tsv",
#                            header = T)

#==============================================================================#
# Get significant pathways
#==============================================================================#

output_tumors <- getSignificance(tumors_res, group = "SPP1_surv_extremes", fit = "Wilcoxon")
output_tumors <- output_tumors %>% filter(FDR < 0.01) %>% rownames_to_column(var = "pathways") %>% arrange(FDR)

write.table(output_tumors, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/SPP1top_bottom_GSEA_sig_tumor.tsv",
            quote = F, sep = "\t", row.names = F)

tumors_res2 <- tumors@meta.data %>%
  dplyr::select(c("orig.ident", "binary_response", colnames(res_tumors)))
tumors_res2 <- tumors_res2 %>% filter(binary_response %in% c("CR_SD", "PD"))
output_tumors2 <- getSignificance(tumors_res2, group = "binary_response", fit = "Wilcoxon")
output_tumors2 <- output_tumors2 %>% filter(FDR < 0.01) %>% rownames_to_column(var = "pathways") %>% arrange(FDR)

write.table(output_tumors2, "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/CSRD_PD_GSEA_sig_tumor.tsv",
            quote = F, sep = "\t", row.names = F)

#==============================================================================#
# Plot
#==============================================================================#

# Barplot of selected pathways for SPP1 high vs. low
path1 <- grep("FATTY_ACID_METABOLISM", colnames(res_tumors), value = T)
path2 <- grep("_SPHINGOLIPID_METABOLISM", colnames(res_tumors), value = T)
path3 <- grep("CERAMIDE", colnames(res_tumors), value = T)
path4 <- grep("_PHOSPHATID", colnames(res_tumors), value = T)

paths <- c(path1, path2, path3, path4)

# Results to a df with cols median.top, median.bottom, and pathways
tumors_res2 <- tumors_res %>%
  pivot_longer(cols = colnames(res_tumors),
               values_to = "val",
               names_to = "pathways") %>%
  filter(pathways %in% paths) %>%
  group_by(SPP1_surv_extremes, pathways) %>%
  dplyr::summarize(median = median(val)) %>%
  ungroup() %>%
  pivot_wider(names_from = SPP1_surv_extremes,
              values_from = median,
              names_prefix = "median.")

delta_plot <- tumors_res2 %>%
  filter(pathways %in% paths) %>% # other_pathways
  mutate(delta = median.bottom - median.top,
         sign = sign(delta),
         signstr = if_else(sign == 1, "SPP1high", "SPP1low")) %>%
  ggplot(aes(x = delta, y = reorder(pathways, delta), fill = signstr)) +
  geom_bar(stat = "identity") +
  #geom_col(width = 0.85) +
  scale_fill_manual(values = c("orangered1", "royalblue3")) +
  theme_classic() +
  #manuscript_theme +
  ylab("") +
  xlab(expression(Delta ~ "median enrichment score")) +
  #coord_flip() +
  #facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
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


# Genes
mm_genes <- c("Cept1", "Cers2", "Chka", "Decr1", "Elovl5", "Kdsr", "Lipc",
              "Mboat2", "Mboat7", "Pcyt1a", "Pnpla8", "Smpd2", "Sptlc2", "Tecr")

# Converting mouse gene names to human
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")

convert_mouse_to_human <- function(gene_list){
  gene_names <- as.data.frame(matrix(nrow = length(gene_list),
                                     ncol = 2))
  colnames(gene_names) <- c("mouse", "human")
  rownames(gene_names) <- gene_list
  gene_names$mouse <- gene_list
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      
      if(length(human_genes)==0){
        gene_names[gene, "human"] <- NA
      } else if (length(human_genes)>1){
        #  human_genes <- paste0(human_genes, collapse = ", ")
        bind_df <- data.frame("mouse" = rep(gene, times = length(human_genes)),
                              "human" = human_genes)
        gene_names <- rbind(gene_names, bind_df)
      } else {
        gene_names[gene, "human"] <- human_genes
      }
    }
  }
  return(gene_names)
}

gene_names <- convert_mouse_to_human(mm_genes)

VlnPlot(tumors,
        features = gene_names$human,
        group.by = "SPP1_surv_extremes",
        #split.by = "SPP1_surv_extremes",
        #split.plot = T,
        pt.size = 0)
