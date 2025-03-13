#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/23/2024
# Description: Gene Set Enrichment Analysis for 13384 tumor scRNAseq data,
# updated escape version
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
library(GSVA)
library(MatrixGenerics)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import data
#==============================================================================#

# The final object
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_immune_fibroblast_reclustered_DoubletFinder.rds")

top <- c("109", "265", "181", "301", "223", "141")
bottom <- c("275", "228", "237", "234", "224", "129", "248", "185", "146")

immune_fibro$SPP1_surv_extremes <- ifelse(immune_fibro$UPN %in% top, "top",
                                          ifelse(immune_fibro$UPN %in% bottom, "bottom", NA))

seurat_object <- immune_fibro

#==============================================================================#
# Modified GSVA functions
#==============================================================================#

.GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  }
  egc <- gene.sets
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

.cntEval <- function(obj, 
                     assay = "RNA", 
                     type = "counts") {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

.split_data.matrix <- function(matrix, chunk.size=1000) {
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1
  
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks-1) {  #make last two chunks of equal size
      left <- ncols-(i-1)*chunk.size
      max <- min+round(left/2)-1
    } else {
      max <- min(i*chunk.size, ncols)
    }
    split.data[[i]] <- matrix[,min:max]
    min <- max+1    #for next chunk
  }
  return(split.data)
}

.gsva.setup <- function(data, egc) {
  params.to.use <- gsvaParam(exprData = data,
                             geneSets = egc,
                             kcdf = "Poisson")
  return(params.to.use)
}

performNormalization <- function(sc.data,
                                 enrichment.data = NULL,
                                 assay = "escape",
                                 gene.sets = NULL,
                                 make.positive = FALSE,
                                 scale.factor = NULL,
                                 groups = NULL) {
    if(!is.null(assay)) {
      if(is_seurat_object(sc.data)) {
        assay.present <- assay %in% Assays(sc.data)
      } else if (is_se_object(sc.data)) {
        assay.present <- assay %in% assays(sc.data)
      }
    } else {
      assay.present <- FALSE
    }
    
    
    if(is_seurat_or_se_object(sc.data) & !is.null(assay) & assay.present) {
      enriched <- .pull.Enrich(sc.data, assay)
    } else {
      enriched <- enrichment.data
    }
    
    if(!is.null(scale.factor) & length(scale.factor) != dim(sc.data)[2]) {
      stop("If using a vector as a scale factor, please ensure the length matches the number of cells.")
    }
    
    #Getting the gene sets that passed filters
    egc <- .GS.check(gene.sets)
    names(egc) <- str_replace_all(names(egc), "_", "-")
    egc <- egc[names(egc) %in% colnames(enriched)]
    
    #Isolating the number of genes per cell expressed
    if(is.null(groups)){
      chunks <- dim(enriched)[[1]]
    }else{
      chunks <- min(groups, dim(enriched)[[1]])
    }
    
    if (is.null(scale.factor)) {
      cnts <- .cntEval(sc.data, assay = "RNA", type = "counts")
      print("Calculating features per cell...")
      egc.sizes <- lapply(egc, function(x){
        scales<-unname(Matrix::colSums(cnts[which(rownames(cnts) %in% x),]!=0))
        scales[scales==0] <- 1
        scales
      })
      egc.sizes <- split_rows(do.call(cbind,egc.sizes), chunk.size=chunks)
      rm(cnts)
    }else{
      egc.sizes <- split_vector(scale.factor, chunk.size=chunks)
    }
    enriched <- split_rows(enriched, chunk.size=chunks)
    
    print("Normalizing enrichment scores per cell...")
    #Dividing the enrichment score by number of genes expressed
    
    enriched <- mapply(function(scores, scales){
      scores/scales
    }, enriched, egc.sizes, SIMPLIFY = FALSE)
    
    enriched <- do.call(rbind, enriched)
    if(make.positive){
      enriched <- apply(enriched, 2, function(x){
        x+max(0, -min(x))
      })
    }
    if(is_seurat_or_se_object(sc.data)) {
      if(is.null(assay)) {
        assay <- "escape"
      }
      sc.data <- .adding.Enrich(sc.data, enriched, paste0(assay, "_normalized"))
      return(sc.data)
    } else {
      return(enriched)
  }
}

escape.matrix <- function(input.data, 
                          gene.sets = NULL, 
                          method = "ssGSEA", 
                          groups = 1000, 
                          min.size = 5,
                          normalize = FALSE,
                          make.positive = FALSE,
                          BPPARAM = SerialParam(),
                          ...) {
  egc <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  egc.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
  if (!is.null(min.size)){
    remove <- unname(which(egc.size < min.size))
    if(length(remove) > 0) {
      egc <- egc[-remove]
      egc.size <- egc.size[-remove]
      if(length(egc) == 0) {
        stop("No gene sets passed the minimum length - please reconsider the 'min.size' parameter")
      }
    }
  }
  
  scores <- list()
  splits <- seq(1, ncol(cnts), by=groups)
  print(paste('Using sets of', groups, 'cells. Running', 
              length(splits), 'times.'))
  split.data <- .split_data.matrix(matrix=cnts, chunk.size=groups)
  
  for (i in seq_along(splits)) {
    last <- min(ncol(cnts), i+groups-1)
    if(method == "GSVA") {
      parameters <- .gsva.setup(split.data[[i]], egc)
    } else if (method == "ssGSEA") {
      parameters <- .ssGSEA.setup(split.data[[i]], egc)
    }
    if(method %in% c("ssGSEA", "GSVA")) {
      a <- suppressWarnings(gsva(param = parameters, 
                                 verbose = FALSE,
                                 BPPARAM = BPPARAM))
    }
    scores[[i]] <- a
  }
  
  # Unable to calculate scores for some pathways for some sets
  complete_pathways <- lapply(scores, function(xx){
    rownames(xx)
  })
  
  complete_pathways <- Reduce(intersect, complete_pathways)
  length(complete_pathways)
  
  complete_GS_CANONICAL <- GS_CANONICAL[complete_pathways]
  
  scores_complete <- lapply(scores, function(xx){
    xx[complete_pathways,]
  })
  
  scores_complete <- do.call(cbind, scores_complete)
  output <- t(as.matrix(scores_complete))
  
  gene.sets_complete <- gene.sets[scores_complete]
  
  #Normalize based on dropout
  #if(normalize) {
  #  output <- performNormalization(sc.data = input.data,
  #                                 enrichment.data = output,
  #                                 #assay = NULL,
  #                                 gene.sets = scores_complete,
  #                                 make.positive = make.positive,
  #                                 groups = groups)
  #}
  return(output)
}

runEscape <- function(input.data, 
                      gene.sets = NULL, 
                      method = "ssGSEA", 
                      groups = 1000, 
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      BPPARAM = BiocParallel::SerialParam(),
                      ...) {
  #.checkSingleObject(input.data)
  enrichment <- escape.matrix(input.data = input.data,
                              gene.sets = gene.sets,
                              method = method,
                              groups = groups,
                              min.size = min.size,
                              BPPARAM = BPPARAM)
  
  input.data <- .adding.Enrich(input.data, enrichment, new.assay.name)
  return(input.data)
}

#==============================================================================#
# Run GSVA
#==============================================================================#

# C2 = curated gene sets,
# H = Hallmark
GS <- getGeneSets(species = "Homo sapiens", library = c("C2", "H"))
GS_CANONICAL <- GS[grep("KEGG|REACTOME|BIOCARTA|HALLMARK", names(GS),
                        ignore.case = TRUE)]

#seurat_object <- subset(seurat_object, subset = celltype %in% c("M1", "M2", "M3"))
escape_res <- runEscape(input.data = seurat_object,
                           method = "GSVA", 
                           new.assay.name = "escapeGSVA",
                           gene.sets = GS_CANONICAL, 
                           groups = 500,
                           min.size = 5)

escape_res <- performNormalization(sc.data = escape_res, 
                                   assay = "escapeGSVA", 
                                   gene.sets = complete_GS_CANONICAL, 
                                   scale.factor = escape_res$nFeature_RNA,
                                   groups = 500)

saveRDS(escape_res, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_C2_H_updated.rds")
#q(save="no")

#==============================================================================#
# Plot results
#==============================================================================#

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_GSEA_C2_H_updated.rds")

complete_pathways <- gsub("_", "-", complete_pathways)

# Plotting
rownames(tumors@assays$escapeGSVA_normalized)
VlnPlot(tumors,
        group.by = "SPP1_surv_extremes",
        features = c("BIOCARTA-MSP-PATHWAY"),
        assay = "escapeGSVA_normalized",
        pt.size = 0)

# Differential enrichment
tumors_subset <- subset(tumors, subset = SPP1_surv_extremes %in% c("top", "bottom"))
myeloid <- subset(tumors_subset, subset = celltype %in% c("N1", paste0("M", seq(1, 9))))

Idents(myeloid) <- myeloid$SPP1_surv_extremes
diff_enrich <- FindAllMarkers(myeloid,
                              assay = "escapeGSVA_normalized", 
                              min.pct = 0.1,
                              logfc.threshold = 2)

head(diff_enrich)
dim(diff_enrich)
diff_enrich_pathways <- unique(gsub("\\.1", "", rownames(diff_enrich)))

# Geyser plot
geyserEnrichment(myeloid, 
                 assay = "escapeGSVA_normalized",
                 gene.set = diff_enrich_pathways)

# Heatmap of celltype specific enriched pathways
Idents(myeloid) <- myeloid$celltype
celltype_enrich <- FindAllMarkers(myeloid,
                                  assay = "escapeGSVA_normalized", 
                                  min.pct = 0.1,
                                  logfc.threshold = 2)

heatmapEnrichment(myeloid, 
                  group.by = "celltype",
                  assay = "escapeGSVA_normalized",
                  gene.set.use = rownames(myeloid@assays$escapeGSVA_normalized@data)[1:12],
                  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)

# TODO
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

####

selected_pathways_subset3 <- selected_pathways %>% filter(Type %in% c("Phagocytosis",
                                                                      "Hypoxia"))
delta_plot_3 <- output_tumors %>%
  filter(pathways %in% selected_pathways_subset3$pathways) %>% # other_pathways
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

delta_plot_3

pdf(file = "/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/GSEA_SPP1highlow_phagocytosis_hypoxia.pdf",
    width=7, height=2)

delta_plot_3

dev.off()
