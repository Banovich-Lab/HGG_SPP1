#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2022/02/22
# Description: Integrating tumor samples
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)
library(gridExtra)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/CART/CART_colors_themes.R")
source("/home/hnatri/CART/basic_integration.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Integrate RNA for all batches
#==============================================================================#

# Processing was done according to CART_Product_leuk_CSF_Tumor_preprocess.R

# The batch_list_filtered object is a list containing all 19 batches
# The batches have been QC filtered, subset to contain only leuk PBMC and Product, 
# RNA has been normalized and scaled using SCTransform
# Protein has been normalized and scaled using NormalizeData(normalization.method = "CLR") and ScaleData()
# "/scratch/hnatri/CART/batch_list_filtered_productonly_210812.rds"

# This object is in the lab storage at /labs/banovich/BCTCSF/Heini
batch_list_filtered <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/batch_list_filtered_all_sampletypes_220707.rds")
# Adding newly sequenced samples
new_tumor_csf_samples <- readRDS("/labs/banovich/BCTCSF/Heini/tumor_csf_list_filtered_220707.rds")
names(new_tumor_csf_samples)
new_tumor_csf <- lapply(names(new_tumor_csf_samples), function(xx){
    seurat_data <- new_tumor_csf_samples[[xx]]
    seurat_data$FID_GEXFB <- gsub("_filtered", "", xx)
    
    seurat_data
})

names(new_tumor_csf) <- names(new_tumor_csf_samples)
length(new_tumor_csf)

# All tumor samples in the list are 13384
new_tumor_csf <- lapply(new_tumor_csf, function(xx){
    if(unique(xx$Sample_Type) %in% c("FDT", "Tumor")){
        return(xx)
    } else {
        return(NULL)
    }
})
new_tumor_csf <- new_tumor_csf %>% purrr::discard(is.null)
length(new_tumor_csf) # 27 + 10 previously sequenced

# Samples sequenced in June 2023
june_tumor_list <- readRDS("/labs/banovich/BCTCSF/Heini/new_tumors_filtered_230609.rds")

# Samples sequenced in August 2023
aug_tumor_list <- readRDS("/labs/banovich/BCTCSF/Heini/new_tumors_filtered_230810.rds")

# Renaming cells
batch_list_filtered <- lapply(batch_list_filtered, function(x){
  renamed.assay <- RenameCells(x, new.names = paste0(unique(x$Batch), "_", colnames(x)))
  renamed.assay@meta.data$IRB <- "13384"
  
  renamed.assay
})

new_tumor_csf <- lapply(new_tumor_csf, function(x){
    renamed.assay <- RenameCells(x,
                                 new.names = paste0(unique(x$FID_GEXFB), "_", colnames(x)))
    
    renamed.assay
})

june_tumor_list <- lapply(june_tumor_list, function(x){
  renamed.assay <- RenameCells(x,
                               new.names = paste0(unique(x$Batch), "_", colnames(x)))
  renamed.assay@meta.data$IRB <- "13384"
  
  renamed.assay
})

aug_tumor_list <- lapply(aug_tumor_list, function(x){
  renamed.assay <- RenameCells(x,
                               new.names = paste0(unique(x$Batch), "_", colnames(x)))
  renamed.assay@meta.data$IRB <- "13384"
  
  renamed.assay
})

# Combining all samples
batch_list_filtered <- c(batch_list_filtered, new_tumor_csf, june_tumor_list, aug_tumor_list)

batch_list_filtered <- lapply(batch_list_filtered, function(xx){
    xx@meta.data$Sample_Type <- gsub("Tumor", "FDT", xx@meta.data$Sample_Type)
    #xx@meta.data$IRB[is.na(xx@meta.data$IRB)] <- "13384"
    
    xx
})

# Tumor only, all batches
# Some samples had too few cell cycle genes expressed, can't use scores in
# scaling
# Samples without cell cycle scores are all 13384
# F04739 = Tumor
# F04742 = Tumor
# F04719 = Tumor
# F04720 = Tumor
# F04714 = TF

tumor_list <- list()
for (i in 1:length(batch_list_filtered)) {
  if (!("FDT" %in% batch_list_filtered[[i]]@meta.data$Sample_Type)){
    next
  } else {
    message(names(batch_list_filtered)[i])
    #Idents(batch_list_filtered[[i]]) = batch_list_filtered[[i]]$Manufacture
    tumor_list[[i]] = subset(batch_list_filtered[[i]], subset = Sample_Type %in% c("FDT", "Tumor"))
    # rerun SCTransform because we're working with a subset
    DefaultAssay(tumor_list[[i]]) = "RNA"
    if(names(batch_list_filtered)[i] %in% c("F04739_filtered",
                                            "F04742_filtered",
                                            "F04719_filtered",
                                            "F04720_filtered",
                                            "F04714_filtered",
                                            c(names(june_tumor_list)))){
        tumor_list[[i]] = SCTransform(tumor_list[[i]],
                                      method = "glmGamPoi",
                                      vars.to.regress = c("percent.mt"),
                                      verbose = F)
    } else{
        tumor_list[[i]] = SCTransform(tumor_list[[i]],
                                     method = "glmGamPoi",
                                     vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                     verbose = F)
    }
}
}

names(tumor_list) <- names(batch_list_filtered)
lapply(tumor_list, function(xx){unique(xx$UPN)})
unique(tumor_list[[30]]$UPN)

# Removing UPN 142
#tumor_list[[30]] <- NULL
tumor_list <- tumor_list %>% purrr::discard(is.null)

#saveRDS(tumor_list, "/scratch/hnatri/CART/tumor_list.rds")
#q(save="no")

tumor_list <- readRDS("/scratch/hnatri/CART/tumor_list.rds")

# Removing UPN 142
# F04774
#tumor_list[["F04774_filtered"]] <- NULL
#tumor_list <- tumor_list %>% purrr::discard(is.null)
#message(names(tumor_list))

# Adding UPN 109 samples
upn109_samples <- list("28360" = "/tgen_labs/banovich/BCTCSF/UPN109/From_CoH/Sherri_11162022/190405_190425_combined_scRNA_2/28360_count",
                       "28361" = "/tgen_labs/banovich/BCTCSF/UPN109/From_CoH/Sherri_11162022/190405_190425_combined_scRNA_2/28361_count")
                       #"28362" = "/labs/banovich/BCTCSF/UPN109/From_CoH/Sherri_11162022/190405_190425_combined_scRNA_2/28362_count")

upn109_data_list <- lapply(upn109_samples, function(xx){
    s10x = Read10X(paste0(xx, "/outs/filtered_feature_bc_matrix/"))
    s10x = CreateSeuratObject(s10x)
    s10x = PercentageFeatureSet(s10x, pattern = "^MT-", col.name = "percent.mt")
    
    s10x
})

# Filtering UPN109 data
bt_merge <- merge(x = upn109_data_list[[1]], y = upn109_data_list[2:length(upn109_data_list)])
#bt_merge <- PercentageFeatureSet(object = bt_merge, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(bt_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Make plots
# subset = nFeature_RNA > 500 & nCount_RNA > 5000 & percent.mt < 10 & nCount_Protein < 10000)
smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nFeature_RNA)
abline(h = 500, v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nFeature_RNA))
abline(h = log(500), v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nCount_RNA)
abline(h = 1000, v = 10)
text(1.5,1200, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nCount_RNA))
abline(h = log(1000), v = 10)
text(1.5,700, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))

# Only the first sample
smoothScatter(upn109_data_list[[1]]@meta.data$percent.mt, upn109_data_list[[1]]@meta.data$nFeature_RNA)
abline(h = 500, v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

# Only the second sample
smoothScatter(upn109_data_list[[2]]@meta.data$percent.mt, upn109_data_list[[2]]@meta.data$nFeature_RNA)
abline(h = 500, v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

upn109_data_list_filtered <- upn109_data_list

# Filtering for features, counts, and proportion of MT reads
# Some outliers in nCounts_Protein >10,000 in the old batches
for (i in 1:length(upn109_data_list_filtered)) { 
    upn109_data_list_filtered[[i]] = subset(upn109_data_list_filtered[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
}

# SCTransform normalization of RNA counts
for (i in 1:length(upn109_data_list_filtered)) {
    DefaultAssay(upn109_data_list_filtered[[i]]) = "RNA"
    upn109_data_list_filtered[[i]] = SCTransform(upn109_data_list_filtered[[i]],
                                               method = "glmGamPoi",
                                               vars.to.regress = c("percent.mt"),
                                               verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

for(i in 1:length(upn109_data_list_filtered)){
    message(names(upn109_data_list_filtered)[i])
    DefaultAssay(upn109_data_list_filtered[[i]]) <- "SCT"
    upn109_data_list_filtered[[i]] <- CellCycleScoring(upn109_data_list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(upn109_data_list_filtered)) {
    DefaultAssay(upn109_data_list_filtered[[i]]) = "RNA"
    upn109_data_list_filtered[[i]] = SCTransform(upn109_data_list_filtered[[i]],
                                               method = "glmGamPoi",
                                               vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                               verbose = T)
}

upn109_data_list_filtered[[1]]@meta.data$UPN <- 109
upn109_data_list_filtered[[1]]@meta.data$Tumor_Sample_Type <- "pretreatment"
upn109_data_list_filtered[[2]]@meta.data$UPN <- 109
upn109_data_list_filtered[[2]]@meta.data$Tumor_Sample_Type <- "posttreatment"

#saveRDS(upn109_data_list_filtered, "/scratch/hnatri/CART/upn109_data_list_filtered.rds")
# upn109_data_list_filtered <- readRDS("/scratch/hnatri/CART/upn109_data_list_filtered.rds")

# Adding to the list
#tumor_list <- c(tumor_list, upn109_data_list_filtered[[1]])
tumor_list <- c(tumor_list, upn109_data_list_filtered)
tumor_list <- tumor_list %>% purrr::discard(is.null)

names(tumor_list)[43] <- "UPN109_pretreatment"
names(tumor_list)[44] <- "UPN109_posttreatment"

saveRDS(tumor_list, "/scratch/hnatri/CART/tumor_list_all_samples_08242023.rds")

# Running SoupX with soupx_for_hashed.R
