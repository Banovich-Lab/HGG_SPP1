#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 8/23/2023
# Description: SoupX on 13384 tumor samples
#==============================================================================

library(Seurat)
library(DropletUtils)
library(SoupX)
library(tidyverse)
library(ggplot2)
library(googlesheets4)
library(dplyr)

#==============================================================================
# Environment
#==============================================================================

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
source("/home/hnatri/Utilities/utilities.R")

#==============================================================================
# Import data
#==============================================================================

# Metadata
gs4_deauth()
mastersheet  <- gs4_get("https://docs.google.com/spreadsheets/d/1hvbygedvsEAU4whYDQV1kjFoXPdCCNkK0v_K9zQRQPs/edit?usp=sharing")
cellranger_paths <- as.data.frame(read_sheet(mastersheet, sheet = "Batch CellRanger paths"))
mastersheet <- as.data.frame(read_sheet(mastersheet, sheet = "Master_Sheet_SC"))

# Filtered Seurat list
tumor_list <- readRDS("/scratch/hnatri/CART/tumor_list_all_samples_08242023.rds")
names(tumor_list) <- gsub("_filtered", "", names(tumor_list))
names(tumor_list) <- gsub("UPN109_", "UPN109", names(tumor_list))

for(i in c("UPN109pretreatment", "UPN109posttreatment")){
  tumor_list[[i]]$Batch <- i
  tumor_list[[i]] <- RenameCells(tumor_list[[i]],
                                 new.names = paste0(i, "_", colnames(tumor_list[[i]])))
}

tumor_list <- lapply(tumor_list, function(xx){
  xx@meta.data$Sample_ID <- paste0(xx@meta.data$UPN, "_", xx@meta.data$Batch)
  
  xx
})

demultiplexed_batches <- lapply(tumor_list, function(xx){
  demultiplexed <- SplitObject(xx, split.by = "Sample_ID")
  
  return(demultiplexed)
})

demultiplexed_batches <- purrr::flatten(demultiplexed_batches)

# Saving as 10x
for (i in 1:length(demultiplexed_batches)) {
  print(paste0("Converting sample ", names(demultiplexed_batches[i])))
  obj.sub <- demultiplexed_batches[[i]] 
  
  DropletUtils::write10xCounts(path = paste0("/scratch/hnatri/CART/SoupX/demultiplexed_", names(demultiplexed_batches[i])), 
                               x = obj.sub[["RNA"]]@data, 
                               barcodes = colnames(obj.sub[["RNA"]]@data), # cell names
                               gene.id = rownames(obj.sub[["RNA"]]@data), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}

#==============================================================================
# Running SoupX
#==============================================================================

tumor_list_SoupX <- sapply(names(demultiplexed_batches), function(i){
  print(i)
  
  if(ncol(demultiplexed_batches[[i]])<80){
    return(demultiplexed_batches[[i]])
  }
  # Read in count and droplet data
  # Converted demultiplexed counts
  message("Reading converted demultiplexed counts")
  d10x_toc <- Read10X(paste0("/scratch/hnatri/CART/SoupX/demultiplexed_", i))
  
  batch_id <- unique(sapply(strsplit(colnames(demultiplexed_batches[[i]]), "_"), `[`, 1))
  
  # Need to read in batch specific empty droplet file
  message("Reading empty droplet data")
  batch_path <- cellranger_paths %>% filter(Batch_name_2 == batch_id) %>% dplyr::pull(CellRanger_output)
  d10x_tod <- Read10X(paste0(batch_path, "/outs/raw_feature_bc_matrix/"))
  
  if(length(names(d10x_tod))>1){
    d10x_tod_gex <- d10x_tod[[1]]
  } else {
    d10x_tod_gex <- d10x_tod
  }
  
  colnames(d10x_tod_gex) <- paste0(i, "_", colnames(d10x_tod_gex))
  
  # Fixing cell names for some
  if(length(grep("CAR", colnames(d10x_tod_gex)))>1){
    colnames(d10x_tod_gex) <- gsub("-FB", "-FB_", colnames(d10x_tod_gex))
    cellnames <- colnames(d10x_toc)
    cellnames <- sapply(strsplit(cellnames,"_"), `[`, 2)
    cellnames <- paste0(i, "_", cellnames)
    
    colnames(d10x_toc) <- cellnames
  }
  
  if(i %in% c("UPN109_pretreatment", "UPN109_posttreatment")){
    colnames(d10x_toc) <- paste0(i, "_", colnames(d10x_toc))
  }
  
  # Some batches only have features that are expressed in at least one cell;
  # need to fix feature order
  rownames(d10x_toc) <- gsub("_", "-", rownames(d10x_toc))
  rownames(d10x_tod_gex) <- gsub("_", "-", rownames(d10x_tod_gex))
  d10x_tod_gex <- d10x_tod_gex[rownames(d10x_toc),]
  
  # Run SoupX
  sc <- SoupChannel(d10x_tod_gex, d10x_toc, calcSoupProfile = FALSE) 
  
  sc <- estimateSoup(sc)
  toc_seu <- CreateSeuratObject(d10x_toc)
  toc_seu <- SCTransform(toc_seu, vst.flavor = "v2")
  if(ncol(toc_seu)<50){
    toc_seu <- RunPCA(toc_seu, npcs = ncol(toc_seu)-1)
  } else{
    toc_seu <- RunPCA(toc_seu)
  }
  toc_seu <- RunUMAP(toc_seu, dims = 1:15)
  toc_seu <- FindNeighbors(toc_seu, dims = 1:15)
  toc_seu <- FindClusters(toc_seu, resolution = 1)

  ## Add meta data to soupX object
  sc <- setClusters(sc, setNames(toc_seu$seurat_clusters, rownames(toc_seu@meta.data)))
  ## Estimate contamination (automated method)
  message(paste0("Getting autoEstCont for: ", i))
  sc <- autoEstCont(sc, tfidfMin = 0.6, soupQuantile = 0.7, forceAccept = TRUE, doPlot = F) 
  out <- adjustCounts(sc)
  
  # Create Seurat object using corrected data
  d10x_seu <- CreateSeuratObject(out, assay = "SoupX_RNA")
  d10x_seu[["RNA"]] <- toc_seu@assays[["RNA"]]
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_SoupX_RNA", assay = "SoupX_RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_SoupX_RNA", assay = "SoupX_RNA")
  
  d10x_seu
})
names(tumor_list_SoupX) <- names(demultiplexed_batches)

saveRDS(tumor_list_SoupX, "/scratch/hnatri/CART/SoupX/tumor_list_demultiplexed_SoupX.rds")

#==============================================================================
# Normalizing for integration
#==============================================================================

tumor_list_SoupX <- readRDS("/scratch/hnatri/CART/SoupX/tumor_list_demultiplexed_SoupX.rds")

# Merging and splitting by batch
tumor_merged <- merge(x = tumor_list_SoupX[[1]], y = tumor_list_SoupX[2:length(tumor_list_SoupX)])

tumor_merged$Batch <- sapply(strsplit(colnames(tumor_merged), "_"), `[`, 1)
tumor_batch_list <- SplitObject(tumor_merged, split.by = "Batch")

for (i in 1:length(tumor_batch_list)) {
  DefaultAssay(tumor_batch_list[[i]]) = "RNA"
  tumor_batch_list[[i]] = PercentageFeatureSet(tumor_batch_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  tumor_batch_list[[i]] = SCTransform(tumor_batch_list[[i]],
                                      method = "glmGamPoi",
                                      vars.to.regress = c("percent.mt"),
                                      verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# F04739_filtered
# F04742_filtered
# F04719_filtered
# F04720_filtered
# F04714_filtered
# don't have enough expressed cell cycle genes to run cell cycle scoring
tumor_batch_list_tmp <- tumor_batch_list

Batch3 <- tumor_batch_listX[["Batch3"]]
Batch5 <- tumor_batch_listX[["Batch5"]]
Batch6 <- tumor_batch_listX[["Batch6"]]

F04739 <- tumor_batch_list[["F04739"]]
F04742 <- tumor_batch_list[["F04742"]]
F04719 <- tumor_batch_list[["F04719"]]
F04720 <- tumor_batch_list[["F04720"]]
F04714 <- tumor_batch_list[["F04714"]]

F05998 <- tumor_batch_list[["F05998-GEX_F05999-CAR_F06001-FB"]]

tumor_batch_list[["Batch3"]] <- NULL
tumor_batch_list[["Batch5"]] <- NULL
tumor_batch_list[["Batch6"]] <- NULL

tumor_batch_list[["F04739"]] <- NULL
tumor_batch_list[["F04742"]] <- NULL
tumor_batch_list[["F04719"]] <- NULL
tumor_batch_list[["F04720"]] <- NULL
tumor_batch_list[["F04714"]] <- NULL

tumor_batch_list[["44-1"]] <- NULL

for(i in 1:length(tumor_batch_list)){
  message(i)
  message(names(tumor_batch_list)[i])
  DefaultAssay(tumor_batch_list[[i]]) <- "SCT"
  tumor_batch_list[[i]] <- CellCycleScoring(tumor_batch_list[[i]],
                                            s.features = s.genes,
                                            g2m.features = g2m.genes,
                                            set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(tumor_batch_list)) {
  DefaultAssay(tumor_batch_list[[i]]) = "RNA"
  tumor_batch_list[[i]] = SCTransform(tumor_batch_list[[i]],
                                      method = "glmGamPoi",
                                      vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                      verbose = T)
}

tumor_batch_list[["F04739"]] <- F04739
tumor_batch_list[["F04742"]] <- F04742
tumor_batch_list[["F04719"]] <- F04719
tumor_batch_list[["F04720"]] <- F04720
tumor_batch_list[["F04714"]] <- F04714

tumor_batch_list[["Batch3"]] <- Batch3
tumor_batch_list[["Batch5"]] <- Batch5
tumor_batch_list[["Batch6"]] <- Batch6

tumor_batch_list[["44-1"]] <- F05998

names(tumor_batch_list)

saveRDS(tumor_batch_list, "/scratch/hnatri/CART/SoupX/tumor_list_demultiplexed_SoupX_normalized.rds")

#==============================================================================
# Integration
#==============================================================================

tumor_batch_list <- readRDS("/scratch/hnatri/CART/SoupX/tumor_list_demultiplexed_SoupX_normalized.rds")

# Excluding UPN 142
# F04774
tumor_batch_list[["F04774"]] <- NULL
tumor_batch_list[["UPN109posttreatment"]] <- NULL
tumor_batch_list <- tumor_batch_list %>% purrr::discard(is.null)

# Set nfeatures to 1000 to avoid error at IntegrateData
features <- SelectIntegrationFeatures(object.list = tumor_batch_list,
                                      nfeatures = 1000)

# Running PCA for rPCA integration
for (i in 1:length(tumor_batch_list)) {
  message(i)
  tumor_batch_list[[i]] <- RunPCA(tumor_batch_list[[i]],
                                  npcs = 30, # changing npcs to avoid an error
                                  features = features,
                                  approx = F,
                                  verbose = F)
}

tumor_batch_list <- PrepSCTIntegration(object.list = tumor_batch_list,
                                       anchor.features = features,
                                       verbose = F)

# Identify anchors and integrate the datasets based on RNA
# Using k=20 neighbors to find anchors (default = 5)?
# Using rPCA to avoid "problem too large" error with CCA
anchors <- FindIntegrationAnchors(object.list = tumor_batch_list,
                                  normalization.method = "SCT", 
                                  anchor.features = features,
                                  reference = seq(1, length(tumor_batch_list), by=3), # using reference batches to avoid an error
                                  k.anchor = 20,
                                  dims = 1:30,
                                  reduction = "rpca")
# Default value for k.weight causes an error here (not enough cells for some 
# samples)
tumor_integrated <- IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT", 
                                  new.assay.name = "integrated_sct",
                                  k.weight = 30,
                                  verbose = T)

# No need to run ScaleData if you've used SCT integration
tumor_integrated <- RunPCA(tumor_integrated,
                           reduction.name = "integrated_sct_pca",
                           verbose = F)
pcs <- get_pcs(tumor_integrated, reduction_name = "integrated_sct_pca")
message(pcs)
tumor_integrated <- RunUMAP(tumor_integrated,
                            reduction = "integrated_sct_pca",
                            reduction.name = "integrated_sct_umap",
                            dims = 1:min(pcs),
                            return.model = TRUE)
tumor_integrated <- FindNeighbors(tumor_integrated,
                                  reduction = "integrated_sct_pca",
                                  dims = 1:min(pcs),
                                  graph.name = c("integrated_sct_nn",
                                                 "integrated_sct_snn"))
# resolution 0.2
tumor_integrated <- FindClusters(tumor_integrated,
                                 resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                 graph.name = "integrated_sct_nn")

DimPlot(tumor_integrated,
        group.by = "integrated_sct_nn_res.0.3",
        reduction = "integrated_sct_umap",
        label = T) + coord_fixed(ratio=1)
DimPlot(tumor_integrated,
        group.by = "integrated_sct_nn_res.0.3",
        split.by = "integrated_sct_nn_res.0.3",
        reduction = "integrated_sct_umap",
        ncol = 4) & coord_fixed(ratio=1)

saveRDS(tumor_integrated, "/labs/banovich/BCTCSF/Heini/tumor_integrated_UPN142_soupX_nn.rds")

