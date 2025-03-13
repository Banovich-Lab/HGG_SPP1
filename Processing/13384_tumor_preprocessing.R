#==============================================================================#
# Original author(s) : Heini M. Natri hnatri@tgen.org
# Date: 2023/6/9
# Description: Preprocessing of newly sequenced 13384 tumors
#==============================================================================#

#==============================================================================#
# Load libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)
library(stringr)
library(readxl)
library(glmGamPoi)
library(gridExtra)
library(googlesheets4)

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Batches 42-44
#==============================================================================#

# Preventing trying to getcredentials
gs4_deauth()

batches_metadata <- gs4_get("https://docs.google.com/spreadsheets/d/1hvbygedvsEAU4whYDQV1kjFoXPdCCNkK0v_K9zQRQPs/edit?usp=sharing")

sc_metadata <- read_sheet(batches_metadata, sheet = "Master_Sheet_SC")
new_tumors_metadata <- sc_metadata[which(sc_metadata$Batch %in% c("42-1", "42-2", "43-1", "43-2", "44-1", "44-2")),]

# Processed samples, output locations
processed_samples <- gs4_get("https://docs.google.com/spreadsheets/d/1DjRQiumC60K0KwDVnq8T2vcXgeZBhOK6ZFMfho7pqM8/edit?usp=sharing")

processed_samples <- read_sheet(processed_samples, sheet = "Sheet1")
new_tumors_seq <- processed_samples %>% filter(`FASTQ Directory/Flowcell` %in% sc_metadata$Flowcell_ID)

# CellRanger output: /tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE
# e.g. F05974-GEX_F05980-CAR_F05989-FB
data_loc <- "/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE/"

# Creating a list of objects
new_tumors_list <- lapply(setdiff(unique(new_tumors_metadata$FID_GEX_CAR_FB), c(NA)), function(xx){
  Read10X(paste0(data_loc, xx, "/outs/filtered_feature_bc_matrix/"))
})

# Simplify antibody names
# Fixing antibody names
for (i in 1:length(new_tumors_list)) {
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(new_tumors_list[[i]]$`Antibody Capture`)), ' ')
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "mouse_")
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "rat_")
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "human_")
}

# Cell hashing antibody names
hash_antibodies = c("TotalSeqC0251_Hashtag1", 
                    "TotalSeqC0252_Hashtag2",  
                    "TotalSeqC0253_Hashtag3", 
                    "TotalSeqC0254_Hashtag4",
                    "TotalSeqC0256_Hashtag6",
                    "TotalSeqC0257_Hashtag7", 
                    "TotalSeqC0258_Hashtag8",
                    "TotalSeqC0259_Hashtag9")

# Iterating through the beet list to create a batch list and a filtered object
# list
# NOTE: no CITEseq, but the protein assay was created
batch_list <- list()
batch_list_filtered <- list()
for (i in 1:length(new_tumors_list)) {
  #i <- 2
  message(names(new_tumors_list)[i])
  # Set up the Seurat object
  # Splitting out the gene expression and the antibody matrices
  GEX = new_tumors_list[[i]][[1]]
  AB = new_tumors_list[[i]][[2]]
  # Creating the object with the expression matrix
  batch_list[[i]] = CreateSeuratObject(counts = GEX, min.cells=1)
  batch_list[[i]] = PercentageFeatureSet(batch_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  
  hash = AB[rownames(AB) %in% hash_antibodies,]
  citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
  
  # Creating the protein and hash assays
  batch_list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
  batch_list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
}

names(batch_list) <- names(new_tumors_list)

# Normalizing and filtering doublets
for (i in names(batch_list)) {
  # Keeping antibodies that were used in the panel
  message(i)
  meta.data <- new_tumors_metadata[which(new_tumors_metadata$FID_GEX_CAR_FB==i),]
  meta.data$CellHashing_Ab <- gsub("_", "-", meta.data$CellHashing_Ab)
  counts <- GetAssayData(batch_list[[i]], assay = "Hash")
  batch_list[[i]][["Hash"]] = CreateAssayObject(counts = counts[which(rownames(counts) %in% meta.data$CellHashing_Ab),])
  
  # Normalizing and scaling
  batch_list[[i]] = NormalizeData(batch_list[[i]], assay = "Protein", normalization.method = "CLR")
  batch_list[[i]] = ScaleData(batch_list[[i]], assay = "Protein")
  batch_list[[i]] = NormalizeData(batch_list[[i]], assay = "Hash", normalization.method = "CLR")
  batch_list[[i]] = ScaleData(batch_list[[i]], assay = "Hash")
  
  #batch.list[[i]]@assays$Hash[c("TotalSeqC0254-Hashtag4"),]
  rowSums(batch_list[[i]]@assays$Hash)
  colSums(batch_list[[i]]@assays$Hash)
  #zero_hto <- colnames(batch.list[[i]])[which(colSums(batch.list[[i]]@assays$Hash)<2)]
  #batch.list[[i]] <- subset(batch.list[[i]], cells = zero_hto, invert = T)
  
  # Assign single cells back to their sample origins
  # Demultiplexing based on the hashing antibodies
  # Singlets kept based on "hash classification global"
  batch_list[[i]] = HTODemux(batch_list[[i]], assay = "Hash", positive.quantile = 0.99, verbose = T)
  
  # Add sample metadata
  #meta.data = read_excel("/tgen_labs/banovich/BCTCSF/Stephanie/Batches_metadata_forR.xlsx", sheet = i)
  #meta.data = meta.data %>% drop_na(UPN)
  batch_list[[i]]@meta.data$UPN = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$UPN))
  
  batch_list[[i]]@meta.data$IRB = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = meta.data$IRB)
  
  batch_list[[i]]@meta.data$Sample_Type = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Sample_Type))
  
  batch_list[[i]]@meta.data$Cycle = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$Cycle))
  
  batch_list[[i]]@meta.data$Day = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$Day))
  
  batch_list[[i]]@meta.data$Manufacture = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Manufacture))
  
  batch_list[[i]]$Cycle_Day = paste("Cycle", batch_list[[i]]$Cycle, "_Day", batch_list[[i]]$Day, sep = "")
  
  batch_list[[i]]@meta.data$Batch = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                        from = meta.data$CellHashing_Ab,
                                                        to = as.character(meta.data$Batch))
}

for (i in 1:length(batch_list)) { 
  # Keeping singlets only
  Idents(batch_list[[i]]) = batch_list[[i]]$Hash_classification.global
  batch_list_filtered[[i]] = subset(batch_list[[i]], idents = "Singlet")
}

# Name list elements
names(batch_list_filtered) <- names(new_tumors_list)

# Smooth-scatter plot of MT reads and RNA counts
batch.list2 <- batch_list_filtered
bt_merge <- merge(x = batch.list2[[1]], y = batch.list2[2:length(batch.list2)])
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

# Filtering for features, counts, and proportion of MT reads
for (i in 1:length(batch_list_filtered)) { 
  batch_list_filtered[[i]] = subset(batch_list_filtered[[i]],
                                    subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
}

for (i in 1:length(batch_list_filtered)) {
  # Keeping all samples
  # Scale protein/antibody data (must re-scale after subsetting)
  batch_list_filtered[[i]] = ScaleData(batch_list_filtered[[i]], assay = "Protein")
  batch_list_filtered[[i]] = ScaleData(batch_list_filtered[[i]], assay = "Hash")
}
names(batch_list_filtered)

# SCTransform normalization of RNA counts
for (i in 1:length(batch_list_filtered)) {
  DefaultAssay(batch_list_filtered[[i]]) = "RNA"
  batch_list_filtered[[i]] = SCTransform(batch_list_filtered[[i]],
                                        method = "glmGamPoi",
                                        vars.to.regress = c("percent.mt"),
                                        verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Unable to calculate cell cycle scores for the first library with 166 cells
for (i in 2:length(batch_list_filtered)){
  message(names(batch_list_filtered)[i])
  DefaultAssay(batch_list_filtered[[i]]) <- "SCT"
  batch_list_filtered[[i]] <- CellCycleScoring(batch_list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 2:length(batch_list_filtered)) {
  DefaultAssay(batch_list_filtered[[i]]) = "RNA"
  batch_list_filtered[[i]] = SCTransform(batch_list_filtered[[i]],
                                              method = "glmGamPoi",
                                              vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                              verbose = T)
}

saveRDS(batch_list_filtered, "/tgen_labs/banovich/BCTCSF/Heini/new_tumors_filtered_230609.rds")

#==============================================================================#
# Batches 45-47
#==============================================================================#

# Preventing trying to getcredentials
gs4_deauth()

#url <- 'docs.google.com/spreadsheets/d/1rCrE8KAelXHX_MHxi48119ZLmJhs29QZb1eGcc11mg0/edit#gid=1097420791'
batches_metadata <- gs4_get("https://docs.google.com/spreadsheets/d/1hvbygedvsEAU4whYDQV1kjFoXPdCCNkK0v_K9zQRQPs/edit?usp=sharing")
sheet_names(batches_metadata)

sc_metadata <- read_sheet(batches_metadata, sheet = "Master_Sheet_SC")
unique(sc_metadata$Batch)
new_tumors_metadata <- sc_metadata[which(sc_metadata$Batch %in% c("45", "46", "47")),]

# Processed samples, output locations
processed_samples <- gs4_get("https://docs.google.com/spreadsheets/d/1DjRQiumC60K0KwDVnq8T2vcXgeZBhOK6ZFMfho7pqM8/edit?usp=sharing")
sheet_names(processed_samples)

processed_samples <- read_sheet(processed_samples, sheet = "Sheet1")
new_tumors_seq <- processed_samples %>% filter(`FASTQ Directory/Flowcell` %in% sc_metadata$Flowcell_ID)

# CellRanger output
data_loc <- "/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE/"

# Creating a list of objects
new_tumors_list <- lapply(setdiff(unique(new_tumors_metadata$FID_GEX_CAR_FB), c(NA)), function(xx){
  Read10X(paste0(data_loc, xx, "/outs/filtered_feature_bc_matrix/"))
})

names(new_tumors_list) <- setdiff(unique(new_tumors_metadata$FID_GEX_CAR_FB), c(NA))

# Simplify antibody names
# Fixing antibody names
for (i in 1:length(new_tumors_list)) {
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(new_tumors_list[[i]]$`Antibody Capture`)), ' ')
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "mouse_")
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "rat_")
  rownames(new_tumors_list[[i]]$`Antibody Capture`) = str_remove(rownames(new_tumors_list[[i]]$`Antibody Capture`), "human_")
}

# Cell hashing antibody names
# The first batch only has a small number of cells, especially for hashtag 8.
# Including it in HTODemux produces an error.
hash_antibodies_1 = c(#"TotalSeqC0258_Hashtag8", 
                      "TotalSeqC0262_Hashtag12",  
                      "TotalSeqC0263_Hashtag13", 
                      "TotalSeqC0264_Hashtag14",
                      "TotalSeqC0265_Hashtag15",
                      "TotalSeqC0266_Hashtag16")

hash_antibodies_2 = c("TotalSeqC0258_Hashtag8", 
                      "TotalSeqC0262_Hashtag12",  
                      "TotalSeqC0263_Hashtag13", 
                      "TotalSeqC0264_Hashtag14",
                      "TotalSeqC0265_Hashtag15",
                      "TotalSeqC0266_Hashtag16")

# Iterating through the beet list to create a batch list and a filtered object
# list
# NOTE: no CITEseq, but the protein assay was created
batch_list <- list()
batch_list_filtered <- list()
for (i in 1:length(new_tumors_list)) {
  #i <- 2
  message(names(new_tumors_list)[i])
  # Set up the Seurat object
  # Splitting out the gene expression and the antibody matrices
  GEX = new_tumors_list[[i]][[1]]
  AB = new_tumors_list[[i]][[2]]
  # Creating the object with the expression matrix
  batch_list[[i]] = CreateSeuratObject(counts = GEX, min.cells=1)
  batch_list[[i]] = PercentageFeatureSet(batch_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  
  if (1==1){
    hash_antibodies <- hash_antibodies_1
  } else {
    hash_antibodies <- hash_antibodies_2
  }
  
  hash = AB[rownames(AB) %in% hash_antibodies,]
  citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
  
  # Creating the protein and hash assays
  batch_list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
  batch_list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
}

names(batch_list) <- names(new_tumors_list)

# Normalizing and filtering doublets
for (i in names(batch_list)) {
  # Keeping antibodies that were used in the panel
  message(i)
  meta.data <- new_tumors_metadata[which(new_tumors_metadata$FID_GEX_CAR_FB==i),]
  meta.data$CellHashing_Ab <- gsub("_", "-", meta.data$CellHashing_Ab)
  batch_list[[i]][["Hash"]] = CreateAssayObject(counts = batch_list[[i]]@assays$Hash[which(rownames(batch_list[[i]]@assays$Hash) %in% meta.data$CellHashing_Ab),])
  
  # Normalizing and scaling
  batch_list[[i]] = NormalizeData(batch_list[[i]], assay = "Protein", normalization.method = "CLR")
  batch_list[[i]] = ScaleData(batch_list[[i]], assay = "Protein")
  batch_list[[i]] = NormalizeData(batch_list[[i]], assay = "Hash", normalization.method = "CLR")
  batch_list[[i]] = ScaleData(batch_list[[i]], assay = "Hash")
  
  # Assign single cells back to their sample origins
  batch_list[[i]] = HTODemux(batch_list[[i]], assay = "Hash", positive.quantile = 0.99, verbose = T)
  
  # Add sample metadata
  batch_list[[i]]@meta.data$UPN = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$UPN))
  
  batch_list[[i]]@meta.data$IRB = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = meta.data$IRB)
  
  batch_list[[i]]@meta.data$Sample_Type = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Sample_Type))
  
  batch_list[[i]]@meta.data$Cycle = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$Cycle))
  
  batch_list[[i]]@meta.data$Day = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$Day))
  
  batch_list[[i]]@meta.data$Manufacture = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Manufacture))
  
  batch_list[[i]]$Cycle_Day = paste("Cycle", batch_list[[i]]$Cycle, "_Day", batch_list[[i]]$Day, sep = "")
  
  batch_list[[i]]@meta.data$Batch = plyr::mapvalues(x = batch_list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$Batch))
}

for (i in 1:length(batch_list)) { 
  # Keeping singlets only
  Idents(batch_list[[i]]) = batch_list[[i]]$Hash_classification.global
  batch_list_filtered[[i]] = subset(batch_list[[i]], idents = "Singlet")
}

# Name list elements
names(batch_list_filtered) <- names(new_tumors_list)

# Smooth-scatter plot of MT reads and RNA counts
batch.list2 <- batch_list_filtered
bt_merge <- merge(x = batch.list2[[1]], y = batch.list2[2:length(batch.list2)])
VlnPlot(bt_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Make plots
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

# Filtering for features, counts, and proportion of MT reads
for (i in 1:length(batch_list_filtered)) { 
  batch_list_filtered[[i]] = subset(batch_list_filtered[[i]],
                                    subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
}

for (i in 1:length(batch_list_filtered)) {
  # Keeping all samples
  # Scale protein/antibody data (must re-scale after subsetting)
  batch_list_filtered[[i]] = ScaleData(batch_list_filtered[[i]], assay = "Protein")
  batch_list_filtered[[i]] = ScaleData(batch_list_filtered[[i]], assay = "Hash")
}
names(batch_list_filtered)

# SCTransform normalization of RNA counts
for (i in 1:length(batch_list_filtered)) {
  DefaultAssay(batch_list_filtered[[i]]) = "RNA"
  batch_list_filtered[[i]] = SCTransform(batch_list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt"),
                                         verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Unable to calculate cell cycle scores for the first library with 166 cells
for (i in 1:length(batch_list_filtered)){
  message(names(batch_list_filtered)[i])
  DefaultAssay(batch_list_filtered[[i]]) <- "SCT"
  batch_list_filtered[[i]] <- CellCycleScoring(batch_list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(batch_list_filtered)) {
  DefaultAssay(batch_list_filtered[[i]]) = "RNA"
  batch_list_filtered[[i]] = SCTransform(batch_list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                         verbose = T)
}

saveRDS(batch_list_filtered, "/tgen_labs/banovich/BCTCSF/Heini/new_tumors_filtered_230810.rds")

