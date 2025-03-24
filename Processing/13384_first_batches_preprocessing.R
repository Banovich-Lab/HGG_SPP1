#==============================================================================#
# Original author(s) : Stephanie L. Yahn,
#                      Heini M. Natri hnatri@tgen.org
# Date: 11/30/2021
# Description: Preprocessing of single cell sequencing data for the CAR T
# projects
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
# Batch 1 preprocessing
# Processed separately because it relies on Demuxlet for demultiplexing
#==============================================================================#

# Read in output from CellRanger (2 data types: Gene Expression and Antibody Capture)
# This reads the data in as a list, because it's multimodal
beet1 = Read10X(data.dir = "/tgen_labs/banovich/BCTCSF/Outs/TGen_Ref_F02323-GEX_F02325-AB/outs/filtered_feature_bc_matrix/")

# Simplify antibody names
# 189 antibodies in this pool
# 3 different antibody panels, some instances where the number differs, e.g. multiple
# CD4s with different antibody name
rownames(beet1$`Antibody Capture`)
rownames(beet1$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(beet1$`Antibody Capture`)), ' ')
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "mouse_")
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "rat_")
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "human_")

# Set up the Seurat object
batch1 = CreateSeuratObject(counts = beet1[["Gene Expression"]])
# Calculating the % of MT reads
batch1 = PercentageFeatureSet(batch1, pattern = "^MT-", col.name = "percent.mt")
# Adding a protein assay based on the antibody matrix. Underscores get replaced
# with dashes.
batch1[["Protein"]] = CreateAssayObject(beet1[["Antibody Capture"]][, colnames(x = batch1)])

# Add demuxlet metadata to assign patient identities to each cell.
# read in the output from demuxlet
demuxlet = fread(file ="/tgen_labs/banovich/BCTCSF/Stephanie/Demuxlet/Batch1/1_TGen_Ref_F02323-GEX_F02325-AB_4xexomes.best")

# Obtain the patient ID and the assignment from the BEST column. IDs are based
# on our exome sequencing data.

# Assignment: doublet or singlet
assignment = sapply(demuxlet$BEST, function(x) {
  strsplit(x,"-")[[1]][[1]]
})

# Patient ID
patient_id = sapply(demuxlet$BEST, function(x) {
  strsplit(x,"-")[[1]][[2]]
})

# Adding the new metadata to the object
CellsMeta = batch1@meta.data
CellsMeta["Exome_Sample_Name"] = patient_id
CellsMeta["Demultiplex_Assignment"] = assignment
batch1 = AddMetaData(batch1, CellsMeta)

# Add sample metadata
# UPN = unique patient number
# Cycle = which infusion
# Day = how many days after the infusion
# Using the exome sample names to map over the new metadata features
b1.meta.data = read_excel("/tgen_labs/banovich/BCTCSF/Stephanie/Batches_metadata_forR.xlsx", sheet = 1)
b1.meta.data = b1.meta.data %>% drop_na(UPN)

batch1@meta.data$UPN = plyr::mapvalues(x = batch1@meta.data$Exome_Sample_Name,
                                       from = b1.meta.data$Exome_Sample_Name,
                                       to = as.character(b1.meta.data$UPN))

batch1@meta.data$Sample_Type = plyr::mapvalues(x = batch1@meta.data$UPN,
                                               from = b1.meta.data$UPN,
                                               to = as.character(b1.meta.data$Sample_Type))

batch1@meta.data$Cycle = plyr::mapvalues(x = batch1@meta.data$UPN,
                                         from = b1.meta.data$UPN,
                                         to = as.character(b1.meta.data$Cycle))

batch1@meta.data$Day = plyr::mapvalues(x = batch1@meta.data$UPN,
                                       from = b1.meta.data$UPN,
                                       to = as.character(b1.meta.data$Day))

batch1@meta.data$Manufacture = plyr::mapvalues(x = batch1@meta.data$UPN,
                                               from = b1.meta.data$UPN,
                                               to = as.character(b1.meta.data$Manufacture))

# Used later to keep the leukapheresis PBMCs
batch1$Cycle_Day = paste("Cycle", batch1$Cycle, "_Day", batch1$Day, sep = "")

# Adding batch
batch1@meta.data["Batch"] = "Batch1"

# keep only singlets
Idents(batch1) = batch1$Demultiplex_Assignment
batch1_filtered = subset(batch1, idents = "SNG")

# Saving tumor sample count matrices for GEO
batch1_filtered_tumor <- subset(batch1_filtered, subset = Sample_Type == "Tumor" )
unique(batch1_filtered_tumor$UPN)

counts_243 <- subset(batch1_filtered_tumor, subset = UPN == "243")
counts_243 <- LayerData(counts_243, assay = "RNA", layer = "counts")
write.csv(counts_243, "/scratch/hnatri/counts_UPN243.csv", quote = F, row.names = T)

counts_185 <- subset(batch1_filtered_tumor, subset = UPN == "185")
counts_185 <- LayerData(counts_185, assay = "RNA", layer = "counts")
write.csv(counts_185, "/scratch/hnatri/counts_UPN185.csv", quote = F, row.names = T)

# Remove unwanted cells/keep wanted cells
# Keeping high quality, singlets, sample types that we want
# subset = nFeature_RNA > 500 & percent.mt < 25 & nCount_Protein < 10000
VlnPlot(batch1_filtered, features = c("nFeature_RNA", "nCount_RNA", "nFeature_Protein", "nCount_Protein", "percent.mt"), ncol = 3)
hist(batch1_filtered$nCount_Protein)
plot(ecdf(batch1_filtered$nCount_RNA))

batch1_filtered = subset(batch1_filtered, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10 & nCount_Protein < 10000)

# Currently has leukapheresis PBMCs and PBMCs matched to CSF
table(batch1_filtered$Manufacture)
table(batch1_filtered$Sample_Type)

# Normalize and scale protein data
# due to the unspecific binding background signal, log-normalization doesn't work well for CITEseq protein data
# instead, Seurat recommends centered log-ratio (CLR) normalization computed independently for each feature
batch1_filtered = NormalizeData(batch1_filtered, assay = "Protein", normalization.method = "CLR")
batch1_filtered = ScaleData(batch1_filtered, assay = "Protein")

rm(CellsMeta, demuxlet, b1.meta.data)

#==============================================================================#
# Batches 2-19 preprocessing
#==============================================================================#

# Using the hashtag oligos instead of demuxlet

# Read in output from CellRanger (2 data types: Gene Expression and Antibody Capture)
beet2 = Read10X("/tgen_labs/banovich/BCTCSF/Outs/F02397-GEX_F02398-AB/outs/filtered_feature_bc_matrix/") ##
beet3 = Read10X("/tgen_labs/banovich/BCTCSF/Outs/F02574-GEX_F02575-AB/outs/filtered_feature_bc_matrix/")
beet4 = Read10X("/tgen_labs/banovich/BCTCSF/Outs/F02588-GEX_F02589-AB/outs/filtered_feature_bc_matrix/")

beet5 = Read10X("/tgen_labs/banovich/BCTCSF/Outs/IL13OP/IL13OP_F02660-GEX_F02661-AB/outs/filtered_feature_bc_matrix/")
beet6 = Read10X("/tgen_labs/banovich/BCTCSF/Outs/IL13OP/IL13OP_F02677-GEX_F02678-AB/outs/filtered_feature_bc_matrix/")

data_loc = "/tgen_labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE/"
beet7 = Read10X(paste0(data_loc, "IL13OP_F02736-GEX_F02755-FB/outs/filtered_feature_bc_matrix/"))
beet8 = Read10X(paste0(data_loc, "IL13OP_F02981-GEX_F02982-FB/outs/filtered_feature_bc_matrix/"))
beet9 = Read10X(paste0(data_loc, "IL13OP_F02985-GEX_F02986-FB/outs/filtered_feature_bc_matrix/"))
beet10 = Read10X(paste0(data_loc, "IL13OP_F03246-GEX_F03247-FB/outs/filtered_feature_bc_matrix/"))
beet11 = Read10X(paste0(data_loc, "IL13OP_F03254-GEX_F03256-FB/outs/filtered_feature_bc_matrix/"))
beet12 = Read10X(paste0(data_loc, "IL13OP_F03253-GEX_F03255-FB/outs/filtered_feature_bc_matrix/"))
beet13 = Read10X(paste0(data_loc, "IL13OP_F03263-GEX_F03265-FB/outs/filtered_feature_bc_matrix/"))
beet14 = Read10X(paste0(data_loc, "IL13OP_F03264-GEX_F03266-FB/F03264-GEX_F03266-FB/outs/filtered_feature_bc_matrix/"))
beet15 = Read10X(paste0(data_loc, "IL13OP_F03285-GEX_F03286-FB/outs/filtered_feature_bc_matrix/"))
beet16 = Read10X(paste0(data_loc, "IL13OP_F03287-GEX_F03288-FB/outs/filtered_feature_bc_matrix/"))
beet17 = Read10X(paste0(data_loc, "IL13OP_F03289-GEX_F03290-FB/F03289-GEX_F03290-FB/outs/filtered_feature_bc_matrix/"))
beet18 = Read10X(paste0(data_loc, "IL13OP_F03291-GEX_F03292-FB/F03291-GEX_F03292-FB/outs/filtered_feature_bc_matrix/"))
beet19 = Read10X(paste0(data_loc, "IL13OP_F03293-GEX_F03294-FB/outs/filtered_feature_bc_matrix/"))

# Simplify antibody names
# Creating a list of objects
bt.list = ls(pattern="beet")
bt.list = str_sort(bt.list, numeric = TRUE)
bt.list = do.call("list", mget(bt.list))

saveRDS(bt.list, "/scratch/hnatri/CART/bt.list.rds")

# Looping through all the objects
for (i in 1:length(bt.list)) {
  rownames(bt.list[[i]]$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(bt.list[[i]]$`Antibody Capture`)), ' ')
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "mouse_")
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "rat_")
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "human_")
}

# Within the antibody capture, there are CITE-seq antibodies and the hashing
# antibodies, need to separate
# Cell hashing antibody names
hash_antibodies = c("TotalSeqC0251_Hashtag1", 
                    "TotalSeqC0252_Hashtag2",  
                    "TotalSeqC0253_Hashtag3", 
                    "TotalSeqC0254_Hashtag4",
                    "TotalSeqC0256_Hashtag6",
                    "TotalSeqC0257_Hashtag7", 
                    "TotalSeqC0258_Hashtag8",
                    "TotalSeqC0259_Hashtag9")

# Not all CITEseq antibodies were used for CITEseq batches 14, 17, and 18, 
# but cellranger was run as if they were, so need to filter the antibodies for these batches
# to remove false positive counts
# all_abs sheet (from Google Drive > feature_barcoding)
all_abs = read.csv(file = "/tgen_labs/banovich/BCTCSF/FeatureBarcoding/8plex_feature_reference - BCTCSF_FR.csv")
all_abs$ab_name = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', all_abs$name), ' ')
all_abs$ab_name = as.character(all_abs$ab_name)
all_abs$ab_name = str_remove(all_abs$ab_name, "mouse_")
all_abs$ab_name = str_remove(all_abs$ab_name, "rat_")
all_abs$ab_name = str_remove(all_abs$ab_name, "human_")

# panel_137 sheet (from Google Drive > feature_barcoding)
panel_137 = read.csv(file = "/tgen_labs/banovich/BCTCSF/FeatureBarcoding/TS-C human panel_137_Antibodies (#399905).xlsx - BarcodeList.csv")
# all_abs sheet has ab names in the correct format to match matrix, panel_137 does not
panel_137 = subset(all_abs, id %in% panel_137$DNA_ID)

# Filter antibodies in batches 14, 17, and 18
keep_prots = c(panel_137$ab_name, hash_antibodies)
bt.list$beet14$`Antibody Capture` = bt.list$beet14$`Antibody Capture`[rownames(bt.list$beet14$`Antibody Capture`) %in% keep_prots,]
bt.list$beet17$`Antibody Capture` = bt.list$beet17$`Antibody Capture`[rownames(bt.list$beet17$`Antibody Capture`) %in% keep_prots,]
bt.list$beet18$`Antibody Capture` = bt.list$beet18$`Antibody Capture`[rownames(bt.list$beet18$`Antibody Capture`) %in% keep_prots,]

# Iterating through the beet list to create a batch list and a filtered object
# list
batch.list = list()
batch.list_filtered = list()
for (i in 2:length(bt.list)) { #excluding Batch 1
  #i <- 2
  message(names(bt.list)[i])
  # Set up the Seurat object
  # Splitting out the gene expression and the antibody matrices
  GEX = bt.list[[i]][[1]]
  AB = bt.list[[i]][[2]]
  # Creating the object with the expression matrix
  batch.list[[i]] = CreateSeuratObject(counts = GEX)
  batch.list[[i]] = PercentageFeatureSet(batch.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  
  hash = AB[rownames(AB) %in% hash_antibodies,]
  citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
  
  # Creating the protein and hash assays
  batch.list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
  batch.list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
}


# Normalizing and filtering doublets
for (i in 2:length(batch.list)) { 
  # Normalizing and scaling
  batch.list[[i]] = NormalizeData(batch.list[[i]], assay = "Protein", normalization.method = "CLR")
  batch.list[[i]] = ScaleData(batch.list[[i]], assay = "Protein")
  batch.list[[i]] = NormalizeData(batch.list[[i]], assay = "Hash", normalization.method = "CLR")
  batch.list[[i]] = ScaleData(batch.list[[i]], assay = "Hash")
  
  # Assign single cells back to their sample origins
  # Demultiplexing based on the hashing antibodies
  # Singlets kept based on "hash classification global"
  batch.list[[i]] = HTODemux(batch.list[[i]], assay = "Hash", positive.quantile = 0.99, verbose = F)
  
  # Add sample metadata
  meta.data = read_excel("/tgen_labs/banovich/BCTCSF/Stephanie/Batches_metadata_forR.xlsx", sheet = i)
  meta.data = meta.data %>% drop_na(UPN)
  batch.list[[i]]@meta.data$UPN = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$UPN))
  
  batch.list[[i]]@meta.data$Sample_Type = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Sample_Type))
  
  batch.list[[i]]@meta.data$Cycle = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$Cycle))
  
  batch.list[[i]]@meta.data$Day = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$Day))
  batch.list[[i]]@meta.data$Manufacture = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                          from = meta.data$CellHashing_Ab,
                                                          to = as.character(meta.data$Manufacture))
  
  batch.list[[i]]$Cycle_Day = paste("Cycle", batch.list[[i]]$Cycle, "_Day", batch.list[[i]]$Day, sep = "")
  
  batch.list[[i]]@meta.data["Batch"] = paste("Batch", i, sep = "")
}

table(batch.list[[2]]$hash.ID)

# Numbers of singlets and doublets
batch.list[[1]] = batch1

saveRDS(batch.list, "/scratch/hnatri/CART/batch.list.rds")

bt_merge <- merge(x = batch.list[[1]], y = batch.list[2:length(batch.list)])

for (i in 2:length(batch.list)) { 
  # Keeping singlets only
  Idents(batch.list[[i]]) = batch.list[[i]]$Hash_classification.global
  # VlnPlot(batch.list[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_Protein", "nFeature_Protein", "nFeature_Hash"), pt.size = 0)
  batch.list_filtered[[i]] = subset(batch.list[[i]], idents = "Singlet")
}

# Smooth-scatter plot of MT reads and RNA counts
batch.list2 <- batch.list_filtered
batch.list2[[1]] = batch1
bt_merge <- merge(x = batch.list2[[1]], y = batch.list2[2:length(batch.list2)])
VlnPlot(bt_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# A few cells in batches 11 and 15 get misassigned as TotalSeqC0259-Hashtag9.
# These are removed later
head(bt_merge@meta.data[which(bt_merge@meta.data$UPN == "TotalSeqC0259-Hashtag9"),])
unique(bt_merge@meta.data[which(bt_merge@meta.data$UPN == "TotalSeqC0259-Hashtag9"),]$Batch)
dim(bt_merge@meta.data[which(bt_merge@meta.data$UPN == "TotalSeqC0259-Hashtag9"),])

# Make plots
# subset = nFeature_RNA > 500 & nCount_RNA > 5000 & percent.mt < 10 & nCount_Protein < 10000)
smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nFeature_RNA)
abline(h = 500, v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nCount_RNA)
abline(h = 1000, v = 10)
text(1.5,1200, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nCount_RNA))
abline(h = log(1000), v = 10)
text(1.5,700, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))
#dev.off()

# Filtering for features, counts, and proportion of MT reads
# Some outliers in nCounts_Protein >10,000
for (i in 2:length(batch.list_filtered)) { 
  batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10 & nCount_Protein < 10000)
}

for (i in 2:length(batch.list_filtered)) {
  # Scale protein/antibody data (must re-scale after subsetting)
  batch.list_filtered[[i]] = ScaleData(batch.list_filtered[[i]], assay = "Protein")
  batch.list_filtered[[i]] = ScaleData(batch.list_filtered[[i]], assay = "Hash")
}

names(batch.list_filtered)

# Remove cells erroneously identified as “TotalSeqC0259-Hashtag9” in batch11
# and batch15 as these batches did not include this hashing antibody
batch.list_filtered[[11]] = subset(batch.list_filtered[[11]], cells = row.names(batch.list_filtered[[11]]@meta.data[!(batch.list_filtered[[11]]@meta.data$UPN %in% "TotalSeqC0259-Hashtag9"), ]))
batch.list_filtered[[15]] = subset(batch.list_filtered[[15]], cells = row.names(batch.list_filtered[[15]]@meta.data[!(batch.list_filtered[[15]]@meta.data$UPN %in% "TotalSeqC0259-Hashtag9"), ]))
# Scale protein/antibody data (must re-scale after subsetting)
batch.list_filtered[[11]] = ScaleData(batch.list_filtered[[11]], assay = "Protein")
batch.list_filtered[[11]] = ScaleData(batch.list_filtered[[11]], assay = "Hash")
batch.list_filtered[[15]] = ScaleData(batch.list_filtered[[15]], assay = "Protein")
batch.list_filtered[[15]] = ScaleData(batch.list_filtered[[15]], assay = "Hash")

# Add batch 1 to lists
batch.list[[1]] = batch1
batch.list_filtered[[1]] = batch1_filtered

# Name list elements
batch_names = paste("Batch", seq(1,19), sep = "")
names(batch.list) = batch_names
batch_filt_names = paste(batch_names, "_filtered", sep = "")
names(batch.list_filtered) = batch_filt_names

# # cells in each batch
lapply(names(batch.list_filtered), function(xx){
  message(xx, " ", nrow(batch.list_filtered[[xx]]@meta.data))
})

# Remove Protein assay from non-CITEseq batches
nonCITE = c("Batch7_filtered", "Batch8_filtered", "Batch9_filtered", "Batch10_filtered", "Batch11_filtered", 
            "Batch12_filtered", "Batch13_filtered", "Batch15_filtered", "Batch16_filtered", "Batch19_filtered")
for (i in nonCITE) {
  message(i)
  batch.list_filtered[[i]]$Protein = NULL
}

# SCTransform normalization of RNA counts (replaces NormalizeData(),
# ScaleData(), and FindVariableFeatures()). SCTransform also supports using
# glmGamPoi package which substantially improves the speed of the learning
# procedure
for (i in 1:length(batch.list_filtered)) {
  DefaultAssay(batch.list_filtered[[i]]) = "RNA"
  batch.list_filtered[[i]] = SCTransform(batch.list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt"),
                                         verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

for(i in 1:length(batch.list_filtered)){
  message(names(batch.list_filtered)[i])
  DefaultAssay(batch.list_filtered[[i]]) <- "SCT"
  batch.list_filtered[[i]] <- CellCycleScoring(batch.list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(batch.list_filtered)) {
  DefaultAssay(batch.list_filtered[[i]]) = "RNA"
  batch.list_filtered[[i]] = SCTransform(batch.list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                         verbose = T)
}

saveRDS(batch.list_filtered, file = "/labs/banovich/BCTCSF/Heini/batch_list_filtered_all_sampletypes.rds")

#==============================================================================#
# Batches 21-39
#==============================================================================#

# Preventing trying to getcredentials
gs4_deauth()

#url <- 'docs.google.com/spreadsheets/d/1rCrE8KAelXHX_MHxi48119ZLmJhs29QZb1eGcc11mg0/edit#gid=1097420791'
new_batches_metadata <- gs4_get("https://docs.google.com/spreadsheets/d/1hvbygedvsEAU4whYDQV1kjFoXPdCCNkK0v_K9zQRQPs/edit?usp=sharing")
sheet_names(new_batches_metadata)

sc_metadata <- read_sheet(new_batches_metadata, sheet = "Master_Sheet_SC")
sc_metadata <- as.data.frame(sc_metadata)
sc_metadata <- sc_metadata[which(sc_metadata$Batch %in% c("21", "22", "23", "24", "25", "26", "27", "28", "37", "38", "39")),]
sc_metadata$FID_GEX_IL13OP <- paste0(sc_metadata$FID_GEXFB, "-GEX_", sc_metadata$FID_IL13, "-IL13")

# Processed samples, output locations
processed_samples <- gs4_get("https://docs.google.com/spreadsheets/d/1DjRQiumC60K0KwDVnq8T2vcXgeZBhOK6ZFMfho7pqM8/edit?usp=sharing")
sheet_names(processed_samples)

processed_samples <- read_sheet(processed_samples, sheet = "Sheet1")

# Product and PBMCs were run with cell hashing
# Tumor and CSF were run solo
# Batches 27, 28, and 39 only have tumor and CSF samples
data_loc <- "/labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE/"
bbeet21 <- Read10X(paste0(data_loc, "F04758-GEX_F04759-FB/outs/filtered_feature_bc_matrix/"))
bbeet22 <- Read10X(paste0(data_loc, "F04760-GEX_F04761-FB/outs/filtered_feature_bc_matrix/"))
bbeet23 <- Read10X(paste0(data_loc, "F04805-GEX_F04804-FB/outs/filtered_feature_bc_matrix/"))
bbeet24 <- Read10X(paste0(data_loc, "F04763-GEX_F04762-FB/outs/filtered_feature_bc_matrix/"))
bbeet25 <- Read10X(paste0(data_loc, "F04771-GEX_F04770-FB/outs/filtered_feature_bc_matrix/"))
bbeet26 <- Read10X(paste0(data_loc, "F04779-GEX_F04778-FB/outs/filtered_feature_bc_matrix/"))
bbeet37 <- Read10X(paste0(data_loc, "F04855-GEX_F04856-FB/outs/filtered_feature_bc_matrix/"))
bbeet38 <- Read10X(paste0(data_loc, "F04860-GEX_F04861-FB/outs/filtered_feature_bc_matrix/"))

# Simplify antibody names
# Creating a list of objects
bt.list = ls(pattern="bbeet")
bt.list = str_sort(bt.list, numeric = TRUE)
bt.list = do.call("list", mget(bt.list))

# NOTE: no CITEseq for new batches, but the protein assay was created
# Fixing antibody names
for (i in 1:length(bt.list)) {
    rownames(bt.list[[i]]$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(bt.list[[i]]$`Antibody Capture`)), ' ')
    rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "mouse_")
    rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "rat_")
    rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "human_")
}

# Iterating through the beet list to create a batch list and a filtered object
# list
new_batch.list = list()
new_batch.list_filtered = list()
for (i in 1:length(bt.list)) { #excluding Batch 1
    #i <- 2
    message(names(bt.list)[i])
    # Set up the Seurat object
    # Splitting out the gene expression and the antibody matrices
    GEX = bt.list[[i]][[1]]
    AB = bt.list[[i]][[2]]
    # Creating the object with the expression matrix
    new_batch.list[[i]] = CreateSeuratObject(counts = GEX)
    new_batch.list[[i]] = PercentageFeatureSet(new_batch.list[[i]], pattern = "^MT-", col.name = "percent.mt")
    
    hash = AB[rownames(AB) %in% hash_antibodies,]
    citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
    
    # Creating the protein and hash assays
    new_batch.list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
    new_batch.list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
}

# Normalizing and filtering doublets
for (i in 1:length(new_batch.list)) {
    # Keeping antibodies that were used in the panel
    batch_num <- as.numeric(as.character(gsub("bbeet", "", names(bt.list)[[i]])))
    message(batch_num)
    meta.data <- sc_metadata[which(sc_metadata$Batch==batch_num),]
    meta.data$CellHashing_Ab <- gsub("_", "-", meta.data$CellHashing_Ab)
    new_batch.list[[i]][["Hash"]] = CreateAssayObject(counts = new_batch.list[[i]]@assays$Hash[which(rownames(new_batch.list[[i]]@assays$Hash) %in% meta.data$CellHashing_Ab),])
    
    # Normalizing and scaling
    new_batch.list[[i]] = NormalizeData(new_batch.list[[i]], assay = "Protein", normalization.method = "CLR")
    new_batch.list[[i]] = ScaleData(new_batch.list[[i]], assay = "Protein")
    new_batch.list[[i]] = NormalizeData(new_batch.list[[i]], assay = "Hash", normalization.method = "CLR")
    new_batch.list[[i]] = ScaleData(new_batch.list[[i]], assay = "Hash")
    
    rowSums(new_batch.list[[i]]@assays$Hash)
    colSums(new_batch.list[[i]]@assays$Hash)

    # Assign single cells back to their sample origins
    # Demultiplexing based on the hashing antibodies
    # Singlets kept based on "hash classification global"
    new_batch.list[[i]] = HTODemux(new_batch.list[[i]], assay = "Hash", positive.quantile = 0.99, verbose = F)
    
    # Add sample metadata
    new_batch.list[[i]]@meta.data$UPN = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$UPN))
    
    new_batch.list[[i]]@meta.data$IRB = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                            from = meta.data$CellHashing_Ab,
                                                            to = meta.data$IRB)
    
    new_batch.list[[i]]@meta.data$Sample_Type = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                            from = meta.data$CellHashing_Ab,
                                                            to = as.character(meta.data$Sample_Type))
    
    new_batch.list[[i]]@meta.data$Cycle = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                      from = meta.data$CellHashing_Ab,
                                                      to = as.character(meta.data$Cycle))
    
    new_batch.list[[i]]@meta.data$Day = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                    from = meta.data$CellHashing_Ab,
                                                    to = as.character(meta.data$Day))
    
    new_batch.list[[i]]@meta.data$Manufacture = plyr::mapvalues(x = new_batch.list[[i]]@meta.data$hash.ID,
                                                            from = meta.data$CellHashing_Ab,
                                                            to = as.character(meta.data$Manufacture))
    
    new_batch.list[[i]]$Cycle_Day = paste("Cycle", new_batch.list[[i]]$Cycle, "_Day", new_batch.list[[i]]$Day, sep = "")
    
    new_batch.list[[i]]@meta.data["Batch"] = paste("Batch", batch_num, sep = "")
    
}

for (i in 1:length(new_batch.list)) { 
    # Keeping singlets only
    Idents(new_batch.list[[i]]) = new_batch.list[[i]]$Hash_classification.global
    # VlnPlot(batch.list[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_Protein", "nFeature_Protein", "nFeature_Hash"), pt.size = 0)
    new_batch.list_filtered[[i]] = subset(new_batch.list[[i]], idents = "Singlet")
}

# Name list elements
batch_names = gsub("bbeet", "Batch", names(bt.list))
names(new_batch.list) = batch_names
batch_filt_names = paste(batch_names, "_filtered", sep = "")
names(new_batch.list_filtered) = batch_filt_names

# Smooth-scatter plot of MT reads and RNA counts
batch.list2 <- new_batch.list_filtered
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
# Some outliers in nCounts_Protein >10,000 in the old batches
for (i in 1:length(new_batch.list_filtered)) { 
    new_batch.list_filtered[[i]] = subset(new_batch.list_filtered[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10 & nCount_Protein < 10000)
}

batch.list_filtered_merge <- merge(x = new_batch.list_filtered[[1]], y = new_batch.list_filtered[2:length(new_batch.list_filtered)])

for (i in 1:length(new_batch.list_filtered)) {
    # Keeping all samples
    
    # remove CSF and Tumor samples
    #batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], 
    #                                  cells = row.names(batch.list_filtered[[i]]@meta.data[!(batch.list_filtered[[i]]@meta.data$Sample_Type %in% c("Tumor", "CSF")), ]))
    
    # keep only leuk_PBMC and Product samples
    #Idents(batch.list_filtered[[i]]) = batch.list_filtered[[i]]$Cycle_Day
    #batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], idents = "CycleNA_DayNA")
    
    # Scale protein/antibody data (must re-scale after subsetting)
    new_batch.list_filtered[[i]] = ScaleData(new_batch.list_filtered[[i]], assay = "Protein")
    new_batch.list_filtered[[i]] = ScaleData(new_batch.list_filtered[[i]], assay = "Hash")
}

# # cells in each batch
lapply(names(new_batch.list_filtered), function(xx){
    message(xx, " ", nrow(new_batch.list_filtered[[xx]]@meta.data))
})

# SCTransform normalization of RNA counts
for (i in 1:length(new_batch.list_filtered)) {
    DefaultAssay(new_batch.list_filtered[[i]]) = "RNA"
    new_batch.list_filtered[[i]] = SCTransform(new_batch.list_filtered[[i]],
                                               method = "glmGamPoi",
                                               vars.to.regress = c("percent.mt"),
                                               verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

for(i in 1:length(new_batch.list_filtered)){
    message(names(new_batch.list_filtered)[i])
    DefaultAssay(new_batch.list_filtered[[i]]) <- "SCT"
    new_batch.list_filtered[[i]] <- CellCycleScoring(new_batch.list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(new_batch.list_filtered)) {
    DefaultAssay(new_batch.list_filtered[[i]]) = "RNA"
    new_batch.list_filtered[[i]] = SCTransform(new_batch.list_filtered[[i]],
                                           method = "glmGamPoi",
                                           vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                           verbose = T)
}

saveRDS(new_batch.list_filtered, "/labs/banovich/BCTCSF/Heini/new_batch_list_filtered_all_sampletypes.rds")

old_new_batches.list_filtered <- c(batch.list_filtered, new_batch.list_filtered)
old_new_batches.list_filtered_merge <- merge(x = old_new_batches.list_filtered[[1]], y = old_new_batches.list_filtered[2:length(old_new_batches.list_filtered)])

#==============================================================================#
# Adding tumor and CSF samples
#==============================================================================#

# Tumor and CSF samples were run solo, without cell hashing
# Batches 27, 28, and 39 only have tumor and CSF samples
head(sc_metadata)
head(processed_samples)
tumor_csf_samples <- processed_samples[which(grepl("BCTCSF", processed_samples$`Isolon Directory`) & processed_samples$Type == "GEX 5' v2" & is.na(processed_samples$`Antibody panel`)),]
dim(tumor_csf_samples)
sc_metadata_tumor_csf <- sc_metadata[which(sc_metadata$FID_GEXFB %in% tumor_csf_samples$Sample),]
dim(sc_metadata_tumor_csf)

# Creating Seurat objects
tumor_csf.list <- lapply(tumor_csf_samples$Sample, function(i){
    # i <- "F04745"
    message(i)
    sample_10x_data <- Read10X(paste0(tumor_csf_samples[which(tumor_csf_samples$Sample==i),]$`Isolon Directory`, "/outs/filtered_feature_bc_matrix/"))
    seurat_object <- CreateSeuratObject(counts = sample_10x_data)
    seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
})
names(tumor_csf.list) <- tumor_csf_samples$Sample

# Adding metadata
for (i in names(tumor_csf.list)) {
    # Keeping antibodies that were used in the panel
    #fid <- names(tumor_csf.list)[[i]]
    message(i)
    meta.data <- sc_metadata[which(sc_metadata$FID_GEX_FB==i),]

    # Add sample metadata
    tumor_csf.list[[i]]@meta.data$UPN = meta.data$UPN
    tumor_csf.list[[i]]@meta.data$Sample_Type = meta.data$Sample_Type
    tumor_csf.list[[i]]@meta.data$Cycle = meta.data$Cycle
    tumor_csf.list[[i]]@meta.data$Day = meta.data$Day
    tumor_csf.list[[i]]@meta.data$Manufacture = meta.data$Manufacture
    tumor_csf.list[[i]]@meta.data$IRB = meta.data$IRB
    tumor_csf.list[[i]]@meta.data$Batch = paste0("Batch", meta.data$Batch)
    tumor_csf.list[[i]]$Cycle_Day = paste("Cycle", tumor_csf.list[[i]]$Cycle, "_Day", tumor_csf.list[[i]]$Day, sep = "")
}

bt_merge <- merge(x = tumor_csf.list[[1]], y = tumor_csf.list[2:length(tumor_csf.list)])
#bt_merge <- PercentageFeatureSet(object = bt_merge, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(bt_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Make density plots
# subset = nFeature_RNA > 500 & nCount_RNA > 5000 & percent.mt < 10 & nCount_Protein < 10000)
smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nFeature_RNA, ylab = "nFeature_RNA", xlab = "% MT")
abline(h = 500, v = 10)
text(15,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nFeature_RNA), ylab = "log(nFeature_RNA)", xlab = "% MT")
abline(h = log(500), v = 10)
text(15,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nCount_RNA, ylab = "nCount_RNA", xlab = "% MT")
abline(h = 1000, v = 10)
text(15,1200, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nCount_RNA), ylab = "log(nCount_RNA)", xlab = "% MT")
abline(h = log(1000), v = 10)
text(15,7, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))

# Filtering for features, counts, and proportion of MT reads
tumor_csf.list_filtered <- list()
for (i in 1:length(tumor_csf.list)) { 
    tumor_csf.list_filtered[[i]] = subset(tumor_csf.list[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
}

tumor_csf_filt_names <- paste(names(tumor_csf.list), "_filtered", sep = "")
names(tumor_csf.list_filtered) <- tumor_csf_filt_names

# # cells per sample
lapply(names(tumor_csf.list_filtered), function(xx){
    message(xx, " ", nrow(tumor_csf.list_filtered[[xx]]@meta.data))
})

# SCTransform normalization of RNA count
for (i in 1:length(tumor_csf.list_filtered)) {
    DefaultAssay(tumor_csf.list_filtered[[i]]) = "RNA"
    tumor_csf.list_filtered[[i]] = SCTransform(tumor_csf.list_filtered[[i]],
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
F04739_filtered <- tumor_csf.list_filtered[["F04739_filtered"]]
F04742_filtered <- tumor_csf.list_filtered[["F04742_filtered"]]
F04719_filtered <- tumor_csf.list_filtered[["F04719_filtered"]]
F04720_filtered <- tumor_csf.list_filtered[["F04720_filtered"]]
F04714_filtered <- tumor_csf.list_filtered[["F04714_filtered"]]
tumor_csf.list_filtered[["F04739_filtered"]] <- NULL
tumor_csf.list_filtered[["F04742_filtered"]] <- NULL
tumor_csf.list_filtered[["F04719_filtered"]] <- NULL
tumor_csf.list_filtered[["F04720_filtered"]] <- NULL
tumor_csf.list_filtered[["F04714_filtered"]] <- NULL
for(i in 1:length(tumor_csf.list_filtered)){
    message(i)
    message(names(tumor_csf.list_filtered)[i])
    DefaultAssay(tumor_csf.list_filtered[[i]]) <- "SCT"
    tumor_csf.list_filtered[[i]] <- CellCycleScoring(tumor_csf.list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(tumor_csf.list_filtered)) {
    DefaultAssay(tumor_csf.list_filtered[[i]]) = "RNA"
    tumor_csf.list_filtered[[i]] = SCTransform(tumor_csf.list_filtered[[i]],
                                           method = "glmGamPoi",
                                           vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                           verbose = T)
}

tumor_csf.list_filtered[["F04739_filtered"]] <- F04739_filtered
tumor_csf.list_filtered[["F04742_filtered"]] <- F04742_filtered
tumor_csf.list_filtered[["F04719_filtered"]] <- F04719_filtered
tumor_csf.list_filtered[["F04720_filtered"]] <- F04720_filtered
tumor_csf.list_filtered[["F04714_filtered"]] <- F04714_filtered

saveRDS(tumor_csf.list_filtered, "/labs/banovich/BCTCSF/Heini/tumor_csf_list_filtered.rds")
