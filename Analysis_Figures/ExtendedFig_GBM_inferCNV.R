#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2023/07/13
# Description: infercnv analysis on 13384 CAR T trial tumor samples
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

.libPaths("/home/hnatri/R/4.1_libs")

library(Seurat)
library(infercnv)

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Prepare files for infercnv
#==============================================================================#

tumors <- readRDS("/labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7.rds")
Idents(tumors) <- tumors$celltype

cellns <- as.data.frame(table(tumors$UPN, tumors$celltype))

# Need to rename cells to avoid ones starting with numbers
tumors <- RenameCells(object = tumors, add.cell.id = "tumor")
tumors$UPN_celltype <- paste0(tumors$UPN, "_", tumors$celltype)

#saveRDS(tumors, "/labs/banovich/BCTCSF/Heini/tumor_integrated_UPN109pre_soupX_snn_metadata.rds")

#==============================================================================#
# infercnv analysis
#==============================================================================#

# A function for creating the inputs for infercnv for each cluster
make_data <- function(celltype){
  message("Creating inputs for celltype ", celltype)
  #tumor_subset <- subset(tumors, subset = cluster == cluster)
  selected <- WhichCells(tumors, idents = celltype)
  tumor_subset <- subset(tumors, cells = selected)
  
  #if(cluster==5){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = c(301, 315, 239), invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  #if(cluster==10){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = c(131, 239, 266, 248, 228, 350), invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  #if(cluster==12){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = 303, invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  
  # Create a cellAnnotations file
  cellAnnotations <- as.data.frame(cbind(rownames(tumor_subset@meta.data),
                                         tumor_subset@meta.data$UPN_celltype))
  
  colnames(cellAnnotations) <- NULL
  cellAnnotations[,1] <- gsub("-", "_", cellAnnotations[,1])
  rownames(cellAnnotations) <- cellAnnotations[,1]
  cellAnnotations[,1] <- NULL
  write.table(cellAnnotations, file=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_cellAnnotations.txt"),
              sep= "\t", quote = F)
  
  # Extract out the count matrix
  singleCell.counts.matrix <- GetAssayData(tumor_subset, slot="counts", assay = "RNA")
  
  # Change hyphen to underscore (infercnv will give erorrs if there's hyphen in cellbarcode)
  colnames(singleCell.counts.matrix) <- gsub("-", "_", colnames(singleCell.counts.matrix))
  
  write.table(round(singleCell.counts.matrix, digits=3), 
              file=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_singleCell.counts.matrix"), quote=F, sep="\t")
}

# A function for running infercnv for each cluster
run_infercnv <- function(celltype){
  message("Running infercnv for celltype ", celltype)
  # Create the infercnv object
  # Note: need to save all the related files in the working directory
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_singleCell.counts.matrix"),
                                       annotations_file=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_cellAnnotations.txt"),
                                       delim="\t",
                                       gene_order_file="/home/hnatri/CART/Tumors/gene_ordering_file.txt", # this file has geneID, chr, start, end
                                       ref_group_names=NULL,
                                       min_max_counts_per_cell=c(100,Inf))
  
  # Perform infercnv operations to reveal CNV signal
  infercnv_obj <- infercnv::run(infercnv_obj,
                                min_cells_per_gene = 10, #default is 3
                                cutoff=0.1,  # expressed genes, use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_HMMi3_samplelevel_output/"),  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # group cells by UPN
                                denoise=T,
                                num_threads=8,
                                analysis_mode='samples', # c("samples", "subclusters", "cells")
                                no_plot=T,
                                HMM=T,
                                HMM_type="i3", # 3 states (neutral, depletion, duplication)
                                hclust_method='ward.D2',
                                sd_amplifier=3,  # sets midpoint for logistic
                                noise_logistic=T, # turns gradient filtering on
                                #num_ref_groups=2, # no normal cells defined
                                tumor_subcluster_partition_method='random_trees')
  #k_obs_groups = 5) # split observation groups into 5 to find malignant + normal cells
  
  # Saving the infercnv result object
  message("Saving the infercnv result object")
  saveRDS(infercnv_obj, file = paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_HMMi3_samplelevel.rds"))
  
  # Plotting
  message("Plotting")
  setwd(paste0("/scratch/hnatri/CART/infercnv/"))
  plot_cnv(infercnv_obj,
           out_dir=".",
           obs_title="Observations (Cells)",
           ref_title="References (Cells)",
           cluster_by_groups=T,
           x.center=1,
           x.range="auto",
           hclust_method='ward.D2',
           # color_safe_pal=TRUE, #using a color blindness safe palette
           output_filename=paste0("tumors_", celltype, "_HMMi3_samplelevel_res"),
           output_format="pdf",
           #k_obs_groups = 5,
           dynamic_resize=0)
  
  message("Adding results to the Seurat object")
  
  # Adjust cell barcode in seurat object to match with inferncv (inferncv uses "_", not "-")
  selected <- WhichCells(tumors, idents = celltype)
  tumor_subset <- subset(tumors, cells = selected)
  #if(cluster==5){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = c(301, 315, 239), invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  #if(cluster==10){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = c(131, 239, 266, 248, 228, 350), invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  #if(cluster==12){
  #  Idents(tumor_subset) <- tumor_subset$UPN
  #  selected <- WhichCells(tumor_subset, idents = 303, invert = T)
  #  tumor_subset <- subset(tumor_subset, cells = selected)
  #}
  rownames(tumor_subset@meta.data) <- gsub("-", "_", rownames(tumor_subset@meta.data))
  colnames(tumor_subset@assays$RNA@counts) <- rownames(tumor_subset@meta.data)
  colnames(tumor_subset@assays$RNA@data) <- rownames(tumor_subset@meta.data)
  
  # Add HMM info into Seurat object 
  tumor_subset <- infercnv::add_to_seurat(infercnv_output_path=paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_HMMi3_samplelevel_output/"),
                                          seurat_obj=tumor_subset, 
                                          top_n=10)
  
  # Adjust cell barcodes to match with other assays in the Seurat object
  #substring(rownames(tumors@meta.data), 24, 24) <- "-"
  colnames(tumor_subset@assays$RNA@counts) <- rownames(tumor_subset@meta.data)
  colnames(tumor_subset@assays$RNA@data) <- rownames(tumor_subset@meta.data)
  tumor_subset[["integrated_sct_umap"]] <- RenameCells(
    object = tumor_subset[["integrated_sct_umap"]],
    new.names = colnames(tumor_subset@assays$RNA@counts))
  
  saveRDS(tumor_subset, file = paste0("/scratch/hnatri/CART/infercnv/tumors_", celltype, "_infercnv_HMMi3_samplelevel_seurat.rds"))
  
  message("Finnished running for celltype ", celltype)
}

# Running for all clusters
for(celltype in unique(tumors$celltype)){
  make_data(celltype)
  run_infercnv(celltype)
  message("Done running for ", celltype)
}
