#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/01/2024
# Description: ProjecTILs T cell annotation for 13384 tumor scRNAseq
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(ProjecTILs)
library(Seurat)
#library(SignatuR)
#library(UCell)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/CART/CART_plot_functions.R")
source("/home/hnatri/CART/13384_Tumors/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import Seurat objects
#==============================================================================#

# Objects with metadata
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/tumor_immune_fibroblast_reclustered.rds")

# Subsetting lymphoid and myeloid compartments
unique(immune_fibro$celltype)
lymphoid <- subset(immune_fibro, subset = celltype %in% paste0("L", seq(1, 10)))
myeloid <- subset(immune_fibro, subset = celltype %in% paste0("L", seq(1, 11)))

query <- lymphoid
DefaultAssay(query) <- "RNA"

unique(query$binary_response)

#==============================================================================#
# Run ProjecTILs
#==============================================================================#

# Loading the reference atlas
ref <- load.reference.map()

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
p0 <- DimPlot(ref, label = T, cols = refCols) + coord_fixed(ratio=1)

# Expression of some marker genes across reference subtypes
markers <- c("Cd4", "Cd8a", "Ccr7", "Tcf7", "Pdcd1", "Havcr2", "Tox", "Izumo1r",
             "Cxcr6", "Xcl1", "Gzmb", "Gzmk", "Ifng", "Foxp3")
VlnPlot(ref, features = markers, stack = T, flip = T, assay = "RNA")

# Run Projection algorithm
query.projected <- Run.ProjecTILs(query, ref = ref)

# Visualize projection
p1 <- plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5) + coord_fixed(ratio=1)

# Plot the predicted composition of the query in terms of reference T cell subtypes
p2 <- plot.statepred.composition(ref, query.projected, metric = "Percent") + NoLegend()

# Compare gene expression
plot.states.radar(ref, query = query.projected)

# Compare gene programs
programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)
query.projected <- AddModuleScore_UCell(query.projected, features = programs, assay = "RNA",
                                        name = NULL)

plot.states.radar(ref, query = query.projected, meta4radar = names(programs))

# Compare cell states across conditions
#product_tnmem <- subset(query.projected, subset = Manufacture == "Tnmem")
#product_tcm <- subset(query.projected, subset = Manufacture == "TCM")
CR_SD_PR <- subset(query.projected, subset = binary_outcome == "CR/SD/PR")
PD <- subset(query.projected, subset = binary_outcome == "PD")

#plot.states.radar(ref, query = list(Tnmem = product_tnmem, TCM = product_tcm))
plot.states.radar(ref, query = list(CR_SD_PR = CR_SD_PR, PD = PD))

# The ProjecTILs.classifier function applies reference-projection but does not
# alter the current embeddings
querydata <- ProjecTILs.classifier(query = query, ref = ref)

palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(ref$functional.cluster), "NA")
p3 <- DimPlot(querydata, group.by = "functional.cluster", cols = palette) + coord_fixed(ratio=1)

(p1 + p2 + p3)

#saveRDS(querydata, "/scratch/hnatri/CART/leukPBMC_Tcells_projecTILs.rds")

# Create barplot
palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(ref$functional.cluster), "NA")
querydata$functional.cluster <- factor(querydata$functional.cluster, levels = levels(ref$functional.cluster))
# seurat_object, plot_var, group_var, group_levels, plot_colors, var_names, legend_title
barplot <- create_barplot(querydata,
                          "functional.cluster",
                          "cluster",
                          as.character(c(0, seq(1, max(unique(as.numeric(as.character(querydata$cluster))))))),
                          palette,
                          c("% cells", "Cluster"),
                          "T cell state")



