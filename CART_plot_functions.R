#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/02/2022
# Description: CAR T plotting scripts
#==============================================================================#

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(tidyr)

#==============================================================================
# Barplot
#==============================================================================

# Inputs
# seurat_object = Seurat object
# group_var = e.g. CD3_status
# group1 = e.g. High
# group2 = e.g. Low
# plot_var = e.g. celltype
# plot_colors = A named vector of colors to use in the plot, corresponding to plot_var (e.g. cell type colors)
# var_names = Used as axis titles, c("group2", "group1")
# legend_title = Legend title, ("" for no title)
create_clusterpropplot <- function(seurat_object,
                                   group_var,
                                   group1,
                                   group2,
                                   plot_var,
                                   plot_colors,
                                   var_names,
                                   legend_title){
  prop_table <- as.data.frame(table(seurat_object@meta.data[,plot_var], as.character(seurat_object@meta.data[,group_var])))
  colnames(prop_table) <- c(plot_var, group_var, "Freq")
  prop_table <- spread(prop_table, plot_var, Freq)
  # Converting to percentage
  prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
  prop_table <- gather(prop_table, plot_var, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)
  prop_table <- spread(prop_table, group_var, Freq)
  colnames(prop_table) <- gsub("plot_var", plot_var, colnames(prop_table))
  
  celltype_prop_scatter_1 <- ggplot(prop_table, aes_string(x = group1, y = group2, color = plot_var, label = plot_var)) +
    geom_point() + 
    scale_color_manual(name = plot_var, values = plot_colors) + 
    theme_bw() +
    xlab(var_names[1]) + 
    ylab(var_names[2]) + 
    ylim(0, max(prop_table[,2:3])+2) + 
    xlim(0, max(prop_table[,2:3])+2) +
    geom_text_repel(size = 2.5) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "gray88")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))
  
  celltype_prop_scatter_1
  
  celltype_prop_scatter_1 <- celltype_prop_scatter_1 + coord_fixed(ratio = 1)
  celltype_prop_scatter_1
}

# Inputs
# seurat_object = Seurat object
# plot_var = e.g. cluster
# group_var = e.g. response (cell proportions are plotted for each group)
# group_levels = a vector for ordering the grouping variable levels
# plot_levels = a vector for ordering the plot variable levels
# plot_colors = A named vector of colors to use in the plot, corresponding to plot_var (e.g. cluster colors)
# var_names = Used as axis titles, c("plot_var_name", "group_var_name")
# legend_title = Legend title, ("" for no title)
create_barplot <- function(seurat_object, plot_var, group_var, group_levels, plot_levels, plot_colors, var_names, legend_title){
    prop_table = as.data.frame(table(seurat_object@meta.data[,plot_var], as.character(seurat_object@meta.data[,group_var])))
    colnames(prop_table) = c("plot_var", "group_var", "Freq")
    prop_table = spread(prop_table, plot_var, Freq)
    # Converting to percetange
    prop_table[,2:length(prop_table)] = (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
    prop_table = gather(prop_table, plot_var, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)
    
    prop_table$group_var <- factor(prop_table$group_var, levels = group_levels)
    
    if(!is.null(plot_levels)){
        prop_table$plot_var <- factor(prop_table$plot_var, levels = plot_levels)
    }
    
    barplot = ggplot(prop_table, aes(x=group_var, y = Freq, fill = plot_var)) +
        geom_bar(stat="identity", position='stack', width = 0.8) +
        scale_fill_manual(name = legend_title, values = plot_colors) +
        xlab(var_names[2]) +
        ylab(var_names[1]) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    barplot
}

# Inputs:
# seurat_object = Seurat object with all features normalized and scaled
# plot_features = a vector a features to plot
# group_var = e.g. cluster
# group_colors = e.g. leukPBMC_cluster_col
# column_title = plot title
create_dotplot_heatmap <- function(seurat_object,
                                   plot_features,
                                   group_var,
                                   group_colors,
                                   column_title,
                                   km=5,
                                   row.order = NULL){
    p <- DotPlot(object = seurat_object, features = plot_features, group.by = group_var)

    # Reproducing the dotplot using ComplexHeatmap
    df <- p$data
    
    # Removing NaN values
    if(nrow(df[which(is.nan(df$avg.exp.scaled)),])>0){
      df <- df[-which(is.nan(df$avg.exp.scaled)),]
    }

    # (Scaled) expression levels
    exp_mat <- df %>% 
        dplyr::select(-pct.exp, -avg.exp) %>%  
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame()
    
    row.names(exp_mat) <- exp_mat$features.plot
    exp_mat <- exp_mat[,-1] %>% as.matrix()
    
    # The percentage of cells expressing a feature
    percent_mat <- df %>% 
        dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame() 
    
    row.names(percent_mat) <- percent_mat$features.plot  
    percent_mat <- percent_mat[,-1] %>% as.matrix()
    
    # Any value that is greater than 2 will be mapped to yellow
    #col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
    col_fun = circlize::colorRamp2(c(-1, 2), c("gray98", "tomato3"))
    
    # Adding annotation colors
    cluster_anno <- colnames(exp_mat)
    dotplot_col <- list("cluster" = group_colors)
    column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
        cluster = cluster_anno,
        col = dotplot_col,
        na_col = "grey"
    )
    
    layer_fun = function(j, i, x, y, w, h, fill){
        grid.rect(x = x, y = y, width = w, height = h, 
                  gp = gpar(col = NA, fill = NA))
        grid.points(x = x, y = y, 
                    gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
                    size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
                    pch = 16
        )
    }
    
    
    lgd_list = list(
        Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
                legend_gp = gpar(col = "black")))
    
    # Row order
    if (!is.null(row.order)){
      # row.order
      row.order <- intersect(row.order, rownames(exp_mat))
      exp_mat[row.order,]
      
      exp_mat <- exp_mat[match(row.order, rownames(exp_mat)),]
      
      hp <- ComplexHeatmap::Heatmap(exp_mat, # t(exp_mat)
                                    heatmap_legend_param=list(title="Scaled expression"),
                                    column_title = column_title, 
                                    col = col_fun,
                                    rect_gp = gpar(type = "none"),
                                    layer_fun = layer_fun,
                                    row_names_gp = gpar(fontsize = 7),
                                    column_names_gp = gpar(fontsize = 7),
                                    cluster_rows = F,
                                    column_km = km,
                                    top_annotation = column_ha, # top_annotation / left_annotation
                                    border = "black")
    } else {
    
    # NAs?
    hp <- ComplexHeatmap::Heatmap(exp_mat, # t(exp_mat)
                  heatmap_legend_param=list(title="Scaled expression"),
                  column_title = column_title, 
                  col = col_fun,
                  rect_gp = gpar(type = "none"),
                  layer_fun = layer_fun,
                  row_names_gp = gpar(fontsize = 7),
                  column_names_gp = gpar(fontsize = 7),
                  #row_km = 4,
                  column_km = km,
                  top_annotation = column_ha, # top_annotation / left_annotation
                  border = "black")
    }
    
    dotplot <- draw(hp, annotation_legend_list = lgd_list)
    
    dotplot
}


create_dotplot_heatmap_horizontal <- function(seurat_object, plot_features, group_var, group_colors, column_title, km=5, col.order = NULL){
    p <- DotPlot(object = seurat_object, features = plot_features, group.by = group_var)
    
    # Reproducing the dotplot using ComplexHeatmap
    df <- p$data
    
    # (Scaled) expression levels
    exp_mat <- df %>% 
        dplyr::select(-pct.exp, -avg.exp) %>%  
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame()
    
    row.names(exp_mat) <- exp_mat$features.plot
    exp_mat <- exp_mat[,-1] %>% as.matrix()
    #exp_mat <- exp_mat[sort(plot_features),]
    
    # The percentage of cells expressing a feature
    percent_mat <- df %>% 
        dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame() 
    
    row.names(percent_mat) <- percent_mat$features.plot  
    percent_mat <- percent_mat[,-1] %>% as.matrix()
    
    # Any value that is greater than 2 will be mapped to yellow
    #col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])
    col_fun = circlize::colorRamp2(c(-1, 2), c("gray98", "tomato3"))
    
    # Adding annotation colors
    cluster_anno <- colnames(exp_mat)
    dotplot_col <- list("cluster" = group_colors)
    column_ha <- rowAnnotation( # HeatmapAnnotation / rowAnnotation
        cluster = cluster_anno,
        col = dotplot_col,
        na_col = "grey"
    )
    
    layer_fun = function(j, i, x, y, w, h, fill){
        grid.rect(x = x, y = y, width = w, height = h, 
                  gp = gpar(col = NA, fill = NA))
        grid.points(x = x, y = y, 
                    gp = gpar(col = col_fun(pindex(t(exp_mat), i, j))), # t(exp_mat)
                    size = pindex(t(percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
                    pch = 16
        )
    }
    
    
    lgd_list = list(
        Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
                legend_gp = gpar(col = "black")))
    
    exp_mat <- exp_mat[complete.cases(exp_mat),]
    
    hp <- ComplexHeatmap::Heatmap(t(exp_mat), # t(exp_mat)
                  heatmap_legend_param=list(title="Scaled expression"),
                  column_title = column_title, 
                  col = col_fun,
                  rect_gp = gpar(type = "none"),
                  layer_fun = layer_fun,
                  row_names_gp = gpar(fontsize = 7),
                  column_names_gp = gpar(fontsize = 7),
                  row_km = km,
                  #column_km = 4,
                  left_annotation = column_ha, # top_annotation / left_annotation
                  border = "black")
    
    dotplot <- draw(hp, annotation_legend_list = lgd_list)
    
    dotplot
}

