#' @title custom umap
#' @description a powerful and flexible function for plotting seurat object's dimensionality reduction (e.g. umap), supporting faceting, custom labels and legends.
#'
#' @param seurat_obj seurat object.
#' @param embedding dimensionality reduction method name, string, default is 'umap'.
#' @param feature metadata column name for coloring points.
#' @param split_by (optional) metadata column name for faceting.
#' @param label_feature (optional) metadata column name for displaying labels on the plot. If NULL, the value of `feature` parameter will be used.
#' @param legend_labels (optional) vector of labels for the legend.
#' @param col_vectors (optional) vector of colors. If NULL, rainbow colors will be automatically generated.
#' @param levels (optional) factor levels of the `feature` column, used to control the order of the legend.
#' @param showLabels logical value, whether to display cluster labels on the plot.
#' @param pointSize size of points on the plot.
#' @param labelSize font size of labels.
#' @param showLabelBackground logical value, whether to display a background for labels (requires `ggrepel` package).
#' @param title main title of the plot.
#' @param ncol number of columns when using `split_by`.
#' @param highlight_clusters (optional) vector of cluster names/numbers to highlight with dashed outlines.
#' @param outline_color color of the dashed outline, default is "red".
#' @param outline_size thickness of the dashed outline, default is 1.
#' @param outline_alpha transparency of the dashed outline, default is 0.8.
#'
#' @return a ggplot object.
#'
#' @import ggplot2
#' @import Seurat
#' @import dplyr
#' @importFrom utils packageVersion
#'
plot_umap_custom <- function(seurat_obj,
                             embedding = 'umap',
                             feature = 'Cell_Type',
                             split_by = NULL,
                             label_feature = NULL,
                             legend_labels = NULL,
                             col_vectors = NULL,
                             levels = NULL,
                             showLabels = TRUE,
                             pointSize = 0.05,
                             labelSize = 10,
                             l_margin = 15,
                             showLabelBackground = FALSE,
                             title = "UMAP",
                             highlight_clusters = NULL,
                             outline_color = "red",
                             outline_size = 1,
                             outline_alpha = 0.8) {
  
  # check dependencies
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function.")
  }
  
  # data preparation
  # get dimensionality reduction coordinates
  coords <- Seurat::Embeddings(seurat_obj, reduction = embedding)
  
  # build base data frame
  plot_data <- data.frame(
    DIM_1 = coords[, 1],
    DIM_2 = coords[, 2],
    feature_col = seurat_obj@meta.data[[feature]]
  )
  
  # determine column names for plotting and labeling
  dim_names <- colnames(coords)[1:2]
  
  # add faceting variable
  if (!is.null(split_by)) {
    plot_data$split_col <- seurat_obj@meta.data[[split_by]]
  }
  
  # set label column, if not specified, use the feature column
  label_feature <- label_feature %||% feature
  plot_data$label_col <- seurat_obj@meta.data[[label_feature]]
  
  # set factor levels
  if (!is.null(levels)) {
    plot_data$feature_col <- factor(plot_data$feature_col, levels = levels)
  }
  
  # color and legend preparation
  # set color vector
  num_colors <- length(unique(plot_data$feature_col))
  col_vectors <- col_vectors %||% rainbow(num_colors)
  
  # ggplot plotting
  # create base layer
  p <- ggplot(plot_data, aes(x = DIM_1, y = DIM_2, color = feature_col)) +
    geom_point(size = pointSize, alpha = 0.8) +
    labs(x = dim_names[1], y = dim_names[2], title = title)
  
  # add faceting
  if (!is.null(split_by)) {
    p <- p + facet_wrap(~split_col)
  }
  
  # custom color and legend
  p <- p +
    scale_color_manual(values = col_vectors, labels = legend_labels, name = "") +
    guides(color = guide_legend(override.aes = list(size = 5),
                                label.theme = element_text(margin = margin(l = -l_margin, unit = 'pt'))))
  
  # add dashed outlines for highlighted clusters
  if (!is.null(highlight_clusters)) {
    if (!requireNamespace("ggforce", quietly = TRUE)) {
      warning("ggforce package is required for cluster outlines. Installing ggforce or use alternative method.")
      # 使用替代方法：凸包 (convex hull)
      for (cluster in highlight_clusters) {
        cluster_data <- plot_data[plot_data$feature_col == cluster | plot_data$label_col == cluster, ]
        if (nrow(cluster_data) > 2) {
          # 计算凸包
          hull_indices <- chull(cluster_data$DIM_1, cluster_data$DIM_2)
          hull_data <- cluster_data[hull_indices, ]
          
          if (!is.null(split_by)) {
            # 如果有分面，需要为每个分面单独计算凸包
            split_groups <- unique(cluster_data$split_col)
            for (group in split_groups) {
              group_data <- cluster_data[cluster_data$split_col == group, ]
              if (nrow(group_data) > 2) {
                hull_indices <- chull(group_data$DIM_1, group_data$DIM_2)
                hull_data <- group_data[hull_indices, ]
                p <- p + geom_polygon(data = hull_data,
                                      aes(x = DIM_1, y = DIM_2),
                                      fill = NA, 
                                      color = outline_color,
                                      linetype = "dashed",
                                      size = outline_size,
                                      alpha = outline_alpha,
                                      inherit.aes = FALSE)
              }
            }
          } else {
            p <- p + geom_polygon(data = hull_data,
                                  aes(x = DIM_1, y = DIM_2),
                                  fill = NA, 
                                  color = outline_color,
                                  linetype = "dashed",
                                  size = outline_size,
                                  alpha = outline_alpha,
                                  inherit.aes = FALSE)
          }
        }
      }
    } else {
      # 使用ggforce包的geom_mark_ellipse (推荐方法)
      for (cluster in highlight_clusters) {
        cluster_data <- plot_data[plot_data$feature_col == cluster | plot_data$label_col == cluster, ]
        if (nrow(cluster_data) > 0) {
          p <- p + ggforce::geom_mark_ellipse(
            data = cluster_data,
            aes(x = DIM_1, y = DIM_2),
            color = outline_color,
            fill = NA,
            linetype = "dashed",
            size = outline_size,
            alpha = outline_alpha,
            inherit.aes = FALSE,
            expand = unit(2, "mm")
          )
        }
      }
    }
  }
  
  # theme setting
  p <- p +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = 'black', linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = 'right',
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"), # for facet titles
      strip.background = element_rect(fill = "grey90", color = "black")
    )
  
  # add labels
  if (showLabels) {
    # calculate label center positions
    group_vars <- if (!is.null(split_by)) c("split_col", "label_col") else "label_col"
    
    label_data <- plot_data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        center_x = mean(DIM_1),
        center_y = mean(DIM_2),
        .groups = 'drop'
      )
    
    # choose label type based on whether to show background
    if (showLabelBackground && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        aes(x = center_x, y = center_y, label = label_col),
        color = "black", size = labelSize / 2.5, inherit.aes = FALSE,
        bg.colour = "white", bg.r = 0.1, min.segment.length = 0
      )
    } else {
      p <- p + geom_text(
        data = label_data,
        aes(x = center_x, y = center_y, label = label_col),
        color = "black", size = labelSize / 2.5, inherit.aes = FALSE
      )
    }
  }
  
  return(p)
}


#' @title create umap with numbered labels and combined legend
#' @description a convenient wrapper function that automatically processes data (converts Cell_Type to numbers) and calls `plot_umap_custom` to plot.
#'
#' @param seurat_obj seurat object, must contain a metadata column named `Cell_Type`.
#' @param group_by metadata column name for grouping and coloring, default is "Cell_Type".
#' @param split_by (optional) metadata column name for faceting.
#' @param col_vectors (optional) vector of colors.
#' @param title main title of the plot.
#' @param highlight_clusters (optional) vector of cluster names to highlight with dashed outlines.
#' @param outline_color color of the dashed outline, default is "red".
#' @param outline_size thickness of the dashed outline, default is 1.
#' @param outline_alpha transparency of the dashed outline, default is 0.8.
#' @param ... other parameters to be passed to `plot_umap_custom` (e.g. `pointSize`, `labelSize`).
#'
#' @return a ggplot object.
#'
create_numbered_umap <- function(seurat_obj,
                                 group_by = "Cell_Type",  # 新增参数，默认为Cell_Type
                                 split_by = NULL,
                                 col_vectors = NULL,
                                 title = "UMAP (colored by cell types)",
                                 highlight_clusters = NULL,
                                 outline_color = "red",
                                 outline_size = 1,
                                 outline_alpha = 0.8,
                                 ...) {
  
  # data preparation
  # ensure the grouping variable is a factor type to ensure consistent order
  if (!is.factor(seurat_obj@meta.data[[group_by]])) {
    seurat_obj@meta.data[[group_by]] <- factor(seurat_obj@meta.data[[group_by]])
  }
  
  # create numbered labels `cls`
  seurat_obj$cls <- as.character(as.numeric(seurat_obj@meta.data[[group_by]]))
  
  # create combined legend labels
  # use forcats::fct_inorder to ensure original factor order of the grouping variable
  legend_df <- seurat_obj@meta.data %>%
    dplyr::select(!!rlang::sym(group_by), cls) %>%  # 使用rlang::sym来处理变量名
    dplyr::distinct() %>%
    dplyr::mutate(cls_num = as.numeric(cls)) %>%
    dplyr::arrange(cls_num)
  
  # 动态创建图例标签
  legend_labels <- paste0(
    # use sprintf to format numbers and ensure alignment
    sprintf("%2s", legend_df$cls), "  ", legend_df[[group_by]]
  )
  
  # 动态更新标题
  if (title == "UMAP (colored by cell types)" && group_by != "Cell_Type") {
    title <- paste0("UMAP (colored by ", group_by, ")")
  }
  
  # 如果指定了highlight_clusters，需要将原始名称转换为对应的数字标签
  highlight_numbers <- NULL
  if (!is.null(highlight_clusters)) {
    # 创建名称到数字的映射
    name_to_number <- setNames(legend_df$cls, legend_df[[group_by]])
    
    # 转换highlight_clusters中的名称为对应的数字
    highlight_numbers <- character(0)
    for (cluster in highlight_clusters) {
      if (cluster %in% names(name_to_number)) {
        # 如果是原始名称，转换为数字
        highlight_numbers <- c(highlight_numbers, name_to_number[cluster])
      } else if (cluster %in% legend_df$cls) {
        # 如果已经是数字，直接使用
        highlight_numbers <- c(highlight_numbers, cluster)
      } else {
        warning(paste("Cluster", cluster, "not found in", group_by, "column"))
      }
    }
  }
  
  # call core function to plot
  plot <- plot_umap_custom(
    seurat_obj = seurat_obj,
    embedding = 'umap',
    feature = 'cls',      # use numbered labels to color
    label_feature = 'cls',# display numbered labels on the plot
    split_by = split_by,
    legend_labels = legend_labels,
    col_vectors = col_vectors,
    levels = as.character(sort(as.numeric(unique(seurat_obj$cls)))), # ensure factor levels are ordered
    title = title,
    highlight_clusters = highlight_numbers,
    outline_color = outline_color,
    outline_size = outline_size,
    outline_alpha = outline_alpha,
    ...                   # pass additional parameters
  )
  
  return(plot)
}