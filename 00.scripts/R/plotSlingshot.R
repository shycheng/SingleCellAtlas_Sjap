plot_slingshot <- function(sds, show_arrow, curve_range = c(0, 1)) {
  library(slingshot)
  library(RColorBrewer)
  library(ggpubr)
  library(dplyr)
  
  test <- data.frame(reducedDims(sds))
  if (length(sds@lineages) > 1) {
    pdt_agg1 = apply(slingPseudotime(sds)[,c(1:length(sds@lineages))], 1, function(t) mean(t, na.rm = T))
  } else {
    pdt_agg1 <- slingPseudotime(sds)
  }
  test$pseudotime <- pdt_agg1
  
  p <- ggscatter(test, x = 'UMAP_1', y = 'UMAP_2', color = 'pseudotime', size = 0.8) + 
    gradient_color(palette = rev(brewer.pal(11, 'Spectral')[-6])) +
    theme_pubr(border = T)
  
  curves <- slingCurves(sds, as.df = TRUE)
  arrow_style <- arrow(type = "closed", length = unit(0.1, "inches"))
  
  if (show_arrow) {
    # 为每个lineage分别处理curve范围
    curves_filtered <- curves %>%
      arrange(Order) %>%
      group_by(Lineage) %>%
      mutate(
        # 计算每个lineage内的相对位置 (0到1之间)
        relative_pos = (row_number() - 1) / (n() - 1),
        # 根据curve_range参数过滤
        keep = relative_pos >= curve_range[1] & relative_pos <= curve_range[2]
      ) %>%
      filter(keep) %>%
      ungroup()
    
    p <- p + geom_path(data = curves_filtered,
                       aes(group = Lineage), show.legend = T, arrow = arrow_style) 
  }
  p
}


plot_trajectory_on_background <- function(sce, 
                                          full_seurat_obj,
                                          lineage_idx = 1,
                                          point_size = 0.5,
                                          show_arrow = TRUE,
                                          arrow_size = 0.1,
                                          line_width = 1.0,
                                          curve_range = c(0,1)) {
  
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(SingleCellExperiment)
  
  # 1. 准备绘图数据框
  # 从完整的Seurat对象中获取所有细胞的UMAP坐标
  all_coords <- as.data.frame(Embeddings(full_seurat_obj, reduction = "umap"))
  all_coords$cell_id <- rownames(all_coords)
  
  # 从sce对象中获取轨迹细胞的伪时间
  traj_cells <- colnames(sce)
  pdt <- slingPseudotime(sce)
  
  # 检查并选择正确的lineage
  if (ncol(pdt) > 1 && lineage_idx > ncol(pdt)) {
    warning("lineage_idx is out of bounds. Using the first lineage.")
    lineage_idx <- 1
  }
  
  pdt_df <- data.frame(
    cell_id = traj_cells,
    pseudotime = pdt[, lineage_idx]
  )
  
  # 2. 合并数据框
  # 将所有细胞的坐标与轨迹细胞的伪时间合并
  # 非轨迹细胞的 pseudotime 值将为 NA
  plot_df <- left_join(all_coords, pdt_df, by = "cell_id")
  
  # 3. 创建ggplot对象
  # 首先绘制所有点，NA值（背景细胞）将被设置为灰色
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = pseudotime), size = point_size) +
    # 使用 scale_color_gradientn 来更好地控制颜色，并将NA值设为灰色
    scale_color_gradientn(
      colors = rev(brewer.pal(11, 'Spectral')[-6]),
      na.value = "grey80", # <-- 这里是关键！
      name = "Pseudotime"
    ) +
    theme_classic() +
    labs(title = paste("Lineage", lineage_idx, "on Full UMAP")) +
    coord_equal() # 保持UMAP的纵横比
  
  if (show_arrow) {
    # 4. 添加轨迹曲线
    curves< - slingCurves(sce, as.df = TRUE)
    if ("Lineage" %in% colnames(curves)) {
      curves_selected <- curves[curves$Lineage == lineage_idx, ] 
    } else {
      curves_selected <- curves
    }
    curves_selected <- curves_selected %>% 
      arrange(Order) %>% 
      mutate(
        relative_pos = (row_number() - 1) / (n() - 1),
        keep = relative_pos >= curve_range[1] & relative_pos <= curve_range[2]
      ) %>% 
      filter(keep)
    if (nrow(curves_selected) > 0) {
      arrow_style <- arrow(type = "closed", length = unit(arrow_size, "inches"))
      
      p <- p + geom_path(data = curves_selected ,
                         aes(x = UMAP_1, y = UMAP_2), # group aesthetic is not needed here
                         color = "black", size = line_width,
                         arrow = arrow_style)
    }
  }
  return(p)
}

# 多基因版本的plotSmoothers函数
plotSmoothers_multiGenes <- function(models, counts, genes, 
                                     nPoints = 100, lwd = 2, size = 2/3,
                                     xlab = "Pseudotime", ylab = "Log(expression + 1)",
                                     border = TRUE, alpha = 1, sample = 1,
                                     pointCol = NULL, curvesCols = NULL,
                                     plotLineages = TRUE, lineagesToPlot = NULL,
                                     same_plot = FALSE, lineage = 1,
                                     ncol = NULL, scales = "free_y") {
  
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)
  library(RColorBrewer)
  
  # 参数检查
  if (is.null(genes)) {
    stop("请指定要绘制的基因列表")
  }
  
  # 确保基因存在
  available_genes <- intersect(genes, rownames(counts))
  if (length(available_genes) == 0) {
    stop("没有找到指定的基因")
  }
  if (length(available_genes) < length(genes)) {
    warning(paste("以下基因未找到:", setdiff(genes, available_genes)))
  }
  genes <- available_genes
  
  # 设置列数
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(length(genes)))
  }
  
  # 从models对象中提取基础信息
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  
  # 检查是否有conditions
  conditions <- suppressWarnings(!is.null(models$tradeSeq$conditions))
  if (conditions) {
    conditionsData <- colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditionsData)
  }
  
  if (same_plot) {
    # 在同一个图中显示多个基因（指定lineage）
    
    if (lineage > nCurves) {
      stop(paste("指定的lineage", lineage, "超出范围，最大为", nCurves))
    }
    
    # 构建时间变量
    lcol <- timeAll <- rep(0, nrow(dm))
    if (!conditions) {
      for (ii in seq_len(nrow(dm))) {
        if (dm[ii, paste0("l", lineage)] == 1) {
          timeAll[ii] <- dm[ii, paste0("t", lineage)]
          lcol[ii] <- lineage
        }
      }
    } else {
      # 处理有conditions的情况
      for (kk in seq_len(nConditions)) {
        for (ii in seq_len(nrow(dm))) {
          if (dm[ii, paste0("l", lineage, "_", kk)] == 1) {
            timeAll[ii] <- dm[ii, paste0("t", lineage)]
            lcol[ii] <- paste0(lineage, "_", kk)
          }
        }
      }
    }
    
    # 设置基因颜色
    if (is.null(curvesCols)) {
      curvesCols <- RColorBrewer::brewer.pal(max(3, length(genes)), "Set1")[1:length(genes)]
    }
    
    # 为每个基因收集数据
    allPlotData <- data.frame()
    allPredData <- data.frame()
    
    for (i in seq_along(genes)) {
      gene <- genes[i]
      
      # 检查基因是否在models中
      if (!gene %in% names(models)) {
        warning(paste("基因", gene, "不在models对象中，跳过"))
        next
      }
      
      id <- which(names(models) %in% gene)
      y <- unname(counts[names(models),][id,])
      beta <- betaMat[id,]
      
      # 原始数据点
      df <- data.frame("time" = timeAll,
                       "gene_count" = y,
                       "lineage" = as.character(lcol),
                       "gene" = gene)
      
      # 筛选该lineage的细胞
      valid_cells <- which(timeAll > 0)  # 只保留分配到该lineage的细胞
      df <- df[valid_cells, ]
      
      # 采样
      if (sample < 1) {
        rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
        df <- df[rows, ]
      }
      
      allPlotData <- rbind(allPlotData, df)
      
      # 预测平滑曲线
      if (!conditions) {
        pred_df <- tradeSeq:::.getPredictRangeDf(dm, lineage, nPoints = nPoints)
        Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = pred_df, pseudotime = pseudotime)
      } else {
        # 为每个condition计算预测
        for (kk in seq_len(nConditions)) {
          pred_df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = lineage, conditionId = kk, nPoints = nPoints)
          Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = pred_df, pseudotime = pseudotime, conditions = conditionsData)
          
          yhat <- c(exp(t(Xdf %*% t(beta)) + pred_df$offset))
          
          pred_data <- data.frame("time" = pred_df[, paste0("t", lineage)],
                                  "gene_count" = yhat,
                                  "gene" = gene,
                                  "condition" = paste0("condition_", kk))
          allPredData <- rbind(allPredData, pred_data)
        }
        next  # 跳到下一个基因
      }
      
      yhat <- c(exp(t(Xdf %*% t(beta)) + pred_df$offset))
      
      pred_data <- data.frame("time" = pred_df[, paste0("t", lineage)],
                              "gene_count" = yhat,
                              "gene" = gene)
      allPredData <- rbind(allPredData, pred_data)
    }
    
    # 创建图形
    p <- ggplot()
    
    # 添加原始数据点
    if (nrow(allPlotData) > 0) {
      p <- p + geom_point(data = allPlotData,
                          aes(x = time, y = log1p(gene_count)),
                          color = "lightgray", size = size, alpha = alpha)
    }
    
    # 添加平滑曲线
    if (nrow(allPredData) > 0) {
      if (!conditions) {
        p <- p + geom_line(data = allPredData,
                           aes(x = time, y = log1p(gene_count), color = gene),
                           size = lwd)
      } else {
        p <- p + geom_line(data = allPredData,
                           aes(x = time, y = log1p(gene_count), color = gene, linetype = condition),
                           size = lwd)
      }
      
      p <- p + scale_color_manual(values = curvesCols, name = "Gene")
    }
    
    p <- p +
      labs(title = paste("Gene Expression in Lineage", lineage),
           x = xlab, y = ylab) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    return(p)
    
  } else {
    # 分别绘制每个基因
    
    plotList <- list()
    
    for (gene in genes) {
      # 使用原始的plotSmoothers函数为每个基因创建图
      tryCatch({
        p <- plotSmoothers(models = models,
                           counts = counts,
                           gene = gene,
                           nPoints = nPoints,
                           lwd = lwd,
                           size = size,
                           xlab = xlab,
                           ylab = ylab,
                           border = border,
                           alpha = alpha,
                           sample = sample,
                           pointCol = pointCol,
                           curvesCols = curvesCols,
                           plotLineages = plotLineages,
                           lineagesToPlot = lineagesToPlot) +
          ggtitle(gene) +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
        
        # 如果是多基因，移除图例以节省空间
        if (length(genes) > 1) {
          p <- p + theme(legend.position = "none")
        }
        
        plotList[[gene]] <- p
        
      }, error = function(e) {
        warning(paste("绘制基因", gene, "时出错:", e$message))
      })
    }
    
    # 组合图形
    if (length(plotList) == 1) {
      finalPlot <- plotList[[1]] + theme(legend.position = "right")
    } else if (length(plotList) > 1) {
      # 使用patchwork组合
      finalPlot <- wrap_plots(plotList, ncol = ncol)
      
      # 添加共同图例
      if (plotLineages && nCurves > 1) {
        # 创建一个带图例的图
        legendPlot <- plotList[[1]] + theme(legend.position = "bottom")
        
        # 提取图例
        library(cowplot)
        legend <- get_legend(legendPlot)
        
        # 组合图形和图例
        finalPlot <- finalPlot / legend + plot_layout(heights = c(1, 0.1))
      }
    } else {
      stop("没有成功绘制任何基因")
    }
    
    return(finalPlot)
  }
}
plotSmoothers_multiGenes_enhanced <- function(models, counts, genes, 
                                              nPoints = 100, lwd = 2, size = 2/3,
                                              xlab = "Pseudotime", ylab = "Log(expression + 1)",
                                              border = TRUE, alpha = 1, sample = 1,
                                              pointCol = NULL, curvesCols = NULL,
                                              plotLineages = TRUE, lineagesToPlot = NULL,
                                              same_plot = FALSE, lineage = 1,
                                              ncol = NULL, scales = "free_y",
                                              # 新增参数
                                              color_by_celltype = FALSE,
                                              celltype_column = "Cell_Type",
                                              celltype_colors = NULL,
                                              show_celltype_legend = TRUE,
                                              point_alpha = 0.6) {
  
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)
  library(RColorBrewer)
  
  # 参数检查
  if (is.null(genes)) {
    stop("请指定要绘制的基因列表")
  }
  
  # 确保基因存在
  available_genes <- intersect(genes, rownames(counts))
  if (length(available_genes) == 0) {
    stop("没有找到指定的基因")
  }
  if (length(available_genes) < length(genes)) {
    warning(paste("以下基因未找到:", setdiff(genes, available_genes)))
  }
  genes <- available_genes
  
  # 设置列数
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(length(genes)))
  }
  
  # 从models对象中提取基础信息
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  
  # 检查细胞类型信息
  if (color_by_celltype) {
    if (celltype_column %in% colnames(colData(models))) {
      celltype_data <- colData(models)[[celltype_column]]
    } else {
      warning(paste("列", celltype_column, "不存在于colData中，将使用默认着色"))
      color_by_celltype <- FALSE
      celltype_data <- NULL
    }
    
    # 设置细胞类型颜色
    if (color_by_celltype && is.null(celltype_colors)) {
      unique_types <- unique(celltype_data)
      n_types <- length(unique_types)
      if (n_types <= 8) {
        celltype_colors <- RColorBrewer::brewer.pal(max(3, n_types), "Set2")[1:n_types]
      } else {
        celltype_colors <- rainbow(n_types)
      }
      names(celltype_colors) <- unique_types
    }
  }
  
  # 检查是否有conditions
  conditions <- suppressWarnings(!is.null(models$tradeSeq$conditions))
  if (conditions) {
    conditionsData <- colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditionsData)
  }
  
  if (same_plot) {
    # 在同一个图中显示多个基因（指定lineage）
    
    if (lineage > nCurves) {
      stop(paste("指定的lineage", lineage, "超出范围，最大为", nCurves))
    }
    
    # 构建时间变量
    lcol <- timeAll <- rep(0, nrow(dm))
    if (!conditions) {
      for (ii in seq_len(nrow(dm))) {
        if (dm[ii, paste0("l", lineage)] == 1) {
          timeAll[ii] <- dm[ii, paste0("t", lineage)]
          lcol[ii] <- lineage
        }
      }
    } else {
      # 处理有conditions的情况
      for (kk in seq_len(nConditions)) {
        for (ii in seq_len(nrow(dm))) {
          if (dm[ii, paste0("l", lineage, "_", kk)] == 1) {
            timeAll[ii] <- dm[ii, paste0("t", lineage)]
            lcol[ii] <- paste0(lineage, "_", kk)
          }
        }
      }
    }
    
    # 设置基因颜色（用于曲线）
    if (is.null(curvesCols)) {
      curvesCols <- RColorBrewer::brewer.pal(max(3, length(genes)), "Set1")[1:length(genes)]
    }
    
    # 为每个基因收集数据
    allPlotData <- data.frame()
    allPredData <- data.frame()
    
    for (i in seq_along(genes)) {
      gene <- genes[i]
      
      # 检查基因是否在models中
      if (!gene %in% names(models)) {
        warning(paste("基因", gene, "不在models对象中，跳过"))
        next
      }
      
      id <- which(names(models) %in% gene)
      y <- unname(counts[names(models),][id,])
      beta <- betaMat[id,]
      
      # 原始数据点
      df <- data.frame("time" = timeAll,
                       "gene_count" = y,
                       "lineage" = as.character(lcol),
                       "gene" = gene)
      
      # 添加细胞类型信息
      if (color_by_celltype) {
        df$celltype <- celltype_data
      }
      
      # 筛选该lineage的细胞
      valid_cells <- which(timeAll > 0)  # 只保留分配到该lineage的细胞
      df <- df[valid_cells, ]
      
      # 采样
      if (sample < 1) {
        rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
        df <- df[rows, ]
      }
      
      allPlotData <- rbind(allPlotData, df)
      
      # 预测平滑曲线
      if (!conditions) {
        pred_df <- tradeSeq:::.getPredictRangeDf(dm, lineage, nPoints = nPoints)
        Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = pred_df, pseudotime = pseudotime)
      } else {
        # 为每个condition计算预测
        for (kk in seq_len(nConditions)) {
          pred_df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = lineage, conditionId = kk, nPoints = nPoints)
          Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = pred_df, pseudotime = pseudotime, conditions = conditionsData)
          
          yhat <- c(exp(t(Xdf %*% t(beta)) + pred_df$offset))
          
          pred_data <- data.frame("time" = pred_df[, paste0("t", lineage)],
                                  "gene_count" = yhat,
                                  "gene" = gene,
                                  "condition" = paste0("condition_", kk))
          allPredData <- rbind(allPredData, pred_data)
        }
        next  # 跳到下一个基因
      }
      
      yhat <- c(exp(t(Xdf %*% t(beta)) + pred_df$offset))
      
      pred_data <- data.frame("time" = pred_df[, paste0("t", lineage)],
                              "gene_count" = yhat,
                              "gene" = gene)
      allPredData <- rbind(allPredData, pred_data)
    }
    
    # 创建图形
    p <- ggplot()
    
    # 添加原始数据点
    if (nrow(allPlotData) > 0) {
      if (color_by_celltype) {
        # 按细胞类型着色
        p <- p + geom_point(data = allPlotData,
                            aes(x = time, y = log1p(gene_count), color = celltype),
                            size = size, alpha = point_alpha) +
          scale_color_manual(values = celltype_colors, name = "Cell Type")
      } else {
        # 默认灰色
        p <- p + geom_point(data = allPlotData,
                            aes(x = time, y = log1p(gene_count)),
                            color = "lightgray", size = size, alpha = point_alpha)
      }
    }
    
    # 添加平滑曲线
    if (nrow(allPredData) > 0) {
      if (!conditions) {
        if (color_by_celltype) {
          # 使用新的颜色映射避免与细胞类型冲突
          p <- p + 
            new_scale_color() +
            geom_line(data = allPredData,
                      aes(x = time, y = log1p(gene_count), color = gene),
                      size = lwd) +
            scale_color_manual(values = curvesCols, name = "Gene")
        } else {
          p <- p + geom_line(data = allPredData,
                             aes(x = time, y = log1p(gene_count), color = gene),
                             size = lwd) +
            scale_color_manual(values = curvesCols, name = "Gene")
        }
      } else {
        p <- p + geom_line(data = allPredData,
                           aes(x = time, y = log1p(gene_count), color = gene, linetype = condition),
                           size = lwd) +
          scale_color_manual(values = curvesCols, name = "Gene")
      }
    }
    
    p <- p +
      labs(title = paste("Gene Expression in Lineage", lineage),
           x = xlab, y = ylab) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    # 调整图例位置
    if (color_by_celltype && show_celltype_legend) {
      p <- p + theme(legend.position = "right")
    }
    
    return(p)
    
  } else {
    # 分别绘制每个基因
    
    plotList <- list()
    
    for (gene in genes) {
      # 修改原始plotSmoothers调用，支持细胞类型着色
      tryCatch({
        if (color_by_celltype) {
          # 使用细胞类型着色
          p <- plotSmoothers(models = models,
                             counts = counts,
                             gene = gene,
                             nPoints = nPoints,
                             lwd = lwd,
                             size = size,
                             xlab = xlab,
                             ylab = ylab,
                             border = border,
                             alpha = point_alpha,
                             sample = sample,
                             pointCol = celltype_column,  # 使用细胞类型列
                             curvesCols = curvesCols,
                             plotLineages = plotLineages,
                             lineagesToPlot = lineagesToPlot)
          
          # 自定义颜色
          if (!is.null(celltype_colors)) {
            p <- p + scale_color_manual(values = celltype_colors, name = "Cell Type")
          }
        } else {
          # 使用默认着色
          p <- p + plotSmoothers(models = models,
                                 counts = counts,
                                 gene = gene,
                                 nPoints = nPoints,
                                 lwd = lwd,
                                 size = size,
                                 xlab = xlab,
                                 ylab = ylab,
                                 border = border,
                                 alpha = alpha,
                                 sample = sample,
                                 pointCol = pointCol,
                                 curvesCols = curvesCols,
                                 plotLineages = plotLineages,
                                 lineagesToPlot = lineagesToPlot)
        }
        
        p <- p + 
          ggtitle(gene) +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
        
        # 如果是多基因且不显示细胞类型图例，移除图例以节省空间
        if (length(genes) > 1 && (!color_by_celltype || !show_celltype_legend)) {
          p <- p + theme(legend.position = "none")
        }
        
        plotList[[gene]] <- p
        
      }, error = function(e) {
        warning(paste("绘制基因", gene, "时出错:", e$message))
      })
    }
    
    # 组合图形
    if (length(plotList) == 1) {
      finalPlot <- plotList[[1]] + theme(legend.position = "right")
    } else if (length(plotList) > 1) {
      # 使用patchwork组合
      finalPlot <- wrap_plots(plotList, ncol = ncol)
      
      # 添加共同图例
      if (color_by_celltype && show_celltype_legend) {
        # 创建一个带细胞类型图例的图
        legendPlot <- plotList[[1]] + theme(legend.position = "bottom")
        
        # 提取图例
        library(cowplot)
        legend <- get_legend(legendPlot)
        
        # 组合图形和图例
        finalPlot <- finalPlot / legend + plot_layout(heights = c(1, 0.15))
      } else if (plotLineages && nCurves > 1) {
        # 添加lineage图例
        legendPlot <- plotList[[1]] + theme(legend.position = "bottom")
        library(cowplot)
        legend <- get_legend(legendPlot)
        finalPlot <- finalPlot / legend + plot_layout(heights = c(1, 0.1))
      }
    } else {
      stop("没有成功绘制任何基因")
    }
    
    return(finalPlot)
  }
}

# 如果需要处理双重颜色映射，安装ggnewscale包
# install.packages("ggnewscale")
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  warning("推荐安装ggnewscale包以支持双重颜色映射: install.packages('ggnewscale')")
  new_scale_color <- function() return(NULL)
} else {
  library(ggnewscale)
}


