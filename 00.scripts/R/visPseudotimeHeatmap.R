smoothMatrix <- function(to_plot){
  pt.matrix <- to_plot
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- rownames(to_plot)
  colnames(pt.matrix) <- colnames(to_plot)
  pt.matrix
}

getGroupMat <- function(
    sce_obj,
    counts = counts,
    Pseudotime = "slingPseudotime_1",
    select_genes = NULL,
    time_reverse = FALSE
){
  
  ptime <- c(sce_obj[[Pseudotime]])
  # 1.1 get cells in that lineage
  lineage_cells <- colnames(sce_obj)[!is.na(ptime)]
  # 1.2 remove values for cells not in the lineage
  ptime <- ptime[!is.na(ptime)]
  names(ptime) <- lineage_cells
  
  breaks <- quantile(ptime, probs = seq(0, 1, length.out = 101))
  # 使用cut函数创建bins
  ptime_bins <- cut(ptime, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  groupList <- split(lineage_cells,ptime_bins)
  breaks <- seq(0, 100, 1)
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  
  Y <- log2(as.matrix(counts[,lineage_cells]) + 1)
  RNA_mat <- Y[,lineage_cells]
  
  mat_group <- matrix(0, nrow = nrow(RNA_mat), ncol = length(groupList))
  
  for(z in seq_along(groupList)){
    
    #Check Cells In Group
    cellsGroupz <- groupList[[z]]
    idx <- BiocGenerics::which(colnames(RNA_mat) %in% cellsGroupz)
    
    #If In Group RowSums
    if(length(idx) > 0){
      mat_group[,z] <- mat_group[,z] + Matrix::rowSums(RNA_mat[,idx,drop=FALSE])
    }
  }
  colnames(mat_group) <- names(groupList)  
  rownames(mat_group) <- rownames(RNA_mat)
  
  if (!is.null(select_genes)) {
    mat_group <- mat_group[select_genes,]
  }
  
  if (time_reverse == TRUE) {
    mat_group <- mat_group[,rev(names(groupList))]
  }
  mat_group
}

getBinCellTypes <- function(sce_obj, Pseudotime = "slingPseudotime_1", 
                            celltype_col = "NewCellType", n_bins = 100, 
                            time_reverse = FALSE) {
  
  ptime <- c(sce_obj[[Pseudotime]])
  celltype <- c(sce_obj[[celltype_col]])
  
  # 获取有效细胞并排序
  valid_idx <- !is.na(ptime)
  lineage_cells <- colnames(sce_obj)[valid_idx]
  ptime <- ptime[valid_idx]
  celltype <- celltype[valid_idx]
  
  # **按拟时间排序**
  ptime_order <- order(ptime)
  ptime <- ptime[ptime_order]
  celltype <- celltype[ptime_order]
  lineage_cells <- lineage_cells[ptime_order]
  
  # 创建bins
  breaks <- quantile(ptime, probs = seq(0, 1, length.out = n_bins + 1))
  ptime_bins <- cut(ptime, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  bin_celltypes <- sapply(1:n_bins, function(i) {
    cells_in_bin <- lineage_cells[ptime_bins == i]
    if (length(cells_in_bin) == 0) return("Unknown")
    
    celltype_counts <- table(celltype[cells_in_bin])
    names(celltype_counts)[which.max(celltype_counts)]
  })
  
  # 如果需要时间反转
  if (time_reverse == TRUE) {
    bin_celltypes <- rev(bin_celltypes)
  }
  
  return(bin_celltypes)
}


getBinCellTypes_fixed <- function(sce_obj, Pseudotime = "slingPseudotime_1", 
                                  celltype_col = "fineCell_Type", n_bins = 100, 
                                  time_reverse = FALSE, debug = FALSE) {
  
  ptime <- c(sce_obj[[Pseudotime]])
  celltype <- c(sce_obj[[celltype_col]])
  
  # 获取有效数据
  valid_idx <- !is.na(ptime) & !is.na(celltype)
  ptime <- ptime[valid_idx]
  celltype <- celltype[valid_idx]
  
  if(debug) {
    cat("有效细胞数:", length(ptime), "\n")
    cat("拟时间范围:", range(ptime), "\n")
  }
  
  # 排序
  ptime_order <- order(ptime)
  ptime <- ptime[ptime_order]
  celltype <- celltype[ptime_order]
  
  # 使用等间距的拟时间分割
  ptime_min <- min(ptime)
  ptime_max <- max(ptime)
  time_points <- seq(ptime_min, ptime_max, length.out = n_bins + 1)
  
  bin_celltypes <- c()
  
  for(i in 1:n_bins) {
    # 找到在当前时间区间内的细胞
    in_bin <- ptime >= time_points[i] & ptime < time_points[i + 1]
    if(i == n_bins) {
      in_bin <- ptime >= time_points[i] & ptime <= time_points[i + 1]
    }
    
    cells_in_bin <- celltype[in_bin]
    
    if(debug && i <= 10) {
      cat("Bin", i, "时间范围: [", round(time_points[i], 2), ",", 
          round(time_points[i+1], 2), "], 细胞数:", sum(in_bin), "\n")
    }
    
    if(length(cells_in_bin) == 0) {
      # 方法1：使用最近邻插值
      # 找到最近的非空bin
      # if(i == 1) {
      #   # 如果是第一个bin，向前寻找
      #   for(j in 2:n_bins) {
      #     next_in_bin <- ptime >= time_points[j] & ptime < time_points[j + 1]
      #     if(j == n_bins) next_in_bin <- ptime >= time_points[j] & ptime <= time_points[j + 1]
      #     if(sum(next_in_bin) > 0) {
      #       next_cells <- celltype[next_in_bin]
      #       next_dominant <- names(table(next_cells))[which.max(table(next_cells))]
      #       bin_celltypes[i] <- next_dominant
      #       break
      #     }
      #   }
      # } else {
        # 使用前一个bin的结果
        bin_celltypes[i] <- bin_celltypes[i-1]
      #}
      
      # 方法2（备选）：使用周围区间的细胞进行插值
      # 扩大搜索范围到周围的区间
      # expand_range <- 2  # 扩大到前后2个区间
      # expanded_min <- max(1, i - expand_range)
      # expanded_max <- min(n_bins, i + expand_range)
      # 
      # expanded_in_bin <- ptime >= time_points[expanded_min] & 
      #                    ptime <= time_points[expanded_max + 1]
      # expanded_cells <- celltype[expanded_in_bin]
      # 
      # if(length(expanded_cells) > 0) {
      #   expanded_dominant <- names(table(expanded_cells))[which.max(table(expanded_cells))]
      #   bin_celltypes[i] <- expanded_dominant
      # } else {
      #   # 如果扩大范围仍然没有细胞，使用全局最常见的细胞类型
      #   global_dominant <- names(table(celltype))[which.max(table(celltype))]
      #   bin_celltypes[i] <- global_dominant
      # }
      
      # if(debug) {
      #   cat("  空bin，使用插值方法赋值为:", bin_celltypes[i], "\n")
      # }
      
    } else {
      # 找到最多的细胞类型
      celltype_counts <- table(cells_in_bin)
      dominant_type <- names(celltype_counts)[which.max(celltype_counts)]
      bin_celltypes[i] <- dominant_type
    }
  }
  
  if (time_reverse == TRUE) {
    bin_celltypes <- rev(bin_celltypes)
  }
  
  if(debug) {
    cat("\n结果统计:\n")
    print(table(bin_celltypes))
    cat("空bin数量:", sum(is.na(bin_celltypes)), "\n")
  }
  
  return(bin_celltypes)
}
