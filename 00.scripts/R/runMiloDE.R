# runMiloDE.R - Milo Differential Expression Analysis
# =============================================================================
# Functions for neighborhood-based differential expression analysis using Milo
# =============================================================================

#' Run Milo Differential Expression Analysis
#'
#' @description Performs differential expression analysis using Milo on 
#'   single-cell data with neighborhood-based approach.
#'
#' @param sce SingleCellExperiment object
#' @param condition_id Column name in colData for condition comparison
#' @param sample_id Column name in colData for sample identification
#' @param k Number of neighbors for graph construction (default: 20)
#' @param output_file Optional path to save results RDS file
#'
#' @return List containing stat_auc, sce, and sce_milo objects
#'
#' @examples
#' \dontrun{}
#' source("R/config.R")
#' sce <- readRDS(get_path("sce_mouseGastrulation.rds", "rds"))
#' results <- run_milo_analysis(sce, condition_id = "tomato")
#' \}
#'
#' @export
run_milo_analysis <- function(sce, 
                               condition_id = "tomato",
                               sample_id = "sample",
                               k = 20,
                               output_file = NULL) {
  
  # Load required packages
  suppressPackageStartupMessages({
    library(miloDE)
    library(miloR)
    library(scuttle)
    library(scran)
    library(uwot)
    library(ggplot2)
    library(viridis)
    library(ggpubr)
  })
  
  # Add logcounts if not present
  if (!"logcounts" %in% assayNames(sce)) {
    message("Adding logcounts...")
    sce <- logNormCounts(sce)
    sce@assays@data$logcounts <- as(sce@assays@data$logcounts, "CsparseMatrix")
  }
  
  # Add UMAPs if not present
  if (!"UMAP" %in% reducedDimNames(sce)) {
    message("Computing UMAP...")
    set.seed(32)
    umaps <- as.data.frame(uwot::umap(reducedDim(sce, "pca.corrected")))
    reducedDim(sce, "UMAP") <- umaps
  }
  
  # Assign neighborhoods
  message("Assigning neighborhoods...")
  set.seed(32)
  sce_milo <- assign_neighbourhoods(
    sce, k = k, order = 2,
    filtering = TRUE, reducedDim_name = "pca.corrected", verbose = FALSE
  )
  
  # Get neighborhood statistics
  nhoods_sce <- nhoods(sce_milo)
  nhood_stat_ct <- data.frame(
    Nhood = 1:ncol(nhoods_sce),
    Nhood_center = colnames(nhoods_sce)
  )
  nhood_stat_ct <- miloR::annotateNhoods(
    sce_milo, nhood_stat_ct,
    coldata_col = "celltype.mapped"
  )
  
  # Calculate AUC per neighborhood
  message("Calculating AUC per neighborhood...")
  stat_auc <- suppressWarnings(
    calc_AUC_per_neighbourhood(sce_milo, condition_id = condition_id)
  )
  
  # Prepare results
  milo_res <- list(
    stat_auc = stat_auc,
    sce = sce,
    sce_milo = sce_milo
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    message("Saving results to: ", output_file)
    saveRDS(milo_res, output_file)
  }
  
  message("✓ Milo analysis complete")
  return(milo_res)
}


# =============================================================================
# Example Usage (commented out - uncomment to run)
# =============================================================================

# # Load configuration
# source("R/config.R")
# 
# # Load data
# sce <- readRDS(get_path("sce_mouseGastrulation.rds", "rds"))
# EmbryoCelltypeColours <- readRDS(get_path("EmbryoCelltypeColours.rds", "rds"))
# 
# # Downsample to selected cell types
# cts <- c("Spinal cord", "Haematoendothelial progenitors", 
#          "Endothelium", "Blood progenitors 1", "Blood progenitors 2")
# sce <- sce[, sce$celltype.mapped %in% cts]
# 
# # Rename for brevity
# sce$celltype.mapped[sce$celltype.mapped == "Haematoendothelial progenitors"] <- "Haem. prog-s."
# 
# # Remove tomato row
# sce <- sce[!rownames(sce) == "tomato-td", ]
# 
# # Filter samples
# sce <- sce[, colData(sce)$sample %in% c(1, 4)]
# 
# # Update tomato field
# sce$tomato <- sapply(sce$tomato, function(x) ifelse(isTRUE(x), "Tal1_KO", "WT"))
# 
# # Subset to highly variable genes
# library(scran)
# dec.sce <- modelGeneVar(sce)
# hvg.genes <- getTopHVGs(dec.sce, n = 3000)
# sce <- sce[hvg.genes, ]
# rowdata <- as.data.frame(rowData(sce))
# rownames(sce) <- rowdata$SYMBOL
# 
# # Run Milo analysis
# milo_res <- run_milo_analysis(
#   sce,
#   condition_id = "tomato",
#   output_file = get_path("milo_res.rds", "rds")
# )

# de_stat = de_test_neighbourhoods(sce_milo ,
#                                  sample_id = "sample",
#                                  design = ~tomato,
#                                  covariates = c("tomato"),
#                                  subset_nhoods = stat_auc$Nhood[!is.na(stat_auc$auc)],
#                                  output_type = "SCE",
#                                  plot_summary_stat = TRUE,
#                                  layout = "UMAP", 
#                                  verbose = T)










