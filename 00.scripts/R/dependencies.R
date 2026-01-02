# dependencies.R - Centralized Package Loading
# =============================================================================
# Load all required packages once to avoid repeated library() calls
# Source this at the beginning of R Markdown files
# =============================================================================

#' @description Loads all required packages for SingleCellAtlas_Sjap analysis
#' @details Suppresses startup messages for cleaner output

# Suppress startup messages
suppressPackageStartupMessages({
    # --- Core Data Analysis ---
    library(Seurat)
    library(tidyverse)
    library(dplyr)
    library(magrittr)

    # --- Single-Cell Specific ---
    library(SingleCellExperiment)
    library(slingshot)
    library(tradeSeq)
    library(scater)

    # --- Plotting Libraries ---
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    library(RColorBrewer)
    library(ggsci)
    library(viridis)

    # --- Gene Enrichment ---
    library(clusterProfiler)

    # --- Utility Libraries ---
    library(reshape2)
    library(ggalluvial)
    library(data.table)
    library(Matrix)

    # --- Workflow Helpers ---
    library(SeuratWrappers)
})

# Check for optional packages
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    message("Note: ggnewscale not installed. Some advanced plotting features may not work.")
    message("Install with: install.packages('ggnewscale')")
}

if (!requireNamespace("cowplot", quietly = TRUE)) {
    message("Note: cowplot not installed. Some legend features may not work.")
    message("Install with: install.packages('cowplot')")
}

message("✓ All core packages loaded successfully")
