# Velocity.R - RNA Velocity Analysis Functions
# =============================================================================
# Functions to prepare Seurat data for scVelo RNA velocity analysis
# =============================================================================

#' Export Seurat object data for scVelo analysis
#'
#' @description Exports spliced, unspliced, and ambiguous count matrices along
#'   with cell metadata and UMAP coordinates for use with Python scVelo.
#'
#' @param seobj A Seurat object with UMAP reduction computed
#' @param filename Base name for output files
#' @param base_dir Output directory (will create subdirectory)
#' @param spliced_data Pre-loaded spliced count matrix (from loom files)
#' @param unspliced_data Pre-loaded unspliced count matrix
#' @param ambiguous_data Pre-loaded ambiguous count matrix
#'
#' @return Invisibly returns the output directory path
#'
#' @details Creates a subdirectory named `{filename}_loom` containing:
#'   - `{filename}.spliced.txt` - Spliced counts
#'   - `{filename}.unspliced.txt` - Unspliced counts
#'   - `{filename}.ambiguous.txt` - Ambiguous counts
#'   - `{filename}.celltype.csv` - Cell barcodes, clusters, and types
#'   - `{filename}.genes.csv` - Gene names
#'   - `{filename}.umap.csv` - UMAP coordinates
#'
#' @examples
#' \dontrun{
#' # Load pre-computed velocity matrices
#' spliced <- readRDS("spliced.dat_loom.rds")
#' unspliced <- readRDS("unspliced.dat_loom.rds")
#' ambiguous <- readRDS("ambiguous.dat_loom.rds")
#'
#' # Export for scVelo
#' export_for_scvelo(
#'   seob, "my_analysis", "~/output",
#'   spliced, unspliced, ambiguous
#' )
#' }
#'
#' @export
export_for_scvelo <- function(seobj, filename, base_dir,
                              spliced_data, unspliced_data, ambiguous_data) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')")
  }

  # Create output directory
  output_path <- file.path(base_dir, paste0(filename, "_loom"))
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  message("Output directory: ", output_path)

  # Find common cells
  cells <- intersect(colnames(seobj), colnames(spliced_data))
  message("Total cells: ", length(cells))

  if (length(cells) == 0) {
    stop("No common cells found between Seurat object and velocity matrices")
  }

  # Subset matrices
  spliced_subset <- spliced_data[, cells]
  unspliced_subset <- unspliced_data[, cells]
  ambiguous_subset <- ambiguous_data[, cells]

  # Write count matrices
  message("Writing count matrices...")

  data.table::fwrite(
    as.data.frame(spliced_subset),
    file.path(output_path, paste0(filename, ".spliced.txt")),
    sep = "\t", row.names = TRUE, col.names = TRUE
  )

  data.table::fwrite(
    as.data.frame(unspliced_subset),
    file.path(output_path, paste0(filename, ".unspliced.txt")),
    sep = "\t", row.names = TRUE, col.names = TRUE
  )

  data.table::fwrite(
    as.data.frame(ambiguous_subset),
    file.path(output_path, paste0(filename, ".ambiguous.txt")),
    sep = "\t", row.names = TRUE, col.names = TRUE
  )

  # Extract and write metadata
  seobj$barcode <- colnames(seobj)

  if ("Cell_Type" %in% colnames(seobj@meta.data)) {
    celltype <- seobj@meta.data[cells, c("barcode", "seurat_clusters", "Cell_Type")]
  } else {
    celltype <- seobj@meta.data[cells, c("barcode", "seurat_clusters")]
  }

  write.csv(celltype,
    file.path(output_path, paste0(filename, ".celltype.csv")),
    row.names = FALSE
  )

  # Write gene names
  genes <- rownames(spliced_subset)
  write.csv(data.frame(gene = genes),
    file.path(output_path, paste0(filename, ".genes.csv")),
    row.names = FALSE
  )

  # Write UMAP coordinates
  if ("umap" %in% names(seobj@reductions)) {
    umap_coords <- seobj@reductions$umap@cell.embeddings[cells, ]
    write.csv(umap_coords,
      file.path(output_path, paste0(filename, ".umap.csv")),
      row.names = TRUE
    )
  } else {
    warning("UMAP reduction not found in Seurat object")
  }

  # Cleanup
  gc(verbose = FALSE)

  message("Export complete: ", output_path)
  invisible(output_path)
}


# =============================================================================
# Legacy function name (deprecated, use export_for_scvelo instead)
# =============================================================================
save_filesFORscVelo <- function(seobj, filename, baseDIR) {
  .Deprecated("export_for_scvelo")
  warning("This function requires global variables: spliced.dat, unspliced.dat, ambigous.dat")
  warning("Consider using export_for_scvelo() with explicit data arguments instead.")

  # Backwards compatibility - expects global variables
  if (!exists("spliced.dat") || !exists("unspliced.dat") || !exists("ambigous.dat")) {
    stop("Global velocity matrices not found. Load them first or use export_for_scvelo().")
  }

  export_for_scvelo(
    seobj, filename, baseDIR,
    get("spliced.dat", envir = .GlobalEnv),
    get("unspliced.dat", envir = .GlobalEnv),
    get("ambigous.dat", envir = .GlobalEnv)
  )
}
