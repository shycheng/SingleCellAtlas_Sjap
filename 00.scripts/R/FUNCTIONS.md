# Utility Functions Reference

This document provides a quick reference for all utility functions available in `00.scripts/R/`.

## Quick Start

```r
# Load configuration first
source("00.scripts/R/config.R")

# Then load desired utility scripts
source("00.scripts/R/plotFunctions.R")
source("00.scripts/R/plot_umap.R")
source("00.scripts/R/plotSlingshot.R")
```

---

## config.R

Central configuration for paths, colors, and analysis parameters.

| Object | Description |
|--------|-------------|
| `CONFIG$data_dir` | Default data directory |
| `CONFIG$figure_dir` | Default figure output directory |
| `CONFIG$colors$cell_types` | Standard cell type color palette (21 colors) |
| `CONFIG$colors$time_points` | Time point colors (7 samples) |
| `CONFIG$cell_type_levels` | Ordered factor levels for cell types |
| `load_gene_info()` | Load gene ID to name mapping |
| `get_path(filename, type)` | Build full path for data/figure/rds files |

---

## plotFunctions.R

General plotting utilities for Seurat objects.

| Function | Description |
|----------|-------------|
| `stacked_vln_plot(obj, features)` | Stacked violin plots for multiple genes |
| `barplot_per(Seob)` | Percentage bar plot by sample/time |
| `barplot_cluster(seob, umap_plot)` | Cluster composition bar plot |
| `barplot_celltype(Seob, umap_plot, n)` | Cell type composition bar plot |
| `getFeaturePlots(Seobj, id)` | Feature plot with gene annotation |
| `plot_alluvium(seob, umapplot, ...)` | Alluvial diagram for cell proportions |
| `simFP(id, seob)` | Simple feature plot |
| `multi_simFP(ids, seob)` | Multiple feature plots in grid |
| `plot_umap(seob, output_file)` | Export UMAP to PDF |

---

## plot_umap.R

Custom UMAP plotting with advanced features.

| Function | Description |
|----------|-------------|
| `plot_umap_custom(seurat_obj, ...)` | Flexible UMAP with faceting, labels, highlights |
| `create_numbered_umap(seurat_obj, ...)` | UMAP with numbered cluster labels |

**Key Parameters:**
- `embedding` - Reduction to use (default: "umap")
- `feature` - Metadata column to color by
- `split_by` - Faceting variable
- `highlight_clusters` - Clusters to outline

---

## plotSlingshot.R

Trajectory analysis visualization.

| Function | Description |
|----------|-------------|
| `plot_slingshot(sds, show_arrow)` | Basic slingshot trajectory plot |
| `plot_trajectory_on_background(sce, full_seurat_obj, ...)` | Trajectory on full UMAP background |
| `plotSmoothers_multiGenes(models, counts, genes, ...)` | Multi-gene expression along pseudotime |

---

## visPseudotimeHeatmap.R

Pseudotime heatmap visualization.

| Function | Description |
|----------|-------------|
| `smoothMatrix(to_plot)` | Smooth and scale expression matrix |
| `getGroupMat(sce_obj, counts, ...)` | Get binned expression by pseudotime |
| `getBinCellTypes(sce_obj, ...)` | Get dominant cell type per pseudotime bin |
| `getBinCellTypes_fixed(sce_obj, ...)` | Fixed-interval version with interpolation |

---

## Velocity.R

RNA velocity export for scVelo.

| Function | Description |
|----------|-------------|
| `export_for_scvelo(seobj, filename, base_dir, ...)` | Export Seurat data for scVelo analysis |

**Output files:**
- `{name}.spliced.txt` - Spliced counts
- `{name}.unspliced.txt` - Unspliced counts  
- `{name}.celltype.csv` - Cell annotations
- `{name}.umap.csv` - UMAP coordinates

---

## runMiloDE.R

Differential expression with Milo neighborhoods.

> **Note:** This is a script template, not a function library. Modify directly for your analysis.

---

## Color Palettes

Use colors from `CONFIG$colors`:

```r
source("R/config.R")

# Cell type colors (named vector)
CONFIG$colors$cell_types

# Time point colors
CONFIG$colors$time_points

# For ggplot2
scale_color_manual(values = CONFIG$colors$cell_types)
```
