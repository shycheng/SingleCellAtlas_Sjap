# SingleCellAtlas_Sjap

Code for the paper: **"Dynamic single-cell transcriptomics reveals lsamp-guided neural network formation in male *S. japonicum* driving female reproduction"**

R Markdown and Jupyter Notebook files can be used to replicate figures and analysis results in the paper.

## Project Structure

```
SingleCellAtlas_Sjap/
├── 00.scripts/
│   ├── *.Rmd                  # R Markdown figure scripts
│   ├── R/                     # Utility R functions
│   └── jupyter/               # Python notebooks
└── README.md
```

### R Markdown Files (`00.scripts/`)

| File | Description |
|------|-------------|
| `Fig1.Rmd` | Figure 1 - Main atlas visualization |
| `Fig2.Rmd` | Figure 2 - Analysis and comparisons |
| `Fig3.Rmd` | Figure 3 - Additional analysis |
| `FigS1-S6.Rmd` | Supplementary Figures S1-S6 |
| `FigS7_S9_S11.Rmd` | Supplementary Figures S7, S9, S11 |
| `FigS11-S15.Rmd` | Supplementary Figures S11-S15 |
| `FigS23.Rmd` | Supplementary Figure S23 |

### Jupyter Notebooks (`00.scripts/jupyter/`)

| File | Description |
|------|-------------|
| `Fig2b.ipynb` | Figure 2B generation |
| `FigS12-S13.ipynb` | Supplementary Figures S12-S13 |

### Utility R Scripts (`00.scripts/R/`)

| Script | Description |
|--------|-------------|
| `plotFunctions.R` | General plotting utilities (violin plots, bar plots, UMAP) |
| `plotSlingshot.R` | Trajectory analysis visualization with Slingshot |
| `plot_umap.R` | UMAP plotting functions |
| `Velocity.R` | RNA velocity analysis functions for scVelo integration |
| `runMiloDE.R` | Differential expression analysis with Milo |
| `visPseudotimeHeatmap.R` | Pseudotime heatmap visualization |

## Usage

1. Clone this repository
2. Open relevant `.Rmd` files in RStudio or `.ipynb` notebooks in JupyterLab
3. Run code cells to reproduce figures and analyses

## Requirements

### R
- Seurat
- tidyverse
- ggplot2, ggpubr, patchwork
- RColorBrewer, ggsci
- slingshot, tradeSeq
- clusterProfiler

### Python
- scanpy
- scvelo
- matplotlib

## Contact

For questions or issues, please contact shawyunchan@gmail.com
