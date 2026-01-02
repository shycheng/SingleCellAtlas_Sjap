# SingleCellAtlas_Sjap

Code for the paper: **"Dynamic single-cell transcriptomics reveals lsamp-guided neural network formation in male *S. japonicum* driving female reproduction"**

R Markdown and Jupyter Notebook files can be used to replicate figures and analysis results in the paper.

## Quick Start

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/SingleCellAtlas_Sjap.git
cd SingleCellAtlas_Sjap

# (Optional) Install Python dependencies
pip install -r requirements.txt
```

```r
# In R, load the configuration first
source("00.scripts/R/config.R")

# Then load utility functions as needed
source("00.scripts/R/plotFunctions.R")
```

## Project Structure

```
SingleCellAtlas_Sjap/
├── 00.scripts/
│   ├── *.Rmd                  # R Markdown figure scripts
│   ├── R/                     # Utility R functions
│   │   ├── config.R           # ⭐ Centralized configuration
│   │   ├── plotFunctions.R    # General plotting utilities
│   │   ├── plot_umap.R        # Custom UMAP functions
│   │   ├── plotSlingshot.R    # Trajectory visualization
│   │   ├── Velocity.R         # scVelo export functions
│   │   ├── visPseudotimeHeatmap.R
│   │   └── FUNCTIONS.md       # Function reference
│   └── jupyter/               # Python notebooks
├── .gitignore
├── requirements.txt           # Python dependencies
└── README.md
```

### R Markdown Files (`00.scripts/`)

| File | Description |
|------|-------------|
| `Fig1.Rmd` | Figure 1 - Main atlas visualization |
| `Fig2.Rmd` | Figure 2 - Sub-clustering & trajectory analysis of germline/vitelline cells |
| `Fig3.Rmd` | Figure 3 - Sub-clustering & trajectory analysis of neuron cells |
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

See [FUNCTIONS.md](00.scripts/R/FUNCTIONS.md) for detailed function documentation.

| Script | Description |
|--------|-------------|
| `config.R` | **Start here** - Paths, colors, and parameters |
| `plotFunctions.R` | General plotting (violin, bar, UMAP) |
| `plotSlingshot.R` | Trajectory analysis with Slingshot |
| `plot_umap.R` | Custom UMAP plotting with faceting |
| `Velocity.R` | RNA velocity export for scVelo |
| `visPseudotimeHeatmap.R` | Pseudotime heatmap visualization |

## Configuration

Edit `00.scripts/R/config.R` to set your local paths, or use environment variables:

```r
# Option 1: Edit config.R directly
CONFIG$data_dir <- "/path/to/your/data"

# Option 2: Set environment variables
Sys.setenv(SCATLAS_DATA_DIR = "/path/to/your/data")
Sys.setenv(SCATLAS_FIG_DIR = "/path/to/figures")
```

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
