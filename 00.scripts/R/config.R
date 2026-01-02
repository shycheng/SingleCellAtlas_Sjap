# Configuration for SingleCellAtlas_Sjap
# =============================================================================
# Edit paths according to your system, or set environment variables:
#   SCATLAS_DATA_DIR, SCATLAS_FIG_DIR, SCATLAS_RDS_DIR
# =============================================================================

# --- Directories ---
CONFIG <- list(
    data_dir = Sys.getenv("SCATLAS_DATA_DIR", unset = "~/data"),
    figure_dir = Sys.getenv("SCATLAS_FIG_DIR", unset = "~/02.figures"),
    rds_dir = Sys.getenv("SCATLAS_RDS_DIR", unset = "~/01.rds_files"),
    genome_dir = Sys.getenv("SCATLAS_GENOME_DIR", unset = "~/genome/Schitosoma_jp")
)

# --- Gene Annotation File ---
CONFIG$gene_info_file <- file.path(
    CONFIG$genome_dir,
    "Schistosoma_japonicum/SJ_GODB/sj_geneInfo_raw.tsv"
)

# --- Analysis Parameters ---
CONFIG$analysis <- list(
    min_cells = 3,
    min_features = 200,
    mt_threshold = 30,
    nfeature_min = 500,
    variable_features = 3000,
    pca_dims = 100,
    cluster_resolution = 2
)

# --- Sample Information ---
CONFIG$samples <- c("14dpi", "18dpi_F", "18dpi_M", "22dpi_F", "22dpi_M", "26dpi_F", "26dpi_M")

CONFIG$time_levels <- c(
    "26dpi_F", "22dpi_F", "18dpi_F", "14dpi",
    "18dpi_M", "22dpi_M", "26dpi_M"
)

# --- Color Palettes ---
CONFIG$colors <- list(
    # Main cell type palette (21 colors)
    cell_types = c(
        "#7FC97F", "#BEAED4", # Neoblast, Germline
        "#FEEDDE", "#FDBE85", "#FD8D3C", "#D94701", # Neural (Oranges)
        "#D95F02", "#7570B3", # Neural continued
        "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", # Tegument (Blues)
        "#E7298A", "#1F78B4", "#B2DF8A", # Muscle, Flame, Gut
        "#33A02C", "#FB9A99", "#E31A1C", "#FFFF99", # Others
        "grey" # Unknown
    ),

    # Time point palette (7 colors for samples)
    time_points = c(
        "#B2182B", "#EF8A62", "#FDDBC7", "yellow",
        "#D1E5F0", "#67A9CF", "#2166AC"
    ),

    # Alternative qualitative palette
    cb_palette = c(
        "#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
        "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27",
        "#e07233", "#ff523f", "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635",
        "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551", "#1a918f", "#ff66fc", "#2927c4",
        "#7149af", "#57e559", "#8e3af4", "#f9a270", "#22547f", "#db5e92", "#edd05e",
        "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f", "#7b34c1", "#0cf29a",
        "#d80fc1", "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5", "#925bea",
        "#63ff4f"
    ),

    # Heatmap color scales
    heatmaps = list(
        white_red = c(
            "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
            "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
        ),
        blue_yellow_red = c(
            "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
            "#FEE090", "#FDAE61", "#F46D43", "#D73027"
        ),
        yellow_green_purple = c(
            "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
            "#2C728E", "#3B528B", "#472D7B", "#440154"
        )
    )
)

# --- Cell Type Definitions ---
CONFIG$cell_type_levels <- c(
    "Neoblast",
    "Germline stem cell",
    "Neural stem cells",
    "Neural precursor cell",
    "Nonciliated neuron",
    "Ciliated neuron",
    "KK7+ Neuron",
    "Parenchyma",
    "Tegument progenitor",
    "Tegument progeny 1",
    "Tegument progeny 2",
    "Syncytial 1",
    "Syncytial 2",
    "Muscle",
    "Flame cells",
    "Gut",
    "Vitellocyte",
    "Mehlis gland",
    "Male gametes",
    "Esophageal gland",
    "Unknown"
)

# Name the cell type colors
names(CONFIG$colors$cell_types) <- CONFIG$cell_type_levels

# --- Helper Function to Load Gene Info ---
load_gene_info <- function(config = CONFIG) {
    if (!file.exists(config$gene_info_file)) {
        warning("Gene info file not found: ", config$gene_info_file)
        return(NULL)
    }

    id2name <- read.delim(config$gene_info_file, sep = "\t", header = FALSE)
    colnames(id2name) <- c("id", "Gene_Name", "Note")
    rownames(id2name) <- id2name$id
    id2name$Gene_Name <- as.character(id2name$Gene_Name)

    # Truncate long gene names

    idx <- which(nchar(as.character(id2name$Gene_Name)) > 13)
    id2name$Gene_Name[idx] <- as.character(rownames(id2name[idx, ]))

    # Clean note field
    id2name$Note <- gsub("\\(.*\\)", "", id2name$Note)
    id2name$id <- gsub("_", "-", id2name$id)
    rownames(id2name) <- id2name$id

    return(id2name)
}

# --- Utility: Get Full Path ---
get_path <- function(filename, type = c("data", "figure", "rds")) {
    type <- match.arg(type)
    base_dir <- switch(type,
        data = CONFIG$data_dir,
        figure = CONFIG$figure_dir,
        rds = CONFIG$rds_dir
    )
    file.path(base_dir, filename)
}

message("Config loaded. Use CONFIG$... to access settings.")
