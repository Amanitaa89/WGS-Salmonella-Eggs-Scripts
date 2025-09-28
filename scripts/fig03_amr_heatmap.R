# scripts/fig03_amr_heatmap.R
# AST heatmap (S/I/R) with letters on tiles, saved as PDF (vector).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
  # for a stable color palette when annotating classes
  library(scales)
})

# ---- paths (edit if needed) ----
ast_path     <- file.path("data","ast_matrix.csv")           # REQUIRED
meta_path    <- file.path("data","isolate_metadata.csv")     # optional (Serovar)
abx_map_path <- file.path("data","antibiotic_classes.csv")   # optional (Class)
out_dir      <- file.path("outputs","figures")
out_file     <- file.path(out_dir, "Fig03_AMR_heatmap.pdf")  # <-- PDF output

if (!file.exists(ast_path)) stop("Missing data/ast_matrix.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load AST matrix ----
ast_raw <- read_csv(ast_path, show_col_types = FALSE)

if (names(ast_raw)[1] != "isolate_id")
  stop("First column in ast_matrix.csv must be named 'isolate_id'.")

isolate_ids <- ast_raw$isolate_id
ast_mat_chr <- as.data.frame(ast_raw[,-1, drop=FALSE])

# keep only S/I/R (others -> NA)
valid_vals <- c("S","I","R")
ast_mat_chr[] <- lapply(ast_mat_chr, function(v) ifelse(v %in% valid_vals, v, NA))

ast_mat <- as.matrix(ast_mat_chr)
rownames(ast_mat) <- isolate_ids

# ---- optional row annotation: serovar ----
row_anno <- NULL
SEROVAR_COLORS <- c(
  "Enteritidis"="#1f78b4","Typhimurium"="#33a02c","Heidelberg"="#e31a1c",
  "Infantis"="#ff7f00","Minnesota"="#6a3d9a","Montevideo"="#b15928",
  "Bareilly"="#a6cee3","Cerro"="#b2df8a","Concord"="#fb9a99",
  "Wien"="#fdbf6f","Lagos"="#cab2d6","Other"="#7f7f7f","Unknown"="#7f7f7f"
)
if (file.exists(meta_path)) {
  meta <- read_csv(meta_path, show_col_types = FALSE)
  if (all(c("isolate_id","serovar") %in% names(meta))) {
    sero <- meta$serovar[match(rownames(ast_mat), meta$isolate_id)]
    sero[is.na(sero)] <- "Unknown"
    sero[!sero %in% names(SEROVAR_COLORS)] <- "Other"
    row_anno <- HeatmapAnnotation(
      df = data.frame(Serovar = factor(sero, levels = names(SEROVAR_COLORS))),
      col = list(Serovar = SEROVAR_COLORS),
      which = "row",
      show_annotation_name = TRUE,
      annotation_name_side = "top"
    )
  }
}

# ---- optional top annotation: antibiotic class ----
top_anno <- NULL
if (file.exists(abx_map_path)) {
  abx_map <- read_csv(abx_map_path, show_col_types = FALSE)
  if (all(c("antibiotic","class") %in% names(abx_map))) {
    abx_cols <- colnames(ast_mat)
    classes <- abx_map$class[match(abx_cols, abx_map$antibiotic)]
    classes[is.na(classes)] <- "Other"
    classes <- factor(classes)
    cls_lvls <- levels(classes)
    cls_cols <- setNames(scales::hue_pal()(length(cls_lvls)), cls_lvls)
    top_anno <- HeatmapAnnotation(Class = classes,
                                  col = list(Class = cls_cols),
                                  annotation_name_side = "left")
  }
}

# ---- heatmap object ----
AMR_COLORS <- c("S"="#2ca02c","I"="#ff7f0e","R"="#d62728")  # green/orange/red

ht <- Heatmap(
  ast_mat,
  name               = "AST",
  col                = AMR_COLORS,
  right_annotation   = row_anno,
  top_annotation     = top_anno,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = grid::gpar(fontsize = 10),
  column_names_gp    = grid::gpar(fontsize = 8, rot = 90, just = "right"),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  rect_gp            = grid::gpar(col = "grey80", lwd = 0.5),

  # show S/I/R letters in cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- ast_mat[i, j]
    if (!is.na(lab) && nzchar(lab)) grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 8))
  },

  heatmap_legend_param = list(
    at     = c("S","I","R"),
    labels = c("Susceptible","Intermediate","Resistant"),
    border = TRUE,
    legend_gp = grid::gpar(col = "black")
  )
)

# ---- save as PDF (vector) ----
# width/height in inches; adjust to fit your layout
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_file, width = 7, height = 5)
} else {
  grDevices::pdf(out_file, width = 7, height = 5, useDingbats = FALSE)
}
draw(ht)
dev.off()
message("Saved: ", out_file)
