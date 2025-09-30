# scripts/fig07_vfgs_by_function_heatmap.R
# VFG presence/absence heatmap grouped by functional categories (don’t need any other scripts).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
  library(scales)
})

# -------- paths (edit if needed) --------
vfg_path   <- file.path("data","vfg_presence.csv")      # REQUIRED: matrix; first col = isolate_id
map_path   <- file.path("data","vfg_class_map.csv")     # REQUIRED: mapping; cols: gene,class
meta_path  <- file.path("data","isolate_metadata.csv")  # optional: isolate_id, serovar
out_dir    <- file.path("outputs","figures")
out_file   <- file.path(out_dir,"Fig07_VFGs_by_function_heatmap.pdf")

if (!file.exists(vfg_path)) stop("Missing data/vfg_presence.csv")
if (!file.exists(map_path)) stop("Missing data/vfg_class_map.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------- load presence/absence matrix --------
vdf <- read_csv(vfg_path, show_col_types = FALSE)
if (names(vdf)[1] != "isolate_id")
  stop("First column in vfg_presence.csv must be 'isolate_id'.")

isolate_ids <- vdf$isolate_id
mat_df <- as.data.frame(vdf[,-1, drop=FALSE])

# accept 0/1 or Present/Absent
mat_df[] <- lapply(mat_df, function(x){
  x <- tolower(as.character(x))
  ifelse(x %in% c("1","present","yes","true"), 1,
  ifelse(x %in% c("0","absent","no","false"), 0, NA))
})
mat <- as.matrix(mat_df)
rownames(mat) <- isolate_ids

# -------- mapping gene -> functional class --------
map <- read_csv(map_path, show_col_types = FALSE)  # cols: gene,class
req <- c("gene","class")
if (!all(req %in% names(map))) stop("vfg_class_map.csv must have columns: gene,class")

# keep genes that exist in the matrix; warn on missing
genes_in_mat <- colnames(mat)
keep <- map$gene %in% genes_in_mat
if (any(!keep)) message("Warning: ", sum(!keep), " mapped gene(s) not found in matrix and will be ignored.")
map <- map[keep, , drop = FALSE]

# reorder columns by class, then by appearance in map
# (you can control the group order by ordering rows in vfg_class_map.csv, or set category_order below)
category_order <- unique(map$class)  # or set manually: c("Adhesins/Outer membrane","F4 fimbriae",...)
map$class <- factor(map$class, levels = category_order)
map <- map[order(map$class), ]
mat <- mat[, map$gene, drop = FALSE]

# top annotation strip = functional class (colored)
cls_lvls <- levels(map$class)
cls_cols <- setNames(hue_pal()(length(cls_lvls)), cls_lvls)
top_anno <- HeatmapAnnotation(
  Class = factor(map$class, levels = cls_lvls),
  col   = list(Class = cls_cols),
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)

# optional row annotation: serovar
row_anno <- NULL
if (file.exists(meta_path)) {
  meta <- read_csv(meta_path, show_col_types = FALSE)
  if (all(c("isolate_id","serovar") %in% names(meta))) {
    SEROVAR_COLORS <- c(
      "Enteritidis"="#1f78b4","Typhimurium"="#33a02c","Heidelberg"="#e31a1c",
      "Infantis"="#ff7f00","Minnesota"="#6a3d9a","Montevideo"="#b15928",
      "Bareilly"="#a6cee3","Cerro"="#b2df8a","Concord"="#fb9a99",
      "Wien"="#fdbf6f","Lagos"="#cab2d6","Other"="#7f7f7f","Unknown"="#7f7f7f"
    )
    sero <- meta$serovar[match(rownames(mat), meta$isolate_id)]
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

# column splits -> one panel per functional class
col_split <- factor(map$class, levels = cls_lvls)

# presence/absence colors (journal-friendly)
PRESENCE_COLORS <- c("0" = "white", "1" = "#3b6fb6")  # deep steel blue

ht <- Heatmap(
  mat,
  name               = "Presence",
  col                = PRESENCE_COLORS,
  top_annotation     = top_anno,
  right_annotation   = row_anno,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = grid::gpar(fontsize = 8.5),
  column_names_gp    = grid::gpar(fontsize = 8, rot = 90, just = "right"),
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,                   # we already ordered columns by class
  column_split       = col_split,               # <-- groups/panels
  column_gap         = unit(2, "mm"),
  border             = FALSE,
  rect_gp            = grid::gpar(col = "grey85", lwd = 0.4),
  heatmap_legend_param = list(
    at = c(0,1), labels = c("Absent","Present"),
    border = TRUE, legend_gp = grid::gpar(col = "black")
  ),
  column_title = "Virulence Factors Presence/Absence by Functional Group",
  column_title_gp = grid::gpar(fontface = "bold")
)

# draw & save (PDF vector)
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_file, width = 11, height = 7.5)  # wide to fit many groups
} else {
  grDevices::pdf(out_file, width = 11, height = 7.5, useDingbats = FALSE)
}
draw(ht, heatmap_legend_side = "right")
dev.off()
message("Saved: ", out_file)
