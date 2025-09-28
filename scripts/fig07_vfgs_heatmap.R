# scripts/fig07_vfgs_heatmap.R
# Virulence factor genes (VFGs) presence/absence heatmap with optional annotations.
# Output: PDF vector file in outputs/figures/.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
  library(scales)
})

# ---------- paths (edit if needed) ----------
vfg_path    <- file.path("data","vfg_presence.csv")        # REQUIRED
meta_path   <- file.path("data","isolate_metadata.csv")    # optional (adds Serovar on rows)
vfg_map_path<- file.path("data","vfg_class_map.csv")       # optional (Gene -> Class/Subsystem)
out_dir     <- file.path("outputs","figures")
out_file    <- file.path(out_dir,"Fig07_VFGs_heatmap.pdf")

if (!file.exists(vfg_path)) stop("Missing data/vfg_presence.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- load & clean ----------
# Accepts 0/1 or Present/Absent; first column must be isolate_id
vfg_raw <- read_csv(vfg_path, show_col_types = FALSE)

if (names(vfg_raw)[1] != "isolate_id") {
  stop("First column in vfg_presence.csv must be 'isolate_id'.")
}

isolate_ids <- vfg_raw$isolate_id
mat_df <- as.data.frame(vfg_raw[,-1, drop=FALSE])

mat_df[] <- lapply(mat_df, function(x){
  x <- tolower(as.character(x))
  ifelse(x %in% c("1","present","yes","true"), 1,
  ifelse(x %in% c("0","absent","no","false"), 0, NA))
})
mat <- as.matrix(mat_df)
rownames(mat) <- isolate_ids

# ---------- optional row annotation: Serovar ----------
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

# ---------- optional top annotation: VFG class/subsystem ----------
top_anno <- NULL
if (file.exists(vfg_map_path)) {
  vmap <- read_csv(vfg_map_path, show_col_types = FALSE)
  # expects columns: gene,class  (gene names must match column names in vfg_presence.csv)
  if (all(c("gene","class") %in% names(vmap))) {
    genes   <- colnames(mat)
    classes <- vmap$class[match(genes, vmap$gene)]
    classes[is.na(classes)] <- "Other"
    classes <- factor(classes)
    cls_lvls <- levels(classes)
    cls_cols <- setNames(hue_pal()(length(cls_lvls)), cls_lvls)
    top_anno <- HeatmapAnnotation(
      Class = classes,
      col   = list(Class = cls_cols),
      annotation_name_side = "left"
    )
  }
}

# ---------- heatmap ----------
PRESENCE_COLORS <- c("0" = "white", "1" = "black")  # print-friendly; change to blue if preferred

ht <- Heatmap(
  mat,
  name               = "Presence",
  col                = PRESENCE_COLORS,
  right_annotation   = row_anno,
  top_annotation     = top_anno,
  show_row_names     = TRUE,
  show_column_names  = FALSE,   # many genes → labels off (toggle to TRUE if few)
  row_names_gp       = grid::gpar(fontsize = 9),
  column_names_gp    = grid::gpar(fontsize = 7, rot = 90, just = "right"),
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  column_title       = "Virulence factor genes (presence/absence)",
  rect_gp            = grid::gpar(col = "grey85", lwd = 0.4),
  heatmap_legend_param = list(
    at = c(0,1), labels = c("Absent","Present"),
    border = TRUE, legend_gp = grid::gpar(col = "black")
  )
)

# ---------- save as PDF ----------
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_file, width = 8, height = 6)
} else {
  grDevices::pdf(out_file, width = 8, height = 6, useDingbats = FALSE)
}
draw(ht)
dev.off()
message("Saved: ", out_file)
