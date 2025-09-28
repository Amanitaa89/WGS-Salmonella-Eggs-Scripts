# scripts/fig04_args_heatmap.R
# Heatmap of ARG presence/absence, grouped by antibiotic class, saved as PDF.

suppressPackageStartupMessages({
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(scales)
})

# ---- paths ----
file_path <- file.path("data","Egg_ARGs.csv")      # your presence/absence CSV
out_dir   <- file.path("outputs","figures")
out_file  <- file.path(out_dir, "Fig04_ARGs_heatmap.pdf")

if (!file.exists(file_path)) stop("Missing data/Egg_ARGs.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- read & preprocess ----
raw <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)

antibiotic_class <- as.character(raw[1, 2:ncol(raw)])
genes            <- as.character(raw[2, 2:ncol(raw)])
samples          <- as.character(raw[3:nrow(raw), 1])

pa_mat <- raw[3:nrow(raw), 2:ncol(raw)]
colnames(pa_mat) <- genes
rownames(pa_mat) <- samples

mat <- apply(pa_mat, c(1,2), function(x){
  if (tolower(x)=="present") return(1)
  if (tolower(x)=="absent")  return(0)
  NA_real_
})
mat <- matrix(mat,
              nrow    = nrow(pa_mat),
              ncol    = ncol(pa_mat),
              dimnames = dimnames(pa_mat))

# ---- build class annotation ----
classes_unique <- unique(antibiotic_class)
class_colors   <- setNames(
  hue_pal()(length(classes_unique)),
  classes_unique
)

col_ha <- HeatmapAnnotation(
  Class = antibiotic_class,
  col   = list(Class = class_colors),
  annotation_name_side = "left"
)

# ---- draw heatmap ----
ht <- Heatmap(
  mat,
  name               = "Presence",
  col                = c("0" = "white", "1" = "steelblue"),
  top_annotation     = col_ha,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = grid::gpar(fontsize = 10),
  column_names_gp    = grid::gpar(fontsize = 8, rot = 90, just = "right"),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  column_title       = "ARGs by Antibiotic Class",
  rect_gp            = grid::gpar(col = "grey80", lwd = 0.5),
  heatmap_legend_param = list(
    at        = c(0, 1),
    labels    = c("Absent", "Present"),
    border    = TRUE,
    legend_gp = grid::gpar(col = "black")
  )
)

# ---- save as PDF ----
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_file, width = 7, height = 5)
} else {
  grDevices::pdf(out_file, width = 7, height = 5, useDingbats = FALSE)
}
draw(ht)
dev.off()
message("Saved: ", out_file)
