# scripts/fig05_plasmid_replicons_heatmap.R
# Heatmap of plasmid replicon presence/absence per isolate, saved as PDF.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
  library(scales)
})

# ---- paths (edit if needed) ----
# Option A: matrix format (rows = isolates, cols = replicons, values 0/1 or Present/Absent)
mat_path  <- file.path("data","plasmid_replicons_matrix.csv")   # preferred
# Option B: long format (isolate_id, replicon, present) -> we’ll pivot
long_path <- file.path("data","plasmid_replicons_long.csv")     # alternative
meta_path <- file.path("data","isolate_metadata.csv")           # optional (serovar)
out_dir   <- file.path("outputs","figures")
out_file  <- file.path(out_dir,"Fig05_Plasmid_replicons_heatmap.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load presence/absence matrix ----
load_matrix <- function() {
  if (file.exists(mat_path)) {
    mraw <- read_csv(mat_path, show_col_types = FALSE)
    if (names(mraw)[1] != "isolate_id")
      stop("First column in plasmid_replicons_matrix.csv must be 'isolate_id'.")
    isolates <- mraw$isolate_id
    m <- as.data.frame(mraw[,-1, drop=FALSE])
    # accept 0/1 or Present/Absent
    m[] <- lapply(m, function(x){
      ifelse(tolower(as.character(x)) %in% c("1","present","yes","true"), 1,
      ifelse(tolower(as.character(x)) %in% c("0","absent","no","false"), 0, NA))
    })
    m <- as.matrix(m); rownames(m) <- isolates
    return(m)
  } else if (file.exists(long_path)) {
    lraw <- read_csv(long_path, show_col_types = FALSE)
    req  <- c("isolate_id","replicon","present")
    if (!all(req %in% names(lraw)))
      stop("plasmid_replicons_long.csv must have columns: isolate_id, replicon, present")
    lraw <- lraw %>%
      mutate(present = ifelse(tolower(as.character(present)) %in%
                                c("1","present","yes","true"), 1,
                              ifelse(tolower(as.character(present)) %in%
                                       c("0","absent","no","false"), 0, NA)))
    isolates  <- unique(lraw$isolate_id)
    replicons <- unique(lraw$replicon)
    m <- matrix(0, nrow = length(isolates), ncol = length(replicons),
                dimnames = list(isolates, replicons))
    for (i in seq_len(nrow(lraw))) {
      m[ as.character(lraw$isolate_id[i]),
         as.character(lraw$replicon[i]) ] <- lraw$present[i]
    }
    return(m)
  } else {
    stop("Provide either data/plasmid_replicons_matrix.csv or data/plasmid_replicons_long.csv")
  }
}

mat <- load_matrix()

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

# ---- heatmap (presence=1, absence=0) ----
PRESENCE_COLORS <- c("0" = "white", "1" = "steelblue")

ht <- Heatmap(
  mat,
  name               = "Presence",
  col                = PRESENCE_COLORS,
  right_annotation   = row_anno,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = grid::gpar(fontsize = 9),
  column_names_gp    = grid::gpar(fontsize = 9, rot = 90, just = "right"),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  column_title       = "Plasmid replicons by isolate",
  rect_gp            = grid::gpar(col = "grey85", lwd = 0.5),
  heatmap_legend_param = list(
    at     = c(0,1),
    labels = c("Absent","Present"),
    border = TRUE,
    legend_gp = grid::gpar(col = "black")
  )
)

# ---- save as PDF (vector) ----
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_file, width = 7, height = 5)
} else {
  grDevices::pdf(out_file, width = 7, height = 5, useDingbats = FALSE)
}
draw(ht)
dev.off()
message("Saved: ", out_file)
