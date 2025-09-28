# scripts/fig02_serovar_barplot.R
# Barplot of Salmonella serovar counts with labels on top of bars
# Output: PDF vector at outputs/figures/Fig02_serovar_barplot.pdf

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ---------------- paths ----------------
data_path <- file.path("data","serovar_counts.csv")   # CSV with columns: serovar,count
out_dir   <- file.path("outputs","figures")
out_file  <- file.path(out_dir,"Fig02_serovar_barplot.pdf")

if (!file.exists(data_path)) stop("Missing data/serovar_counts.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- load data ----------------
df <- read_csv(data_path, show_col_types = FALSE)

stopifnot(all(c("serovar","count") %in% names(df)))

# order serovars by descending count
df <- df %>% mutate(serovar = reorder(serovar, -count))

# ---------------- plot ----------------
p <- ggplot(df, aes(x = serovar, y = count)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = count), vjust = -0.2, size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Number of isolates",
    title = "Distribution of Salmonella serovars among isolates"
  )

# ---------------- save ----------------
ggsave(out_file, plot = p, width = 8, height = 4, device = cairo_pdf)
message("Saved: ", out_file)
