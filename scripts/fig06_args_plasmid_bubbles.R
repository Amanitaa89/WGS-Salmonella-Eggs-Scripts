# scripts/fig06_args_plasmid_bubbles.R
# Bubble plot showing associations between ARGs and plasmid replicons across serovars.
# x = serovar (or any group), y = ARG gene, bubble size = count of isolates,
# color = plasmid replicon. Output is a PDF (vector).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---------------- paths (edit if needed) ----------------
# Long-format input with ONE row per isolate–gene–replicon observation
# Required columns: isolate_id, serovar, gene, replicon, present (1/0 or yes/no)
in_path  <- file.path("data","args_plasmid_long.csv")
out_dir  <- file.path("outputs","figures")
out_file <- file.path(out_dir,"Fig06_ARGs_plasmid_bubbles.pdf")

if (!file.exists(in_path)) stop("Missing data/args_plasmid_long.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- load & clean ----------------
dat <- read_csv(in_path, show_col_types = FALSE)

req <- c("isolate_id","serovar","gene","replicon","present")
if (!all(req %in% names(dat))) {
  stop("Input must contain columns: isolate_id, serovar, gene, replicon, present")
}

# Normalize 'present' to 0/1
dat <- dat %>%
  mutate(present = tolower(as.character(present)),
         present = ifelse(present %in% c("1","yes","present","true"), 1,
                   ifelse(present %in% c("0","no","absent","false"), 0, NA)))

# Keep only present=1 rows
dat1 <- dat %>% filter(present == 1)

# Count co-occurrences per (serovar, gene, replicon)
counts <- dat1 %>%
  count(serovar, gene, replicon, name = "n") %>%
  arrange(desc(n))

# --------------- plot ----------------
# Color palette for replicons (auto, color-blind friendly)
replicon_lvls <- sort(unique(counts$replicon))
replicon_cols <- setNames(hue_pal()(length(replicon_lvls)), replicon_lvls)

p <- ggplot(counts, aes(x = serovar, y = gene)) +
  geom_point(aes(size = n, color = replicon), alpha = 0.85) +
  scale_size_area(max_size = 10, breaks = pretty(counts$n, n = 4), name = "Isolate count") +
  scale_color_manual(values = replicon_cols, name = "Plasmid replicon") +
  labs(
    title = "ARG–plasmid co-occurrence by serovar (bubble plot)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# --------------- save ----------------
pdf(out_file, width = 9, height = 6, useDingbats = FALSE)
print(p)
dev.off()
message("Saved: ", out_file)
