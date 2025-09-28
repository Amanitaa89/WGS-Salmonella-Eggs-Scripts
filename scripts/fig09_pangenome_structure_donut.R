# scripts/fig09_pangenome_structure_donut.R
# Donut chart of pangenome categories (Core, Soft core, Shell, Cloud), saved as PDF.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---------- paths ----------
in_path  <- file.path("data", "pangenome_categories.csv")
out_dir  <- file.path("outputs","figures")
out_file <- file.path(out_dir, "Fig09_Pangenome_structure_donut.pdf")

if (!file.exists(in_path)) stop("Missing data/pangenome_categories.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- read & normalize ----------
# Accept either:
#   category,percent
# or category,count
raw <- read_csv(in_path, show_col_types = FALSE)

nm <- tolower(names(raw))
names(raw) <- nm

if (!"category" %in% names(raw)) stop("Input must have a 'category' column.")
if (!("percent" %in% names(raw) | "count" %in% names(raw))) {
  stop("Input must have either a 'percent' or 'count' column.")
}

df <- raw %>%
  mutate(category = tolower(trimws(category))) %>%
  mutate(category = recode(category,
                           "softcore"="soft core",
                           "soft_core"="soft core")) %>%
  group_by(category) %>% summarise(
    percent = if ("percent" %in% names(raw))
                sum(as.numeric(percent), na.rm=TRUE)
              else NA_real_,
    count   = if ("count" %in% names(raw))
                sum(as.numeric(count), na.rm=TRUE)
              else NA_real_,
    .groups = "drop"
  )

if (is.na(df$percent[1])) {
  total <- sum(df$count, na.rm=TRUE)
  if (total == 0) stop("Counts sum to zero.")
  df <- df %>% mutate(percent = 100 * count / total)
} else {
  total <- if ("count" %in% names(df)) sum(df$count, na.rm=TRUE) else NA_real_
}

# order & palette (consistent, color-blind friendly)
lvl <- c("core","soft core","shell","cloud")
pal <- c("core"="#2171b5",       # dark blue
         "soft core"="#6baed6",  # light blue
         "shell"="#fdbf6f",      # orange
         "cloud"="#b3b3b3")      # grey

df <- df %>%
  mutate(category = factor(category, levels = lvl)) %>%
  arrange(category) %>%
  mutate(pct_lab = round(percent, 1),
         ypos = cumsum(percent) - percent/2)

# ---------- donut ----------
p <- ggplot(df, aes(x = 2, y = percent, fill = category)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  # inner white ring to create the "donut"
  geom_col(data = data.frame(x=1, percent=100),
           aes(x = x, y = percent), inherit.aes = FALSE,
           fill = "white", color = NA, width = 1) +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = pal, name = "Category") +
  geom_text(aes(y = ypos, label = paste0(pct_lab, "%")), color="black", size = 3.5) +
  theme_void(base_size = 11) +
  labs(title = "Pangenome structure") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# optional center label: total gene clusters if counts were provided
if (!is.na(total)) {
  p <- p + annotate("text", x = 1, y = 0, label = paste0("n = ", total),
                    size = 4, fontface = "bold")
}

# ---------- save ----------
pdf(out_file, width = 6.5, height = 5, useDingbats = FALSE)
print(p)
dev.off()
message("Saved: ", out_file)
