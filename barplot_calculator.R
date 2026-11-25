library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(protti)

data=read_protti("C:/Users/ni9566/Desktop/Large_Drug_Screen_official/Official_with_WOx/WOX_LDS_full_mods_final_filtering.txt")
cell_count=read_protti("C:/Users/ni9566/Desktop/Large_Drug_Screen_official/cell_counts_reordered.txt")
cell_count <- cell_count[,-59 ]

# First 95 columns in data are peptide intensity columns
intensity_cols <- names(data)[1:95]

# condition_levels should match the column order in cell_count
condition_levels <- colnames(cell_count)

# Count detected peptides per condition
peptide_counts <- data %>%
  summarise(across(all_of(intensity_cols), ~ sum(!is.na(.x)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "condition",
    values_to = "n_peptides"
  )

# Convert cell_count from wide to long format
cell_counts <- cell_count %>%
  summarise(across(all_of(condition_levels), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "condition",
    values_to = "cell_count"
  )

# Combine peptide counts and cell counts
plot_data <- peptide_counts %>%
  left_join(cell_counts, by = "condition") %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    condition_base = str_remove(as.character(condition), "_\\d+$"),
    x_pos = as.numeric(condition)
  ) %>%
  arrange(condition)

# Add spacing after each condition ending in "_1"
gap_after <- grepl("_1$", condition_levels)

extra_gap <- cumsum(c(FALSE, gap_after[-length(gap_after)]))

plot_data <- plot_data %>%
  mutate(
    x_pos_gap = x_pos + extra_gap
  )

# Create one x axis label per condition, centered across its replicates
condition_labels <- plot_data %>%
  group_by(condition_base) %>%
  summarise(
    x_label_pos = mean(x_pos_gap),
    .groups = "drop"
  ) %>%
  mutate(
    condition_base = toupper(condition_base)
  )

# Scale cell count to the peptide count axis
scale_factor <- max(plot_data$n_peptides, na.rm = TRUE) /
  max(plot_data$cell_count, na.rm = TRUE)

# Plot peptide counts as bars and cell counts as a line
p = ggplot(plot_data, aes(x = x_pos_gap)) +
  geom_col(
    aes(y = n_peptides),
    fill = "turquoise3"
  ) +
  geom_line(
    aes(y = cell_count * scale_factor, group = 1),
    linewidth = 0.2
  ) +
  geom_point(
    aes(y = cell_count * scale_factor),
    size = 1.2
  ) +
  scale_x_continuous(
    breaks = plot_data$x_pos_gap,
    labels = rep("", nrow(plot_data)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = "Number of detected peptides",
    sec.axis = sec_axis(
      trans = ~ .x / scale_factor,
      name = "Estimated cell count"
    )
  ) +
  geom_text(
    data = condition_labels,
    aes(
      x = x_label_pos,
      y = -0.07 * max(plot_data$n_peptides, na.rm = TRUE),
      label = condition_base
    ),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.5,
    hjust = 1,
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 80, l = 5.5)
  )


p

ggsave(
  filename = "HLA_I_barplot.pdf",
  plot = p,
  width = 6,
  height = 5,
  units = "in"
)
getwd()
