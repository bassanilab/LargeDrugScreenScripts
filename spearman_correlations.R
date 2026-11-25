library(protti)
setwd("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/Spearmancorrelationplots/")
data=read_protti("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/Spearmancorrelationplots/HLAII_official_fusedMOx_Carbamidomethyl_removednonbinders_normalized.txt")
cell_ordering=read_protti("cell_counts_reordered.txt")

data <- data %>%
  select(any_of(names(cell_ordering)), everything())


sts_cor <- data %>%
  select(starts_with("dms3")) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

sts_cor

conditions <- unique(
  sub("_\\d+$", "", names(cell_ordering))
)

conditions

# Compute all pairwise Spearman correlations for each condition
cor_list <- lapply(conditions, function(cond) {
  
  cols <- names(data)[startsWith(names(data), paste0(cond, "_"))]
  
  cor_mat <- cor(
    data[, cols, drop = FALSE],
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  # Keep each pairwise correlation only once
  cor_mat[upper.tri(cor_mat)]
})

names(cor_list) <- conditions

# Pad shorter conditions with NA if necessary
max_length <- max(lengths(cor_list))

cor_list <- lapply(cor_list, function(x) {
  length(x) <- max_length
  x
})

# Create dataframe
spearman_conditions <- as.data.frame(cor_list)

spearman_conditions


library(dplyr)
library(tidyr)
library(ggplot2)

spearman_long <- spearman_conditions %>%
  pivot_longer(
    cols = everything(),
    names_to = "condition",
    values_to = "spearman"
  ) %>%
  filter(!is.na(spearman))

median_df <- spearman_long %>%
  group_by(condition) %>%
  summarise(
    median_spearman = median(spearman),
    .groups = "drop"
  )

# Extract conditions in the same order as cell_ordering
conditions <- unique(
  sub("_\\d+$", "", names(cell_ordering))
)

# Apply that ordering
spearman_long$condition <- factor(
  spearman_long$condition,
  levels = conditions
)

median_df$condition <- factor(
  median_df$condition,
  levels = conditions
)

# Plot
p <- ggplot() +
  geom_col(
    data = median_df,
    aes(
      x = condition,
      y = median_spearman
    ),
    width = 0.7,
    fill = "orange"
  ) +
  geom_jitter(
    data = spearman_long,
    aes(
      x = condition,
      y = spearman
    ),
    width = 0.12,
    height = 0,
    size = 2
  ) +
  scale_y_reverse() +
  labs(
    x = "Condition",
    y = "Spearman correlation"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

p

ggsave(
  filename = "HLAII_spearman_correlation.pdf",
  plot = p,
  width = 6,
  height = 5,
  units = "in"
)


