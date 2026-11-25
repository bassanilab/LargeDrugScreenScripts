library(protti)
library(limma)


#setwd("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/GSEA_enrichments/")
data=read_protti("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/GSEA_enrichments/class1_official_with_GO.txt")


# ============================================================
# Immunopeptidomics preprocessing
# ============================================================

# User defined three letter condition code
condition_code <- "vlm"


# ============================================================
# 1. Identify annotation and sample columns
# ============================================================

annotation_cols <- c(
  "t_t_pg_genes",
  "t_t_eg_modified_sequence",
  "c_gobp_slim_name",
  "c_kegg_name",
  "c_gocc_slim_name"
)

# Experimental condition columns
condition_cols <- grep(
  paste0("^", condition_code),
  colnames(data),
  value = TRUE,
  ignore.case = TRUE
)

# Control columns beginning with dms3 or dms4
control_cols <- grep(
  "^dms[34]",
  colnames(data),
  value = TRUE,
  ignore.case = TRUE
)

sample_cols <- c(condition_cols, control_cols)


# Keep only relevant columns
data_reduced <- data[
  ,
  c(annotation_cols, sample_cols),
  drop = FALSE
]


# ============================================================
# 2. Filter peptides based on replicate coverage
# ============================================================

# Required number of observed replicates
condition_required <- ceiling(length(condition_cols) * 2 / 3)
control_required   <- ceiling(length(control_cols) * 2 / 3)

# Number of observed values per peptide
condition_observed <- rowSums(
  !is.na(data_reduced[, condition_cols, drop = FALSE])
)

control_observed <- rowSums(
  !is.na(data_reduced[, control_cols, drop = FALSE])
)

# Keep peptide if sufficient observations are present in
# either the condition OR control group
keep_peptide <-
  condition_observed >= condition_required |
  control_observed >= control_required

data_filtered <- data_reduced[keep_peptide, , drop = FALSE]


# ============================================================
# 3. Median normalization
# ============================================================

# Calculate median of each sample
sample_medians <- sapply(
  data_filtered[, sample_cols, drop = FALSE],
  median,
  na.rm = TRUE
)

# Target median shared by all samples
target_median <- median(sample_medians, na.rm = TRUE)

# Shift every sample so that its median equals target_median
data_normalized <- data_filtered

for (col in sample_cols) {
  
  current_median <- median(
    data_normalized[[col]],
    na.rm = TRUE
  )
  
  data_normalized[[col]] <-
    data_normalized[[col]] +
    (target_median - current_median)
}


# ============================================================
# 4. Left shifted normal distribution imputation
# ============================================================

set.seed(123)

data_imputed <- data_normalized

for (col in sample_cols) {
  
  values <- data_normalized[[col]]
  
  observed_values <- values[!is.na(values)]
  
  # Distribution of observed peptide intensities
  observed_mean <- mean(observed_values)
  observed_sd   <- sd(observed_values)
  
  # Imputation distribution
  imputation_mean <- observed_mean - 1.8 * observed_sd
  imputation_sd   <- 0.3 * observed_sd
  
  missing <- is.na(values)
  
  # Generate a separate random value for every missing peptide
  data_imputed[[col]][missing] <- rnorm(
    sum(missing),
    mean = imputation_mean,
    sd = imputation_sd
  )
}


# ============================================================
# Final processed dataframe
# ============================================================

data_processed <- data_imputed


# ============================================================
# 5. Gene set enrichment analysis with limma
# ============================================================




# ============================================================
# 5.1 Expression matrix
# ============================================================

# Samples used for the comparison
analysis_cols <- c(condition_cols, control_cols)

# Create peptide abundance matrix
expr_matrix <- as.matrix(
  data_processed[, analysis_cols, drop = FALSE]
)

storage.mode(expr_matrix) <- "numeric"

# Give every peptide row a unique identifier
rownames(expr_matrix) <- paste0("peptide_", seq_len(nrow(expr_matrix)))


# ============================================================
# 5.2 Experimental design
# ============================================================

group <- factor(
  c(
    rep("condition", length(condition_cols)),
    rep("control", length(control_cols))
  ),
  levels = c("control", "condition")
)

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)

contrast_matrix <- makeContrasts(
  condition_vs_control = condition - control,
  levels = design
)


# ============================================================
# 5.3 Function to convert annotation column into gene sets
# ============================================================

make_gene_sets <- function(annotation_vector,
                           row_ids,
                           separator = ";") {
  
  gene_sets <- list()
  
  for (i in seq_along(annotation_vector)) {
    
    annotation <- annotation_vector[i]
    
    # Ignore missing or empty annotations
    if (
      is.na(annotation) ||
      annotation == "" ||
      annotation == "NA"
    ) {
      next
    }
    
    # Split multiple annotations belonging to the same peptide
    terms <- unlist(
      strsplit(
        annotation,
        split = separator,
        fixed = TRUE
      )
    )
    
    terms <- trimws(terms)
    
    terms <- terms[
      terms != "" &
        !is.na(terms)
    ]
    
    # Add peptide to each annotation set
    for (term in terms) {
      
      gene_sets[[term]] <- c(
        gene_sets[[term]],
        row_ids[i]
      )
    }
  }
  
  # Remove duplicated entries
  gene_sets <- lapply(
    gene_sets,
    unique
  )
  
  return(gene_sets)
}


# ============================================================
# 5.4 Build BP, KEGG and CC gene sets
# ============================================================

bp_sets <- make_gene_sets(
  annotation_vector = data_processed$c_gobp_slim_name,
  row_ids = rownames(expr_matrix)
)

kegg_sets <- make_gene_sets(
  annotation_vector = data_processed$c_kegg_name,
  row_ids = rownames(expr_matrix)
)

cc_sets <- make_gene_sets(
  annotation_vector = data_processed$c_gocc_slim_name,
  row_ids = rownames(expr_matrix)
)


# ============================================================
# 5.5 Remove very small gene sets
# ============================================================

minimum_set_size <- 5


bp_sets <- bp_sets[
  lengths(bp_sets) >= minimum_set_size
]

kegg_sets <- kegg_sets[
  lengths(kegg_sets) >= minimum_set_size
]

cc_sets <- cc_sets[
  lengths(cc_sets) >= minimum_set_size
]


# ============================================================
# 5.6 Convert peptide names to matrix indices
# ============================================================

convert_sets_to_indices <- function(gene_sets, expression_matrix) {
  
  ids <- rownames(expression_matrix)
  
  lapply(
    gene_sets,
    function(x) {
      match(x, ids)
    }
  )
}


bp_index <- convert_sets_to_indices(
  bp_sets,
  expr_matrix
)

kegg_index <- convert_sets_to_indices(
  kegg_sets,
  expr_matrix
)

cc_index <- convert_sets_to_indices(
  cc_sets,
  expr_matrix
)


# ============================================================
# 5.7 CAMERA enrichment analysis
# ============================================================

camera_bp <- camera(
  y = expr_matrix,
  index = bp_index,
  design = design,
  contrast = contrast_matrix[, "condition_vs_control"]
)

camera_kegg <- camera(
  y = expr_matrix,
  index = kegg_index,
  design = design,
  contrast = contrast_matrix[, "condition_vs_control"]
)

camera_cc <- camera(
  y = expr_matrix,
  index = cc_index,
  design = design,
  contrast = contrast_matrix[, "condition_vs_control"]
)


# ============================================================
# 5.8 Convert results into tables
# ============================================================

camera_bp$GeneSet <- rownames(camera_bp)
camera_bp$Database <- "GO_BP"

camera_kegg$GeneSet <- rownames(camera_kegg)
camera_kegg$Database <- "KEGG"

camera_cc$GeneSet <- rownames(camera_cc)
camera_cc$Database <- "GO_CC"


# Combine all enrichment results
enrichment_results <- rbind(
  camera_bp,
  camera_kegg,
  camera_cc
)

rownames(enrichment_results) <- NULL


# Reorder columns
enrichment_results <- enrichment_results[
  ,
  c(
    "Database",
    "GeneSet",
    "NGenes",
    "Direction",
    "PValue",
    "FDR"
  )
]


# Sort by FDR
enrichment_results <- enrichment_results[
  order(enrichment_results$FDR),
]


# ============================================================
# Final enrichment table
# ============================================================

enrichment_results

# ============================================================
# 6. Lollipop plot of top enriched and depleted terms
# ============================================================

library(ggplot2)
library(dplyr)


# ============================================================
# 6.1 Select top 10 enriched and depleted terms
# ============================================================

top_enriched <- enrichment_results %>%
  filter(Direction == "Up") %>%
  arrange(PValue) %>%
  slice_head(n = 10)


top_depleted <- enrichment_results %>%
  filter(Direction == "Down") %>%
  arrange(PValue) %>%
  slice_head(n = 10)


# Combine results
plot_data <- bind_rows(
  top_enriched,
  top_depleted
)


# ============================================================
# 6.2 Calculate signed significance score
# ============================================================

plot_data <- plot_data %>%
  mutate(
    
    log10_pvalue = -log10(PValue),
    
    signed_log10_pvalue = ifelse(
      Direction == "Up",
      log10_pvalue,
      -log10_pvalue
    ),
    
    # Add database name to distinguish terms if required
    plot_label = paste0(
      GeneSet,
      " [",
      Database,
      "]"
    )
  )


# ============================================================
# 6.3 Order terms according to significance
# ============================================================

plot_data$plot_label <- factor(
  plot_data$plot_label,
  levels = plot_data$plot_label[
    order(plot_data$signed_log10_pvalue)
  ]
)


# ============================================================
# 6.4 Create lollipop plot
# ============================================================

lollipop_plot <- ggplot(
  plot_data,
  aes(
    x = signed_log10_pvalue,
    y = plot_label
  )
) +
  
  # Lollipop stems
  geom_segment(
    aes(
      x = 0,
      xend = signed_log10_pvalue,
      y = plot_label,
      yend = plot_label,
      colour = Direction
    ),
    linewidth = 0.8
  ) +
  
  # Lollipop heads
  geom_point(
    aes(
      colour = Direction,
      size = NGenes
    )
  ) +
  
  # Centre line
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  
  # Labels
  labs(
    title = paste0(
      toupper(condition_code),
      " pathway enrichment"
    ),
    
    subtitle = paste0(
      toupper(condition_code),
      " versus DMS control"
    ),
    
    x = expression(
      "Signed " ~ -log[10](PValue)
    ),
    
    y = NULL,
    
    colour = "Direction",
    
    size = "Peptides"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    legend.position = "right"
  )


# Display plot
lollipop_plot

# ============================================================
# 6.5 Save lollipop plot as PDF
# ============================================================

pdf_filename <- paste0(
  toupper(condition_code),
  "_pathway_enrichment_lollipop.pdf"
)

ggsave(
  filename = pdf_filename,
  plot = lollipop_plot,
  device = "pdf",
  width = 10,
  height = 8,
  units = "in"
)
