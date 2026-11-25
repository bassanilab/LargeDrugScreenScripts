library(protti)

setwd("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/")

data1=read_protti("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/PTM_probing_search/HLA_class1_PTM_probing_AUGUST_official_highcell.txt")
data2=read_protti("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/PTM_probing_search/HLA_class1_PTM_probing_AUGUST_official_lowcell.txt")

library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

#data1=data1[-9]


modified_sequence_col <- "t_t_t_eg_modified_sequence"
reproducibility_fraction <- 0.60
minimum_ptm_peptidoforms <- 50

forced_phospho_ptms <- c(
  "S[Phospho (STY)]",
  "T[Phospho (STY)]",
  "Y[Phospho (STY)]"
)

N1=79
N2=15

# ==============================================================================
# Helper functions
# ==============================================================================

validate_input <- function(data, N, data_name) {
  
  if (!is.data.frame(data)) {
    stop(data_name, " must be a data.frame.")
  }
  
  if (
    length(N) != 1L ||
    !is.numeric(N) ||
    is.na(N) ||
    !is.finite(N) ||
    N != as.integer(N) ||
    N < 1L ||
    N > ncol(data)
  ) {
    stop("Invalid quantitative column count for ", data_name, ".")
  }
  
  if (!(modified_sequence_col %in% colnames(data))) {
    stop(
      "Column '",
      modified_sequence_col,
      "' was not found in ",
      data_name,
      "."
    )
  }
  
  sample_columns <- colnames(data)[seq_len(as.integer(N))]
  
  if (any(is.na(sample_columns)) || any(sample_columns == "")) {
    stop("Quantitative sample columns must have valid names.")
  }
  
  if (anyDuplicated(sample_columns)) {
    stop("Quantitative sample column names must be unique.")
  }
  
  if (any(!str_detect(sample_columns, "_"))) {
    stop(
      "All quantitative sample columns must follow condition_replicate in ",
      data_name,
      "."
    )
  }
  
  invisible(TRUE)
}


extract_ptms <- function(sequences) {
  
  sequences <- as.character(sequences)
  sequences[is.na(sequences)] <- ""
  
  # Recognizes residue linked modifications:
  # M[Oxidation]
  # S[Phospho (STY)]
  #
  # Also recognizes modifications without a preceding residue:
  # [Acetyl (Protein N-term)]
  # [Carboxymethyl (Any N-term)]
  
  pattern <- paste0(
    "(?:",
    "[A-Za-z]\\[[^\\]]+\\]",
    "|",
    "(?<![A-Za-z])\\[[^\\]]+\\]",
    ")"
  )
  
  str_extract_all(
    sequences,
    regex(pattern)
  )
}


jaccard_distance <- function(set_a, set_b) {
  
  set_a <- unique(set_a)
  set_b <- unique(set_b)
  
  union_set <- union(set_a, set_b)
  
  if (length(union_set) == 0L) {
    return(NA_real_)
  }
  
  1 - length(intersect(set_a, set_b)) / length(union_set)
}


build_detection_matrix <- function(data, sample_columns) {
  
  sequences <- as.character(
    data[[modified_sequence_col]]
  )
  
  keep <- !is.na(sequences) & trimws(sequences) != ""
  
  sequences <- sequences[keep]
  
  quantitative_data <- data[
    keep,
    sample_columns,
    drop = FALSE
  ]
  
  if (length(sequences) == 0L) {
    stop("No valid modified peptide sequences were found.")
  }
  
  detection <- as.matrix(
    data.frame(
      lapply(
        quantitative_data,
        function(x) {
          !is.na(x) & x != 0
        }
      ),
      check.names = FALSE
    )
  )
  
  storage.mode(detection) <- "integer"
  
  # Fast aggregation of duplicate modified sequences
  aggregated <- rowsum(
    detection,
    group = sequences,
    reorder = FALSE
  )
  
  aggregated > 0
}


build_ptm_map <- function(sequences) {
  
  sequences <- unique(
    as.character(sequences)
  )
  
  sequences <- sequences[
    !is.na(sequences) &
      trimws(sequences) != ""
  ]
  
  if (length(sequences) == 0L) {
    return(list())
  }
  
  ptm_lists <- extract_ptms(sequences)
  n_ptms <- lengths(ptm_lists)
  
  if (sum(n_ptms) == 0L) {
    return(list())
  }
  
  pairs <- data.frame(
    PTM = unlist(
      ptm_lists,
      use.names = FALSE
    ),
    sequence = rep(
      sequences,
      n_ptms
    ),
    stringsAsFactors = FALSE
  )
  
  pairs <- unique(pairs)
  
  lapply(
    split(
      pairs$sequence,
      pairs$PTM
    ),
    unique
  )
}


get_ptm_count_df <- function(ptm_map) {
  
  if (length(ptm_map) == 0L) {
    
    return(
      data.frame(
        PTM = character(0),
        total_unique_peptidoforms = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  counts <- lengths(ptm_map)
  
  result <- data.frame(
    PTM = names(counts),
    total_unique_peptidoforms = as.integer(counts),
    stringsAsFactors = FALSE
  )
  
  result <- result[
    order(
      result$total_unique_peptidoforms,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  rownames(result) <- NULL
  
  result
}


get_ptm_count <- function(ptm, count_df) {
  
  index <- match(
    ptm,
    count_df$PTM
  )
  
  if (is.na(index)) {
    return(0L)
  }
  
  as.integer(
    count_df$total_unique_peptidoforms[index]
  )
}


# ==============================================================================
# Process one dataset
# ==============================================================================

process_dataset <- function(
    data,
    N,
    data_name,
    control_conditions
) {
  
  validate_input(
    data,
    N,
    data_name
  )
  
  N <- as.integer(N)
  
  sample_columns <- colnames(data)[seq_len(N)]
  
  cat(
    "\n================ ",
    data_name,
    " ================\n",
    sep = ""
  )
  
  cat("\nSample columns used:\n")
  print(sample_columns)
  
  
  # Convert only quantitative sample columns to numeric
  
  for (sample_column in sample_columns) {
    
    original <- data[[sample_column]]
    
    numeric_values <- suppressWarnings(
      as.numeric(
        as.character(original)
      )
    )
    
    bad_values <- (
      !is.na(original) &
        trimws(as.character(original)) != "" &
        is.na(numeric_values)
    )
    
    if (any(bad_values)) {
      
      warning(
        "Non numeric values converted to NA in ",
        data_name,
        ": ",
        sample_column
      )
    }
    
    data[[sample_column]] <- numeric_values
  }
  
  
  # Parse condition and replicate names
  
  original_conditions <- str_to_lower(
    str_extract(
      sample_columns,
      "^[^_]+"
    )
  )
  
  conditions <- original_conditions
  
  conditions[
    original_conditions %in%
      str_to_lower(control_conditions)
  ] <- "dms"
  
  sample_info <- data.frame(
    sample = sample_columns,
    original_condition = original_conditions,
    condition = conditions,
    replicate = str_replace(
      sample_columns,
      "^[^_]+_",
      ""
    ),
    stringsAsFactors = FALSE
  )
  
  cat("\nCondition structure:\n")
  print(sample_info)
  
  if (!("dms" %in% conditions)) {
    stop(
      "No DMS control was identified in ",
      data_name,
      "."
    )
  }
  
  condition_order <- unique(conditions)
  
  perturbation_conditions <- condition_order[
    condition_order != "dms"
  ]
  
  if (length(perturbation_conditions) == 0L) {
    stop(
      "No perturbation conditions were found in ",
      data_name,
      "."
    )
  }
  
  
  # Build unique peptidoform detection matrix
  
  cat("\nBuilding peptidoform detection matrix...\n")
  
  detection_matrix <- build_detection_matrix(
    data,
    sample_columns
  )
  
  cat("Detection matrix complete.\n")
  
  
  # Apply 60 percent reproducibility criterion
  
  reproducible_sets <- setNames(
    vector(
      "list",
      length(condition_order)
    ),
    condition_order
  )
  
  reproducibility_summary <- data.frame(
    condition = character(0),
    replicates = integer(0),
    minimum_required = integer(0),
    reproducible_peptidoforms = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (condition_name in condition_order) {
    
    condition_samples <- sample_info$sample[
      sample_info$condition == condition_name
    ]
    
    number_of_replicates <- length(
      condition_samples
    )
    
    threshold <- if (condition_name == "snt") 1L else ceiling(reproducibility_fraction * number_of_replicates)
    
    detection_counts <- rowSums(
      detection_matrix[
        ,
        condition_samples,
        drop = FALSE
      ]
    )
    
    reproducible_sequences <- rownames(
      detection_matrix
    )[
      detection_counts >= threshold
    ]
    
    reproducible_sets[[condition_name]] <-
      reproducible_sequences
    
    reproducibility_summary <- rbind(
      reproducibility_summary,
      data.frame(
        condition = condition_name,
        replicates = number_of_replicates,
        minimum_required = threshold,
        reproducible_peptidoforms =
          length(reproducible_sequences),
        stringsAsFactors = FALSE
      )
    )
  }
  
  cat(
    "\nReproducibly detected peptides per condition:\n"
  )
  
  print(
    reproducibility_summary
  )
  
  
  # Build PTM to peptidoform map from complete dataset
  
  cat("\nBuilding PTM map...\n")
  
  ptm_map <- build_ptm_map(
    data[[modified_sequence_col]]
  )
  
  cat("PTM map complete.\n")
  
  total_ptm_counts <- get_ptm_count_df(
    ptm_map
  )
  
  cat("\nPTM counts before filtering:\n")
  print(total_ptm_counts)
  
  
  # PTM filtering
  
  retained_ptms <- unique(
    c(
      total_ptm_counts$PTM[
        total_ptm_counts$total_unique_peptidoforms >=
          minimum_ptm_peptidoforms
      ],
      forced_phospho_ptms
    )
  )
  
  removed_ptms <- setdiff(
    total_ptm_counts$PTM,
    retained_ptms
  )
  
  cat("\nPTMs retained:\n")
  print(retained_ptms)
  
  cat("\nPTMs removed:\n")
  
  if (length(removed_ptms) == 0L) {
    cat("None\n")
  } else {
    print(removed_ptms)
  }
  
  
  # Jaccard distances
  
  cat("\nCalculating Jaccard distances...\n")
  
  jaccard_matrix <- matrix(
    NA_real_,
    nrow = length(retained_ptms),
    ncol = length(perturbation_conditions),
    dimnames = list(
      retained_ptms,
      perturbation_conditions
    )
  )
  
  for (ptm in retained_ptms) {
    
    ptm_sequences <- ptm_map[[ptm]]
    
    if (is.null(ptm_sequences)) {
      ptm_sequences <- character(0)
    }
    
    control_set <- intersect(
      reproducible_sets[["dms"]],
      ptm_sequences
    )
    
    for (condition_name in perturbation_conditions) {
      
      treatment_set <- intersect(
        reproducible_sets[[condition_name]],
        ptm_sequences
      )
      
      jaccard_matrix[
        ptm,
        condition_name
      ] <- jaccard_distance(
        control_set,
        treatment_set
      )
    }
  }
  
  cat("Jaccard calculation complete.\n")
  
  
  list(
    sample_columns = sample_columns,
    sample_info = sample_info,
    reproducible_sets = reproducible_sets,
    reproducibility_summary = reproducibility_summary,
    detection_matrix = detection_matrix,
    ptm_map = ptm_map,
    total_ptm_counts = total_ptm_counts,
    retained_ptms = retained_ptms,
    removed_ptms = removed_ptms,
    jaccard_unaligned = jaccard_matrix
  )
}


# ==============================================================================
# Run analysis
# ==============================================================================

results_data1 <- process_dataset(
  data = data1,
  N = N1,
  data_name = "data1",
  control_conditions = c(
    "dms3",
    "dms4"
  )
)

results_data2 <- process_dataset(
  data = data2,
  N = N2,
  data_name = "data2",
  control_conditions = "dms2"
)


# Union of retained PTMs

union_ptms <- union(
  results_data1$retained_ptms,
  results_data2$retained_ptms
)

counts_data1 <- vapply(
  union_ptms,
  get_ptm_count,
  integer(1),
  count_df = results_data1$total_ptm_counts
)

counts_data2 <- vapply(
  union_ptms,
  get_ptm_count,
  integer(1),
  count_df = results_data2$total_ptm_counts
)


# Order PTMs by prevalence without clustering

ptm_order <- order(
  pmax(
    counts_data1,
    counts_data2
  ),
  counts_data1 + counts_data2,
  decreasing = TRUE
)

all_retained_ptms <- union_ptms[
  ptm_order
]


align_matrix <- function(matrix_in, ptms) {
  
  result <- matrix(
    NA_real_,
    nrow = length(ptms),
    ncol = ncol(matrix_in),
    dimnames = list(
      ptms,
      colnames(matrix_in)
    )
  )
  
  matching_ptms <- intersect(
    ptms,
    rownames(matrix_in)
  )
  
  if (length(matching_ptms) > 0L) {
    
    result[
      matching_ptms,
      colnames(matrix_in)
    ] <- matrix_in[
      matching_ptms,
      ,
      drop = FALSE
    ]
  }
  
  result
}


jaccard_data1 <- align_matrix(
  results_data1$jaccard_unaligned,
  all_retained_ptms
)

jaccard_data2 <- align_matrix(
  results_data2$jaccard_unaligned,
  all_retained_ptms
)


jaccard_df_data1 <- data.frame(
  PTM = rownames(jaccard_data1),
  jaccard_data1,
  check.names = FALSE,
  row.names = NULL
)

jaccard_df_data2 <- data.frame(
  PTM = rownames(jaccard_data2),
  jaccard_data2,
  check.names = FALSE,
  row.names = NULL
)


total_mod_counts_df_data1 <- data.frame(
  PTM = all_retained_ptms,
  total_unique_peptidoforms = vapply(
    all_retained_ptms,
    get_ptm_count,
    integer(1),
    count_df = results_data1$total_ptm_counts
  ),
  stringsAsFactors = FALSE
)

total_mod_counts_df_data2 <- data.frame(
  PTM = all_retained_ptms,
  total_unique_peptidoforms = vapply(
    all_retained_ptms,
    get_ptm_count,
    integer(1),
    count_df = results_data2$total_ptm_counts
  ),
  stringsAsFactors = FALSE
)


results_data1$jaccard <- jaccard_data1
results_data2$jaccard <- jaccard_data2


# ==============================================================================
# Shared color scale
# ==============================================================================

finite_values <- c(
  jaccard_data1[
    is.finite(jaccard_data1)
  ],
  jaccard_data2[
    is.finite(jaccard_data2)
  ]
)

if (length(finite_values) == 0L) {
  stop(
    "No finite Jaccard distances were calculated."
  )
}

maximum_observed_jaccard <- max(
  finite_values
)

cat(
  "\nMaximum observed Jaccard distance: ",
  maximum_observed_jaccard,
  "\n",
  sep = ""
)


if (maximum_observed_jaccard > 0) {
  
  legend_breaks <- seq(
    0,
    1,
    length.out = 5
  )
  
  jaccard_col_fun <- colorRamp2(
    legend_breaks,
    c(
      "white",
      "yellow",
      "orange",
      "red",
      "darkred"
    )
  )
  
} else {
  
  legend_breaks <- 0
  
  jaccard_col_fun <- colorRamp2(
    c(
      0,
      .Machine$double.eps
    ),
    c(
      "white",
      "white"
    )
  )
}


cell_labels <- function(matrix_in) {
  
  function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    value <- matrix_in[
      i,
      j
    ]
    
    if (is.na(value)) {
      return(
        invisible(NULL)
      )
    }
    
    text_color <- if (
      maximum_observed_jaccard > 0 &&
      value / maximum_observed_jaccard >= 0.65
    ) {
      "white"
    } else {
      "black"
    }
    
    grid.text(
      sprintf(
        "%.2f",
        value
      ),
      x = x,
      y = y,
      gp = gpar(
        col = text_color,
        fontsize = 8
      )
    )
  }
}


# ==============================================================================
# Heatmap and barplots
# ==============================================================================

bar_values_data1 <-
  total_mod_counts_df_data1$total_unique_peptidoforms

bar_values_data2 <-
  total_mod_counts_df_data2$total_unique_peptidoforms

maximum_bar_count <- max(
  c(
    bar_values_data1,
    bar_values_data2,
    1
  )
)


# Reserve enough width for full modification names

ptm_name_width <- max_text_width(
  rownames(jaccard_data2),
  gp = gpar(
    fontsize = 8
  )
) + unit(
  4,
  "mm"
)


heatmap_data2 <- Heatmap(
  jaccard_data2,
  
  name = "Jaccard distance",
  
  col = jaccard_col_fun,
  
  na_col = "lightgrey",
  
  cluster_rows = FALSE,
  
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  
  row_names_side = "left",
  
  row_labels = rownames(jaccard_data2),
  
  row_names_gp = gpar(
    fontsize = 8
  ),
  
  row_names_max_width =
    ptm_name_width,
  
  show_column_names = TRUE,
  
  column_names_rot = 45,
  
  column_names_gp = gpar(
    fontsize = 9
  ),
  
  column_title =
    "Low cell input",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  cell_fun =
    cell_labels(
      jaccard_data2
    ),
  
  heatmap_legend_param = list(
    title =
      "Jaccard distance",
    
    at =
      legend_breaks,
    
    labels =
      sprintf(
        "%.2f",
        legend_breaks
      )
  ),
  
  show_heatmap_legend = TRUE
)


bar_annotation_data2 <- rowAnnotation(
  
  `Total peptidoforms` =
    anno_barplot(
      
      bar_values_data2,
      
      baseline = 0,
      
      ylim = c(
        0,
        maximum_bar_count
      ),
      
      border = FALSE,
      
      gp = gpar(
        fill = "grey40"
      ),
      
      axis = TRUE
    ),
  
  width = unit(
    28,
    "mm"
  )
)


heatmap_data1 <- Heatmap(
  jaccard_data1,
  
  name =
    "Jaccard distance high",
  
  col =
    jaccard_col_fun,
  
  na_col =
    "lightgrey",
  
  cluster_rows =
    FALSE,
  
  cluster_columns =
    FALSE,
  
  show_row_names =
    FALSE,
  
  show_column_names =
    TRUE,
  
  column_names_rot =
    45,
  
  column_names_gp = gpar(
    fontsize = 9
  ),
  
  column_title =
    "High cell input",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  cell_fun =
    cell_labels(
      jaccard_data1
    ),
  
  show_heatmap_legend =
    FALSE
)


bar_annotation_data1 <- rowAnnotation(
  
  `Total peptidoforms` =
    anno_barplot(
      
      bar_values_data1,
      
      baseline = 0,
      
      ylim = c(
        0,
        maximum_bar_count
      ),
      
      border = FALSE,
      
      gp = gpar(
        fill = "grey40"
      ),
      
      axis = TRUE
    ),
  
  width = unit(
    28,
    "mm"
  )
)


combined_heatmap <-
  heatmap_data2 +
  bar_annotation_data2 +
  heatmap_data1 +
  bar_annotation_data1


# ==============================================================================
# Summary statistics
# ==============================================================================

all_jaccard_values <- c(
  as.vector(
    jaccard_data1
  ),
  as.vector(
    jaccard_data2
  )
)

all_jaccard_values <- all_jaccard_values[
  is.finite(
    all_jaccard_values
  )
]

median_jaccard <- median(
  all_jaccard_values
)

top_decile_jaccard <- quantile(
  all_jaccard_values,
  probs = 0.90,
  names = FALSE
)

cat(
  "\nMedian Jaccard distance: ",
  median_jaccard,
  "\n",
  sep = ""
)

cat(
  "Top decile Jaccard threshold: ",
  top_decile_jaccard,
  "\n",
  sep = ""
)

cat(
  "Number of finite Jaccard values: ",
  length(all_jaccard_values),
  "\n",
  sep = ""
)


# ==============================================================================
# Display and save figure
# ==============================================================================

# Close any previously open graphics devices.
# This also releases PDF files that may still be locked by R.

graphics.off()


# Display heatmap in RStudio Plots pane

draw(
  combined_heatmap,
  heatmap_legend_side = "right",
  merge_legends = TRUE,
  padding = unit(
    c(
      5,
      5,
      5,
      5
    ),
    "mm"
  )
)


# Save heatmap to PDF

pdf(
  "heatmap_jaccard_side_by_side_HLAI_50_limit.pdf",
  width = 16,
  height = 10
)

tryCatch(
  {
    draw(
      combined_heatmap,
      heatmap_legend_side = "right",
      merge_legends = TRUE,
      padding = unit(
        c(
          5,
          5,
          5,
          5
        ),
        "mm"
      )
    )
  },
  finally = {
    dev.off()
  }
)





