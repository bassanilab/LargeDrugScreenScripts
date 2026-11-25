library(protti)
setwd("//FILERCH/data/PEPMS/Common/Common_beta/Laboratory/Manuscripts/Nicolas_papers/LDS/Method_and_Results/Heatmaps")
df1=read_protti("low_cells_normalized_zscored_HLAII.txt")
df2=read_protti("high_cells_normalized_zscored_HLAII.txt")
###############################################################################
# Two independently clustered peptide heatmaps displayed side by side
###############################################################################


###############################################################################
# User settings
###############################################################################

# df1 intensity columns are columns 1 through N1
N1 <- 15

# df2 intensity columns are columns 1 through N2
N2 <- 81

# Peptide sequence column
peptide_column <- "t_t_t_eg_modified_sequence"

# Panel titles
dataset_title1 <- "Dataset 1"
dataset_title2 <- "Dataset 2"

# FALSE keeps individual replicates and groups them by condition
# TRUE averages replicates and displays one column per condition
collapse_replicates1 <- FALSE
collapse_replicates2 <- FALSE

# Missing value replacement used only for clustering
#
# Options:
# "zero"
# "row_mean"
clustering_na_fill <- "zero"

# Common colour scale for both heatmaps
z_limit <- 3

# Missing cells remain visible with this colour
missing_value_colour <- "#BDBDBD"

# Peptide names are generally unreadable with approximately 20,000 rows
show_peptide_names <- FALSE

# Rows are still clustered when the dendrogram is hidden
show_row_dendrogram <- TRUE

# Panel width options:
#
# "proportional" gives more width to the dataset with more columns
# "equal" gives both panels equal width
panel_width_mode <- "proportional"

# Output files
output_pdf <- "two_peptide_heatmaps.pdf"
output_row_order1 <- "df1_peptide_row_order.csv"
output_row_order2 <- "df2_peptide_row_order.csv"

# PDF dimensions in inches
pdf_width <- 16
pdf_height <- 10


###############################################################################
# Check required packages
###############################################################################

required_packages <- c(
  "ComplexHeatmap",
  "circlize",
  "scico",
  "fastcluster"
)

installed_packages <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (any(!installed_packages)) {
  
  stop(
    "The following packages are missing: ",
    paste(
      required_packages[!installed_packages],
      collapse = ", "
    ),
    call. = FALSE
  )
}


###############################################################################
# Check that df1 and df2 exist
###############################################################################

if (!exists("df1")) {
  stop("The object 'df1' was not found.", call. = FALSE)
}

if (!exists("df2")) {
  stop("The object 'df2' was not found.", call. = FALSE)
}


###############################################################################
# Create the shared colour function
###############################################################################

if (
  length(z_limit) != 1L ||
  is.na(z_limit) ||
  !is.finite(z_limit) ||
  z_limit <= 0
) {
  stop("'z_limit' must be one positive number.", call. = FALSE)
}

heatmap_palette <- scico::scico(
  n = 101,
  palette = "vik"
)

colour_function <- circlize::colorRamp2(
  seq(
    from = -z_limit,
    to = z_limit,
    length.out = length(heatmap_palette)
  ),
  heatmap_palette
)


###############################################################################
# Function that prepares and clusters one dataset
###############################################################################

prepare_peptide_heatmap <- function(
    data,
    N,
    dataset_label,
    heatmap_name,
    peptide_column,
    collapse_replicates,
    clustering_na_fill,
    colour_function,
    missing_value_colour,
    show_peptide_names,
    show_row_dendrogram
) {
  
  dataset_error <- function(...) {
    
    stop(
      "[",
      dataset_label,
      "] ",
      ...,
      call. = FALSE
    )
  }
  
  dataset_warning <- function(...) {
    
    warning(
      "[",
      dataset_label,
      "] ",
      ...,
      call. = FALSE
    )
  }
  
  
  ###########################################################################
  # Validate the dataset
  ###########################################################################
  
  if (!is.data.frame(data)) {
    dataset_error("The object must be a data frame.")
  }
  
  if (nrow(data) < 2L) {
    dataset_error(
      "At least two peptide rows are required for clustering."
    )
  }
  
  if (
    length(N) != 1L ||
    is.na(N) ||
    !is.numeric(N) ||
    N < 1L ||
    N > ncol(data)
  ) {
    dataset_error(
      "'N' must be a single integer between 1 and ncol(data)."
    )
  }
  
  N <- as.integer(N)
  
  if (!peptide_column %in% names(data)) {
    dataset_error(
      "The peptide column was not found: ",
      peptide_column
    )
  }
  
  intensity_columns <- seq_len(N)
  
  if (peptide_column %in% names(data)[intensity_columns]) {
    dataset_error(
      "The peptide sequence column is included among columns 1 through N. ",
      "N must include intensity columns only."
    )
  }
  
  numeric_columns <- vapply(
    data[, intensity_columns, drop = FALSE],
    is.numeric,
    logical(1)
  )
  
  if (any(!numeric_columns)) {
    
    dataset_error(
      "The following intensity columns are not numeric: ",
      paste(
        names(data)[intensity_columns][!numeric_columns],
        collapse = ", "
      )
    )
  }
  
  
  ###########################################################################
  # Construct the intensity matrix
  ###########################################################################
  
  intensity_matrix <- as.matrix(
    data[, intensity_columns, drop = FALSE]
  )
  
  storage.mode(intensity_matrix) <- "double"
  
  # Convert NaN values to standard NA values
  intensity_matrix[is.nan(intensity_matrix)] <- NA_real_
  
  if (any(is.infinite(intensity_matrix))) {
    
    dataset_warning(
      "Infinite intensity values were converted to NA."
    )
    
    intensity_matrix[
      is.infinite(intensity_matrix)
    ] <- NA_real_
  }
  
  sample_names <- colnames(intensity_matrix)
  
  if (
    is.null(sample_names) ||
    any(is.na(sample_names)) ||
    any(!nzchar(sample_names))
  ) {
    dataset_error(
      "Every intensity column must have a nonempty column name."
    )
  }
  
  
  ###########################################################################
  # Extract peptide labels
  ###########################################################################
  
  peptide_labels <- as.character(
    data[[peptide_column]]
  )
  
  missing_peptide_labels <- (
    is.na(peptide_labels) |
      !nzchar(peptide_labels)
  )
  
  peptide_labels[missing_peptide_labels] <- paste0(
    "unnamed_peptide_",
    which(missing_peptide_labels)
  )
  
  # Internal row names are made unique
  rownames(intensity_matrix) <- sprintf(
    "%s_peptide_row_%06d",
    heatmap_name,
    seq_len(nrow(intensity_matrix))
  )
  
  
  ###########################################################################
  # Extract condition and replicate names
  ###########################################################################
  
  has_separator <- grepl(
    "_",
    sample_names,
    fixed = TRUE
  )
  
  if (any(!has_separator)) {
    
    dataset_error(
      "Every intensity column must contain an underscore. ",
      "Columns without an underscore: ",
      paste(
        sample_names[!has_separator],
        collapse = ", "
      )
    )
  }
  
  # Everything before the first underscore is the condition
  condition <- sub(
    "_.*$",
    "",
    sample_names
  )
  
  # Everything after the first underscore is the replicate
  replicate <- sub(
    "^[^_]*_",
    "",
    sample_names
  )
  
  if (any(!nzchar(condition))) {
    dataset_error(
      "At least one column has an empty condition name."
    )
  }
  
  if (any(!nzchar(replicate))) {
    dataset_error(
      "At least one column has an empty replicate name."
    )
  }
  
  # Preserve the order in which conditions first appear
  condition_levels <- unique(condition)
  
  condition_factor <- factor(
    condition,
    levels = condition_levels
  )
  
  # Place all replicates from the same condition together
  column_order <- order(
    as.integer(condition_factor),
    seq_along(condition_factor)
  )
  
  intensity_matrix <- intensity_matrix[
    ,
    column_order,
    drop = FALSE
  ]
  
  condition_factor <- factor(
    condition[column_order],
    levels = condition_levels
  )
  
  replicate <- replicate[column_order]
  sample_names <- sample_names[column_order]
  
  
  ###########################################################################
  # Optionally average replicates within each condition
  ###########################################################################
  
  if (collapse_replicates) {
    
    condition_character <- as.character(
      condition_factor
    )
    
    heatmap_matrix <- vapply(
      condition_levels,
      FUN = function(current_condition) {
        
        current_columns <- (
          condition_character == current_condition
        )
        
        condition_mean <- rowMeans(
          intensity_matrix[
            ,
            current_columns,
            drop = FALSE
          ],
          na.rm = TRUE
        )
        
        # All missing replicates produce NaN
        # Keep these cells missing in the displayed heatmap
        condition_mean[
          is.nan(condition_mean)
        ] <- NA_real_
        
        condition_mean
      },
      FUN.VALUE = numeric(
        nrow(intensity_matrix)
      )
    )
    
    heatmap_matrix <- as.matrix(
      heatmap_matrix
    )
    
    rownames(heatmap_matrix) <- rownames(
      intensity_matrix
    )
    
    colnames(heatmap_matrix) <- condition_levels
    
    column_split_value <- NULL
    displayed_column_labels <- condition_levels
    
  } else {
    
    heatmap_matrix <- intensity_matrix
    
    column_split_value <- condition_factor
    displayed_column_labels <- replicate
  }
  
  
  ###########################################################################
  # Prepare a complete matrix for clustering
  ###########################################################################
  
  # Missing values remain NA in heatmap_matrix.
  # They are replaced only in clustering_matrix.
  
  clustering_matrix <- heatmap_matrix
  
  if (clustering_na_fill == "zero") {
    
    clustering_matrix[
      is.na(clustering_matrix)
    ] <- 0
    
  } else if (clustering_na_fill == "row_mean") {
    
    peptide_means <- rowMeans(
      clustering_matrix,
      na.rm = TRUE
    )
    
    # Rows containing only missing values receive zero
    peptide_means[
      !is.finite(peptide_means)
    ] <- 0
    
    for (
      column_index in seq_len(
        ncol(clustering_matrix)
      )
    ) {
      
      missing_rows <- is.na(
        clustering_matrix[, column_index]
      )
      
      clustering_matrix[
        missing_rows,
        column_index
      ] <- peptide_means[missing_rows]
    }
    
  } else {
    
    dataset_error(
      "'clustering_na_fill' must be either ",
      "'zero' or 'row_mean'."
    )
  }
  
  storage.mode(clustering_matrix) <- "double"
  
  
  ###########################################################################
  # Independent hierarchical clustering
  ###########################################################################
  
  row_hclust <- fastcluster::hclust.vector(
    clustering_matrix,
    method = "ward",
    metric = "euclidean"
  )
  
  
  ###########################################################################
  # Construct the heatmap
  ###########################################################################
  
  heatmap_object <- ComplexHeatmap::Heatmap(
    
    heatmap_matrix,
    
    name = heatmap_name,
    
    col = colour_function,
    
    na_col = missing_value_colour,
    
    show_heatmap_legend = FALSE,
    
    # Independent row clustering
    cluster_rows = row_hclust,
    row_dend_reorder = FALSE,
    show_row_dend = show_row_dendrogram,
    row_dend_width = grid::unit(
      14,
      "mm"
    ),
    
    # Columns are grouped but not clustered
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    
    column_split = column_split_value,
    cluster_column_slices = FALSE,
    column_gap = grid::unit(
      2,
      "mm"
    ),
    
    # Replicate or condition labels
    column_labels = displayed_column_labels,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(
      fontsize = 7
    ),
    column_names_max_height = grid::unit(
      35,
      "mm"
    ),
    
    # Condition headings
    column_title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold",
      fill = "#F2F2F2",
      col = "black"
    ),
    
    # Peptide labels
    row_labels = peptide_labels,
    show_row_names = show_peptide_names,
    row_names_gp = grid::gpar(
      fontsize = 4
    ),
    row_names_max_width = grid::unit(
      60,
      "mm"
    ),
    
    # No cell borders in the dense matrix
    rect_gp = grid::gpar(
      col = NA
    ),
    
    # Rasterize the heatmap body
    use_raster = TRUE,
    raster_device = "png",
    raster_quality = 4,
    
    border = FALSE
  )
  
  
  ###########################################################################
  # Prepare clustered row order table
  ###########################################################################
  
  clustered_row_order <- row_hclust$order
  
  row_order_table <- data.frame(
    
    heatmap_position = seq_along(
      clustered_row_order
    ),
    
    original_df_row = clustered_row_order,
    
    peptide = peptide_labels[
      clustered_row_order
    ],
    
    stringsAsFactors = FALSE
  )
  
  
  ###########################################################################
  # Return all required objects
  ###########################################################################
  
  list(
    
    heatmap = heatmap_object,
    
    row_hclust = row_hclust,
    
    row_order_table = row_order_table,
    
    number_of_rows = nrow(
      heatmap_matrix
    ),
    
    number_of_columns = ncol(
      heatmap_matrix
    )
  )
}


###############################################################################
# Compute df1 and df2 separately
###############################################################################

result1 <- prepare_peptide_heatmap(
  
  data = df1,
  
  N = N1,
  
  dataset_label = dataset_title1,
  
  heatmap_name = "df1_z_score",
  
  peptide_column = peptide_column,
  
  collapse_replicates = collapse_replicates1,
  
  clustering_na_fill = clustering_na_fill,
  
  colour_function = colour_function,
  
  missing_value_colour = missing_value_colour,
  
  show_peptide_names = show_peptide_names,
  
  show_row_dendrogram = show_row_dendrogram
)


result2 <- prepare_peptide_heatmap(
  
  data = df2,
  
  N = N2,
  
  dataset_label = dataset_title2,
  
  heatmap_name = "df2_z_score",
  
  peptide_column = peptide_column,
  
  collapse_replicates = collapse_replicates2,
  
  clustering_na_fill = clustering_na_fill,
  
  colour_function = colour_function,
  
  missing_value_colour = missing_value_colour,
  
  show_peptide_names = show_peptide_names,
  
  show_row_dendrogram = show_row_dendrogram
)


###############################################################################
# Save the independently clustered peptide orders
###############################################################################

utils::write.csv(
  result1$row_order_table,
  file = output_row_order1,
  row.names = FALSE
)

utils::write.csv(
  result2$row_order_table,
  file = output_row_order2,
  row.names = FALSE
)


###############################################################################
# Create one shared legend
###############################################################################

z_score_legend <- ComplexHeatmap::Legend(
  
  title = "Peptide\nz score",
  
  col_fun = colour_function,
  
  at = c(
    -z_limit,
    0,
    z_limit
  ),
  
  labels = c(
    sprintf(
      "\u2264 %.1f",
      -z_limit
    ),
    "0",
    sprintf(
      "\u2265 %.1f",
      z_limit
    )
  ),
  
  legend_height = grid::unit(
    45,
    "mm"
  ),
  
  title_gp = grid::gpar(
    fontsize = 9,
    fontface = "bold"
  ),
  
  labels_gp = grid::gpar(
    fontsize = 8
  )
)


missing_value_legend <- ComplexHeatmap::Legend(
  
  labels = "Missing",
  
  legend_gp = grid::gpar(
    fill = missing_value_colour,
    col = NA
  ),
  
  labels_gp = grid::gpar(
    fontsize = 8
  )
)


shared_legend <- ComplexHeatmap::packLegend(
  
  z_score_legend,
  
  missing_value_legend,
  
  direction = "vertical",
  
  gap = grid::unit(
    5,
    "mm"
  )
)


###############################################################################
# Determine relative panel widths
###############################################################################

if (panel_width_mode == "proportional") {
  
  panel_widths <- c(
    result1$number_of_columns,
    result2$number_of_columns
  )
  
} else if (panel_width_mode == "equal") {
  
  panel_widths <- c(
    1,
    1
  )
  
} else {
  
  stop(
    "'panel_width_mode' must be either ",
    "'proportional' or 'equal'.",
    call. = FALSE
  )
}


###############################################################################
# Draw both heatmaps side by side
###############################################################################

grDevices::pdf(
  
  file = output_pdf,
  
  width = pdf_width,
  
  height = pdf_height,
  
  useDingbats = FALSE
)


tryCatch(
  
  {
    
    grid::grid.newpage()
    
    comparison_layout <- grid::grid.layout(
      
      nrow = 2,
      
      ncol = 5,
      
      heights = grid::unit.c(
        
        grid::unit(
          8,
          "mm"
        ),
        
        grid::unit(
          1,
          "null"
        )
      ),
      
      widths = grid::unit.c(
        
        grid::unit(
          panel_widths[1],
          "null"
        ),
        
        grid::unit(
          5,
          "mm"
        ),
        
        grid::unit(
          panel_widths[2],
          "null"
        ),
        
        grid::unit(
          5,
          "mm"
        ),
        
        grid::unit(
          20,
          "mm"
        )
      )
    )
    
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout = comparison_layout,
        
        name = "two_heatmap_layout"
      )
    )
    
    
    #########################################################################
    # df1 title
    #########################################################################
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout.pos.row = 1,
        
        layout.pos.col = 1
      )
    )
    
    grid::grid.text(
      
      dataset_title1,
      
      gp = grid::gpar(
        
        fontsize = 12,
        
        fontface = "bold"
      )
    )
    
    grid::popViewport()
    
    
    #########################################################################
    # df2 title
    #########################################################################
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout.pos.row = 1,
        
        layout.pos.col = 3
      )
    )
    
    grid::grid.text(
      
      dataset_title2,
      
      gp = grid::gpar(
        
        fontsize = 12,
        
        fontface = "bold"
      )
    )
    
    grid::popViewport()
    
    
    #########################################################################
    # df1 heatmap
    #########################################################################
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout.pos.row = 2,
        
        layout.pos.col = 1,
        
        clip = "off"
      )
    )
    
    ComplexHeatmap::draw(
      
      result1$heatmap,
      
      newpage = FALSE,
      
      show_heatmap_legend = FALSE,
      
      show_annotation_legend = FALSE,
      
      padding = grid::unit(
        
        c(
          2,
          2,
          2,
          2
        ),
        
        "mm"
      )
    )
    
    grid::popViewport()
    
    
    #########################################################################
    # df2 heatmap
    #########################################################################
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout.pos.row = 2,
        
        layout.pos.col = 3,
        
        clip = "off"
      )
    )
    
    ComplexHeatmap::draw(
      
      result2$heatmap,
      
      newpage = FALSE,
      
      show_heatmap_legend = FALSE,
      
      show_annotation_legend = FALSE,
      
      padding = grid::unit(
        
        c(
          2,
          2,
          2,
          2
        ),
        
        "mm"
      )
    )
    
    grid::popViewport()
    
    
    #########################################################################
    # Shared legend
    #########################################################################
    
    grid::pushViewport(
      
      grid::viewport(
        
        layout.pos.row = 2,
        
        layout.pos.col = 5
      )
    )
    
    ComplexHeatmap::draw(
      
      shared_legend,
      
      x = grid::unit(
        0.5,
        "npc"
      ),
      
      y = grid::unit(
        0.5,
        "npc"
      ),
      
      just = c(
        "center",
        "center"
      )
    )
    
    grid::popViewport()
    
    
    grid::popViewport()
  },
  
  finally = {
    
    grDevices::dev.off()
  }
)


###############################################################################
# Report output locations and dimensions
###############################################################################

message(
  dataset_title1,
  ": ",
  result1$number_of_rows,
  " peptides and ",
  result1$number_of_columns,
  " displayed columns."
)

message(
  dataset_title2,
  ": ",
  result2$number_of_rows,
  " peptides and ",
  result2$number_of_columns,
  " displayed columns."
)

message(
  "Heatmaps written to: ",
  normalizePath(
    output_pdf,
    mustWork = FALSE
  )
)

message(
  "df1 row order written to: ",
  normalizePath(
    output_row_order1,
    mustWork = FALSE
  )
)

message(
  "df2 row order written to: ",
  normalizePath(
    output_row_order2,
    mustWork = FALSE
  )
)

df4=read_protti("low_cells_normalized_HLAII.txt")
df3=read_protti("high_cells_normalized_HLAII.txt")



###############################################################################
# Differential peptide analysis
#
# df3 = high_cells
# df4 = low_cells
#
# Experimental condition examples:
#
# sts_1
# sts_2
# sts_3
#
# jq1_1
# jq1_2
# jq1_3
#
# n14_1
# n14_2
# n14_3
#
# Control examples:
#
# dms2_1
# dms2_2
#
# dms3_1
# dms3_2
#
# dms4_1
# dms4_2
#
###############################################################################


###############################################################################
# User settings
###############################################################################

peptide_column <- "t_t_t_eg_modified_sequence"


# Peptide must be observed in at least this fraction of condition
# replicates OR control replicates

presence_fraction <- 2 / 3


# Missing value imputation parameters

imputation_downshift <- 1.8

imputation_width <- 0.3


# Benjamini Hochberg adjusted significance threshold

alpha <- 0.05


# Seeds for reproducible imputation

seed_df3 <- 31415

seed_df4 <- 27182


###############################################################################
# Main analysis function
###############################################################################

analyse_peptide_screen <- function(
    data,
    dataset_name,
    peptide_column,
    presence_fraction = 2 / 3,
    imputation_downshift = 1.8,
    imputation_width = 0.3,
    alpha = 0.05,
    random_seed = 12345
) {
  
  
  ###########################################################################
  # Check input
  ###########################################################################
  
  if (!is.data.frame(data)) {
    
    stop(
      "[",
      dataset_name,
      "] Input must be a data frame.",
      call. = FALSE
    )
  }
  
  
  if (!peptide_column %in% names(data)) {
    
    stop(
      "[",
      dataset_name,
      "] Peptide column not found: ",
      peptide_column,
      call. = FALSE
    )
  }
  
  
  if (
    presence_fraction <= 0 ||
    presence_fraction > 1
  ) {
    
    stop(
      "[",
      dataset_name,
      "] presence_fraction must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  
  ###########################################################################
  # Identify intensity columns BY NAME
  ###########################################################################
  
  # IMPORTANT:
  #
  # We no longer require columns to already be numeric.
  #
  # This allows columns such as jq1_1 or n14_1 to be included even if
  # R imported them as character columns.
  #
  # Experimental conditions must contain exactly three alphanumeric
  # characters before the first underscore.
  #
  # Examples:
  #
  # sts_1
  # vbl_2
  # jq1_3
  # n14_1
  #
  # Controls may contain:
  #
  # dms_1
  # dms2_1
  # dms3_2
  # dms4_3
  
  
  experimental_name_match <- grepl(
    "^[A-Za-z0-9]{3}_[^_]+$",
    names(data)
  )
  
  
  control_name_match <- grepl(
    "^dms[0-9]*_[^_]+$",
    names(data),
    ignore.case = TRUE
  )
  
  
  intensity_columns <- which(
    experimental_name_match |
      control_name_match
  )
  
  
  if (length(intensity_columns) == 0L) {
    
    stop(
      "[",
      dataset_name,
      "] No intensity columns matching the expected naming pattern were found.",
      call. = FALSE
    )
  }
  
  
  intensity_column_names <- names(data)[
    intensity_columns
  ]
  
  
  ###########################################################################
  # Print every detected intensity column
  ###########################################################################
  
  message(
    "\n[",
    dataset_name,
    "] Detected intensity columns:"
  )
  
  
  message(
    paste(
      intensity_column_names,
      collapse = ", "
    )
  )
  
  
  ###########################################################################
  # Specifically report jq1 and n14 detection
  ###########################################################################
  
  jq1_columns <- grep(
    "^jq1_",
    intensity_column_names,
    value = TRUE,
    ignore.case = TRUE
  )
  
  
  n14_columns <- grep(
    "^n14_",
    intensity_column_names,
    value = TRUE,
    ignore.case = TRUE
  )
  
  
  message(
    "\n[",
    dataset_name,
    "] jq1 columns detected: ",
    ifelse(
      length(jq1_columns) > 0,
      paste(
        jq1_columns,
        collapse = ", "
      ),
      "NONE"
    )
  )
  
  
  message(
    "[",
    dataset_name,
    "] n14 columns detected: ",
    ifelse(
      length(n14_columns) > 0,
      paste(
        n14_columns,
        collapse = ", "
      ),
      "NONE"
    )
  )
  
  
  ###########################################################################
  # Convert intensity columns explicitly to numeric
  ###########################################################################
  
  # This handles columns imported as character or factor.
  #
  # Values such as:
  #
  # ""
  # "NA"
  # "NaN"
  #
  # will become NA.
  
  intensity_dataframe <- data[
    ,
    intensity_columns,
    drop = FALSE
  ]
  
  
  intensity_dataframe[] <- lapply(
    intensity_dataframe,
    function(x) {
      
      suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )
    }
  )
  
  
  ###########################################################################
  # Construct numeric intensity matrix
  ###########################################################################
  
  intensity_matrix <- as.matrix(
    intensity_dataframe
  )
  
  
  storage.mode(
    intensity_matrix
  ) <- "double"
  
  
  colnames(
    intensity_matrix
  ) <- intensity_column_names
  
  
  ###########################################################################
  # Standardize missing values
  ###########################################################################
  
  intensity_matrix[
    is.nan(intensity_matrix)
  ] <- NA_real_
  
  
  intensity_matrix[
    is.infinite(intensity_matrix)
  ] <- NA_real_
  
  
  ###########################################################################
  # Check that intensity columns contain usable numeric data
  ###########################################################################
  
  usable_values_per_column <- colSums(
    !is.na(
      intensity_matrix
    )
  )
  
  
  completely_empty_columns <- names(
    usable_values_per_column
  )[
    usable_values_per_column == 0
  ]
  
  
  if (length(completely_empty_columns) > 0L) {
    
    warning(
      "[",
      dataset_name,
      "] The following detected intensity columns contain no usable ",
      "numeric values after conversion: ",
      paste(
        completely_empty_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  ###########################################################################
  # Extract peptide names
  ###########################################################################
  
  peptide_names <- as.character(
    data[[peptide_column]]
  )
  
  
  missing_peptide_names <- (
    is.na(peptide_names) |
      !nzchar(peptide_names)
  )
  
  
  peptide_names[
    missing_peptide_names
  ] <- paste0(
    "unnamed_peptide_",
    which(
      missing_peptide_names
    )
  )
  
  
  rownames(
    intensity_matrix
  ) <- make.unique(
    peptide_names
  )
  
  
  ###########################################################################
  # Extract condition prefix
  ###########################################################################
  
  # Examples:
  #
  # sts_1  -> sts
  # jq1_1  -> jq1
  # n14_3  -> n14
  # dms2_1 -> dms2
  # dms3_1 -> dms3
  # dms4_1 -> dms4
  
  condition_prefix <- sub(
    "_.*$",
    "",
    intensity_column_names
  )
  
  
  condition_prefix <- tolower(
    condition_prefix
  )
  
  
  ###########################################################################
  # Identify DMS controls
  ###########################################################################
  
  control_mask <- grepl(
    "^dms[0-9]*$",
    condition_prefix,
    ignore.case = TRUE
  )
  
  
  control_columns <- which(
    control_mask
  )
  
  
  if (length(control_columns) == 0L) {
    
    stop(
      "[",
      dataset_name,
      "] No DMS control columns were found.",
      call. = FALSE
    )
  }
  
  
  number_control_replicates <- length(
    control_columns
  )
  
  
  control_column_names <- intensity_column_names[
    control_columns
  ]
  
  
  control_prefixes <- unique(
    condition_prefix[
      control_columns
    ]
  )
  
  
  message(
    "\n[",
    dataset_name,
    "] DMS groups detected: ",
    paste(
      control_prefixes,
      collapse = ", "
    )
  )
  
  
  message(
    "[",
    dataset_name,
    "] DMS columns: ",
    paste(
      control_column_names,
      collapse = ", "
    )
  )
  
  
  message(
    "[",
    dataset_name,
    "] Number of pooled DMS replicates: ",
    number_control_replicates
  )
  
  
  ###########################################################################
  # Identify experimental conditions
  ###########################################################################
  
  experimental_conditions <- unique(
    condition_prefix[
      !control_mask
    ]
  )
  
  
  if (length(experimental_conditions) == 0L) {
    
    stop(
      "[",
      dataset_name,
      "] No experimental conditions were found.",
      call. = FALSE
    )
  }
  
  
  message(
    "\n[",
    dataset_name,
    "] Experimental conditions detected:"
  )
  
  
  message(
    paste(
      experimental_conditions,
      collapse = ", "
    )
  )
  
  
  ###########################################################################
  # Verify jq1 and n14 at the condition level
  ###########################################################################
  
  if ("jq1" %in% experimental_conditions) {
    
    message(
      "[",
      dataset_name,
      "] jq1 successfully included."
    )
  }
  
  
  if ("n14" %in% experimental_conditions) {
    
    message(
      "[",
      dataset_name,
      "] n14 successfully included."
    )
    )


###########################################################################
# Preserve original peptide presence
###########################################################################

# Presence is determined BEFORE imputation.

peptide_present <- !is.na(
  intensity_matrix
)


###########################################################################
# Global distribution for fallback imputation
###########################################################################

all_observed_values <- intensity_matrix[
  !is.na(
    intensity_matrix
  )
]


if (length(all_observed_values) < 2L) {
  
  stop(
    "[",
    dataset_name,
    "] Too few observed intensity values for imputation.",
    call. = FALSE
  )
}


global_mean <- mean(
  all_observed_values
)


global_sd <- stats::sd(
  all_observed_values
)


if (
  !is.finite(global_sd) ||
  global_sd <= 0
) {
  
  stop(
    "[",
    dataset_name,
    "] Could not calculate a valid global standard deviation.",
    call. = FALSE
  )
}


###########################################################################
# Create imputed matrix
###########################################################################

imputed_matrix <- intensity_matrix


set.seed(
  random_seed
)


###########################################################################
# Impute each replicate independently
###########################################################################

for (
  column_index in seq_len(
    ncol(
      imputed_matrix
    )
  )
) {
  
  
  #########################################################################
  # Observed values in current replicate
  #########################################################################
  
  observed_values <- intensity_matrix[
    !is.na(
      intensity_matrix[
        ,
        column_index
      ]
    ),
    column_index
  ]
  
  
  #########################################################################
  # Calculate sample mean and SD
  #########################################################################
  
  if (length(observed_values) >= 2L) {
    
    sample_mean <- mean(
      observed_values
    )
    
    
    sample_sd <- stats::sd(
      observed_values
    )
    
  } else {
    
    sample_mean <- global_mean
    
    sample_sd <- global_sd
  }
  
  
  #########################################################################
  # Fallback if sample distribution cannot be estimated
  #########################################################################
  
  if (!is.finite(sample_mean)) {
    
    sample_mean <- global_mean
  }
  
  
  if (
    !is.finite(sample_sd) ||
    sample_sd <= 0
  ) {
    
    sample_sd <- global_sd
  }
  
  
  #########################################################################
  # Find missing values
  #########################################################################
  
  missing_positions <- which(
    is.na(
      intensity_matrix[
        ,
        column_index
      ]
    )
  )
  
  
  number_missing <- length(
    missing_positions
  )
  
  
  #########################################################################
  # Downshifted normal imputation
  #########################################################################
  
  if (number_missing > 0L) {
    
    
    imputation_mean <- (
      sample_mean -
        imputation_downshift *
        sample_sd
    )
    
    
    imputation_sd <- (
      imputation_width *
        sample_sd
    )
    
    
    imputed_matrix[
      missing_positions,
      column_index
    ] <- stats::rnorm(
      n = number_missing,
      mean = imputation_mean,
      sd = imputation_sd
    )
  }
}


###########################################################################
# Prepare summary result storage
###########################################################################

summary_results <- vector(
  mode = "list",
  length = length(
    experimental_conditions
  )
)


###########################################################################
# Test every experimental condition
###########################################################################

for (
  condition_index in seq_along(
    experimental_conditions
  )
) {
  
  
  current_condition <- experimental_conditions[
    condition_index
  ]
  
  
  message(
    "[",
    dataset_name,
    "] Testing ",
    current_condition,
    " versus pooled DMS"
  )
  
  
  #########################################################################
  # Identify condition replicates
  #########################################################################
  
  condition_columns <- which(
    condition_prefix ==
      current_condition
  )
  
  
  number_condition_replicates <- length(
    condition_columns
  )
  
  
  #########################################################################
  # Required peptide presence
  #########################################################################
  
  required_condition_presence <- ceiling(
    presence_fraction *
      number_condition_replicates
  )
  
  
  required_control_presence <- ceiling(
    presence_fraction *
      number_control_replicates
  )
  
  
  #########################################################################
  # Count observed condition replicates
  #########################################################################
  
  condition_presence_count <- rowSums(
    peptide_present[
      ,
      condition_columns,
      drop = FALSE
    ]
  )
  
  
  #########################################################################
  # Count observed DMS replicates
  #########################################################################
  
  control_presence_count <- rowSums(
    peptide_present[
      ,
      control_columns,
      drop = FALSE
    ]
  )
  
  
  #########################################################################
  # Apply requested OR filter
  #########################################################################
  
  passes_presence_filter <- (
    condition_presence_count >=
      required_condition_presence
  ) |
    (
      control_presence_count >=
        required_control_presence
    )
  
  
  retained_rows <- which(
    passes_presence_filter
  )
  
  
  number_retained <- length(
    retained_rows
  )
  
  
  #########################################################################
  # Require at least two replicates per group
  #########################################################################
  
  if (
    number_condition_replicates < 2L ||
    number_control_replicates < 2L
  ) {
    
    
    summary_results[[
      condition_index
    ]] <- data.frame(
      
      condition =
        current_condition,
      
      condition_replicates =
        number_condition_replicates,
      
      control_replicates =
        number_control_replicates,
      
      minimum_condition_presence =
        required_condition_presence,
      
      minimum_control_presence =
        required_control_presence,
      
      peptides_passing_filter =
        number_retained,
      
      peptides_tested =
        NA_integer_,
      
      significant_total =
        NA_integer_,
      
      significant_positive =
        NA_integer_,
      
      significant_negative =
        NA_integer_,
      
      status =
        "not_tested_less_than_2_replicates",
      
      stringsAsFactors = FALSE
    )
    
    
    next
  }
  
  
  #########################################################################
  # No peptides passed filter
  #########################################################################
  
  if (number_retained == 0L) {
    
    
    summary_results[[
      condition_index
    ]] <- data.frame(
      
      condition =
        current_condition,
      
      condition_replicates =
        number_condition_replicates,
      
      control_replicates =
        number_control_replicates,
      
      minimum_condition_presence =
        required_condition_presence,
      
      minimum_control_presence =
        required_control_presence,
      
      peptides_passing_filter =
        0L,
      
      peptides_tested =
        0L,
      
      significant_total =
        0L,
      
      significant_positive =
        0L,
      
      significant_negative =
        0L,
      
      status =
        "no_peptides_passed_filter",
      
      stringsAsFactors = FALSE
    )
    
    
    next
  }
  
  
  #########################################################################
  # Extract imputed values
  #########################################################################
  
  condition_values <- imputed_matrix[
    retained_rows,
    condition_columns,
    drop = FALSE
  ]
  
  
  control_values <- imputed_matrix[
    retained_rows,
    control_columns,
    drop = FALSE
  ]
  
  
  #########################################################################
  # Calculate group means
  #########################################################################
  
  condition_means <- rowMeans(
    condition_values
  )
  
  
  control_means <- rowMeans(
    control_values
  )
  
  
  #########################################################################
  # Calculate mean difference
  #########################################################################
  
  mean_difference <- (
    condition_means -
      control_means
  )
  
  
  #########################################################################
  # Group variances
  #########################################################################
  
  condition_variances <- apply(
    condition_values,
    1,
    stats::var
  )
  
  
  control_variances <- apply(
    control_values,
    1,
    stats::var
  )
  
  
  #########################################################################
  # Student t test degrees of freedom
  #########################################################################
  
  degrees_freedom <- (
    number_condition_replicates +
      number_control_replicates -
      2
  )
  
  
  #########################################################################
  # Pooled variance
  #########################################################################
  
  pooled_variance <- (
    
    (
      (
        number_condition_replicates -
          1
      ) *
        condition_variances
    ) +
      
      (
        (
          number_control_replicates -
            1
        ) *
          control_variances
      )
    
  ) / degrees_freedom
  
  
  #########################################################################
  # Standard error
  #########################################################################
  
  standard_error <- sqrt(
    
    pooled_variance *
      
      (
        1 /
          number_condition_replicates +
          
          1 /
          number_control_replicates
      )
  )
  
  
  #########################################################################
  # Student t statistic
  #########################################################################
  
  t_statistic <- (
    mean_difference /
      standard_error
  )
  
  
  #########################################################################
  # Two sided p values
  #########################################################################
  
  raw_p_values <- rep(
    NA_real_,
    number_retained
  )
  
  
  valid_test <- (
    is.finite(
      t_statistic
    ) &
      is.finite(
        standard_error
      ) &
      standard_error > 0
  )
  
  
  raw_p_values[
    valid_test
  ] <- 2 * stats::pt(
    q = abs(
      t_statistic[
        valid_test
      ]
    ),
    df = degrees_freedom,
    lower.tail = FALSE
  )
  
  
  #########################################################################
  # Handle zero variance cases
  #########################################################################
  
  zero_standard_error <- (
    is.finite(
      standard_error
    ) &
      standard_error == 0
  )
  
  
  raw_p_values[
    zero_standard_error &
      mean_difference == 0
  ] <- 1
  
  
  raw_p_values[
    zero_standard_error &
      mean_difference != 0
  ] <- 0
  
  
  #########################################################################
  # Benjamini Hochberg correction
  #########################################################################
  
  adjusted_p_values <- stats::p.adjust(
    raw_p_values,
    method = "BH"
  )
  
  
  #########################################################################
  # Significant peptides
  #########################################################################
  
  significant <- (
    !is.na(
      adjusted_p_values
    ) &
      adjusted_p_values < alpha
  )
  
  
  #########################################################################
  # Positive enrichment
  #########################################################################
  
  significant_positive <- (
    significant &
      mean_difference > 0
  )
  
  
  #########################################################################
  # Negative enrichment
  #########################################################################
  
  significant_negative <- (
    significant &
      mean_difference < 0
  )
  
  
  #########################################################################
  # Store results
  #########################################################################
  
  summary_results[[
    condition_index
  ]] <- data.frame(
    
    condition =
      current_condition,
    
    condition_replicates =
      number_condition_replicates,
    
    control_replicates =
      number_control_replicates,
    
    minimum_condition_presence =
      required_condition_presence,
    
    minimum_control_presence =
      required_control_presence,
    
    peptides_passing_filter =
      number_retained,
    
    peptides_tested =
      sum(
        !is.na(
          raw_p_values
        )
      ),
    
    significant_total =
      sum(
        significant
      ),
    
    significant_positive =
      sum(
        significant_positive
      ),
    
    significant_negative =
      sum(
        significant_negative
      ),
    
    status =
      "tested",
    
    stringsAsFactors = FALSE
  )
}


###########################################################################
# Combine results
###########################################################################

summary_dataframe <- do.call(
  rbind,
  summary_results
)


rownames(
  summary_dataframe
) <- NULL


###########################################################################
# Return results
###########################################################################

summary_dataframe
  }
  
  
  ###############################################################################
  # Run high cell analysis
  #
  # df3
  #
  # dms3 and dms4 are pooled as controls
  ###############################################################################
  
  high_cells_stats <- analyse_peptide_screen(
    
    data = df3,
    
    dataset_name = "high_cells",
    
    peptide_column = peptide_column,
    
    presence_fraction = presence_fraction,
    
    imputation_downshift = imputation_downshift,
    
    imputation_width = imputation_width,
    
    alpha = alpha,
    
    random_seed = seed_df3
  )
  
  
  ###############################################################################
  # Run low cell analysis
  #
  # df4
  #
  # dms2 samples are used as controls
  ###############################################################################
  
  low_cells_stats <- analyse_peptide_screen(
    
    data = df4,
    
    dataset_name = "low_cells",
    
    peptide_column = peptide_column,
    
    presence_fraction = presence_fraction,
    
    imputation_downshift = imputation_downshift,
    
    imputation_width = imputation_width,
    
    alpha = alpha,
    
    random_seed = seed_df4
  )
  
  
  ###############################################################################
  # Display results
  ###############################################################################
  
  print(
    high_cells_stats
  )
  
  
  print(
    low_cells_stats
  )
  
  
  ###############################################################################
  # Specifically verify jq1 and n14
  ###############################################################################
  
  high_cells_stats[
    high_cells_stats$condition %in%
      c(
        "jq1",
        "n14"
      ),
  ]
  
  
  low_cells_stats[
    low_cells_stats$condition %in%
      c(
        "jq1",
        "n14"
      ),
  ]
  
  
  ###############################################################################
  # Optional RStudio inspection
  ###############################################################################
  
  View(
    high_cells_stats
  )
  
  
  View(
    low_cells_stats
  )
  
  #############################################################################################################################################################################################################################################
  # Combine Heatmap and statistics
  #############################################################################################################################################################################################################################################
  
  ###############################################################################
  # TWO PEPTIDE HEATMAPS WITH HIERARCHICAL CLUSTERING
  # AND DIVERGING PERCENT SIGNIFICANT PEPTIDE BARPLOTS
  #
  # df1 corresponds to high_cells_stats
  # df2 corresponds to low_cells_stats
  ###############################################################################
  
  
  ###############################################################################
  # USER SETTINGS
  ###############################################################################
  
  N1 <- 15L
  N2 <- 81L
  
  peptide_column <- "t_t_t_eg_modified_sequence"
  
  dataset_title1 <- "High input"
  dataset_title2 <- "Low input"
  
  collapse_replicates1 <- FALSE
  collapse_replicates2 <- FALSE
  
  clustering_na_fill <- "zero"
  
  z_limit <- 3
  
  missing_value_colour <- "#BDBDBD"
  
  show_peptide_names <- FALSE
  
  show_row_dendrogram <- TRUE
  
  row_dendrogram_width_mm <- 24
  
  row_dendrogram_line_width <- 0.35
  
  panel_width_mode <- "proportional"
  
  
  ###############################################################################
  # BARPLOT SETTINGS
  ###############################################################################
  
  # Publication suitable orange
  upregulated_colour <- "#E69F00"
  
  # Turquoise
  downregulated_colour <- "#00A6A6"
  
  # Height of the upper and lower portions of the diverging barplot
  barplot_half_height_mm <- 11
  
  # NULL automatically calculates a common maximum for both heatmaps
  #
  # You can alternatively specify a fixed value such as:
  #
  # barplot_limit_percent <- 20
  
  barplot_limit_percent <- NULL
  
  
  ###############################################################################
  # DEFINE WHICH STATISTICS BELONG TO EACH HEATMAP
  ###############################################################################
  
  stats_for_df1 <- low_cells_stats
  
  stats_for_df2 <- high_cells_stats
  
  
  ###############################################################################
  # OUTPUT SETTINGS
  ###############################################################################
  
  output_pdf <- "two_peptide_heatmaps_with_regulation.pdf"
  
  output_row_order1 <- "df1_peptide_row_order.csv"
  
  output_row_order2 <- "df2_peptide_row_order.csv"
  
  pdf_width <- 18
  
  # Slightly taller than before to accommodate barplots
  pdf_height <- 12
  
  legend_width_mm <- 32
  
  layout_gap_mm <- 5
  
  
  ###############################################################################
  # CHECK REQUIRED PACKAGES
  ###############################################################################
  
  required_packages <- c(
    "ComplexHeatmap",
    "circlize",
    "scico",
    "fastcluster"
  )
  
  installed_packages <- vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  
  if (any(!installed_packages)) {
    
    stop(
      "The following packages are missing: ",
      paste(
        required_packages[!installed_packages],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  ###############################################################################
  # CHECK INPUT DATA
  ###############################################################################
  
  if (!exists("df1")) {
    
    stop(
      "The object 'df1' was not found.",
      call. = FALSE
    )
  }
  
  if (!exists("df2")) {
    
    stop(
      "The object 'df2' was not found.",
      call. = FALSE
    )
  }
  
  if (!exists("high_cells_stats")) {
    
    stop(
      "The object 'high_cells_stats' was not found.",
      call. = FALSE
    )
  }
  
  if (!exists("low_cells_stats")) {
    
    stop(
      "The object 'low_cells_stats' was not found.",
      call. = FALSE
    )
  }
  
  
  ###############################################################################
  # CREATE SHARED HEATMAP COLOUR FUNCTION
  ###############################################################################
  
  heatmap_palette <- scico::scico(
    n = 101,
    palette = "vik"
  )
  
  colour_function <- circlize::colorRamp2(
    seq(
      from = -z_limit,
      to = z_limit,
      length.out = length(heatmap_palette)
    ),
    heatmap_palette
  )
  
  
  ###############################################################################
  # CALCULATE PERCENTAGES FROM T TEST SUMMARY TABLE
  ###############################################################################
  
  calculate_percentages_from_stats <- function(
    stats_dataframe
  ) {
    
    denominator <- stats_dataframe$peptides_tested
    
    percent_positive <- rep(
      NA_real_,
      nrow(stats_dataframe)
    )
    
    percent_negative <- rep(
      NA_real_,
      nrow(stats_dataframe)
    )
    
    usable <- (
      !is.na(denominator) &
        denominator > 0
    )
    
    percent_positive[
      usable
    ] <- 100 *
      stats_dataframe$significant_positive[usable] /
      denominator[usable]
    
    percent_negative[
      usable
    ] <- 100 *
      stats_dataframe$significant_negative[usable] /
      denominator[usable]
    
    data.frame(
      
      condition = tolower(
        as.character(
          stats_dataframe$condition
        )
      ),
      
      percent_positive = percent_positive,
      
      percent_negative = percent_negative,
      
      stringsAsFactors = FALSE
    )
  }
  
  
  ###############################################################################
  # CALCULATE PERCENTAGES FOR BOTH DATASETS
  ###############################################################################
  
  percentage_stats1 <- calculate_percentages_from_stats(
    stats_for_df1
  )
  
  percentage_stats2 <- calculate_percentages_from_stats(
    stats_for_df2
  )
  
  
  ###############################################################################
  # DETERMINE COMMON BARPLOT SCALE
  ###############################################################################
  
  if (is.null(barplot_limit_percent)) {
    
    all_percentages <- c(
      
      percentage_stats1$percent_positive,
      percentage_stats1$percent_negative,
      
      percentage_stats2$percent_positive,
      percentage_stats2$percent_negative
    )
    
    all_percentages <- all_percentages[
      is.finite(all_percentages)
    ]
    
    if (length(all_percentages) == 0L) {
      
      calculated_limit <- 5
      
    } else {
      
      calculated_limit <- max(
        all_percentages
      )
      
      calculated_limit <- ceiling(
        calculated_limit / 5
      ) * 5
      
      calculated_limit <- max(
        calculated_limit,
        5
      )
    }
    
    barplot_limit_percent <- calculated_limit
  }
  
  
  message(
    "Barplot scale: +/- ",
    barplot_limit_percent,
    "%"
  )
  
  
  ###############################################################################
  # MAP CONDITION NAMES TO PERCENTAGES
  ###############################################################################
  
  get_condition_percentages <- function(
    conditions,
    percentage_dataframe,
    dataset_label
  ) {
    
    conditions <- tolower(
      as.character(
        conditions
      )
    )
    
    percent_positive <- numeric(
      length(conditions)
    )
    
    percent_negative <- numeric(
      length(conditions)
    )
    
    
    for (
      condition_index in seq_along(
        conditions
      )
    ) {
      
      current_condition <- conditions[
        condition_index
      ]
      
      
      #########################################################################
      # EVERY DMS CONTROL IS ZERO
      #########################################################################
      
      if (
        grepl(
          "^dms",
          current_condition
        )
      ) {
        
        percent_positive[
          condition_index
        ] <- 0
        
        percent_negative[
          condition_index
        ] <- 0
        
        next
      }
      
      
      #########################################################################
      # FIND CONDITION IN STATISTICS TABLE
      #########################################################################
      
      stats_position <- match(
        current_condition,
        percentage_dataframe$condition
      )
      
      
      if (is.na(stats_position)) {
        
        warning(
          "[",
          dataset_label,
          "] No statistics found for condition ",
          current_condition,
          ". Barplot value set to zero.",
          call. = FALSE
        )
        
        percent_positive[
          condition_index
        ] <- 0
        
        percent_negative[
          condition_index
        ] <- 0
        
        next
      }
      
      
      #########################################################################
      # RETRIEVE PERCENTAGES
      #########################################################################
      
      positive_value <- percentage_dataframe$percent_positive[
        stats_position
      ]
      
      negative_value <- percentage_dataframe$percent_negative[
        stats_position
      ]
      
      
      if (!is.finite(positive_value)) {
        
        positive_value <- 0
      }
      
      if (!is.finite(negative_value)) {
        
        negative_value <- 0
      }
      
      
      percent_positive[
        condition_index
      ] <- positive_value
      
      percent_negative[
        condition_index
      ] <- negative_value
    }
    
    
    data.frame(
      
      condition = conditions,
      
      percent_positive = percent_positive,
      
      percent_negative = percent_negative,
      
      stringsAsFactors = FALSE
    )
  }
  
  
  ###############################################################################
  # FUNCTION TO PREPARE ONE HEATMAP
  ###############################################################################
  
  prepare_peptide_heatmap <- function(
    data,
    N,
    dataset_label,
    heatmap_name,
    peptide_column,
    stats_dataframe,
    collapse_replicates,
    clustering_na_fill,
    colour_function,
    missing_value_colour,
    show_peptide_names,
    show_row_dendrogram,
    row_dendrogram_width_mm,
    row_dendrogram_line_width,
    barplot_limit_percent,
    barplot_half_height_mm,
    upregulated_colour,
    downregulated_colour
  ) {
    
    
    ###########################################################################
    # INTERNAL ERROR FUNCTION
    ###########################################################################
    
    dataset_error <- function(...) {
      
      stop(
        "[",
        dataset_label,
        "] ",
        ...,
        call. = FALSE
      )
    }
    
    
    ###########################################################################
    # IDENTIFY NUMERIC INTENSITY COLUMNS
    ###########################################################################
    
    numeric_column_index <- which(
      vapply(
        data,
        is.numeric,
        logical(1)
      )
    )
    
    
    if (length(numeric_column_index) < N) {
      
      dataset_error(
        "Only ",
        length(numeric_column_index),
        " numeric columns were found, but N = ",
        N,
        "."
      )
    }
    
    
    intensity_columns <- numeric_column_index[
      seq_len(N)
    ]
    
    
    message(
      "[",
      dataset_label,
      "] Heatmap intensity columns:"
    )
    
    message(
      paste(
        names(data)[intensity_columns],
        collapse = ", "
      )
    )
    
    
    ###########################################################################
    # CONSTRUCT INTENSITY MATRIX
    ###########################################################################
    
    intensity_matrix <- as.matrix(
      data[
        ,
        intensity_columns,
        drop = FALSE
      ]
    )
    
    storage.mode(
      intensity_matrix
    ) <- "double"
    
    
    intensity_matrix[
      is.nan(intensity_matrix)
    ] <- NA_real_
    
    
    intensity_matrix[
      is.infinite(intensity_matrix)
    ] <- NA_real_
    
    
    sample_names <- colnames(
      intensity_matrix
    )
    
    
    ###########################################################################
    # PEPTIDE LABELS
    ###########################################################################
    
    peptide_labels <- as.character(
      data[[peptide_column]]
    )
    
    
    missing_peptide_labels <- (
      is.na(peptide_labels) |
        !nzchar(peptide_labels)
    )
    
    
    peptide_labels[
      missing_peptide_labels
    ] <- paste0(
      "unnamed_peptide_",
      which(
        missing_peptide_labels
      )
    )
    
    
    rownames(
      intensity_matrix
    ) <- sprintf(
      "%s_peptide_row_%06d",
      heatmap_name,
      seq_len(
        nrow(
          intensity_matrix
        )
      )
    )
    
    
    ###########################################################################
    # PARSE CONDITION AND REPLICATE
    ###########################################################################
    
    condition <- sub(
      "_.*$",
      "",
      sample_names
    )
    
    replicate <- sub(
      "^[^_]*_",
      "",
      sample_names
    )
    
    
    ###########################################################################
    # CONDITION ORDER
    ###########################################################################
    
    condition_levels <- unique(
      condition
    )
    
    
    condition_factor <- factor(
      condition,
      levels = condition_levels
    )
    
    
    column_order <- order(
      as.integer(
        condition_factor
      ),
      seq_along(
        condition_factor
      )
    )
    
    
    intensity_matrix <- intensity_matrix[
      ,
      column_order,
      drop = FALSE
    ]
    
    
    condition_factor <- factor(
      condition[column_order],
      levels = condition_levels
    )
    
    
    replicate <- replicate[
      column_order
    ]
    
    
    sample_names <- sample_names[
      column_order
    ]
    
    
    ###########################################################################
    # COLLAPSE REPLICATES IF REQUESTED
    ###########################################################################
    
    if (collapse_replicates) {
      
      
      condition_character <- as.character(
        condition_factor
      )
      
      
      heatmap_matrix <- vapply(
        condition_levels,
        FUN = function(
    current_condition
        ) {
          
          current_columns <- (
            condition_character ==
              current_condition
          )
          
          
          condition_mean <- rowMeans(
            intensity_matrix[
              ,
              current_columns,
              drop = FALSE
            ],
            na.rm = TRUE
          )
          
          
          condition_mean[
            is.nan(
              condition_mean
            )
          ] <- NA_real_
          
          
          condition_mean
        },
    FUN.VALUE = numeric(
      nrow(
        intensity_matrix
      )
    )
      )
      
      
      heatmap_matrix <- as.matrix(
        heatmap_matrix
      )
      
      
      rownames(
        heatmap_matrix
      ) <- rownames(
        intensity_matrix
      )
      
      
      colnames(
        heatmap_matrix
      ) <- condition_levels
      
      
      column_split_value <- NULL
      
      
      displayed_column_labels <- condition_levels
      
      
      annotation_conditions <- condition_levels
      
      
    } else {
      
      
      heatmap_matrix <- intensity_matrix
      
      
      column_split_value <- condition_factor
      
      
      displayed_column_labels <- replicate
      
      
      annotation_conditions <- as.character(
        condition_factor
      )
    }
    
    
    ###########################################################################
    # GET PERCENT UP AND DOWN FOR DISPLAYED CONDITIONS
    ###########################################################################
    
    condition_percentages <- get_condition_percentages(
      
      conditions =
        annotation_conditions,
      
      percentage_dataframe =
        stats_dataframe,
      
      dataset_label =
        dataset_label
    )
    
    
    positive_percentages <- condition_percentages$percent_positive
    
    
    # Negative values are deliberately plotted below zero
    
    negative_percentages <- (
      -1 *
        condition_percentages$percent_negative
    )
    
    
    ###########################################################################
    # CREATE BOTTOM DIVERGING BARPLOT
    ###########################################################################
    
    # The two annotation tracks have zero at their common boundary.
    #
    # The first extends upward from zero.
    #
    # The second extends downward from zero.
    #
    # With zero gap they visually form one diverging barplot.
    
    
    upper_axis_points <- unique(
      round(
        c(
          0,
          barplot_limit_percent / 2,
          barplot_limit_percent
        ),
        1
      )
    )
    
    
    lower_axis_points <- unique(
      round(
        c(
          -barplot_limit_percent,
          -barplot_limit_percent / 2,
          0
        ),
        1
      )
    )
    
    
    bottom_bar_annotation <- ComplexHeatmap::HeatmapAnnotation(
      
      
      #########################################################################
      # UPREGULATED PEPTIDES
      #########################################################################
      
      Upregulated = ComplexHeatmap::anno_barplot(
        
        positive_percentages,
        
        baseline = 0,
        
        ylim = c(
          0,
          barplot_limit_percent
        ),
        
        bar_width = 1,
        
        border = FALSE,
        
        gp = grid::gpar(
          fill = upregulated_colour,
          col = NA
        ),
        
        height = grid::unit(
          barplot_half_height_mm,
          "mm"
        ),
        
        axis = TRUE,
        
        axis_param = list(
          
          side = "left",
          
          at = upper_axis_points,
          
          labels = paste0(
            upper_axis_points,
            "%"
          ),
          
          gp = grid::gpar(
            fontsize = 6
          )
        )
      ),
      
      
      #########################################################################
      # DOWNREGULATED PEPTIDES
      #########################################################################
      
      Downregulated = ComplexHeatmap::anno_barplot(
        
        negative_percentages,
        
        baseline = 0,
        
        ylim = c(
          -barplot_limit_percent,
          0
        ),
        
        bar_width = 1,
        
        border = FALSE,
        
        gp = grid::gpar(
          fill = downregulated_colour,
          col = NA
        ),
        
        height = grid::unit(
          barplot_half_height_mm,
          "mm"
        ),
        
        axis = TRUE,
        
        axis_param = list(
          
          side = "left",
          
          at = lower_axis_points,
          
          labels = paste0(
            lower_axis_points,
            "%"
          ),
          
          gp = grid::gpar(
            fontsize = 6
          )
        )
      ),
      
      
      #########################################################################
      # MAKE BOTH TRACKS TOUCH
      #########################################################################
      
      gap = grid::unit(
        0,
        "mm"
      ),
      
      show_annotation_name = FALSE
    )
    
    
    ###########################################################################
    # MATRIX USED ONLY FOR CLUSTERING
    ###########################################################################
    
    clustering_matrix <- heatmap_matrix
    
    
    if (
      clustering_na_fill ==
      "zero"
    ) {
      
      
      clustering_matrix[
        is.na(
          clustering_matrix
        )
      ] <- 0
      
      
    } else if (
      clustering_na_fill ==
      "row_mean"
    ) {
      
      
      peptide_means <- rowMeans(
        clustering_matrix,
        na.rm = TRUE
      )
      
      
      peptide_means[
        !is.finite(
          peptide_means
        )
      ] <- 0
      
      
      for (
        column_index in seq_len(
          ncol(
            clustering_matrix
          )
        )
      ) {
        
        
        missing_rows <- is.na(
          clustering_matrix[
            ,
            column_index
          ]
        )
        
        
        clustering_matrix[
          missing_rows,
          column_index
        ] <- peptide_means[
          missing_rows
        ]
      }
      
      
    } else {
      
      
      dataset_error(
        "'clustering_na_fill' must be 'zero' or 'row_mean'."
      )
    }
    
    
    ###########################################################################
    # HIERARCHICAL CLUSTERING
    ###########################################################################
    
    row_hclust <- fastcluster::hclust.vector(
      
      clustering_matrix,
      
      method = "ward",
      
      metric = "euclidean"
    )
    
    
    ###########################################################################
    # BUILD HEATMAP
    ###########################################################################
    
    heatmap_object <- ComplexHeatmap::Heatmap(
      
      heatmap_matrix,
      
      name = heatmap_name,
      
      col = colour_function,
      
      na_col = missing_value_colour,
      
      show_heatmap_legend = FALSE,
      
      
      #########################################################################
      # ROW CLUSTERING
      #########################################################################
      
      cluster_rows = row_hclust,
      
      row_dend_reorder = FALSE,
      
      show_row_dend = show_row_dendrogram,
      
      row_dend_side = "left",
      
      row_dend_width = grid::unit(
        row_dendrogram_width_mm,
        "mm"
      ),
      
      row_dend_gp = grid::gpar(
        col = "black",
        lwd = row_dendrogram_line_width
      ),
      
      
      #########################################################################
      # COLUMN GROUPING
      #########################################################################
      
      cluster_columns = FALSE,
      
      show_column_dend = FALSE,
      
      column_split = column_split_value,
      
      cluster_column_slices = FALSE,
      
      column_gap = grid::unit(
        2,
        "mm"
      ),
      
      
      #########################################################################
      # COLUMN LABELS
      #########################################################################
      
      column_labels = displayed_column_labels,
      
      show_column_names = TRUE,
      
      column_names_side = "bottom",
      
      column_names_rot = 45,
      
      column_names_centered = TRUE,
      
      column_names_gp = grid::gpar(
        fontsize = 7
      ),
      
      column_names_max_height = grid::unit(
        30,
        "mm"
      ),
      
      
      #########################################################################
      # CONDITION HEADINGS
      #########################################################################
      
      column_title_gp = grid::gpar(
        fontsize = 9,
        fontface = "bold",
        fill = "#F2F2F2",
        col = "black"
      ),
      
      
      #########################################################################
      # BOTTOM BARPLOT
      #########################################################################
      
      bottom_annotation =
        bottom_bar_annotation,
      
      
      #########################################################################
      # ROW LABELS
      #########################################################################
      
      row_labels = peptide_labels,
      
      show_row_names = show_peptide_names,
      
      row_names_gp = grid::gpar(
        fontsize = 4
      ),
      
      
      #########################################################################
      # RASTERIZATION
      #########################################################################
      
      rect_gp = grid::gpar(
        col = NA
      ),
      
      use_raster = TRUE,
      
      raster_device = "png",
      
      raster_quality = 4,
      
      border = FALSE
    )
    
    
    ###########################################################################
    # CLUSTERED ROW ORDER
    ###########################################################################
    
    clustered_row_order <- row_hclust$order
    
    
    row_order_table <- data.frame(
      
      heatmap_position =
        seq_along(
          clustered_row_order
        ),
      
      original_df_row =
        clustered_row_order,
      
      peptide =
        peptide_labels[
          clustered_row_order
        ],
      
      stringsAsFactors = FALSE
    )
    
    
    ###########################################################################
    # RETURN RESULTS
    ###########################################################################
    
    list(
      
      heatmap =
        heatmap_object,
      
      row_hclust =
        row_hclust,
      
      row_order_table =
        row_order_table,
      
      number_of_rows =
        nrow(
          heatmap_matrix
        ),
      
      number_of_columns =
        ncol(
          heatmap_matrix
        ),
      
      condition_percentages =
        condition_percentages
    )
  }
  
  
  ###############################################################################
  # PREPARE HEATMAP 1
  ###############################################################################
  
  result1 <- prepare_peptide_heatmap(
    
    data = df1,
    
    N = N1,
    
    dataset_label = dataset_title1,
    
    heatmap_name = "df1_z_score",
    
    peptide_column = peptide_column,
    
    stats_dataframe = percentage_stats1,
    
    collapse_replicates = collapse_replicates1,
    
    clustering_na_fill = clustering_na_fill,
    
    colour_function = colour_function,
    
    missing_value_colour = missing_value_colour,
    
    show_peptide_names = show_peptide_names,
    
    show_row_dendrogram = show_row_dendrogram,
    
    row_dendrogram_width_mm = row_dendrogram_width_mm,
    
    row_dendrogram_line_width =
      row_dendrogram_line_width,
    
    barplot_limit_percent =
      barplot_limit_percent,
    
    barplot_half_height_mm =
      barplot_half_height_mm,
    
    upregulated_colour =
      upregulated_colour,
    
    downregulated_colour =
      downregulated_colour
  )
  
  
  ###############################################################################
  # PREPARE HEATMAP 2
  ###############################################################################
  
  result2 <- prepare_peptide_heatmap(
    
    data = df2,
    
    N = N2,
    
    dataset_label = dataset_title2,
    
    heatmap_name = "df2_z_score",
    
    peptide_column = peptide_column,
    
    stats_dataframe = percentage_stats2,
    
    collapse_replicates = collapse_replicates2,
    
    clustering_na_fill = clustering_na_fill,
    
    colour_function = colour_function,
    
    missing_value_colour = missing_value_colour,
    
    show_peptide_names = show_peptide_names,
    
    show_row_dendrogram = show_row_dendrogram,
    
    row_dendrogram_width_mm = row_dendrogram_width_mm,
    
    row_dendrogram_line_width =
      row_dendrogram_line_width,
    
    barplot_limit_percent =
      barplot_limit_percent,
    
    barplot_half_height_mm =
      barplot_half_height_mm,
    
    upregulated_colour =
      upregulated_colour,
    
    downregulated_colour =
      downregulated_colour
  )
  
  
  ###############################################################################
  # SAVE ROW ORDERS
  ###############################################################################
  
  utils::write.csv(
    
    result1$row_order_table,
    
    file = output_row_order1,
    
    row.names = FALSE
  )
  
  
  utils::write.csv(
    
    result2$row_order_table,
    
    file = output_row_order2,
    
    row.names = FALSE
  )
  
  
  ###############################################################################
  # HEATMAP Z SCORE LEGEND
  ###############################################################################
  
  z_score_legend <- ComplexHeatmap::Legend(
    
    title = "Peptide\nZ score",
    
    col_fun = colour_function,
    
    at = c(
      -z_limit,
      0,
      z_limit
    ),
    
    labels = c(
      paste0(
        "\u2264 ",
        -z_limit
      ),
      "0",
      paste0(
        "\u2265 ",
        z_limit
      )
    ),
    
    legend_height = grid::unit(
      45,
      "mm"
    ),
    
    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    
    labels_gp = grid::gpar(
      fontsize = 8
    )
  )
  
  
  ###############################################################################
  # MISSING VALUE LEGEND
  ###############################################################################
  
  missing_value_legend <- ComplexHeatmap::Legend(
    
    labels = "Missing",
    
    legend_gp = grid::gpar(
      fill = missing_value_colour,
      col = NA
    ),
    
    labels_gp = grid::gpar(
      fontsize = 8
    )
  )
  
  
  ###############################################################################
  # REGULATION LEGEND
  ###############################################################################
  
  regulation_legend <- ComplexHeatmap::Legend(
    
    title = "Significant\npeptides",
    
    labels = c(
      "Upregulated",
      "Downregulated"
    ),
    
    legend_gp = grid::gpar(
      
      fill = c(
        upregulated_colour,
        downregulated_colour
      ),
      
      col = NA
    ),
    
    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    
    labels_gp = grid::gpar(
      fontsize = 8
    )
  )
  
  
  ###############################################################################
  # COMBINE LEGENDS
  ###############################################################################
  
  shared_legend <- ComplexHeatmap::packLegend(
    
    z_score_legend,
    
    missing_value_legend,
    
    regulation_legend,
    
    direction = "vertical",
    
    gap = grid::unit(
      6,
      "mm"
    )
  )
  
  
  ###############################################################################
  # DETERMINE RELATIVE PANEL WIDTHS
  ###############################################################################
  
  if (
    panel_width_mode ==
    "proportional"
  ) {
    
    panel_widths <- c(
      result1$number_of_columns,
      result2$number_of_columns
    )
    
  } else {
    
    panel_widths <- c(
      1,
      1
    )
  }
  
  
  ###############################################################################
  # CREATE PDF
  ###############################################################################
  
  grDevices::pdf(
    
    file = output_pdf,
    
    width = pdf_width,
    
    height = pdf_height,
    
    useDingbats = FALSE
  )
  
  
  tryCatch(
    
    {
      
      
      grid::grid.newpage()
      
      
      comparison_layout <- grid::grid.layout(
        
        nrow = 2,
        
        ncol = 5,
        
        heights = grid::unit.c(
          
          grid::unit(
            8,
            "mm"
          ),
          
          grid::unit(
            1,
            "null"
          )
        ),
        
        widths = grid::unit.c(
          
          grid::unit(
            legend_width_mm,
            "mm"
          ),
          
          grid::unit(
            layout_gap_mm,
            "mm"
          ),
          
          grid::unit(
            panel_widths[1],
            "null"
          ),
          
          grid::unit(
            layout_gap_mm,
            "mm"
          ),
          
          grid::unit(
            panel_widths[2],
            "null"
          )
        )
      )
      
      
      grid::pushViewport(
        
        grid::viewport(
          layout = comparison_layout
        )
      )
      
      
      #########################################################################
      # DATASET 1 TITLE
      #########################################################################
      
      grid::pushViewport(
        
        grid::viewport(
          layout.pos.row = 1,
          layout.pos.col = 3
        )
      )
      
      
      grid::grid.text(
        
        dataset_title1,
        
        gp = grid::gpar(
          fontsize = 12,
          fontface = "bold"
        )
      )
      
      
      grid::popViewport()
      
      
      #########################################################################
      # DATASET 2 TITLE
      #########################################################################
      
      grid::pushViewport(
        
        grid::viewport(
          layout.pos.row = 1,
          layout.pos.col = 5
        )
      )
      
      
      grid::grid.text(
        
        dataset_title2,
        
        gp = grid::gpar(
          fontsize = 12,
          fontface = "bold"
        )
      )
      
      
      grid::popViewport()
      
      
      #########################################################################
      # SHARED LEGEND
      #########################################################################
      
      grid::pushViewport(
        
        grid::viewport(
          layout.pos.row = 2,
          layout.pos.col = 1
        )
      )
      
      
      ComplexHeatmap::draw(
        
        shared_legend,
        
        x = grid::unit(
          0.5,
          "npc"
        ),
        
        y = grid::unit(
          0.5,
          "npc"
        ),
        
        just = c(
          "center",
          "center"
        )
      )
      
      
      grid::popViewport()
      
      
      #########################################################################
      # DATASET 1
      #########################################################################
      
      grid::pushViewport(
        
        grid::viewport(
          
          layout.pos.row = 2,
          
          layout.pos.col = 3,
          
          clip = "off"
        )
      )
      
      
      ComplexHeatmap::draw(
        
        result1$heatmap,
        
        newpage = FALSE,
        
        show_heatmap_legend = FALSE,
        
        show_annotation_legend = FALSE,
        
        padding = grid::unit(
          c(
            2,
            2,
            2,
            2
          ),
          "mm"
        )
      )
      
      
      grid::popViewport()
      
      
      #########################################################################
      # DATASET 2
      #########################################################################
      
      grid::pushViewport(
        
        grid::viewport(
          
          layout.pos.row = 2,
          
          layout.pos.col = 5,
          
          clip = "off"
        )
      )
      
      
      ComplexHeatmap::draw(
        
        result2$heatmap,
        
        newpage = FALSE,
        
        show_heatmap_legend = FALSE,
        
        show_annotation_legend = FALSE,
        
        padding = grid::unit(
          c(
            2,
            2,
            2,
            2
          ),
          "mm"
        )
      )
      
      
      grid::popViewport()
      
      
      grid::popViewport()
    },
    
    finally = {
      
      grDevices::dev.off()
    }
  )
  
  
  ###############################################################################
  # INSPECT THE PERCENTAGES USED IN THE FIGURE
  ###############################################################################
  
  cat(
    "\nPercentages used for ",
    dataset_title1,
    ":\n",
    sep = ""
  )
  
  print(
    result1$condition_percentages
  )
  
  
  cat(
    "\nPercentages used for ",
    dataset_title2,
    ":\n",
    sep = ""
  )
  
  print(
    result2$condition_percentages
  )
  
  
  ###############################################################################
  # OUTPUT
  ###############################################################################
  
  message(
    "Heatmap written to: ",
    normalizePath(
      output_pdf,
      mustWork = FALSE
    )
  )
  
  