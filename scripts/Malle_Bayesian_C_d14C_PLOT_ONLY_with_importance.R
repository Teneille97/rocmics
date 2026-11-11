
# =============================================================================
# MALLE ERW: BAYESIAN C AND d14C PLOT-ONLY SCRIPT
#
# Purpose:
#   1. Load existing Bayesian prediction CSV files.
#   2. Load corresponding coefficient CSV files.
#   3. Regenerate observed-versus-predicted plots.
#   4. Add relative-importance bar plots below each prediction plot.
#   5. Save plots and relative-importance CSV summaries.
#
# IMPORTANT:
#   - Does NOT fit or refit any Bayesian models.
#   - Does NOT run CmdStanR or projpred.
#   - Uses the saved prediction and coefficient CSV files.
#   - Overwrites existing plots with the same filenames.
#   - Original variable names in the CSV files are retained.
#
# Relative importance:
#   abs(standardised posterior mean coefficient) /
#   sum(abs(standardised posterior mean coefficients)) * 100
#
#   Intercepts are excluded.
#   This represents relative coefficient magnitude, NOT variance explained.
#
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

required_packages <- c(
  "tidyverse",
  "patchwork",
  "here",
  "scales"
)

missing_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)


# =============================================================================
# 2. DIRECTORIES
# =============================================================================

project_dir <- here::here()

input_dir <- file.path(
  project_dir,
  "outputs",
  "bayesian_C_d14C_predictions_withftir"
)

tables_dir <- file.path(
  input_dir,
  "tables"
)

plots_dir <- file.path(
  input_dir,
  "plots"
)

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(tables_dir)) {
  stop("The tables directory does not exist: ", tables_dir)
}

message("Project directory: ", project_dir)
message("Input tables: ", tables_dir)
message("Output plots: ", plots_dir)


# =============================================================================
# 3. VARIABLE DISPLAY LABELS
#
# These labels are used only in figures.
# Original column names remain unchanged in the CSV files.
# =============================================================================

variable_labels <- c(
  
  # ---------------------------------------------------------------------------
  # Aggregate stability
  # ---------------------------------------------------------------------------
  
  "GMD_h2o" =
    "GMD-H2O (mm)",
  
  "D_h2o" =
    "D-H2O (no units)",
  
  "MWD_h2o" =
    "MWD-H2O (mm)",
  
  "GMD_oxalate_over_h2o" =
    "GMD-Oxalate / GMD-H2O (no units)",
  
  "D_oxalate_over_h2o" =
    "D-Oxalate / D-H2O (no units)",
  
  "MWD_oxalate_over_h2o" =
    "MWD-Oxalate / MWD-H2O (no units)",
  
  
  # ---------------------------------------------------------------------------
  # Carbon and radiocarbon
  # ---------------------------------------------------------------------------
  
  "C_Bulk" =
    "Bulk Carbon (%)",
  
  "C_MAOM" =
    "MAOM Carbon (%)",
  
  "C_oPOM" =
    "oPOM Carbon (%)",
  
  "C_fPOM" =
    "fPOM Carbon (%)",
  
  "d14C_2026_Bulk" =
    "Bulk Δ¹⁴C (‰)",
  
  "D14C_2026_Bulk" =
    "Bulk Δ¹⁴C (‰)",
  
  "d14C_2026_MAOM" =
    "MAOM Δ¹⁴C (‰)",
  
  "D14C_2026_MAOM" =
    "MAOM Δ¹⁴C (‰)",
  
  "d14C_2026_oPOM" =
    "oPOM Δ¹⁴C (‰)",
  
  "D14C_2026_oPOM" =
    "oPOM Δ¹⁴C (‰)",
  
  "d14C_2026_fPOM" =
    "fPOM Δ¹⁴C (‰)",
  
  "D14C_2026_fPOM" =
    "fPOM Δ¹⁴C (‰)",
  
  
  # ---------------------------------------------------------------------------
  # CEC and exchangeable cations
  # ---------------------------------------------------------------------------
  
  "CEC_CEC" =
    "CEC",
  
  "CEC_Base_saturation" =
    "Base saturation (%)",
  
  "CEC_H" =
    "Exchangeable H",
  
  "CEC_Ca" =
    "Exchangeable Ca",
  
  "CEC_K" =
    "Exchangeable K",
  
  "CEC_Mg" =
    "Exchangeable Mg",
  
  "CEC_Na" =
    "Exchangeable Na",
  
  
  # ---------------------------------------------------------------------------
  # FTIR
  # ---------------------------------------------------------------------------
  
  "FTIR_delta_aromatic" =
    "Aromatic (no units)",
  
  "FTIR_delta_Carbonate_bend" =
    "Carbonate bend (no units)",
  
  "FTIR_delta_carbonate_stretch" =
    "Carbonate stretch (no units)",
  
  "FTIR_delta_carbonate stretch" =
    "Carbonate stretch (no units)",
  
  "FTIR_delta_carbonyl" =
    "Carbonyl (no units)",
  
  "FTIR_delta_Alphatic_total" =
    "Aliphatic (no units)",
  
  "FTIR_delta_Aliphatic_total" =
    "Aliphatic (no units)",
  
  "FTIR_delta_COO_total" =
    "Carboxylic (no units)",
  
  "FTIR_delta_polysaccharide" =
    "Polysaccharide (no units)",
  
  "FTIR_delta_poorly_crystalline" =
    "Amorphous (no units)",
  
  "FTIR_delta_Silicate" =
    "Silicate (no units)",
  
  "FTIR_delta_total_area" =
    "Total (no units)",
  
  
  # ---------------------------------------------------------------------------
  # Oxalate / DCB ratios
  # ---------------------------------------------------------------------------
  
  "Fe_ox_over_DCB" =
    "Fe-ox / Fe-DCB (no units)",
  
  "Al_ox_over_DCB" =
    "Al-ox / Al-DCB (no units)",
  
  "Mn_ox_over_DCB" =
    "Mn-ox / Mn-DCB (no units)",
  
  "Si_ox_over_DCB" =
    "Si-ox / Si-DCB (no units)",
  
  "Fe_ox_over_Si_ox" =
    "Fe-ox / Si-ox (no units)",
  
  "Al_ox_over_Si_ox" =
    "Al-ox / Si-ox (no units)",
  
  
  # ---------------------------------------------------------------------------
  # Soil chemistry
  # ---------------------------------------------------------------------------
  
  "soil_pH" =
    "Soil pH (no units)",
  
  "soil_EC" =
    "Soil EC (μS cm⁻¹)",
  
  
  # ---------------------------------------------------------------------------
  # Tessier carbonate fraction: acetic acid
  # ---------------------------------------------------------------------------
  
  "Tessier_Acetic_Ca_mg_g" =
    "Ca-Acetic",
  
  "Tessier_Acetic_Na_mg_g" =
    "Na-Acetic",
  
  "Tessier_Acetic_Mg_mg_g" =
    "Mg-Acetic",
  
  "Tessier_Acetic_K_mg_g" =
    "K-Acetic",
  
  
  # ---------------------------------------------------------------------------
  # Tessier oxidizable fraction: H2O2
  # ---------------------------------------------------------------------------
  
  "Tessier_H2O2_Ca_mg_g" =
    "Ca-H2O2",
  
  "Tessier_H2O2_Na_mg_g" =
    "Na-H2O2",
  
  "Tessier_H2O2_Mg_mg_g" =
    "Mg-H2O2",
  
  "Tessier_H2O2_K_mg_g" =
    "K-H2O2",
  
  
  # ---------------------------------------------------------------------------
  # Tessier reducible fraction: hydroxylamine
  # ---------------------------------------------------------------------------
  
  "Tessier_Hydroxylamine_Ca_mg_g" =
    "Ca-Hydroxylamine",
  
  "Tessier_Hydroxylamine_Na_mg_g" =
    "Na-Hydroxylamine",
  
  "Tessier_Hydroxylamine_Mg_mg_g" =
    "Mg-Hydroxylamine",
  
  "Tessier_Hydroxylamine_K_mg_g" =
    "K-Hydroxylamine"
)

# =============================================================================
# 4. VARIABLE LABEL FUNCTIONS
# =============================================================================

get_variable_labels <- function(variable_names) {
  
  variable_names <- as.character(variable_names)
  
  labels <- unname(
    variable_labels[variable_names]
  )
  
  missing_labels <- is.na(labels)
  
  labels[missing_labels] <- variable_names[missing_labels]
  
  labels
}


# Works for both individual names and vectors
get_variable_label <- function(variable_name) {
  
  get_variable_labels(variable_name)
}

# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================

# ---------------------------------------------------------------------------
# 5.1 Identify predictor column
# ---------------------------------------------------------------------------

identify_predictor_column <- function(coef_data) {
  
  possible_columns <- c(
    "term",
    "predictor",
    "variable",
    "parameter",
    "Parameter"
  )
  
  found <- intersect(
    possible_columns,
    names(coef_data)
  )
  
  if (length(found) == 0) {
    stop(
      "Could not identify predictor column in coefficient CSV. ",
      "Expected one of: ",
      paste(possible_columns, collapse = ", ")
    )
  }
  
  found[1]
}


# ---------------------------------------------------------------------------
# 5.2 Identify coefficient column
# ---------------------------------------------------------------------------

identify_estimate_column <- function(coef_data) {
  
  possible_columns <- c(
    "Estimate",
    "estimate",
    "mean",
    "posterior_mean"
  )
  
  found <- intersect(
    possible_columns,
    names(coef_data)
  )
  
  if (length(found) == 0) {
    stop(
      "Could not identify standardised coefficient column. ",
      "Expected 'Estimate'."
    )
  }
  
  found[1]
}


# ---------------------------------------------------------------------------
# 5.3 Calculate relative importance
# ---------------------------------------------------------------------------

calculate_relative_importance <- function(coef_data) {
  
  predictor_column <- identify_predictor_column(coef_data)
  
  estimate_column <- identify_estimate_column(coef_data)
  
  importance_data <- coef_data %>%
    
    transmute(
      predictor = as.character(.data[[predictor_column]]),
      Estimate = as.numeric(.data[[estimate_column]])
    ) %>%
    
    filter(
      !is.na(predictor),
      !is.na(Estimate)
    ) %>%
    
    # Remove intercepts
    filter(
      !tolower(predictor) %in% c(
        "intercept",
        "(intercept)",
        "constant"
      )
    ) %>%
    
    group_by(predictor) %>%
    
    # Handle duplicated predictor terms, if present
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    mutate(
      absolute_estimate = abs(Estimate)
    )
  
  total_absolute_estimate <- sum(
    importance_data$absolute_estimate,
    na.rm = TRUE
  )
  
  if (
    nrow(importance_data) == 0 ||
    !is.finite(total_absolute_estimate) ||
    total_absolute_estimate == 0
  ) {
    return(
      tibble(
        predictor = character(),
        Estimate = numeric(),
        absolute_estimate = numeric(),
        relative_importance = numeric()
      )
    )
  }
  
  importance_data %>%
    
    mutate(
      relative_importance =
        absolute_estimate / total_absolute_estimate * 100
    ) %>%
    
    arrange(
      relative_importance
    )
}


# ---------------------------------------------------------------------------
# 5.4 Identify response from prediction filename
# ---------------------------------------------------------------------------

identify_response_from_filename <- function(file_name) {
  
  response <- basename(file_name)
  
  response <- sub(
    "^predictions_",
    "",
    response
  )
  
  response <- sub(
    "\\.csv$",
    "",
    response
  )
  
  # Remove projpred suffix, if present
  response <- sub(
    "__projpred_.*$",
    "",
    response
  )
  
  response
}


# ---------------------------------------------------------------------------
# 5.5 Read CSV safely
# ---------------------------------------------------------------------------

read_model_csv <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  readr::read_csv(
    file_path,
    show_col_types = FALSE,
    progress = FALSE
  )
}


# =============================================================================
# 6. PLOT THEME
# =============================================================================

theme_bayesian <- function(base_size = 11) {
  
  theme_bw(base_size = base_size) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = base_size - 1,
        hjust = 0.5,
        colour = "grey30"
      ),
      
      axis.title = element_text(
        face = "bold",
        colour = "black"
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom",
      
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey50"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      plot.margin = margin(
        8, 12, 8, 8
      )
    )
}


# =============================================================================
# 7. OBSERVED-VERSUS-PREDICTED PLOT WITH UNCERTAINTY RIBBONS
# =============================================================================

plot_observed_predicted <- function(
    prediction_data,
    response,
    model_id
) {
  
  required_columns <- c(
    "observed",
    "predicted_mean",
    "fitted_lower_95",
    "fitted_upper_95",
    "predictive_lower_95",
    "predictive_upper_95"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(prediction_data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing prediction columns for ",
      model_id,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  prediction_data <- prediction_data %>%
    mutate(
      observed = as.numeric(observed),
      predicted_mean = as.numeric(predicted_mean),
      fitted_lower_95 = as.numeric(fitted_lower_95),
      fitted_upper_95 = as.numeric(fitted_upper_95),
      predictive_lower_95 = as.numeric(predictive_lower_95),
      predictive_upper_95 = as.numeric(predictive_upper_95)
    ) %>%
    filter(
      is.finite(observed),
      is.finite(predicted_mean)
    )
  
  if (nrow(prediction_data) == 0) {
    stop("No valid prediction observations for ", model_id)
  }
  
  # Treatment column is optional
  if (!"Treatment_model" %in% names(prediction_data)) {
    prediction_data$Treatment_model <- "All observations"
  }
  
  prediction_data <- prediction_data %>%
    mutate(
      Treatment_model = as.factor(Treatment_model)
    )
  
  # Common axis limits
  all_values <- c(
    prediction_data$observed,
    prediction_data$predicted_mean,
    prediction_data$fitted_lower_95,
    prediction_data$fitted_upper_95,
    prediction_data$predictive_lower_95,
    prediction_data$predictive_upper_95
  )
  
  all_values <- all_values[
    is.finite(all_values)
  ]
  
  axis_min <- min(all_values, na.rm = TRUE)
  axis_max <- max(all_values, na.rm = TRUE)
  
  axis_padding <- 0.05 * (axis_max - axis_min)
  
  if (!is.finite(axis_padding) || axis_padding == 0) {
    axis_padding <- 0.1
  }
  
  axis_min <- axis_min - axis_padding
  axis_max <- axis_max + axis_padding
  
  # Plot
  p <- ggplot(
    prediction_data,
    aes(
      x = observed,
      y = predicted_mean,
      shape = Treatment_model
    )
  ) +
    
    # Wider, lighter 95% predictive interval
    geom_errorbar(
      aes(
        ymin = predictive_lower_95,
        ymax = predictive_upper_95
      ),
      width = 0,
      colour = "grey65",
      linewidth = 0.4,
      alpha = 0.35
    ) +
    
    # Narrower, darker 95% fitted credible interval
    geom_errorbar(
      aes(
        ymin = fitted_lower_95,
        ymax = fitted_upper_95
      ),
      width = 0,
      colour = "grey25",
      linewidth = 0.8,
      alpha = 0.8
    ) +
    
    # 1:1 reference line
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.8
    ) +
    
    # Predicted values
    geom_point(
      size = 2.8,
      alpha = 0.85,
      stroke = 0.6
    ) +
    
    scale_x_continuous(
      limits = c(axis_min, axis_max)
    ) +
    
    scale_y_continuous(
      limits = c(axis_min, axis_max)
    ) +
    
    coord_equal() +
    
    labs(
      title = get_variable_label(response),
      
      subtitle = paste0(
        "Observed versus predicted values: ",
        model_id
      ),
      
      x = paste0(
        "Observed ",
        get_variable_label(response)
      ),
      
      y = paste0(
        "Predicted ",
        get_variable_label(response)
      ),
      
      shape = "Treatment"
    ) +
    
    theme_bayesian()
  
  p
}


# =============================================================================
# 8. RELATIVE-IMPORTANCE BAR PLOT
# =============================================================================

plot_relative_importance <- function(
    importance_data,
    response
) {
  
  if (nrow(importance_data) == 0) {
    
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No coefficient data available",
          size = 4
        ) +
        theme_void()
    )
  }
  
  plot_data <- importance_data %>%
    
    mutate(
      predictor_label = get_variable_labels(predictor),
      
      predictor_label = factor(
        predictor_label,
        levels = predictor_label
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = relative_importance,
      y = predictor_label
    )
  ) +
    
    geom_col(
      fill = "grey45",
      colour = "black",
      linewidth = 0.25,
      width = 0.7
    ) +
    
    geom_text(
      aes(
        label = sprintf("%.1f%%", relative_importance)
      ),
      hjust = -0.15,
      size = 3.2
    ) +
    
    scale_x_continuous(
      limits = c(
        0,
        max(plot_data$relative_importance, na.rm = TRUE) * 1.22
      ),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    labs(
      title = "Relative importance of predictors",
      
      subtitle = paste0(
        "Absolute standardised posterior mean coefficients: ",
        get_variable_label(response)
      ),
      
      x = "Relative importance (%)",
      
      y = NULL
    ) +
    
    theme_bayesian(base_size = 9) +
    
    theme(
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0
      ),
      
      plot.subtitle = element_text(
        size = 8,
        hjust = 0
      ),
      
      axis.text.y = element_text(
        size = 8
      ),
      
      axis.title.x = element_text(
        size = 9
      ),
      
      legend.position = "none",
      
      panel.grid.major.y = element_blank(),
      
      plot.margin = margin(
        5, 20, 5, 8
      )
    )
}


# =============================================================================
# 9. COMBINE PREDICTION AND IMPORTANCE PLOTS
# =============================================================================

combine_model_plots <- function(
    prediction_plot,
    importance_plot,
    n_predictors
) {
  
  # Allocate more height to the importance panel when there are more predictors
  importance_height <- max(
    1.1,
    0.22 * n_predictors
  )
  
  prediction_plot /
    
    importance_plot +
    
    patchwork::plot_layout(
      heights = c(
        3,
        importance_height
      )
    )
}


# =============================================================================
# 10. FIND SAVED PREDICTION FILES
# =============================================================================

prediction_files <- list.files(
  tables_dir,
  pattern = "^predictions_.*\\.csv$",
  full.names = TRUE
)

# Exclude combined prediction files
prediction_files <- prediction_files[
  !grepl(
    "all_responses|combined",
    basename(prediction_files),
    ignore.case = TRUE
  )
]

if (length(prediction_files) == 0) {
  
  stop(
    "No prediction CSV files found in: ",
    tables_dir
  )
}

message(
  "\nFound ",
  length(prediction_files),
  " prediction CSV files."
)


# =============================================================================
# 11. PROCESS EACH SAVED MODEL
# =============================================================================

plot_summary <- list()
all_importance <- list()

for (i in seq_along(prediction_files)) {
  
  prediction_file <- prediction_files[i]
  
  model_id <- sub(
    "^predictions_",
    "",
    basename(prediction_file)
  )
  
  model_id <- sub(
    "\\.csv$",
    "",
    model_id
  )
  
  response <- identify_response_from_filename(
    prediction_file
  )
  
  coefficient_file <- file.path(
    tables_dir,
    paste0(
      "coefficients_",
      model_id,
      ".csv"
    )
  )
  
  message(
    "\n--------------------------------------------------"
  )
  
  message(
    "Processing model ",
    i,
    " of ",
    length(prediction_files)
  )
  
  message(
    "Model: ",
    model_id
  )
  
  message(
    "Response: ",
    response
  )
  
  # ---------------------------------------------------------------------------
  # Check corresponding coefficient file
  # ---------------------------------------------------------------------------
  
  if (!file.exists(coefficient_file)) {
    
    warning(
      "No matching coefficient file found for ",
      model_id,
      ". Skipping this model."
    )
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      status = "Skipped: coefficient file missing",
      n_observations = NA_integer_,
      n_predictors = NA_integer_
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Read prediction and coefficient data
  # ---------------------------------------------------------------------------
  
  prediction_data <- tryCatch(
    
    read_model_csv(prediction_file),
    
    error = function(e) {
      warning(
        "Could not read prediction file for ",
        model_id,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  coefficient_data <- tryCatch(
    
    read_model_csv(coefficient_file),
    
    error = function(e) {
      warning(
        "Could not read coefficient file for ",
        model_id,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (
    is.null(prediction_data) ||
    is.null(coefficient_data)
  ) {
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      status = "Skipped: file-reading error",
      n_observations = NA_integer_,
      n_predictors = NA_integer_
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Calculate relative importance
  # ---------------------------------------------------------------------------
  
  importance_data <- tryCatch(
    
    calculate_relative_importance(
      coefficient_data
    ),
    
    error = function(e) {
      warning(
        "Could not calculate importance for ",
        model_id,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (
    is.null(importance_data) ||
    nrow(importance_data) == 0
  ) {
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      status = "Skipped: no valid predictor coefficients",
      n_observations = nrow(prediction_data),
      n_predictors = 0L
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Add model identification to importance table
  # ---------------------------------------------------------------------------
  
  importance_data <- importance_data %>%
    mutate(
      model = model_id,
      response = response,
      response_label = get_variable_labels(response),
      predictor_label = get_variable_labels(predictor)
    ) %>%
    
    select(
      model,
      response,
      response_label,
      predictor,
      predictor_label,
      Estimate,
      absolute_estimate,
      relative_importance
    )
  
  all_importance[[model_id]] <- importance_data
  
  
  # ---------------------------------------------------------------------------
  # Save individual relative-importance CSV
  # ---------------------------------------------------------------------------
  
  importance_csv <- file.path(
    tables_dir,
    paste0(
      "relative_importance_",
      model_id,
      ".csv"
    )
  )
  
  readr::write_csv(
    importance_data,
    importance_csv
  )
  
  
  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------
  
  model_plots <- tryCatch({
    
    prediction_plot <- plot_observed_predicted(
      prediction_data = prediction_data,
      response = response,
      model_id = model_id
    )
    
    importance_plot <- plot_relative_importance(
      importance_data = importance_data,
      response = response
    )
    
    combined_plot <- combine_model_plots(
      prediction_plot = prediction_plot,
      importance_plot = importance_plot,
      n_predictors = nrow(importance_data)
    )
    
    combined_plot
    
  }, error = function(e) {
    
    warning(
      "Plot generation failed for ",
      model_id,
      ": ",
      conditionMessage(e)
    )
    
    NULL
  })
  
  
  if (is.null(model_plots)) {
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      status = "Skipped: plotting error",
      n_observations = nrow(prediction_data),
      n_predictors = nrow(importance_data)
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Save PNG
  # ---------------------------------------------------------------------------
  
  output_png <- file.path(
    plots_dir,
    paste0(
      "observed_vs_predicted_",
      model_id,
      ".png"
    )
  )
  
  ggsave(
    filename = output_png,
    plot = model_plots,
    width = 10,
    height = max(
      7,
      5 + 0.22 * nrow(importance_data)
    ),
    dpi = 400,
    bg = "white"
  )
  
  
  # ---------------------------------------------------------------------------
  # Save PDF
  # ---------------------------------------------------------------------------
  
  output_pdf <- file.path(
    plots_dir,
    paste0(
      "observed_vs_predicted_",
      model_id,
      ".pdf"
    )
  )
  
  ggsave(
    filename = output_pdf,
    plot = model_plots,
    width = 10,
    height = max(
      7,
      5 + 0.22 * nrow(importance_data)
    ),
    bg = "white"
  )
  
  
  # ---------------------------------------------------------------------------
  # Record successful model
  # ---------------------------------------------------------------------------
  
  plot_summary[[model_id]] <- tibble(
    model = model_id,
    response = response,
    response_label = get_variable_label(response),
    status = "Successfully plotted",
    n_observations = nrow(prediction_data),
    n_predictors = nrow(importance_data),
    png_file = basename(output_png),
    pdf_file = basename(output_pdf)
  )
  
  
  message(
    "Successfully saved: ",
    basename(output_png)
  )
  
}


# =============================================================================
# 12. SAVE COMBINED RELATIVE-IMPORTANCE TABLE
# =============================================================================

if (length(all_importance) > 0) {
  
  all_importance_df <- bind_rows(
    all_importance
  )
  
  readr::write_csv(
    all_importance_df,
    file.path(
      tables_dir,
      "relative_importance_all_models.csv"
    )
  )
  
  message(
    "\nSaved combined relative-importance table."
  )
  
} else {
  
  warning(
    "No relative-importance tables were generated."
  )
}


# =============================================================================
# 13. SAVE PLOT GENERATION SUMMARY
# =============================================================================

plot_summary_df <- bind_rows(
  plot_summary
)

readr::write_csv(
  plot_summary_df,
  file.path(
    tables_dir,
    "plot_generation_summary.csv"
  )
)


# =============================================================================
# 14. FINAL SUMMARY
# =============================================================================

n_success <- sum(
  plot_summary_df$status == "Successfully plotted",
  na.rm = TRUE
)

n_skipped <- nrow(plot_summary_df) - n_success

message(
  "\n=================================================="
)

message(
  "BAYESIAN PLOT GENERATION COMPLETE"
)

message(
  "=================================================="
)

message(
  "Prediction files found: ",
  length(prediction_files)
)

message(
  "Successfully plotted: ",
  n_success
)

message(
  "Skipped: ",
  n_skipped
)

message(
  "\nPlots saved to:"
)

message(
  plots_dir
)

message(
  "\nTables saved to:"
)

message(
  tables_dir
)

message(
  "\nPlot generation summary:"
)

print(plot_summary_df)

message(
  "\n=================================================="
)
