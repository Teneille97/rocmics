
# =============================================================================
# MALLE ERW: BAYESIAN C AND d14C BEST-MODEL PLOT-ONLY SCRIPT
#
# Purpose:
#   1. Load existing prediction CSV files for best-selected models only.
#   2. Load their corresponding best-model coefficient CSV files.
#   3. Regenerate observed-versus-predicted plots.
#   4. Show 95% fitted credible intervals and predictive intervals.
#   5. Display multiline coefficient equations below prediction plots.
#   6. Add relative-importance bar plots below coefficient equations.
#   7. Save individual and combined relative-importance CSV summaries.
#   8. Save a summary of successfully plotted and skipped models.
#
# IMPORTANT:
#   - Does NOT fit or refit any Bayesian models.
#   - Does NOT run CmdStanR or projpred.
#   - Does NOT perform model selection.
#   - Uses existing saved prediction and coefficient CSV files only.
#   - Processes ONLY files ending in _bestmod.csv.
#   - Ignores all candidate model files.
#   - Overwrites existing plots with the same filenames.
#
# COEFFICIENT EQUATIONS:
#   - Uses standardised posterior mean coefficients.
#   - Displays the five strongest predictors.
#   - Positive coefficients are indicated by (+).
#   - Negative coefficients are indicated by (-).
#   - Includes the intercept where available.
#   - Automatically wraps equations across multiple lines.
#   - Appears below the prediction plot without a border.
#
# RELATIVE IMPORTANCE:
#   abs(standardised posterior mean coefficient) /
#   sum(abs(standardised posterior mean coefficients)) * 100
#
#   Intercepts are excluded.
#   This represents relative coefficient magnitude,
#   NOT variance explained.
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
  "best_model",
  "tables"
)

plots_dir <- file.path(
  input_dir,
  "best_model",
  "plots"
)

if (!dir.exists(tables_dir)) {
  
  stop(
    "The best-model tables directory does not exist: ",
    tables_dir
  )
}

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project directory: ", project_dir)
message("Best-model input tables: ", tables_dir)
message("Best-model output plots: ", plots_dir)


# =============================================================================
# 3. VARIABLE DISPLAY LABELS
# =============================================================================

variable_labels <- c(
  
  # ---------------------------------------------------------------------------
  # Aggregate stability
  # ---------------------------------------------------------------------------
  
  "GMD_h2o" = "GMD-H2O (mm)",
  "D_h2o" = "D-H2O (no units)",
  "MWD_h2o" = "MWD-H2O (mm)",
  
  "GMD_oxalate_over_h2o" =
    "GMD-Oxalate / GMD-H2O (no units)",
  
  "D_oxalate_over_h2o" =
    "D-Oxalate / D-H2O (no units)",
  
  "MWD_oxalate_over_h2o" =
    "MWD-Oxalate / MWD-H2O (no units)",
  
  
  # ---------------------------------------------------------------------------
  # Carbon and radiocarbon
  # ---------------------------------------------------------------------------
  
  "C_Bulk" = "Bulk Carbon (%)",
  "C_MAOM" = "MAOM Carbon (%)",
  "C_oPOM" = "oPOM Carbon (%)",
  "C_fPOM" = "fPOM Carbon (%)",
  
  "d14C_2026_Bulk" = "Bulk Δ¹⁴C (‰)",
  "D14C_2026_Bulk" = "Bulk Δ¹⁴C (‰)",
  
  "d14C_2026_MAOM" = "MAOM Δ¹⁴C (‰)",
  "D14C_2026_MAOM" = "MAOM Δ¹⁴C (‰)",
  
  "d14C_2026_oPOM" = "oPOM Δ¹⁴C (‰)",
  "D14C_2026_oPOM" = "oPOM Δ¹⁴C (‰)",
  
  "d14C_2026_fPOM" = "fPOM Δ¹⁴C (‰)",
  "D14C_2026_fPOM" = "fPOM Δ¹⁴C (‰)",
  
  
  # ---------------------------------------------------------------------------
  # CEC and exchangeable cations
  # ---------------------------------------------------------------------------
  
  "CEC_CEC" = "CEC",
  "CEC_Base_saturation" = "Base saturation (%)",
  "CEC_H" = "Exchangeable H",
  "CEC_Ca" = "Exchangeable Ca",
  "CEC_K" = "Exchangeable K",
  "CEC_Mg" = "Exchangeable Mg",
  "CEC_Na" = "Exchangeable Na",
  
  
  # ---------------------------------------------------------------------------
  # FTIR
  # ---------------------------------------------------------------------------
  
  "FTIR_delta_aromatic" = "Aromatic (no units)",
  
  "FTIR_delta_Carbonate_bend" =
    "Carbonate bend (no units)",
  
  "FTIR_delta_carbonate_stretch" =
    "Carbonate stretch (no units)",
  
  "FTIR_delta_carbonate stretch" =
    "Carbonate stretch (no units)",
  
  "FTIR_delta_carbonyl" = "Carbonyl (no units)",
  
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
  
  "Fe_ox_over_DCB" = "Fe-ox / Fe-DCB (no units)",
  "Al_ox_over_DCB" = "Al-ox / Al-DCB (no units)",
  "Mn_ox_over_DCB" = "Mn-ox / Mn-DCB (no units)",
  "Si_ox_over_DCB" = "Si-ox / Si-DCB (no units)",
  
  "Fe_ox_over_Si_ox" = "Fe-ox / Si-ox (no units)",
  "Al_ox_over_Si_ox" = "Al-ox / Si-ox (no units)",
  
  
  # ---------------------------------------------------------------------------
  # Soil chemistry
  # ---------------------------------------------------------------------------
  
  "soil_pH" = "Soil pH (no units)",
  "soil_EC" = "Soil EC (μS cm⁻¹)",
  
  
  # ---------------------------------------------------------------------------
  # Tessier carbonate fraction: acetic acid
  # ---------------------------------------------------------------------------
  
  "Tessier_Acetic_Ca_mg_g" = "Ca-Acetic",
  "Tessier_Acetic_Na_mg_g" = "Na-Acetic",
  "Tessier_Acetic_Mg_mg_g" = "Mg-Acetic",
  "Tessier_Acetic_K_mg_g" = "K-Acetic",
  
  
  # ---------------------------------------------------------------------------
  # Tessier oxidizable fraction: H2O2
  # ---------------------------------------------------------------------------
  
  "Tessier_H2O2_Ca_mg_g" = "Ca-H2O2",
  "Tessier_H2O2_Na_mg_g" = "Na-H2O2",
  "Tessier_H2O2_Mg_mg_g" = "Mg-H2O2",
  "Tessier_H2O2_K_mg_g" = "K-H2O2",
  
  "Tessier_ H2O2_Na_mg_g" = "Na-H2O2",
  "Tessier_ H2O2_Mg_mg_g" = "Mg-H2O2",
  "Tessier_ H2O2_K_mg_g" = "K-H2O2",
  
  
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
    "K-Hydroxylamine",
  
  "Tessier_ Hydroxylamine _Na_mg_g" =
    "Na-Hydroxylamine",
  
  "Tessier_ Hydroxylamine _Mg_mg_g" =
    "Mg-Hydroxylamine",
  
  "Tessier_ Hydroxylamine _K_mg_g" =
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


get_variable_label <- function(variable_name) {
  
  get_variable_labels(variable_name)
}


# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1 Identify predictor column
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# 5.2 Identify coefficient column
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# 5.3 Identify intercept rows
# -----------------------------------------------------------------------------

is_intercept <- function(x) {
  
  tolower(trimws(as.character(x))) %in% c(
    "intercept",
    "(intercept)",
    "constant",
    "alpha",
    "beta0"
  )
}


# -----------------------------------------------------------------------------
# 5.4 Calculate relative importance
# -----------------------------------------------------------------------------

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
    
    filter(
      !is_intercept(predictor)
    ) %>%
    
    group_by(predictor) %>%
    
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


# -----------------------------------------------------------------------------
# 5.5 Extract intercept
# -----------------------------------------------------------------------------

extract_intercept <- function(coef_data) {
  
  predictor_column <- identify_predictor_column(coef_data)
  
  estimate_column <- identify_estimate_column(coef_data)
  
  intercept_data <- coef_data %>%
    
    transmute(
      predictor = as.character(.data[[predictor_column]]),
      Estimate = as.numeric(.data[[estimate_column]])
    ) %>%
    
    filter(
      is_intercept(predictor),
      is.finite(Estimate)
    )
  
  if (nrow(intercept_data) == 0) {
    return(0)
  }
  
  mean(intercept_data$Estimate)
}


# -----------------------------------------------------------------------------
# 5.6 Build multiline coefficient equation
# -----------------------------------------------------------------------------

build_coefficient_equation <- function(
    coef_data,
    max_predictors = 5,
    predictors_per_line = 3
) {
  
  predictor_column <- identify_predictor_column(coef_data)
  
  estimate_column <- identify_estimate_column(coef_data)
  
  coefficient_data <- coef_data %>%
    
    transmute(
      predictor = as.character(.data[[predictor_column]]),
      Estimate = as.numeric(.data[[estimate_column]])
    ) %>%
    
    filter(
      !is.na(predictor),
      is.finite(Estimate)
    )
  
  intercept <- extract_intercept(coef_data)
  
  predictor_data <- coefficient_data %>%
    
    filter(
      !is_intercept(predictor)
    ) %>%
    
    group_by(predictor) %>%
    
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    arrange(
      desc(abs(Estimate))
    ) %>%
    
    slice_head(
      n = max_predictors
    )
  
  if (nrow(predictor_data) == 0) {
    
    return(
      paste0(
        "\u0177 = ",
        sprintf("%.2f", intercept)
      )
    )
  }
  
  equation_terms <- vapply(
    seq_len(nrow(predictor_data)),
    function(i) {
      
      estimate <- predictor_data$Estimate[i]
      
      predictor_label <- get_variable_label(
        predictor_data$predictor[i]
      )
      
      paste0(
        ifelse(estimate >= 0, "+", "\u2212"),
        " ",
        sprintf("%.2f", abs(estimate)),
        " \u00D7 ",
        predictor_label
      )
    },
    character(1)
  )
  
  # Split equation into groups to prevent long lines.
  
  term_groups <- split(
    equation_terms,
    ceiling(
      seq_along(equation_terms) / predictors_per_line
    )
  )
  
  equation_lines <- vapply(
    seq_along(term_groups),
    function(i) {
      
      paste(
        term_groups[[i]],
        collapse = " "
      )
    },
    character(1)
  )
  
  # Add intercept to final line.
  
  intercept_term <- paste0(
    ifelse(intercept >= 0, "+", "\u2212"),
    " ",
    sprintf("%.2f", abs(intercept))
  )
  
  equation_lines[length(equation_lines)] <- paste(
    equation_lines[length(equation_lines)],
    intercept_term
  )
  
  paste(
    c(
      "\u0177 =",
      equation_lines
    ),
    collapse = "\n"
  )
}


# -----------------------------------------------------------------------------
# 5.7 Identify response from filename
# -----------------------------------------------------------------------------

identify_response_from_filename <- function(file_name) {
  
  response <- basename(file_name)
  
  response <- sub(
    "^predictions_",
    "",
    response
  )
  
  response <- sub(
    "_bestmod\\.csv$",
    "",
    response
  )
  
  response
}


# -----------------------------------------------------------------------------
# 5.8 Read CSV safely
# -----------------------------------------------------------------------------

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
# 7. OBSERVED-VERSUS-PREDICTED PLOT
# =============================================================================

plot_prediction <- function(
    prediction_data,
    response_label,
    model_id
) {
  
  prediction_data <- prediction_data %>%
    
    mutate(
      across(
        all_of(c(
          "observed",
          "predicted_mean",
          "fitted_lower_95",
          "fitted_upper_95",
          "predictive_lower_95",
          "predictive_upper_95"
        )),
        as.numeric
      )
    ) %>%
    
    filter(
      is.finite(observed),
      is.finite(predicted_mean),
      is.finite(fitted_lower_95),
      is.finite(fitted_upper_95),
      is.finite(predictive_lower_95),
      is.finite(predictive_upper_95)
    )
  
  ggplot(
    prediction_data,
    aes(
      x = observed,
      y = predicted_mean
    )
  ) +
    
    # 1. 95% predictive interval
    
    geom_linerange(
      aes(
        ymin = predictive_lower_95,
        ymax = predictive_upper_95
      ),
      colour = "grey55",
      linewidth = 0.7,
      alpha = 0.7
    ) +
    
    # 2. 95% fitted credible interval
    
    geom_linerange(
      aes(
        ymin = fitted_lower_95,
        ymax = fitted_upper_95
      ),
      colour = "steelblue",
      linewidth = 3,
      alpha = 0.35
    ) +
    
    # 3. 1:1 reference line
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.7
    ) +
    
    # 4. Predicted values
    
    geom_point(
      colour = "black",
      fill = "steelblue",
      shape = 21,
      size = 2.5,
      alpha = 0.9
    ) +
    
    labs(
      title = response_label,
      
      subtitle = paste(
        "Best model:",
        model_id
      ),
      
      x = "Observed",
      
      y = "Predicted"
    ) +
    
    coord_equal() +
    
    theme_bayesian(base_size = 11)
}


# =============================================================================
# 8. COEFFICIENT EQUATION PANEL
# =============================================================================

plot_equation <- function(coefficient_equation) {
  
  equation_data <- tibble(
    x = 0,
    y = 1,
    equation = coefficient_equation
  )
  
  ggplot(
    equation_data,
    aes(
      x = x,
      y = y,
      label = equation
    )
  ) +
    
    # Equation heading
    
    annotate(
      "text",
      x = 0,
      y = 1.45,
      label = "Standardised posterior mean coefficients",
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 3.0,
      colour = "grey25"
    ) +
    
    # Multiline equation without a border or background
    
    geom_text(
      hjust = 0,
      vjust = 1,
      size = 3.1,
      lineheight = 1.2,
      colour = "grey20"
    ) +
    
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(-0.1, 1.55),
      clip = "off"
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        5, 15, 5, 15
      )
    )
}


# =============================================================================
# 9. RELATIVE-IMPORTANCE BAR PLOT
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
        label = sprintf(
          "%.1f%%",
          relative_importance
        )
      ),
      hjust = -0.15,
      size = 3.2
    ) +
    
    scale_x_continuous(
      limits = c(
        0,
        max(
          plot_data$relative_importance,
          na.rm = TRUE
        ) * 1.22
      ),
      expand = expansion(
        mult = c(0, 0.02)
      )
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
# 10. COMBINE PREDICTION, EQUATION AND IMPORTANCE PLOTS
# =============================================================================

combine_model_plots <- function(
    prediction_plot,
    equation_plot,
    importance_plot,
    n_predictors
) {
  
  # Dedicated space for multiline coefficient equation.
  
  equation_height <- 0.9
  
  importance_height <- max(
    1.1,
    0.22 * n_predictors
  )
  
  prediction_plot /
    
    equation_plot /
    
    importance_plot +
    
    patchwork::plot_layout(
      heights = c(
        3,
        equation_height,
        importance_height
      )
    )
}


# =============================================================================
# 11. FIND SAVED BEST-MODEL PREDICTION FILES ONLY
# =============================================================================

prediction_files <- list.files(
  path = tables_dir,
  pattern = "^predictions_.*_bestmod\\.csv$",
  full.names = TRUE,
  recursive = FALSE
)

prediction_files <- sort(
  prediction_files
)

if (length(prediction_files) == 0) {
  
  stop(
    paste0(
      "\nNo saved best-model prediction CSV files were found.\n\n",
      "Expected filenames such as:\n",
      "  predictions_C_Bulk_bestmod.csv\n",
      "  predictions_C_MAOM_bestmod.csv\n\n",
      "Directory searched:\n",
      tables_dir,
      "\n\nCheck that the original Bayesian modelling script ",
      "successfully generated the best-model prediction files."
    )
  )
}

message(
  "\nFound ",
  length(prediction_files),
  " saved best-model prediction files."
)

message(
  "\nBest-model files to be processed:"
)

for (f in prediction_files) {
  message("  - ", basename(f))
}


# =============================================================================
# 12. PROCESS EACH SAVED BEST MODEL
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
    "Processing best model ",
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
  # Check matching coefficient file
  # ---------------------------------------------------------------------------
  
  if (!file.exists(coefficient_file)) {
    
    warning(
      "Matching coefficient file not found for ",
      model_id,
      ". Expected: ",
      basename(coefficient_file),
      ". Skipping this model."
    )
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      response_label = get_variable_label(response),
      status = "Skipped: matching best-model coefficient file missing",
      n_observations = NA_integer_,
      n_predictors = NA_integer_,
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Read existing prediction and coefficient files
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
      response_label = get_variable_label(response),
      status = "Skipped: file-reading error",
      n_observations = NA_integer_,
      n_predictors = NA_integer_,
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Check required prediction columns
  # ---------------------------------------------------------------------------
  
  required_prediction_columns <- c(
    "observed",
    "predicted_mean",
    "fitted_lower_95",
    "fitted_upper_95",
    "predictive_lower_95",
    "predictive_upper_95"
  )
  
  missing_prediction_columns <- setdiff(
    required_prediction_columns,
    names(prediction_data)
  )
  
  if (length(missing_prediction_columns) > 0) {
    
    warning(
      "Prediction file for ",
      model_id,
      " is missing required columns: ",
      paste(
        missing_prediction_columns,
        collapse = ", "
      ),
      ". Skipping this model."
    )
    
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      response_label = get_variable_label(response),
      status = "Skipped: required prediction columns missing",
      n_observations = nrow(prediction_data),
      n_predictors = NA_integer_,
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
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
        "Could not calculate relative importance for ",
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
      response_label = get_variable_label(response),
      status = "Skipped: no valid predictor coefficients",
      n_observations = nrow(prediction_data),
      n_predictors = 0L,
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # Build multiline coefficient equation
  # ---------------------------------------------------------------------------
  
  coefficient_equation <- tryCatch(
    
    build_coefficient_equation(
      coef_data = coefficient_data,
      max_predictors = 5,
      predictors_per_line = 3
    ),
    
    error = function(e) {
      
      warning(
        "Could not build coefficient equation for ",
        model_id,
        ": ",
        conditionMessage(e)
      )
      
      "Coefficient equation unavailable"
    }
  )
  
  message(
    "\nCoefficient equation:\n",
    coefficient_equation
  )
  
  
  # ---------------------------------------------------------------------------
  # Add model identification to importance table
  # ---------------------------------------------------------------------------
  
  importance_data <- importance_data %>%
    
    mutate(
      model = model_id,
      response = response,
      response_label = get_variable_label(response),
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
  # Generate prediction, equation and importance plots
  # ---------------------------------------------------------------------------
  
  model_plots <- tryCatch({
    
    prediction_plot <- plot_prediction(
      
      prediction_data = prediction_data,
      
      response_label = get_variable_label(response),
      
      model_id = model_id
    )
    
    equation_plot <- plot_equation(
      coefficient_equation
    )
    
    importance_plot <- plot_relative_importance(
      
      importance_data = importance_data,
      
      response = response
    )
    
    combine_model_plots(
      
      prediction_plot = prediction_plot,
      
      equation_plot = equation_plot,
      
      importance_plot = importance_plot,
      
      n_predictors = nrow(importance_data)
    )
    
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
      response_label = get_variable_label(response),
      status = "Skipped: plotting error",
      n_observations = nrow(prediction_data),
      n_predictors = nrow(importance_data),
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
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
      5.5 + 0.22 * nrow(importance_data)
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
      5.5 + 0.22 * nrow(importance_data)
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
    coefficient_equation = coefficient_equation,
    prediction_file = basename(prediction_file),
    coefficient_file = basename(coefficient_file),
    png_file = basename(output_png),
    pdf_file = basename(output_pdf)
  )
  
  message(
    "Successfully saved: ",
    basename(output_png)
  )
  
}


# =============================================================================
# 13. SAVE COMBINED RELATIVE-IMPORTANCE TABLE
# =============================================================================

if (length(all_importance) > 0) {
  
  all_importance_df <- bind_rows(
    all_importance
  )
  
  readr::write_csv(
    all_importance_df,
    file.path(
      tables_dir,
      "relative_importance_best_models_all_responses.csv"
    )
  )
  
  message(
    "\nSaved combined relative-importance table for best models."
  )
  
} else {
  
  warning(
    "No relative-importance tables were generated."
  )
}


# =============================================================================
# 14. SAVE BEST-MODEL PLOT GENERATION SUMMARY
# =============================================================================

if (length(plot_summary) > 0) {
  
  plot_summary_df <- bind_rows(
    plot_summary
  )
  
} else {
  
  plot_summary_df <- tibble(
    model = character(),
    response = character(),
    response_label = character(),
    status = character(),
    n_observations = integer(),
    n_predictors = integer()
  )
}

readr::write_csv(
  plot_summary_df,
  file.path(
    tables_dir,
    "best_model_plot_generation_summary.csv"
  )
)


# =============================================================================
# 15. FINAL SUMMARY
# =============================================================================

n_success <- sum(
  plot_summary_df$status == "Successfully plotted",
  na.rm = TRUE
)

n_skipped <- sum(
  grepl(
    "^Skipped:",
    plot_summary_df$status
  ),
  na.rm = TRUE
)

message(
  "\n=================================================="
)

message(
  "BEST-MODEL PLOT GENERATION COMPLETE"
)

message(
  "=================================================="
)

message(
  "Saved best-model prediction files found: ",
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
  "\nOnly saved _bestmod.csv model outputs were processed."
)

message(
  "No Bayesian models were fitted or refitted."
)

message(
  "\nCoefficient equations display the five strongest predictors."
)

message(
  "Equations are displayed below prediction plots without borders."
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
  "\nBest-model plot generation summary:"
)

print(plot_summary_df)

message(
  "\n=================================================="
)
