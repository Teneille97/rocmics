# =============================================================================
# MALLE ERW: BAYESIAN C AND d14C BEST-MODEL PLOT-ONLY SCRIPT
#
# Purpose:
#   1. Load existing prediction CSV files for best-selected models only.
#   2. Load their corresponding best-model coefficient CSV files.
#   3. Plot observed vs. posterior predictive replicated values (y_rep) with
#      95% fitted credible intervals (mu), 95% posterior predictive intervals,
#      R², and RMSE.
#   4. Display coefficient equations and relative-importance bar plots
#      directly beneath each observed-vs-replicated plot.
#   5. Save individual and combined relative-importance CSV summaries,
#      2x2 observed-vs-replicated panels (including equations and importance
#      barplots in each cell), and posterior residual plots.
#
# IMPORTANT:
#   - Does NOT fit or refit any Bayesian models.
#   - Does NOT run CmdStanR or projpred.
#   - Does NOT perform model selection.
#   - Uses existing saved prediction and coefficient CSV files only.
#   - Processes ONLY files ending in _bestmod.csv.
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

required_packages <- c(
  "tidyverse",
  "patchwork",
  "here",
  "scales",
  "posterior"
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
# 2. PROJECT PATHS AND DATA PREPARATION
# =============================================================================

project_dir <- normalizePath(
  "C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics",
  winslash = "/",
  mustWork = TRUE
)

setwd(project_dir)

model_script <- file.path(
  project_dir,
  "scripts",
  "Malle_Bayesian_C_d14C_pool_predictions_withftir.R"
)

if (!file.exists(model_script)) {
  stop(
    "Cannot find the original Bayesian model script: ",
    model_script
  )
}

input_file <- file.path(
  project_dir,
  "outputs",
  "integrated_geochemistry_for_bayesian.csv"
)

out_dir <- file.path(
  project_dir,
  "outputs",
  "bayesian_C_d14C_predictions_withftir"
)

table_dir <- file.path(out_dir, "tables")

best_dir <- file.path(out_dir, "best_model")
best_model_dir <- file.path(best_dir, "models")
best_table_dir <- file.path(best_dir, "tables")
best_plot_dir <- file.path(best_dir, "plots")
best_diag_dir <- file.path(best_dir, "diagnostics")

min_predictor_observed_fraction <- 0.70
max_abs_predictor_correlation <- 0.95

for (d in c(
  out_dir, table_dir, best_dir, best_model_dir,
  best_table_dir, best_plot_dir, best_diag_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

model_script_lines <- readLines(model_script, warn = FALSE)

prep_start <- grep(
  "^# 1\\. READ AND FILTER DATA",
  model_script_lines
)[1]

prep_end <- grep(
  "^# 4B\\. Extract predictions",
  model_script_lines
)[1]

if (is.na(prep_start) || is.na(prep_end) || prep_end <= prep_start) {
  stop(
    "Could not identify the data-preparation range in the original model script. ",
    "Expected Sections 1 through 4A, ending before Section 4B."
  )
}

eval(
  parse(text = model_script_lines[seq.int(prep_start, prep_end - 1L)]),
  envir = .GlobalEnv
)

dat_for_residuals <- get("dat", envir = .GlobalEnv)
candidate_predictors_for_residuals <- get(
  "candidate_predictors",
  envir = .GlobalEnv
)
prepare_model_data_for_residuals <- get(
  "prepare_model_data",
  envir = .GlobalEnv
)
safe_name_for_files <- get("safe_name", envir = .GlobalEnv)

stopifnot(
  is.data.frame(dat_for_residuals),
  length(candidate_predictors_for_residuals) > 0L,
  is.function(prepare_model_data_for_residuals),
  is.function(safe_name_for_files)
)

tables_dir <- best_table_dir
plots_dir <- best_plot_dir

if (!dir.exists(tables_dir)) {
  stop("The best-model tables directory does not exist: ", tables_dir)
}

message("Project directory: ", project_dir)
message("Best-model input tables: ", tables_dir)
message("Best-model output plots: ", plots_dir)
message("Residual diagnostics directory: ", best_diag_dir)


# =============================================================================
# 3. VARIABLE DISPLAY LABELS
# =============================================================================

variable_labels <- c(
  # Aggregate stability
  "GMD_h2o" = "GMD-H₂O (mm)",
  "D_h2o" = "D-H₂O",
  "MWD_h2o" = "MWD-H₂O (mm)",
  "GMD_oxalate_over_h2o" = "GMD-Oxalate / GMD-H₂O",
  "D_oxalate_over_h2o" = "D-Oxalate / D-H₂O",
  "MWD_oxalate_over_h2o" = "MWD-Oxalate / MWD-H₂O",
  
  # Carbon and radiocarbon
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
  
  # CEC and exchangeable cations
  "CEC_CEC" = "CEC",
  "CEC_Base_saturation" = "Base saturation (%)",
  "CEC_H" = "Exchangeable H",
  "CEC_Ca" = "Exchangeable Ca",
  "CEC_K" = "Exchangeable K",
  "CEC_Mg" = "Exchangeable Mg",
  "CEC_Na" = "Exchangeable Na",
  "CEC_Mn" = "Exchangeable Mn",
  
  # FTIR
  "FTIR_delta_aromatic" = "Aromatic",
  "FTIR_delta_Carbonate_bend" = "Carbonate bend",
  "FTIR_delta_carbonate_stretch" = "Carbonate stretch",
  "FTIR_delta_carbonate stretch" = "Carbonate stretch",
  "FTIR_delta_carbonyl" = "Carbonyl",
  "FTIR_delta_Alphatic_total" = "Aliphatic",
  "FTIR_delta_Aliphatic_total" = "Aliphatic",
  "FTIR_delta_COO_total" = "Carboxylic",
  "FTIR_delta_polysaccharide" = "Polysaccharide",
  "FTIR_delta_poorly_crystalline" = "Amorphous",
  "FTIR_delta_Silicate" = "Silicate",
  "FTIR_delta_total_area" = "Total",
  
  # Oxalate / DCB ratios
  "Fe_ox_over_DCB" = "Fe-ox / Fe-DCB",
  "Al_ox_over_DCB" = "Al-ox / Al-DCB",
  "Mn_ox_over_DCB" = "Mn-ox / Mn-DCB",
  "Si_ox_over_DCB" = "Si-ox / Si-DCB",
  "Fe_ox_over_Si_ox" = "Fe-ox / Si-ox",
  "Al_ox_over_Si_ox" = "Al-ox / Si-ox",
  
  # Soil chemistry
  "soil_pH" = "Soil pH",
  "soil_EC" = "Soil EC (μS cm⁻¹)",
  
  # Tessier fractions
  "Tessier_Acetic_Ca_mg_g" = "Ca-Acetic",
  "Tessier_Acetic_Na_mg_g" = "Na-Acetic",
  "Tessier_Acetic_Mg_mg_g" = "Mg-Acetic",
  "Tessier_Acetic_K_mg_g" = "K-Acetic",
  "Tessier_H2O2_Ca_mg_g" = "Ca-H₂O₂",
  "Tessier_H2O2_Na_mg_g" = "Na-H₂O₂",
  "Tessier_H2O2_Mg_mg_g" = "Mg-H₂O₂",
  "Tessier_H2O2_K_mg_g" = "K-H₂O₂",
  "Tessier_Hydroxylamine_Ca_mg_g" = "Ca-Hydroxylamine",
  "Tessier_Hydroxylamine_Na_mg_g" = "Na-Hydroxylamine",
  "Tessier_Hydroxylamine_Mg_mg_g" = "Mg-Hydroxylamine",
  "Tessier_Hydroxylamine_K_mg_g" = "K-Hydroxylamine",
  "Tessier_Hydroxylamine_Fe_mg_g" = "Fe-Hydroxylamine",
  "Tessier_Hydroxylamine_Si_mg_g" = "Si-Hydroxylamine",
  
  "pH_Oct2023" = "pH Oct 2023",
  "Al_Oct2023" = "Al Oct 2023",
  "Mn_Oct2023" = "Mn Oct 2023"
)


# =============================================================================
# 4. VARIABLE LABEL FUNCTIONS
# =============================================================================

get_variable_labels <- function(variable_names) {
  variable_names <- as.character(variable_names)
  labels <- unname(variable_labels[variable_names])
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

identify_predictor_column <- function(coef_data) {
  possible_columns <- c(
    "term",
    "predictor",
    "variable",
    "parameter",
    "Parameter"
  )
  found <- intersect(possible_columns, names(coef_data))
  if (length(found) == 0) {
    stop(
      "Could not identify predictor column in coefficient CSV. ",
      "Expected one of: ",
      paste(possible_columns, collapse = ", ")
    )
  }
  found[1]
}

identify_estimate_column <- function(coef_data) {
  possible_columns <- c(
    "Estimate",
    "estimate",
    "mean",
    "posterior_mean"
  )
  found <- intersect(possible_columns, names(coef_data))
  if (length(found) == 0) {
    stop(
      "Could not identify standardised coefficient column. ",
      "Expected one of: ",
      paste(possible_columns, collapse = ", ")
    )
  }
  found[1]
}

is_intercept <- function(x) {
  tolower(trimws(as.character(x))) %in% c(
    "intercept",
    "(intercept)",
    "constant",
    "alpha",
    "beta0"
  )
}

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
      is.finite(Estimate),
      !is_intercept(predictor)
    ) %>%
    group_by(predictor) %>%
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(abs(Estimate) > 1e-8) %>%
    mutate(absolute_estimate = abs(Estimate))
  
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
    arrange(relative_importance)
}

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

build_coefficient_equation <- function(
    coef_data,
    predictors_per_line = 3,
    zero_tolerance = 1e-8
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
    filter(!is_intercept(predictor)) %>%
    group_by(predictor) %>%
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(abs(Estimate) > zero_tolerance) %>%
    arrange(desc(abs(Estimate)))
  
  intercept_term <- sprintf("%.2f", intercept)
  
  if (nrow(predictor_data) == 0) {
    return(paste0("\u0177 = ", intercept_term))
  }
  
  equation_terms <- vapply(
    seq_len(nrow(predictor_data)),
    function(i) {
      estimate <- predictor_data$Estimate[i]
      predictor_label <- get_variable_label(predictor_data$predictor[i])
      
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
  
  term_groups <- split(
    equation_terms,
    ceiling(seq_along(equation_terms) / predictors_per_line)
  )
  
  equation_lines <- vapply(
    term_groups,
    function(group) paste(group, collapse = " "),
    character(1)
  )
  
  equation_lines <- c(
    paste0("\u0177 = ", intercept_term),
    equation_lines
  )
  
  paste(equation_lines, collapse = "\n")
}

identify_response_from_filename <- function(file_name) {
  response <- basename(file_name)
  response <- sub("^predictions_", "", response)
  response <- sub("_bestmod\\.csv$", "", response)
  response
}

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
      plot.margin = margin(8, 12, 8, 8)
    )
}


# =============================================================================
# 7. OBSERVED-VERSUS-REPLICATED PLOT
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
  
  prediction_r2 <- if (
    nrow(prediction_data) >= 2 &&
    is.finite(stats::var(prediction_data$observed)) &&
    stats::var(prediction_data$observed) > 0
  ) {
    1 - sum(
      (prediction_data$observed - prediction_data$predicted_mean)^2
    ) / sum(
      (prediction_data$observed - mean(prediction_data$observed))^2
    )
  } else {
    NA_real_
  }
  
  prediction_rmse <- sqrt(
    mean(
      (prediction_data$observed - prediction_data$predicted_mean)^2
    )
  )
  
  performance_label <- paste0(
    "Best model: ", model_id,
    "\nR² = ",
    ifelse(is.finite(prediction_r2), sprintf("%.2f", prediction_r2), "NA"),
    "   |   RMSE = ",
    ifelse(is.finite(prediction_rmse), sprintf("%.2f", prediction_rmse), "NA")
  )
  
  ggplot(
    prediction_data,
    aes(x = observed, y = predicted_mean)
  ) +
    geom_linerange(
      aes(
        ymin = predictive_lower_95,
        ymax = predictive_upper_95,
        colour = "95% posterior predictive interval (y_rep)"
      ),
      linewidth = 0.9,
      alpha = 0.75
    ) +
    geom_linerange(
      aes(
        ymin = fitted_lower_95,
        ymax = fitted_upper_95,
        colour = "95% credible interval (mu)"
      ),
      linewidth = 1.1,
      alpha = 0.95
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.7
    ) +
    geom_point(
      colour = "black",
      fill = "steelblue",
      shape = 21,
      size = 2.5,
      alpha = 0.9
    ) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "95% posterior predictive interval (y_rep)" = "grey55",
        "95% credible interval (mu)" = "steelblue4"
      )
    ) +
    labs(
      title = response_label,
      subtitle = performance_label,
      x = "Observed values",
      y = expression(
        "Posterior predictive replicated values (" *
          italic(y)[rep] * ")"
      )
    ) +
    coord_equal() +
    theme_bayesian(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    )
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
    aes(x = x, y = y, label = equation)
  ) +
    annotate(
      "text",
      x = 0,
      y = 1.45,
      label = "Standardised posterior mean coefficients",
      hjust = 0.5,
      vjust = 1,
      fontface = "bold",
      size = 3.0,
      colour = "grey25"
    ) +
    geom_text(
      hjust = 0.5,
      vjust = 1,
      size = 3.1,
      lineheight = 1.2,
      colour = "grey20"
    ) +
    coord_cartesian(
      xlim = c(-1, 1),
      ylim = c(-0.1, 1.55),
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(5, 15, 5, 15)
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
    aes(x = relative_importance, y = predictor_label)
  ) +
    geom_col(
      fill = "grey45",
      colour = "black",
      linewidth = 0.25,
      width = 0.45
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
    theme_bayesian(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9, hjust = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5, 15, 5, 8)
    )
}


# =============================================================================
# 10. COMBINE REPLICATED-VS-OBSERVED PLOT WITH IMPORTANCE BARPLOT BENEATH
# =============================================================================

combine_model_plots <- function(
    prediction_plot,
    equation_plot,
    importance_plot,
    n_predictors
) {
  equation_height <- 0.85
  importance_height <- max(1.1, 0.22 * n_predictors)
  
  (prediction_plot / equation_plot / importance_plot) +
    patchwork::plot_layout(
      heights = c(3.8, equation_height, importance_height)
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

prediction_files <- sort(prediction_files)

prediction_files <- prediction_files[
  basename(prediction_files) != "predictions_all_responses_bestmod.csv"
]

if (length(prediction_files) == 0) {
  stop(
    paste0(
      "\nNo saved best-model prediction CSV files were found.\n\n",
      "Expected filenames such as:\n",
      "  predictions_C_Bulk_bestmod.csv\n",
      "  predictions_C_MAOM_bestmod.csv\n\n",
      "Directory searched:\n",
      tables_dir
    )
  )
}

message("\nFound ", length(prediction_files), " saved best-model prediction files.")
for (f in prediction_files) {
  message("  - ", basename(f))
}


# =============================================================================
# 12. PROCESS EACH SAVED BEST MODEL
# =============================================================================

plot_summary <- list()
all_importance <- list()
prediction_panel_plots <- list()

for (i in seq_along(prediction_files)) {
  prediction_file <- prediction_files[i]
  
  model_id <- sub("^predictions_", "", basename(prediction_file))
  model_id <- sub("\\.csv$", "", model_id)
  response <- identify_response_from_filename(prediction_file)
  
  coefficient_file <- file.path(
    tables_dir,
    paste0("coefficients_", model_id, ".csv")
  )
  
  message("\n--------------------------------------------------")
  message("Processing best model ", i, " of ", length(prediction_files))
  message("Model: ", model_id)
  message("Response: ", response)
  
  if (!file.exists(coefficient_file)) {
    warning("Matching coefficient file not found for ", model_id, ". Skipping.")
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
  
  prediction_data <- tryCatch(read_model_csv(prediction_file), error = function(e) NULL)
  coefficient_data <- tryCatch(read_model_csv(coefficient_file), error = function(e) NULL)
  
  if (is.null(prediction_data) || is.null(coefficient_data)) {
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
    warning("Prediction file for ", model_id, " is missing columns. Skipping.")
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
  
  importance_data <- tryCatch(
    calculate_relative_importance(coefficient_data),
    error = function(e) NULL
  )
  
  if (is.null(importance_data) || nrow(importance_data) == 0) {
    warning("No valid non-zero predictor coefficients for ", model_id, ". Skipping.")
    plot_summary[[model_id]] <- tibble(
      model = model_id,
      response = response,
      response_label = get_variable_label(response),
      status = "Skipped: no valid non-zero predictor coefficients",
      n_observations = nrow(prediction_data),
      n_predictors = 0L,
      prediction_file = basename(prediction_file),
      coefficient_file = basename(coefficient_file)
    )
    next
  }
  
  coefficient_equation <- tryCatch(
    build_coefficient_equation(
      coef_data = coefficient_data,
      predictors_per_line = 3,
      zero_tolerance = 1e-8
    ),
    error = function(e) "Coefficient equation unavailable"
  )
  
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
  
  readr::write_csv(
    importance_data,
    file.path(tables_dir, paste0("relative_importance_", model_id, ".csv"))
  )
  
  built_parts <- tryCatch(
    {
      pred_p <- plot_prediction(
        prediction_data = prediction_data,
        response_label = get_variable_label(response),
        model_id = model_id
      )
      
      eq_p <- plot_equation(coefficient_equation)
      
      imp_p <- plot_relative_importance(
        importance_data = importance_data,
        response = response
      )
      
      comb_p <- combine_model_plots(
        prediction_plot = pred_p,
        equation_plot = eq_p,
        importance_plot = imp_p,
        n_predictors = nrow(importance_data)
      )
      
      list(
        prediction_plot = pred_p,
        combined_plot = comb_p
      )
    },
    error = function(e) {
      warning("Plot generation failed for ", model_id, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(built_parts)) {
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
  
  # Store the full combined plot (prediction + equation + importance barplot)
  # so the grouped 2x2 panels include the barplots underneath each response.
  prediction_panel_plots[[response]] <- built_parts$combined_plot
  output_height <- max(7.2, 6.2 + 0.22 * nrow(importance_data))
  
  output_png <- file.path(
    plots_dir,
    paste0("observed_vs_predicted_", model_id, ".png")
  )
  output_pdf <- file.path(
    plots_dir,
    paste0("observed_vs_predicted_", model_id, ".pdf")
  )
  
  ggsave(output_png, built_parts$combined_plot, width = 10, height = output_height, dpi = 400, bg = "white")
  ggsave(output_pdf, built_parts$combined_plot, width = 10, height = output_height, bg = "white")
  
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
  
  message("Successfully saved: ", basename(output_png))
}


# =============================================================================
# 12B. SAVE FACETED 2 x 2 OBSERVED-VS-REPLICATED SUMMARY FIGURES
# =============================================================================

response_key <- function(x) tolower(gsub("[^[:alnum:]]", "", x))

save_prediction_panel <- function(responses, panel_name, output_dir) {
  if (length(prediction_panel_plots) == 0L || is.null(names(prediction_panel_plots))) {
    warning("No prediction plots available for ", panel_name, ".")
    return(invisible(NULL))
  }
  
  available_keys <- response_key(names(prediction_panel_plots))
  requested_keys <- response_key(responses)
  matched <- match(requested_keys, available_keys)
  matched <- matched[!is.na(matched)]
  available_names <- names(prediction_panel_plots)[matched]
  
  if (length(available_names) == 0L) {
    warning(
      "No prediction plots matched requested responses for ", panel_name,
      ". Available: ", paste(names(prediction_panel_plots), collapse = ", ")
    )
    return(invisible(NULL))
  }
  
  panel <- patchwork::wrap_plots(
    prediction_panel_plots[available_names],
    ncol = 2
  )
  
  ggplot2::ggsave(
    file.path(output_dir, paste0("observed_vs_predicted_panel_", panel_name, ".png")),
    panel, width = 18, height = 16, dpi = 400, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0("observed_vs_predicted_panel_", panel_name, ".pdf")),
    panel, width = 18, height = 16, bg = "white"
  )
  
  message(
    "Saved observed-vs-predicted panel ", panel_name, " with ",
    length(available_names), " response plot(s): ",
    paste(available_names, collapse = ", ")
  )
}

save_prediction_panel(c("C_Bulk", "C_fPOM", "C_oPOM", "C_MAOM"), "C_pools", plots_dir)
save_prediction_panel(
  c("d14C_2026_Bulk", "d14C_2026_fPOM", "d14C_2026_oPOM", "d14C_2026_MAOM"),
  "d14C_pools",
  plots_dir
)


# =============================================================================
# 13–15. SAVE COMBINED SUMMARIES
# =============================================================================

if (length(all_importance) > 0) {
  readr::write_csv(
    bind_rows(all_importance),
    file.path(tables_dir, "relative_importance_best_models_all_responses.csv")
  )
}

plot_summary_df <- if (length(plot_summary) > 0) {
  bind_rows(plot_summary)
} else {
  tibble(
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
  file.path(tables_dir, "best_model_plot_generation_summary.csv")
)

message("\n==================================================")
message("BEST-MODEL PLOT GENERATION COMPLETE")
message("Successfully plotted: ", sum(plot_summary_df$status == "Successfully plotted", na.rm = TRUE))
message("==================================================")


# =============================================================================
# 16. POSTERIOR RESIDUALS VS SELECTED PREDICTORS
# =============================================================================

residual_dir <- file.path(best_diag_dir, "posterior_residuals")
dir.create(residual_dir, recursive = TRUE, showWarnings = FALSE)

selected_files <- list.files(
  best_table_dir,
  pattern = "^selected_predictors_.*_bestmod\\.csv$",
  full.names = TRUE
)

selected_files <- selected_files[
  basename(selected_files) != "selected_predictors_all_responses_bestmod.csv"
]
selected_files <- sort(selected_files)

if (length(selected_files) > 0L && is.function(prepare_model_data_for_residuals)) {
  all_residual_summaries <- list()
  residual_plot_count <- 0L
  
  for (selected_file in selected_files) {
    selected_tbl <- readr::read_csv(selected_file, show_col_types = FALSE, progress = FALSE)
    if (nrow(selected_tbl) == 0L || !all(c("response", "predictor") %in% names(selected_tbl))) next
    
    response <- as.character(selected_tbl$response[1])
    if (is.na(response) || !nzchar(response) || response == "all_responses") next
    
    selected_predictors <- unique(as.character(selected_tbl$predictor))
    selected_predictors <- selected_predictors[!is.na(selected_predictors) & nzchar(selected_predictors)]
    if (length(selected_predictors) == 0L) next
    
    model_id <- safe_name_for_files(response)
    draws_file <- file.path(best_model_dir, paste0("draws_", model_id, "_bestmod.rds"))
    if (!file.exists(draws_file)) next
    
    draw_matrix <- posterior::as_draws_matrix(readRDS(draws_file))
    prep <- tryCatch(
      prepare_model_data_for_residuals(
        data = dat_for_residuals,
        response = response,
        candidates = candidate_predictors_for_residuals
      ),
      error = function(e) NULL
    )
    if (is.null(prep)) next
    if (length(setdiff(selected_predictors, prep$predictors)) > 0L) next
    
    mu_cols <- grep("^mu\\[[0-9]+\\]$", colnames(draw_matrix), value = TRUE)
    if (length(mu_cols) != prep$n) next
    mu_cols <- mu_cols[order(as.integer(sub("^mu\\[([0-9]+)\\]$", "\\1", mu_cols)))]
    mu_draws <- draw_matrix[, mu_cols, drop = FALSE]
    
    residual_draws <- sweep(
      mu_draws,
      MARGIN = 2,
      STATS = prep$data$y_std,
      FUN = function(mu, obs) obs - mu
    ) * prep$y_sd
    
    residual_summary <- tibble::tibble(
      response = response,
      sample = prep$data$sample,
      Treatment = as.character(prep$data$Treatment_model),
      residual_mean = colMeans(residual_draws, na.rm = TRUE),
      residual_median = apply(residual_draws, 2, stats::median, na.rm = TRUE),
      residual_lower_95 = apply(residual_draws, 2, stats::quantile, probs = 0.025, na.rm = TRUE, names = FALSE),
      residual_upper_95 = apply(residual_draws, 2, stats::quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
    )
    
    readr::write_csv(
      residual_summary,
      file.path(residual_dir, paste0("posterior_residuals_", model_id, "_bestmod.csv"))
    )
    
    residual_data <- dplyr::bind_cols(
      residual_summary,
      prep$data %>% dplyr::select(dplyr::all_of(selected_predictors))
    )
    
    response_label <- get_variable_label(response)
    
    for (predictor in selected_predictors) {
      plot_data <- residual_data %>%
        dplyr::transmute(
          sample = sample,
          Treatment = Treatment,
          predictor_value = .data[[predictor]],
          residual_mean = residual_mean,
          residual_lower_95 = residual_lower_95,
          residual_upper_95 = residual_upper_95
        ) %>%
        dplyr::filter(
          is.finite(predictor_value),
          is.finite(residual_mean),
          is.finite(residual_lower_95),
          is.finite(residual_upper_95)
        )
      
      if (nrow(plot_data) < 4L) next
      
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = predictor_value, y = residual_mean)) +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = residual_lower_95, ymax = residual_upper_95),
          colour = "grey65", linewidth = 0.7, alpha = 0.8
        ) +
        ggplot2::geom_point(
          ggplot2::aes(fill = Treatment),
          shape = 21, colour = "black", size = 2.7, alpha = 0.9
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
        {
          if (dplyr::n_distinct(plot_data$predictor_value) >= 4L) {
            ggplot2::geom_smooth(
              method = "loess", formula = y ~ x, se = TRUE,
              colour = "firebrick", fill = "firebrick", alpha = 0.12, linewidth = 0.9
            )
          } else {
            NULL
          }
        } +
        ggplot2::labs(
          title = paste(response_label, "- residuals vs", get_variable_label(predictor)),
          subtitle = "Posterior mean residuals; intervals reflect uncertainty in fitted mean (mu). Red line = LOESS smoother.",
          x = paste0(get_variable_label(predictor), " (standardised)"),
          y = paste0("Residual (observed - fitted; ", response, " units)"),
          fill = "Treatment"
        ) +
        theme_bayesian(base_size = 11) +
        ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
      
      ggplot2::ggsave(
        file.path(residual_dir, paste0("residual_", model_id, "_", safe_name_for_files(predictor), ".png")),
        p, width = 7, height = 5.5, dpi = 300, bg = "white"
      )
      residual_plot_count <- residual_plot_count + 1L
    }
    
    all_residual_summaries[[response]] <- residual_summary
  }
  
  if (length(all_residual_summaries) > 0L) {
    readr::write_csv(
      dplyr::bind_rows(all_residual_summaries),
      file.path(residual_dir, "posterior_residuals_all_responses.csv")
    )
  }
}
