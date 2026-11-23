
# =============================================================================
# MALLE ERW: Bayesian prediction of C and d14C using projpred
#
# PURPOSE
#   Fit regularised Bayesian reference models for:
#     C_Bulk, C_MAOM, C_oPOM, C_fPOM
#     d14C_2026_Bulk, d14C_2026_MAOM, d14C_2026_oPOM, d14C_2026_fPOM
#
#   Use projection predictive variable selection to identify the most
#   parsimonious model for each response using:
#
#     - Forward variable selection
#     - PSIS-LOO predictive performance evaluation
#     - One-standard-error parsimony criterion
#
#   Refit the selected predictors using CmdStan.
#
# OUTPUTS
#   - Full reference models and posterior draws
#   - projpred selection objects and performance plots
#   - Selected best-model predictor sets
#   - Best-model posterior coefficient summaries
#   - Best-model predictions
#   - Best-model diagnostics and posterior predictive checks
#   - Observed vs predicted plots with uncertainty ribbons
#   - Relative predictor importance bars underneath each plot
#
# IMPORTANT
#   - Predictors and responses are standardised for fitting.
#   - Missing predictors are median-imputed.
#   - Predictor collinearity is reduced before modelling.
#   - Best-model plots end in "_bestmod.png".
#   - Predictions are in-sample fitted predictions.
#   - projpred uses LOO predictive evaluation for selection.
#   - Relative importance is based on absolute standardised posterior
#     mean coefficients and is NOT variance explained.
# =============================================================================


suppressPackageStartupMessages({
  
  library(tidyverse)
  library(here)
  library(posterior)
  library(bayesplot)
  library(patchwork)
  library(cmdstanr)
  library(projpred)
  
})

set.seed(20261028)


# =============================================================================
# 0. SETTINGS AND PATHS
# =============================================================================

input_file <- here::here(
  "outputs",
  "integrated_geochemistry_for_bayesian.csv"
)

out_dir <- here::here(
  "outputs",
  "bayesian_C_d14C_predictions_withftir"
)

model_dir <- file.path(out_dir, "models")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
diag_dir <- file.path(out_dir, "diagnostics")

selection_dir <- file.path(out_dir, "projpred_selection")
selection_plot_dir <- file.path(selection_dir, "plots")
selection_table_dir <- file.path(selection_dir, "tables")

best_dir <- file.path(out_dir, "best_model")
best_model_dir <- file.path(best_dir, "models")
best_plot_dir <- file.path(best_dir, "plots")
best_table_dir <- file.path(best_dir, "tables")
best_diag_dir <- file.path(best_dir, "diagnostics")

all_dirs <- c(
  out_dir,
  model_dir,
  plot_dir,
  table_dir,
  diag_dir,
  selection_dir,
  selection_plot_dir,
  selection_table_dir,
  best_dir,
  best_model_dir,
  best_plot_dir,
  best_table_dir,
  best_diag_dir
)

for (d in all_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install cmdstanr before running this script.")
}

if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not configured. Check cmdstanr::cmdstan_path().")
}

if (!requireNamespace("projpred", quietly = TRUE)) {
  stop("Install projpred before running this script.")
}


message(
  "CmdStan version: ",
  as.character(cmdstanr::cmdstan_version())
)

message(
  "projpred version: ",
  as.character(utils::packageVersion("projpred"))
)


# MCMC settings

n_chains <- 4L
n_warmup <- 1500L
n_sampling <- 2000L

n_cores <- min(
  n_chains,
  parallel::detectCores(logical = FALSE)
)

adapt_delta <- 0.99
max_treedepth <- 12L


# Predictor filtering

min_predictor_observed_fraction <- 0.70
max_abs_predictor_correlation <- 0.95


# Regularisation

slope_prior_sd <- 0.5


# projpred settings

# Maximum number of terms explored during forward selection.
# The selected model can be substantially smaller.
max_submodel_terms <- 10L

# LOO evaluation settings.
# PSIS-LOO is used rather than refitting a complete reference model
# separately for every fold.

projpred_cv_method <- "LOO"

# Number of posterior draws/clusters used in projection.
projpred_nclusters <- 20L
projpred_ndraws <- 400L

# One-standard-error parsimony criterion.
parsimony_alpha <- 0.32
parsimony_pct <- 0
parsimony_type <- "upper"


# =============================================================================
# 1. READ AND FILTER DATA
# =============================================================================

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

integrated <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  na = c("", "NA", "#VALUE!")
)

if (!"sample" %in% names(integrated)) {
  stop("The integrated dataset must contain a sample column.")
}

if (!"Treatment" %in% names(integrated)) {
  stop("The integrated dataset must contain a Treatment column.")
}


normalise_treatment <- function(x) {
  
  z <- stringr::str_to_lower(
    stringr::str_squish(as.character(x))
  )
  
  z <- stringr::str_replace_all(z, "ü", "u")
  z <- stringr::str_replace_all(z, "[^a-z0-9]+", "")
  
  dplyr::case_when(
    stringr::str_detect(z, "^control$") ~ "Control",
    stringr::str_detect(z, "^lime") ~ "Lime",
    stringr::str_detect(z, "^bolsdorfer") ~ "Bolsdorfer",
    stringr::str_detect(z, "^huhnerberg") ~ "Huhnerberg",
    stringr::str_detect(z, "^eifelgold") ~ "Eifelgold",
    TRUE ~ NA_character_
  )
}


allowed_treatments <- c(
  "Control",
  "Lime_2",
  "Bolsdorfer_50",
  "Huhnerberg_50",
  "Eifelgold_50"
)


dat <- integrated %>%
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(
      as.integer(as.character(sample))
    ),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(
    !is.na(sample),
    !is.na(Treatment_model)
  )


if (nrow(dat) == 0L) {
  stop("No observations remain after treatment filtering.")
}


readr::write_csv(
  dat %>%
    count(Treatment, Treatment_model, name = "n_plots"),
  file.path(table_dir, "treatment_sample_counts.csv")
)


# =============================================================================
# 2. RESPONSE VARIABLES
# =============================================================================

c_responses <- c(
  "C_Bulk",
  "C_MAOM",
  "C_oPOM",
  "C_fPOM"
)

d14c_responses <- c(
  "d14C_2026_Bulk",
  "d14C_2026_MAOM",
  "d14C_2026_oPOM",
  "d14C_2026_fPOM"
)

expected_responses <- c(
  c_responses,
  d14c_responses
)

responses <- intersect(
  expected_responses,
  names(dat)
)

if (length(responses) == 0L) {
  stop("No expected response columns were found.")
}

missing_responses <- setdiff(
  expected_responses,
  names(dat)
)

if (length(missing_responses) > 0L) {
  message(
    "Missing responses (skipped): ",
    paste(missing_responses, collapse = ", ")
  )
}


# =============================================================================
# 3. CANDIDATE GEOCHEMICAL PREDICTORS
# =============================================================================

explicit_geochem_predictors <- c(
  
  # Aggregate stability: water
  
  "GMD_h2o",
  "D_h2o",
  "MWD_h2o",
  
  # Aggregate stability: oxalate / water ratios
  
  "GMD_oxalate_over_h2o",
  "D_oxalate_over_h2o",
  "MWD_oxalate_over_h2o",
  
  # Oxalate / DCB ratios
  
  "Fe_ox_over_DCB",
  "Al_ox_over_DCB",
  "Mn_ox_over_DCB",
  "Si_ox_over_DCB",
  
  # Mineral ratios
  
  "Fe_ox_over_Si_ox",
  "Al_ox_over_Si_ox",
  
  # Tessier exchangeable
  
  "Tessier_Acetic_Ca_mg_g",
  "Tessier_Acetic_K_mg_g",
  "Tessier_Acetic_Na_mg_g",
  "Tessier_Acetic_Mg_mg_g",
  
  # Tessier reducible
  
  "Tessier_Hydroxylamine_Fe_mg_g",
  "Tessier_Hydroxylamine_Al_mg_g",
  "Tessier_Hydroxylamine_Mn_mg_g",
  "Tessier_Hydroxylamine_Si_mg_g",
  
  # February 2026 soil chemistry
  
  "soil_pH",
  "soil_EC",
  
  # October 2023 soil chemistry
  
  "pH_Oct2023",
  "EC_Oct2023",
  "H_Oct2023",
  "Ca_Oct2023",
  "K_Oct2023",
  "Mg_Oct2023",
  "Na_Oct2023",
  "Al_Oct2023",
  "Fe_Oct2023",
  "Mn_Oct2023",
  "CEC_Oct2023",
  "Base_saturation_Oct2023"
)


numeric_cols <- names(dat)[
  vapply(dat, is.numeric, logical(1))
]


pattern_predictors <- numeric_cols[
  stringr::str_starts(numeric_cols, "Tessier_") |
    stringr::str_starts(numeric_cols, "CEC") |
    stringr::str_starts(numeric_cols, "Fe_ox") |
    stringr::str_starts(numeric_cols, "Al_ox") |
    stringr::str_starts(numeric_cols, "Mn_ox") |
    stringr::str_starts(numeric_cols, "Si_ox") |
    stringr::str_starts(numeric_cols, "FTIR_delta_")
]


candidate_predictors <- unique(c(
  intersect(explicit_geochem_predictors, names(dat)),
  intersect(pattern_predictors, names(dat))
))


# Exclude response leakage and metadata.

response_like <- names(dat)[
  stringr::str_detect(
    names(dat),
    "^C_(Bulk|MAOM|oPOM|fPOM)$"
  ) |
    stringr::str_detect(
      names(dat),
      "^d14C"
    )
]


metadata_cols <- c(
  "sample",
  "Treatment",
  "Treatment_model",
  "Basalt_type_exp",
  "Dose_resp_exp",
  "Date",
  "year",
  "depth",
  "fraction",
  "SOM_fraction"
)


candidate_predictors <- setdiff(
  candidate_predictors,
  c(response_like, metadata_cols)
)


candidate_predictors <- candidate_predictors[
  !stringr::str_detect(
    candidate_predictors,
    regex(
      "dataset|dose_resp_exp|basalt_type_exp|treatment|^sample$",
      ignore_case = TRUE
    )
  )
]


candidate_predictors <- candidate_predictors[
  vapply(
    candidate_predictors,
    function(v) is.numeric(dat[[v]]),
    logical(1)
  )
]


if (length(candidate_predictors) == 0L) {
  stop("No numeric geochemical predictors were detected.")
}


# Predictor availability report.

updated_predictors <- c(
  "Fe_ox_over_Si_ox",
  "Al_ox_over_Si_ox",
  "GMD_h2o",
  "D_h2o",
  "MWD_h2o",
  "GMD_oxalate_over_h2o",
  "D_oxalate_over_h2o",
  "MWD_oxalate_over_h2o"
)


readr::write_csv(
  tibble(
    predictor = updated_predictors,
    available_in_dataset = updated_predictors %in% names(dat),
    included_as_candidate = updated_predictors %in% candidate_predictors
  ),
  file.path(table_dir, "updated_predictor_availability_check.csv")
)


ftir_predictors <- candidate_predictors[
  stringr::str_starts(candidate_predictors, "FTIR_delta_")
]

message(
  "\nNumber of candidate predictors: ",
  length(candidate_predictors)
)

message(
  "Number of FTIR predictors: ",
  length(ftir_predictors)
)

print(candidate_predictors)


readr::write_csv(
  tibble(candidate_predictor = candidate_predictors),
  file.path(table_dir, "candidate_predictors.csv")
)


# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

safe_name <- function(x) {
  
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}


# -----------------------------------------------------------------------------
# 4A. Prepare response-specific data
# -----------------------------------------------------------------------------

prepare_model_data <- function(data, response, candidates) {
  
  y_raw <- suppressWarnings(
    as.numeric(data[[response]])
  )
  
  keep_y <- is.finite(y_raw)
  
  work <- data[keep_y, , drop = FALSE]
  y_raw <- y_raw[keep_y]
  
  if (length(y_raw) < 8L) {
    stop(
      response,
      ": fewer than 8 finite response observations."
    )
  }
  
  eligible <- candidates[
    vapply(candidates, function(v) {
      
      x <- suppressWarnings(
        as.numeric(work[[v]])
      )
      
      observed_fraction <- mean(is.finite(x))
      sx <- stats::sd(x, na.rm = TRUE)
      
      is.finite(observed_fraction) &&
        observed_fraction >= min_predictor_observed_fraction &&
        is.finite(sx) &&
        sx > 0
      
    }, logical(1))
  ]
  
  if (length(eligible) == 0L) {
    stop(response, ": no eligible predictors.")
  }
  
  
  # Median imputation and standardisation.
  
  X <- as.data.frame(
    lapply(eligible, function(v) {
      
      x <- suppressWarnings(
        as.numeric(work[[v]])
      )
      
      med <- stats::median(
        x[is.finite(x)],
        na.rm = TRUE
      )
      
      x[!is.finite(x)] <- med
      
      as.numeric(scale(x))
      
    })
  )
  
  names(X) <- eligible
  
  
  valid <- vapply(
    X,
    function(x) {
      all(is.finite(x)) &&
        is.finite(stats::sd(x)) &&
        stats::sd(x) > 0
    },
    logical(1)
  )
  
  X <- X[, valid, drop = FALSE]
  
  
  if (ncol(X) == 0L) {
    stop(response, ": no predictors remain after standardisation.")
  }
  
  
  # Iteratively remove one variable from highly correlated pairs.
  
  dropped_correlated <- character(0)
  
  repeat {
    
    if (ncol(X) < 2L) break
    
    cm <- suppressWarnings(
      stats::cor(X)
    )
    
    diag(cm) <- 0
    
    max_cor <- max(abs(cm), na.rm = TRUE)
    
    if (
      !is.finite(max_cor) ||
      max_cor <= max_abs_predictor_correlation
    ) {
      break
    }
    
    idx <- which(
      abs(cm) == max_cor,
      arr.ind = TRUE
    )[1, ]
    
    drop_name <- colnames(cm)[idx[2]]
    
    dropped_correlated <- c(
      dropped_correlated,
      drop_name
    )
    
    X <- X[
      ,
      setdiff(names(X), drop_name),
      drop = FALSE
    ]
  }
  
  
  # Standardise response.
  
  y_mean <- mean(y_raw)
  y_sd <- stats::sd(y_raw)
  
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop(response, ": response has no variation.")
  }
  
  y_std <- (y_raw - y_mean) / y_sd
  
  
  model_data <- as_tibble(X)
  
  model_data$y_std <- y_std
  model_data$sample <- work$sample
  model_data$Treatment_model <- factor(
    work$Treatment_model
  )
  
  
  list(
    data = model_data,
    y_raw = y_raw,
    y_mean = y_mean,
    y_sd = y_sd,
    predictors = names(X),
    dropped_correlated = dropped_correlated,
    n = nrow(model_data)
  )
}


# -----------------------------------------------------------------------------
# 4B. Extract predictions
# -----------------------------------------------------------------------------

make_predictions <- function(
    draw_matrix,
    model_data,
    prep,
    response) {
  
  mu_cols <- grep(
    "^mu\\[",
    colnames(draw_matrix),
    value = TRUE
  )
  
  rep_cols <- grep(
    "^y_rep\\[",
    colnames(draw_matrix),
    value = TRUE
  )
  
  if (
    length(mu_cols) != nrow(model_data) ||
    length(rep_cols) != nrow(model_data)
  ) {
    stop("Could not extract all mu/y_rep draws for ", response)
  }
  
  
  mu_draws <- draw_matrix[
    ,
    mu_cols,
    drop = FALSE
  ]
  
  rep_draws <- draw_matrix[
    ,
    rep_cols,
    drop = FALSE
  ]
  
  
  pred_mean <- apply(mu_draws, 2, mean)
  pred_median <- apply(mu_draws, 2, median)
  
  pred_lower <- apply(
    mu_draws,
    2,
    stats::quantile,
    probs = 0.025
  )
  
  pred_upper <- apply(
    mu_draws,
    2,
    stats::quantile,
    probs = 0.975
  )
  
  pred_lower_pi <- apply(
    rep_draws,
    2,
    stats::quantile,
    probs = 0.025
  )
  
  pred_upper_pi <- apply(
    rep_draws,
    2,
    stats::quantile,
    probs = 0.975
  )
  
  
  tibble(
    response = response,
    sample = model_data$sample,
    Treatment_model = as.character(
      model_data$Treatment_model
    ),
    observed = prep$y_raw,
    predicted_mean =
      pred_mean * prep$y_sd + prep$y_mean,
    predicted_median =
      pred_median * prep$y_sd + prep$y_mean,
    fitted_lower_95 =
      pred_lower * prep$y_sd + prep$y_mean,
    fitted_upper_95 =
      pred_upper * prep$y_sd + prep$y_mean,
    predictive_lower_95 =
      pred_lower_pi * prep$y_sd + prep$y_mean,
    predictive_upper_95 =
      pred_upper_pi * prep$y_sd + prep$y_mean
  )
}


# -----------------------------------------------------------------------------
# 4C. Coefficient summaries
# -----------------------------------------------------------------------------

summarise_coefficients <- function(
    draw_matrix,
    predictors,
    response,
    prep,
    model_label) {
  
  coef_names <- c(
    "alpha",
    paste0(
      "beta[",
      seq_along(predictors),
      "]"
    )
  )
  
  coef_labels <- c(
    "Intercept",
    predictors
  )
  
  
  purrr::map_dfr(
    seq_along(coef_names),
    function(j) {
      
      vals <- draw_matrix[
        ,
        coef_names[j]
      ]
      
      tibble(
        term = coef_labels[j],
        Estimate = mean(vals),
        Est.Error = stats::sd(vals),
        Q2.5 = unname(
          stats::quantile(vals, 0.025)
        ),
        Q50 = unname(
          stats::quantile(vals, 0.50)
        ),
        Q97.5 = unname(
          stats::quantile(vals, 0.975)
        ),
        prob_positive = mean(vals > 0),
        prob_negative = mean(vals < 0)
      )
      
    }
  ) %>%
    mutate(
      response = response,
      model = model_label,
      n_predictors = length(predictors),
      predictors = paste(
        predictors,
        collapse = " + "
      ),
      estimate_original_response_per_predictor_SD =
        Estimate * prep$y_sd,
      lower_original_response_per_predictor_SD =
        Q2.5 * prep$y_sd,
      upper_original_response_per_predictor_SD =
        Q97.5 * prep$y_sd
    )
}


# -----------------------------------------------------------------------------
# 4D. Relative predictor importance
# -----------------------------------------------------------------------------

calculate_relative_importance <- function(
    coef_tbl,
    response) {
  
  importance_tbl <- coef_tbl %>%
    filter(
      term != "Intercept",
      is.finite(Estimate)
    ) %>%
    transmute(
      Predictor = as.character(term),
      Standardised_coefficient = Estimate
    )
  
  if (nrow(importance_tbl) == 0L) {
    stop("No predictor coefficients available for ", response)
  }
  
  total_absolute <- sum(
    abs(importance_tbl$Standardised_coefficient)
  )
  
  if (
    !is.finite(total_absolute) ||
    total_absolute == 0
  ) {
    stop("Invalid total coefficient magnitude for ", response)
  }
  
  importance_tbl %>%
    mutate(
      Relative_importance =
        abs(Standardised_coefficient) /
        total_absolute * 100
    ) %>%
    arrange(Relative_importance) %>%
    mutate(
      Predictor = factor(
        Predictor,
        levels = Predictor
      )
    )
}


# -----------------------------------------------------------------------------
# 4E. Combined observed-predicted and importance plot
# -----------------------------------------------------------------------------

plot_observed_predicted <- function(
    pred,
    coef_tbl,
    response,
    output_path,
    model_label) {
  
  pdat <- pred %>%
    arrange(observed) %>%
    mutate(row_order = row_number())
  
  
  # Upper panel.
  
  p_obs <- ggplot(
    pdat,
    aes(
      x = observed,
      y = predicted_mean
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = fitted_lower_95,
        ymax = fitted_upper_95
      ),
      fill = "grey70",
      alpha = 0.55
    ) +
    
    geom_errorbar(
      aes(
        ymin = predictive_lower_95,
        ymax = predictive_upper_95
      ),
      width = 0,
      alpha = 0.35,
      colour = "grey40"
    ) +
    
    geom_point(
      aes(shape = Treatment_model),
      size = 2.8,
      alpha = 0.9
    ) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.7,
      colour = "black"
    ) +
    
    coord_equal() +
    
    labs(
      title = paste(
        "Observed vs Bayesian predicted:",
        response
      ),
      subtitle = paste0(
        model_label,
        " | Shaded band = 95% credible interval for fitted mean; ",
        "whiskers = 95% posterior predictive interval"
      ),
      x = paste("Observed", response),
      y = paste("Predicted", response),
      shape = "Treatment"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  
  # Lower panel.
  
  importance_tbl <- calculate_relative_importance(
    coef_tbl,
    response
  )
  
  n_predictors_plot <- nrow(importance_tbl)
  
  p_importance <- ggplot(
    importance_tbl,
    aes(
      x = Relative_importance,
      y = forcats::fct_reorder(
        Predictor,
        Relative_importance
      )
    )
  ) +
    
    geom_col(
      fill = "#2878B5",
      width = 0.72,
      alpha = 0.9
    ) +
    
    geom_text(
      aes(
        label = sprintf(
          "%.1f%%",
          Relative_importance
        )
      ),
      hjust = -0.15,
      size = 3
    ) +
    
    scale_x_continuous(
      limits = c(
        0,
        max(importance_tbl$Relative_importance) * 1.20
      ),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    labs(
      title = "Relative importance of predictors",
      subtitle = paste0(
        "Absolute standardised posterior mean coefficient ",
        "as percentage of total absolute coefficient magnitude"
      ),
      x = "Relative importance (%)",
      y = NULL
    ) +
    
    theme_bw(base_size = 10) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      plot.subtitle = element_text(size = 8.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(
        t = 5,
        r = 25,
        b = 5,
        l = 5
      )
    )
  
  
  # Combine panels.
  
  combined_plot <- p_obs / p_importance +
    
    patchwork::plot_layout(
      heights = c(
        3,
        max(1.1, 0.20 * n_predictors_plot)
      )
    ) +
    
    patchwork::plot_annotation(
      caption = paste0(
        "Relative importance is calculated from absolute posterior ",
        "mean standardised coefficients. Percentages sum to 100% ",
        "within each model and do not represent variance explained."
      ),
      theme = theme(
        plot.caption = element_text(
          size = 8,
          colour = "grey35",
          hjust = 0
        )
      )
    )
  
  
  output_height <- max(
    8,
    6.2 + 0.20 * n_predictors_plot
  )
  
  
  ggsave(
    filename = output_path,
    plot = combined_plot,
    width = 8.5,
    height = output_height,
    dpi = 320,
    limitsize = FALSE
  )
  
  combined_plot
}


# -----------------------------------------------------------------------------
# 4F. In-sample model performance
# -----------------------------------------------------------------------------

calculate_metrics <- function(
    pred,
    response,
    model_label,
    predictors) {
  
  rmse <- sqrt(
    mean(
      (pred$observed - pred$predicted_mean)^2,
      na.rm = TRUE
    )
  )
  
  mae <- mean(
    abs(
      pred$observed - pred$predicted_mean
    ),
    na.rm = TRUE
  )
  
  r2 <- suppressWarnings(
    stats::cor(
      pred$observed,
      pred$predicted_mean,
      use = "complete.obs"
    )^2
  )
  
  coverage_mean <- mean(
    pred$observed >= pred$fitted_lower_95 &
      pred$observed <= pred$fitted_upper_95,
    na.rm = TRUE
  )
  
  coverage_predictive <- mean(
    pred$observed >= pred$predictive_lower_95 &
      pred$observed <= pred$predictive_upper_95,
    na.rm = TRUE
  )
  
  tibble(
    response = response,
    model = model_label,
    n = nrow(pred),
    n_predictors = length(predictors),
    predictors = paste(
      predictors,
      collapse = " + "
    ),
    RMSE_in_sample = rmse,
    MAE_in_sample = mae,
    squared_correlation_in_sample = r2,
    observed_coverage_fitted_mean_95 = coverage_mean,
    observed_coverage_predictive_95 = coverage_predictive
  )
}


# -----------------------------------------------------------------------------
# 4G. CmdStan fitting function
# -----------------------------------------------------------------------------

fit_stan_model <- function(
    stan_model,
    model_data,
    predictors,
    model_id,
    seed_offset) {
  
  X <- as.matrix(
    model_data[
      ,
      predictors,
      drop = FALSE
    ]
  )
  
  storage.mode(X) <- "double"
  
  stan_data <- list(
    N = nrow(X),
    K = ncol(X),
    X = X,
    y = as.numeric(model_data$y_std)
  )
  
  csv_dir <- file.path(
    model_dir,
    paste0("cmdstan_csv_", model_id)
  )
  
  dir.create(
    csv_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  message(
    "\nSampling model: ",
    model_id,
    " | n = ",
    stan_data$N,
    " | K = ",
    stan_data$K
  )
  
  fit <- stan_model$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_cores,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = 20261028 + seed_offset,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 200,
    output_dir = csv_dir
  )
  
  saveRDS(
    fit,
    file.path(
      model_dir,
      paste0("fit_", model_id, ".rds")
    )
  )
  
  draws_array <- fit$draws(
    format = "draws_array"
  )
  
  saveRDS(
    draws_array,
    file.path(
      model_dir,
      paste0("draws_", model_id, ".rds")
    )
  )
  
  list(
    fit = fit,
    draws_array = draws_array,
    draw_matrix = posterior::as_draws_matrix(draws_array)
  )
}


# =============================================================================
# 5. STAN MODEL
# =============================================================================

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, SLOPE_PRIOR_SD);
  sigma ~ exponential(1);

  y ~ normal(alpha + X * beta, sigma);
}

generated quantities {

  vector[N] mu;
  vector[N] y_rep;

  mu = alpha + X * beta;

  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
}
"


stan_code <- sub(
  "SLOPE_PRIOR_SD",
  format(slope_prior_sd, scientific = FALSE),
  stan_code,
  fixed = TRUE
)


stan_file <- file.path(
  model_dir,
  "regularised_multiple_regression.stan"
)

writeLines(
  stan_code,
  stan_file,
  useBytes = TRUE
)


stan_model <- tryCatch(
  
  cmdstanr::cmdstan_model(
    stan_file,
    quiet = TRUE
  ),
  
  error = function(e) {
    stop(
      "Stan compilation failed: ",
      conditionMessage(e)
    )
  }
  
)


# =============================================================================
# 6. STORAGE OBJECTS
# =============================================================================

all_predictions <- list()
all_metrics <- list()
all_coefficients <- list()
all_model_info <- list()
all_draws <- list()

all_best_predictions <- list()
all_best_metrics <- list()
all_best_coefficients <- list()
all_best_info <- list()
all_best_diagnostics <- list()

all_selection_info <- list()


# =============================================================================
# 7. FIT AND SELECT MODELS FOR EACH RESPONSE
# =============================================================================

for (response_index in seq_along(responses)) {
  
  response <- responses[response_index]
  
  message("\n\n============================================================")
  message("RESPONSE: ", response)
  message("============================================================")
  
  
  model_id <- safe_name(response)
  
  
  # ---------------------------------------------------------------------------
  # 7A. Prepare data
  # ---------------------------------------------------------------------------
  
  prep <- tryCatch(
    
    prepare_model_data(
      dat,
      response,
      candidate_predictors
    ),
    
    error = function(e) {
      
      message(
        "Skipping ",
        response,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(prep)) next
  
  
  message("Observations: ", prep$n)
  message("Available predictors: ", length(prep$predictors))
  
  if (length(prep$dropped_correlated) > 0L) {
    message(
      "Dropped for high collinearity: ",
      paste(prep$dropped_correlated, collapse = ", ")
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 7B. Fit or reuse full reference model
  # ---------------------------------------------------------------------------
  
  draws_file <- file.path(
    model_dir,
    paste0("draws_", model_id, ".rds")
  )
  
  if (file.exists(draws_file)) {
    
    message(
      "Reusing existing reference model posterior draws: ",
      response
    )
    
    draws_array <- readRDS(draws_file)
    
    draw_matrix <- posterior::as_draws_matrix(
      draws_array
    )
    
    fit <- NULL
    
  } else {
    
    ref_result <- tryCatch(
      
      fit_stan_model(
        stan_model = stan_model,
        model_data = prep$data,
        predictors = prep$predictors,
        model_id = model_id,
        seed_offset = response_index
      ),
      
      error = function(e) {
        
        message(
          "Reference model failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (is.null(ref_result)) next
    
    fit <- ref_result$fit
    draws_array <- ref_result$draws_array
    draw_matrix <- ref_result$draw_matrix
    
  }
  
  
  all_draws[[response]] <- draws_array
  
  
  # ---------------------------------------------------------------------------
  # 7C. Create custom projpred reference model
  # ---------------------------------------------------------------------------
  
  # This function returns the posterior expected response for each
  # reference-model draw and observation.
  
  ref_fit <- list(
    draw_matrix = as.matrix(draw_matrix),
    predictors = prep$predictors,
    data = as.data.frame(prep$data)
  )
  
  
  ref_predfun <- function(fit, newdata = NULL) {
    
    if (is.null(newdata)) {
      newdata <- fit$data
    }
    
    X_new <- as.matrix(
      newdata[
        ,
        fit$predictors,
        drop = FALSE
      ]
    )
    
    storage.mode(X_new) <- "double"
    
    beta_names <- paste0(
      "beta[",
      seq_along(fit$predictors),
      "]"
    )
    
    beta_draws <- as.matrix(
      fit$draw_matrix[
        ,
        beta_names,
        drop = FALSE
      ]
    )
    
    alpha_draws <- as.numeric(
      fit$draw_matrix[, "alpha"]
    )
    
    eta <- X_new %*% t(beta_draws)
    
    sweep(
      eta,
      MARGIN = 2,
      STATS = alpha_draws,
      FUN = "+"
    )
  }
  
  
  full_formula <- stats::as.formula(
    paste(
      "y_std ~",
      paste(
        prep$predictors,
        collapse = " + "
      )
    )
  )
  
  
  ref_obj <- tryCatch(
    
    projpred::init_refmodel(
      object = ref_fit,
      data = as.data.frame(prep$data),
      formula = full_formula,
      family = stats::gaussian(),
      ref_predfun = ref_predfun,
      dis = as.numeric(draw_matrix[, "sigma"])
    ),
    
    error = function(e) {
      
      message(
        "Could not initialise projpred reference model: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(ref_obj)) next
  
  
  # ---------------------------------------------------------------------------
  # 7D. Run cross-validated projection predictive selection
  # ---------------------------------------------------------------------------
  
  max_terms <- min(
    max_submodel_terms,
    length(prep$predictors)
  )
  
  
  message(
    "\nRunning projpred forward selection with LOO-CV."
  )
  
  message(
    "Maximum explored submodel size: ",
    max_terms
  )
  
  
  vsel <- tryCatch(
    
    projpred::cv_varsel(
      ref_obj,
      method = "forward",
      cv_method = projpred_cv_method,
      nterms_max = max_terms,
      nclusters = min(
        projpred_nclusters,
        nrow(draw_matrix)
      ),
      nclusters_pred = min(
        projpred_nclusters,
        nrow(draw_matrix)
      ),
      seed = 20261028 + response_index,
      verbose = 1
    ),
    
    error = function(e) {
      
      message(
        "projpred CV failed for ",
        response,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (is.null(vsel)) {
    
    message(
      "No best model can be selected for ",
      response,
      ". Continuing to next response."
    )
    
    next
  }
  
  
  # Save selection object.
  
  saveRDS(
    vsel,
    file.path(
      selection_dir,
      paste0("cv_projpred_", model_id, ".rds")
    )
  )
  
  
  # Save selection summary.
  
  capture.output(
    print(summary(vsel)),
    file = file.path(
      selection_table_dir,
      paste0("cv_summary_", model_id, ".txt")
    )
  )
  
  
  # Save predictor ranking.
  
  ranking_result <- projpred::ranking(vsel)
  
  capture.output(
    print(ranking_result),
    file = file.path(
      selection_table_dir,
      paste0("ranking_", model_id, ".txt")
    )
  )
  
  
  # Save CV predictor inclusion proportions.
  
  cv_prop <- tryCatch(
    projpred::cv_proportions(vsel),
    error = function(e) {
      message(
        "Could not calculate CV proportions: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (!is.null(cv_prop)) {
    
    capture.output(
      print(cv_prop),
      file = file.path(
        selection_table_dir,
        paste0("cv_proportions_", model_id, ".txt")
      )
    )
    
    saveRDS(
      cv_prop,
      file.path(
        selection_dir,
        paste0("cv_proportions_", model_id, ".rds")
      )
    )
  }
  
  
  # Plot predictive performance.
  
  p_vsel <- plot(
    vsel,
    stats = "elpd",
    deltas = TRUE
  ) +
    
    ggplot2::ggtitle(
      paste(
        "Cross-validated projpred performance:",
        response
      )
    ) +
    
    ggplot2::theme_bw(base_size = 12)
  
  
  ggplot2::ggsave(
    filename = file.path(
      selection_plot_dir,
      paste0("projpred_cv_performance_", model_id, ".png")
    ),
    plot = p_vsel,
    width = 7.5,
    height = 5.5,
    dpi = 300
  )
  
  
  # ---------------------------------------------------------------------------
  # 7E. Select the most parsimonious model using the 1-SE criterion
  # ---------------------------------------------------------------------------
  
  best_size <- tryCatch(
    
    projpred::suggest_size(
      vsel,
      stat = "elpd",
      alpha = parsimony_alpha,
      pct = parsimony_pct,
      type = parsimony_type,
      baseline = "ref"
    ),
    
    error = function(e) {
      
      message(
        "Could not determine parsimonious model size: ",
        conditionMessage(e)
      )
      
      NA_real_
    }
  )
  
  
  if (
    length(best_size) != 1L ||
    !is.finite(best_size) ||
    is.na(best_size)
  ) {
    
    message(
      "The 1-SE criterion did not identify a model size for ",
      response,
      ". No arbitrary model size will be substituted."
    )
    
    next
  }
  
  
  best_size <- as.integer(best_size)
  
  
  message("\n--------------------------------------------")
  message("PARSIMONIOUS MODEL SELECTED")
  message("--------------------------------------------")
  message("Response: ", response)
  message("Selected number of predictors: ", best_size)
  
  
  # A zero-predictor result represents an intercept-only model.
  # It is not silently replaced with a larger model.
  
  if (best_size < 1L) {
    
    message(
      "The criterion selected the intercept-only model for ",
      response,
      ". No predictor importance plot will be generated."
    )
    
    all_selection_info[[response]] <- tibble(
      response = response,
      selected_n_predictors = 0L,
      selected_predictors = "",
      selection_status = "Intercept-only model"
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # 7F. Project the reference model onto the selected submodel
  # ---------------------------------------------------------------------------
  
  best_projection <- tryCatch(
    
    projpred::project(
      vsel,
      nterms = best_size,
      ndraws = min(
        projpred_ndraws,
        nrow(draw_matrix)
      ),
      seed = 20261028 + 100L * response_index,
      verbose = 0
    ),
    
    error = function(e) {
      
      message(
        "Projection failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (is.null(best_projection)) next
  
  
  selected_predictors <- as.character(
    projpred::predictor_terms(best_projection)
  )
  
  selected_predictors <- intersect(
    selected_predictors,
    prep$predictors
  )
  
  
  if (length(selected_predictors) != best_size) {
    
    stop(
      "Expected ",
      best_size,
      " selected predictors but obtained ",
      length(selected_predictors),
      " for ",
      response,
      ". Check projpred::predictor_terms()."
    )
  }
  
  
  message(
    "\nSelected predictors:\n",
    paste(
      selected_predictors,
      collapse = "\n"
    )
  )
  
  
  # Save projection.
  
  saveRDS(
    best_projection,
    file.path(
      best_model_dir,
      paste0("projection_", model_id, "_bestmod.rds")
    )
  )
  
  
  # Save selected predictors.
  
  selected_terms_tbl <- tibble(
    response = response,
    selected_model_size = best_size,
    selection_order = seq_along(selected_predictors),
    predictor = selected_predictors,
    selection_criterion = "One-standard-error ELPD",
    alpha = parsimony_alpha,
    pct = parsimony_pct,
    type = parsimony_type,
    cv_method = projpred_cv_method
  )
  
  
  readr::write_csv(
    selected_terms_tbl,
    file.path(
      best_table_dir,
      paste0("selected_predictors_", model_id, "_bestmod.csv")
    )
  )
  
  
  all_selection_info[[response]] <- tibble(
    response = response,
    selected_n_predictors = best_size,
    selected_predictors = paste(
      selected_predictors,
      collapse = " + "
    ),
    selection_status = "Successfully selected"
  )
  
  
  # ---------------------------------------------------------------------------
  # 7G. Refit the selected model using CmdStan
  # ---------------------------------------------------------------------------
  
  best_model_id <- paste0(
    model_id,
    "_bestmod"
  )
  
  
  best_fit_result <- tryCatch(
    
    fit_stan_model(
      stan_model = stan_model,
      model_data = prep$data,
      predictors = selected_predictors,
      model_id = best_model_id,
      seed_offset = 10000L + response_index
    ),
    
    error = function(e) {
      
      message(
        "Best-model refitting failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (is.null(best_fit_result)) next
  
  
  best_fit <- best_fit_result$fit
  
  best_draws_array <- best_fit_result$draws_array
  
  best_draw_matrix <- best_fit_result$draw_matrix
  
  
  # Save a dedicated best-model posterior draw object.
  
  saveRDS(
    best_draws_array,
    file.path(
      best_model_dir,
      paste0("draws_", model_id, "_bestmod.rds")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7H. Calculate best-model posterior predictions
  # ---------------------------------------------------------------------------
  
  best_pred <- make_predictions(
    draw_matrix = best_draw_matrix,
    model_data = prep$data,
    prep = prep,
    response = response
  ) %>%
    
    mutate(
      model = "projpred one-SE best model",
      n_predictors = best_size,
      predictors = paste(
        selected_predictors,
        collapse = " + "
      ),
      .before = 1
    )
  
  
  all_best_predictions[[response]] <- best_pred
  
  
  readr::write_csv(
    best_pred,
    file.path(
      best_table_dir,
      paste0("predictions_", model_id, "_bestmod.csv")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7I. Calculate best-model coefficient summaries
  # ---------------------------------------------------------------------------
  
  best_coef <- summarise_coefficients(
    draw_matrix = best_draw_matrix,
    predictors = selected_predictors,
    response = response,
    prep = prep,
    model_label = "projpred one-SE best model"
  )
  
  
  all_best_coefficients[[response]] <- best_coef
  
  
  readr::write_csv(
    best_coef,
    file.path(
      best_table_dir,
      paste0("coefficients_", model_id, "_bestmod.csv")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7J. Generate best-model plot with predictor importance bars
  # ---------------------------------------------------------------------------
  
  best_plot <- plot_observed_predicted(
    pred = best_pred,
    coef_tbl = best_coef,
    response = response,
    output_path = file.path(
      best_plot_dir,
      paste0(
        "observed_vs_predicted_",
        model_id,
        "_bestmod.png"
      )
    ),
    model_label = paste0(
      "projpred one-SE selected model | ",
      best_size,
      " predictors"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7K. Best-model diagnostics
  # ---------------------------------------------------------------------------
  
  best_summary <- best_fit$summary(
    variables = c(
      "alpha",
      "beta",
      "sigma"
    )
  )
  
  
  best_sampler <- as.data.frame(
    best_fit$sampler_diagnostics(
      format = "draws_df"
    )
  )
  
  
  best_diag <- tibble(
    
    response = response,
    
    model = "projpred one-SE best model",
    
    n = prep$n,
    
    n_predictors = best_size,
    
    predictors = paste(
      selected_predictors,
      collapse = " + "
    ),
    
    max_Rhat = if (
      all(is.na(best_summary$rhat))
    ) {
      NA_real_
    } else {
      max(best_summary$rhat, na.rm = TRUE)
    },
    
    min_bulk_ESS = if (
      all(is.na(best_summary$ess_bulk))
    ) {
      NA_real_
    } else {
      min(best_summary$ess_bulk, na.rm = TRUE)
    },
    
    min_tail_ESS = if (
      all(is.na(best_summary$ess_tail))
    ) {
      NA_real_
    } else {
      min(best_summary$ess_tail, na.rm = TRUE)
    },
    
    divergences = sum(
      best_sampler$divergent__ == 1,
      na.rm = TRUE
    ),
    
    treedepth_hits = sum(
      best_sampler$treedepth__ >= max_treedepth,
      na.rm = TRUE
    )
  )
  
  
  all_best_diagnostics[[response]] <- best_diag
  
  
  readr::write_csv(
    best_diag,
    file.path(
      best_diag_dir,
      paste0("diagnostics_", model_id, "_bestmod.csv")
    )
  )
  
  
  capture.output(
    best_fit$summary(),
    file = file.path(
      best_diag_dir,
      paste0("summary_", model_id, "_bestmod.txt")
    )
  )
  
  
  # Posterior predictive check.
  
  best_rep_cols <- grep(
    "^y_rep\\[",
    colnames(best_draw_matrix),
    value = TRUE
  )
  
  
  best_yrep <- best_draw_matrix[
    seq_len(
      min(100L, nrow(best_draw_matrix))
    ),
    best_rep_cols,
    drop = FALSE
  ]
  
  
  best_pp <- bayesplot::ppc_dens_overlay(
    y = as.numeric(prep$data$y_std),
    yrep = best_yrep
  ) +
    
    ggplot2::ggtitle(
      paste(
        "Best-model posterior predictive check:",
        response
      )
    ) +
    
    ggplot2::theme_bw(base_size = 12)
  
  
  ggplot2::ggsave(
    filename = file.path(
      best_diag_dir,
      paste0("pp_check_", model_id, "_bestmod.png")
    ),
    plot = best_pp,
    width = 7,
    height = 5.5,
    dpi = 300
  )
  
  
  # ---------------------------------------------------------------------------
  # 7L. Best-model in-sample performance
  # ---------------------------------------------------------------------------
  
  best_metrics <- calculate_metrics(
    pred = best_pred,
    response = response,
    model_label = "projpred one-SE best model",
    predictors = selected_predictors
  )
  
  
  all_best_metrics[[response]] <- best_metrics
  
  
  readr::write_csv(
    best_metrics,
    file.path(
      best_table_dir,
      paste0("performance_", model_id, "_bestmod.csv")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 7M. Full reference model predictions and coefficients
  # ---------------------------------------------------------------------------
  
  # These are retained for comparison with the selected best model.
  
  ref_coef <- summarise_coefficients(
    draw_matrix = draw_matrix,
    predictors = prep$predictors,
    response = response,
    prep = prep,
    model_label = "full reference model"
  )
  
  
  ref_pred <- make_predictions(
    draw_matrix = draw_matrix,
    model_data = prep$data,
    prep = prep,
    response = response
  ) %>%
    
    mutate(
      model = "full reference model",
      n_predictors = length(prep$predictors),
      predictors = paste(
        prep$predictors,
        collapse = " + "
      ),
      .before = 1
    )
  
  
  all_predictions[[response]] <- ref_pred
  
  all_coefficients[[response]] <- ref_coef
  
  
  all_metrics[[response]] <- calculate_metrics(
    pred = ref_pred,
    response = response,
    model_label = "full reference model",
    predictors = prep$predictors
  )
  
  
  readr::write_csv(
    ref_pred,
    file.path(
      table_dir,
      paste0("predictions_", model_id, ".csv")
    )
  )
  
  
  readr::write_csv(
    ref_coef,
    file.path(
      table_dir,
      paste0("coefficients_", model_id, ".csv")
    )
  )
  
  
  # Full-model plot.
  
  plot_observed_predicted(
    pred = ref_pred,
    coef_tbl = ref_coef,
    response = response,
    output_path = file.path(
      plot_dir,
      paste0(
        "observed_vs_predicted_",
        model_id,
        ".png"
      )
    ),
    model_label = "full reference model"
  )
  
  
  # Full-model performance metrics.
  
  all_model_info[[response]] <- tibble(
    response = response,
    model = "full reference model",
    n = prep$n,
    n_predictors = length(prep$predictors),
    predictors = paste(
      prep$predictors,
      collapse = " + "
    ),
    dropped_correlated = paste(
      prep$dropped_correlated,
      collapse = " + "
    )
  )
  
  
  all_best_info[[response]] <- tibble(
    response = response,
    model = "projpred one-SE best model",
    n = prep$n,
    n_predictors = best_size,
    predictors = paste(
      selected_predictors,
      collapse = " + "
    ),
    dropped_correlated = paste(
      prep$dropped_correlated,
      collapse = " + "
    )
  )
  
  
  message("\nCompleted response: ", response)
  
  message(
    "Selected best model contains ",
    best_size,
    " predictors."
  )
  
  message(
    "Best-model plot saved: ",
    file.path(
      best_plot_dir,
      paste0(
        "observed_vs_predicted_",
        model_id,
        "_bestmod.png"
      )
    )
  )
  
}


# =============================================================================
# 8. SAVE COMBINED OUTPUTS
# =============================================================================

if (length(all_best_predictions) == 0L) {
  
  stop(
    "No parsimonious models were successfully fitted. ",
    "Check the projpred selection results and model diagnostics."
  )
}


# -----------------------------------------------------------------------------
# 8A. Full reference model outputs
# -----------------------------------------------------------------------------

predictions_all <- bind_rows(all_predictions)

metrics_all <- bind_rows(all_metrics)

coefficients_all <- bind_rows(all_coefficients)

model_info_all <- bind_rows(all_model_info)


readr::write_csv(
  predictions_all,
  file.path(
    table_dir,
    "predictions_all_responses.csv"
  )
)


readr::write_csv(
  metrics_all,
  file.path(
    table_dir,
    "model_performance_in_sample.csv"
  )
)


readr::write_csv(
  coefficients_all,
  file.path(
    table_dir,
    "coefficients_all_responses.csv"
  )
)


readr::write_csv(
  model_info_all,
  file.path(
    table_dir,
    "model_predictor_sets.csv"
  )
)


# -----------------------------------------------------------------------------
# 8B. Best-model outputs
# -----------------------------------------------------------------------------

best_predictions_all <- bind_rows(
  all_best_predictions
)

best_metrics_all <- bind_rows(
  all_best_metrics
)

best_coefficients_all <- bind_rows(
  all_best_coefficients
)

best_diagnostics_all <- bind_rows(
  all_best_diagnostics
)

best_model_info_all <- bind_rows(
  all_best_info
)


readr::write_csv(
  best_predictions_all,
  file.path(
    best_table_dir,
    "predictions_all_responses_bestmod.csv"
  )
)


readr::write_csv(
  best_metrics_all,
  file.path(
    best_table_dir,
    "model_performance_all_responses_bestmod.csv"
  )
)


readr::write_csv(
  best_coefficients_all,
  file.path(
    best_table_dir,
    "coefficients_all_responses_bestmod.csv"
  )
)


readr::write_csv(
  best_diagnostics_all,
  file.path(
    best_diag_dir,
    "diagnostics_all_responses_bestmod.csv"
  )
)


readr::write_csv(
  best_model_info_all,
  file.path(
    best_table_dir,
    "model_predictor_sets_bestmod.csv"
  )
)


# -----------------------------------------------------------------------------
# 8C. Combined model-selection summary
# -----------------------------------------------------------------------------

if (length(all_selection_info) > 0L) {
  
  selection_summary_all <- bind_rows(
    all_selection_info
  )
  
  readr::write_csv(
    selection_summary_all,
    file.path(
      best_table_dir,
      "projpred_parsimony_selection_summary.csv"
    )
  )
}


# -----------------------------------------------------------------------------
# 8D. Combined selected predictor table
# -----------------------------------------------------------------------------

selected_predictor_files <- list.files(
  best_table_dir,
  pattern = "^selected_predictors_.*_bestmod\\.csv$",
  full.names = TRUE
)


if (length(selected_predictor_files) > 0L) {
  
  selected_predictors_all <- purrr::map_dfr(
    selected_predictor_files,
    readr::read_csv,
    show_col_types = FALSE
  )
  
  readr::write_csv(
    selected_predictors_all,
    file.path(
      best_table_dir,
      "selected_predictors_all_responses_bestmod.csv"
    )
  )
}


# -----------------------------------------------------------------------------
# 8E. Save posterior draws
# -----------------------------------------------------------------------------

saveRDS(
  all_draws,
  file.path(
    model_dir,
    "all_bayesian_C_d14C_reference_posterior_draws.rds"
  )
)


saveRDS(
  all_best_predictions,
  file.path(
    best_model_dir,
    "all_best_model_predictions.rds"
  )
)


# -----------------------------------------------------------------------------
# 8F. Save complete results bundle
# -----------------------------------------------------------------------------

saveRDS(
  
  list(
    
    full_reference = list(
      predictions = predictions_all,
      metrics = metrics_all,
      coefficients = coefficients_all,
      model_info = model_info_all
    ),
    
    best_models = list(
      predictions = best_predictions_all,
      metrics = best_metrics_all,
      coefficients = best_coefficients_all,
      diagnostics = best_diagnostics_all,
      model_info = best_model_info_all
    ),
    
    selection = all_selection_info,
    
    candidate_predictors = candidate_predictors,
    
    settings = list(
      
      min_predictor_observed_fraction =
        min_predictor_observed_fraction,
      
      max_abs_predictor_correlation =
        max_abs_predictor_correlation,
      
      slope_prior_sd =
        slope_prior_sd,
      
      max_submodel_terms =
        max_submodel_terms,
      
      projpred_cv_method =
        projpred_cv_method,
      
      parsimony_alpha =
        parsimony_alpha,
      
      parsimony_pct =
        parsimony_pct,
      
      parsimony_type =
        parsimony_type,
      
      n_chains =
        n_chains,
      
      n_warmup =
        n_warmup,
      
      n_sampling =
        n_sampling,
      
      adapt_delta =
        adapt_delta,
      
      max_treedepth =
        max_treedepth,
      
      backend =
        "cmdstanr_direct",
      
      selection_method =
        "projpred forward selection with PSIS-LOO evaluation",
      
      selection_criterion =
        "One-standard-error ELPD parsimony criterion"
      
    )
    
  ),
  
  file.path(
    out_dir,
    "bayesian_C_d14C_results_bundle_bestmod.rds"
  )
)


# =============================================================================
# 9. FINAL COMPARISON TABLE
# =============================================================================

if (
  nrow(metrics_all) > 0L &&
  nrow(best_metrics_all) > 0L
) {
  
  full_comparison <- metrics_all %>%
    filter(model == "full reference model") %>%
    select(
      response,
      n,
      full_n_predictors = n_predictors,
      full_RMSE = RMSE_in_sample,
      full_MAE = MAE_in_sample,
      full_R2 = squared_correlation_in_sample
    )
  
  best_comparison <- best_metrics_all %>%
    select(
      response,
      best_n_predictors = n_predictors,
      best_RMSE = RMSE_in_sample,
      best_MAE = MAE_in_sample,
      best_R2 = squared_correlation_in_sample
    )
  
  comparison <- full_comparison %>%
    left_join(
      best_comparison,
      by = "response"
    ) %>%
    mutate(
      predictors_removed =
        full_n_predictors - best_n_predictors,
      
      proportion_predictors_retained =
        best_n_predictors / full_n_predictors,
      
      RMSE_change =
        best_RMSE - full_RMSE,
      
      MAE_change =
        best_MAE - full_MAE,
      
      R2_change =
        best_R2 - full_R2
    )
  
  
  readr::write_csv(
    comparison,
    file.path(
      best_table_dir,
      "full_vs_best_model_comparison.csv"
    )
  )
}


# =============================================================================
# 10. COMPLETION
# =============================================================================

message("\n")
message("============================================================")
message("MALLE ERW BAYESIAN MODELLING COMPLETED")
message("============================================================")

message("\nMain output directory:")
message(out_dir)

message("\nReference model plots:")
message(plot_dir)

message("\nprojpred selection outputs:")
message(selection_dir)

message("\nBest-model plots:")
message(best_plot_dir)

message("\nBest-model tables:")
message(best_table_dir)

message("\nBest-model diagnostics:")
message(best_diag_dir)


message(
  "\nEach successfully selected response now has one best-model plot ",
  "ending in _bestmod.png."
)


message(
  "\nThe best-model plots contain observed versus predicted values, ",
  "95% fitted-mean credible ribbons, 95% posterior predictive intervals, ",
  "and relative predictor-importance bars underneath."
)


message(
  "\nThe selected model size is determined using projpred::suggest_size() ",
  "with the one-standard-error ELPD criterion."
)


message(
  "\nImportant: model performance metrics in the final refitted models ",
  "are in-sample metrics. The projpred LOO results are the relevant ",
  "cross-validated predictive performance estimates for model selection."
)


message(
  "\nRelative predictor importance is based on absolute posterior mean ",
  "standardised coefficients. It is not variance explained or a ",
  "causal contribution."
)


message(
  "\nExisting reference-model posterior draws are reused when present. ",
  "Delete the relevant saved draws and fit files if you change the data, ",
  "predictor preparation or reference-model specification."
)

message("\n============================================================")


# =============================================================================
# 11. POSTERIOR CORRELATIONS BETWEEN BEST-MODEL REGRESSION COEFFICIENTS
# =============================================================================

message("\nChecking posterior correlations between regression coefficients...")

# Settings
coefficient_correlation_threshold <- 0.70

# Storage
all_coefficient_correlations <- list()
all_coefficient_correlations_summary <- list()


# Helper: ensure CSV tables contain only ordinary atomic columns
check_csv_columns <- function(x, object_name) {
  
  invalid <- vapply(
    x,
    function(col) is.list(col) || is.matrix(col),
    logical(1)
  )
  
  if (any(invalid)) {
    stop(
      object_name,
      " contains unsupported columns: ",
      paste(names(x)[invalid], collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(x)
}


# ---------------------------------------------------------------------------
# 11A. Calculate posterior coefficient correlations
# ---------------------------------------------------------------------------

for (response in names(all_best_info)) {
  
  message("\nCoefficient correlation analysis: ", response)
  
  model_id <- safe_name(response)
  
  draws_file <- file.path(
    best_model_dir,
    paste0("draws_", model_id, "_bestmod.rds")
  )
  
  if (!file.exists(draws_file)) {
    message("Posterior draws not found. Skipping.")
    next
  }
  
  # Load posterior draws
  draws_array <- readRDS(draws_file)
  
  draw_matrix <- posterior::as_draws_matrix(draws_array)
  
  # Identify regression coefficients
  beta_names <- grep(
    "^beta\\[",
    colnames(draw_matrix),
    value = TRUE
  )
  
  if (length(beta_names) < 2L) {
    message("Fewer than two predictors: no coefficient correlations.")
    next
  }
  
  beta_draws <- draw_matrix[
    ,
    beta_names,
    drop = FALSE
  ]
  
  # Explicitly retain a numeric matrix
  beta_draws <- as.matrix(beta_draws)
  storage.mode(beta_draws) <- "double"
  
  # Retrieve actual predictor names
  predictor_names <- strsplit(
    as.character(all_best_info[[response]]$predictors[1]),
    " \\+ "
  )[[1]]
  
  if (length(predictor_names) != length(beta_names)) {
    warning(
      "Predictor names do not match coefficient draws for ",
      response
    )
    next
  }
  
  colnames(beta_draws) <- predictor_names
  
  
  # -------------------------------------------------------------------------
  # Posterior correlation matrix and long-format table
  # -------------------------------------------------------------------------
  
  correlation_matrix <- stats::cor(
    beta_draws,
    method = "pearson"
  )
  
  correlation_long <- tibble::tibble(
    response = rep(as.character(response), length(predictor_names)^2),
    Predictor_1 = rep(predictor_names, times = length(predictor_names)),
    Predictor_2 = rep(predictor_names, each = length(predictor_names)),
    Posterior_correlation = as.numeric(as.vector(correlation_matrix))
  )
  
  
  # -------------------------------------------------------------------------
  # Pairwise coefficient comparisons
  # -------------------------------------------------------------------------
  
  predictor_pairs <- utils::combn(
    seq_along(predictor_names),
    2
  )
  
  pairwise_results <- purrr::map_dfr(
    seq_len(ncol(predictor_pairs)),
    function(i) {
      
      j1 <- predictor_pairs[1, i]
      j2 <- predictor_pairs[2, i]
      
      b1 <- as.numeric(beta_draws[, j1])
      b2 <- as.numeric(beta_draws[, j2])
      
      r <- as.numeric(stats::cor(b1, b2))
      
      tibble::tibble(
        response = as.character(response),
        Predictor_1 = as.character(predictor_names[j1]),
        Predictor_2 = as.character(predictor_names[j2]),
        
        Posterior_correlation = r,
        
        Absolute_posterior_correlation = abs(r),
        
        Probability_same_sign = as.numeric(
          mean(sign(b1) == sign(b2))
        ),
        
        Probability_opposite_sign = as.numeric(
          mean(sign(b1) != sign(b2))
        ),
        
        Probability_both_positive = as.numeric(
          mean(b1 > 0 & b2 > 0)
        ),
        
        Probability_both_negative = as.numeric(
          mean(b1 < 0 & b2 < 0)
        ),
        
        Strong_posterior_correlation =
          abs(r) >= coefficient_correlation_threshold
      )
    }
  ) %>%
    arrange(desc(Absolute_posterior_correlation))
  
  
  # Check and save individual response results
  check_csv_columns(pairwise_results, "pairwise_results")
  
  readr::write_csv(
    pairwise_results,
    file.path(
      best_table_dir,
      paste0(
        "posterior_coefficient_correlations_",
        model_id,
        "_bestmod.csv"
      )
    )
  )
  
  
  # Store combined results
  all_coefficient_correlations[[response]] <-
    correlation_long
  
  all_coefficient_correlations_summary[[response]] <-
    pairwise_results
  
  
  # -------------------------------------------------------------------------
  # Heatmap
  # -------------------------------------------------------------------------
  
  p_cor <- ggplot2::ggplot(
    correlation_long,
    ggplot2::aes(
      x = Predictor_1,
      y = Predictor_2,
      fill = Posterior_correlation
    )
  ) +
    
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.4
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", Posterior_correlation)
      ),
      size = 3
    ) +
    
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    
    ggplot2::coord_equal() +
    
    ggplot2::labs(
      title = paste(
        "Posterior coefficient correlations:",
        response
      ),
      subtitle = paste0(
        "Selected Bayesian model | ",
        nrow(pairwise_results),
        " coefficient pairs"
      ),
      x = NULL,
      y = NULL,
      fill = "Posterior\ncorrelation"
    ) +
    
    ggplot2::theme_bw(base_size = 11) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid = ggplot2::element_blank()
    )
  
  
  ggplot2::ggsave(
    filename = file.path(
      best_plot_dir,
      paste0(
        "posterior_coefficient_correlation_",
        model_id,
        "_bestmod.png"
      )
    ),
    plot = p_cor,
    width = max(7, 0.65 * length(predictor_names)),
    height = max(6, 0.65 * length(predictor_names)),
    dpi = 320,
    limitsize = FALSE
  )
  
  
  # Report highly correlated coefficient pairs
  high_cor <- pairwise_results %>%
    filter(Strong_posterior_correlation)
  
  if (nrow(high_cor) > 0L) {
    
    message(
      "Strongly correlated coefficient pairs (|r| >= ",
      coefficient_correlation_threshold,
      "):"
    )
    
    print(high_cor)
    
  } else {
    
    message(
      "No coefficient pairs exceed |r| = ",
      coefficient_correlation_threshold
    )
  }
  
}


# ---------------------------------------------------------------------------
# 11B. Combined correlation tables
# ---------------------------------------------------------------------------

if (length(all_coefficient_correlations_summary) > 0L) {
  
  coefficient_correlations_all <- dplyr::bind_rows(
    all_coefficient_correlations_summary
  )
  
  coefficient_correlation_matrices_all <- dplyr::bind_rows(
    all_coefficient_correlations
  )
  
  # Validate both combined tables
  check_csv_columns(
    coefficient_correlations_all,
    "coefficient_correlations_all"
  )
  
  check_csv_columns(
    coefficient_correlation_matrices_all,
    "coefficient_correlation_matrices_all"
  )
  
  # Export combined pairwise results
  readr::write_csv(
    coefficient_correlations_all,
    file.path(
      best_table_dir,
      "posterior_coefficient_correlations_all_responses_bestmod.csv"
    )
  )
  
  # Export complete long-format correlation matrices
  readr::write_csv(
    coefficient_correlation_matrices_all,
    file.path(
      best_table_dir,
      "posterior_coefficient_correlation_matrices_all_responses_bestmod.csv"
    )
  )
  
  
  # -------------------------------------------------------------------------
  # Summary per response
  # -------------------------------------------------------------------------
  
  correlation_summary <- coefficient_correlations_all %>%
    group_by(response) %>%
    summarise(
      n_coefficient_pairs = n(),
      
      maximum_absolute_correlation =
        max(Absolute_posterior_correlation, na.rm = TRUE),
      
      mean_absolute_correlation =
        mean(Absolute_posterior_correlation, na.rm = TRUE),
      
      n_strongly_correlated_pairs =
        sum(Strong_posterior_correlation),
      
      proportion_strongly_correlated =
        mean(Strong_posterior_correlation),
      
      .groups = "drop"
    )
  
  readr::write_csv(
    correlation_summary,
    file.path(
      best_table_dir,
      "posterior_coefficient_correlation_summary_bestmod.csv"
    )
  )
  
  print(correlation_summary)
}


message("\nPosterior coefficient correlation analysis completed.")
message("Tables saved in: ", best_table_dir)
message("Heatmaps saved in: ", best_plot_dir)


# =============================================================================
# SECTION 12. DIAGNOSE POSTERIOR COEFFICIENT CORRELATIONS
#
# PURPOSE:
#   Determine whether posterior coefficient correlations are associated
#   with correlations between the original predictor variables.
#
#   This section:
#     1. Identifies posterior coefficient pairs with |r| >= 0.50.
#     2. Calculates raw Pearson and Spearman predictor correlations.
#     3. Calculates correlations within the observations available for
#        the corresponding model response.
#     4. Reports sample sizes and missingness.
#     5. Provides recommendations for targeted reduced-model sensitivity tests.
#
# IMPORTANT:
#   - Does NOT modify the original predictor sets.
#   - Does NOT refit Bayesian models.
#   - Does NOT perform projpred selection.
#   - Does NOT remove predictors automatically.
# =============================================================================


# -----------------------------------------------------------------------------
# 12A. SETTINGS
# -----------------------------------------------------------------------------

posterior_correlation_file <- file.path(
  best_table_dir,
  "posterior_coefficient_correlations_all_responses_bestmod.csv"
)

if (!file.exists(posterior_correlation_file)) {
  stop(
    "Cannot find posterior coefficient correlation file: ",
    posterior_correlation_file
  )
}

# Include moderately strong posterior coefficient correlations.
posterior_threshold <- 0.50

# Threshold for flagging substantial raw predictor collinearity.
raw_correlation_threshold <- 0.70


# -----------------------------------------------------------------------------
# 12A.1. ENSURE ORIGINAL MODEL DATA ARE AVAILABLE
# -----------------------------------------------------------------------------

# Recreate dat using exactly the same treatment filtering as the original
# Bayesian modelling script, if dat is not already present.

if (!exists("dat", inherits = TRUE)) {
  
  message("Object 'dat' not found. Reloading the original integrated dataset.")
  
  input_file <- here::here(
    "outputs",
    "integrated_geochemistry_for_bayesian.csv"
  )
  
  if (!file.exists(input_file)) {
    stop(
      "Cannot find integrated dataset: ",
      input_file,
      "\nCheck that your R working directory is the project root."
    )
  }
  
  integrated <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    na = c("", "NA", "#VALUE!")
  )
  
  if (!all(c("sample", "Treatment") %in% names(integrated))) {
    stop("Integrated dataset must contain sample and Treatment columns.")
  }
  
  normalise_treatment <- function(x) {
    z <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
    z <- stringr::str_replace_all(z, "ü", "u")
    z <- stringr::str_replace_all(z, "[^a-z0-9]+", "")
    
    dplyr::case_when(
      stringr::str_detect(z, "^control$") ~ "Control",
      stringr::str_detect(z, "^lime") ~ "Lime",
      stringr::str_detect(z, "^bolsdorfer") ~ "Bolsdorfer",
      stringr::str_detect(z, "^huhnerberg") ~ "Huhnerberg",
      stringr::str_detect(z, "^eifelgold") ~ "Eifelgold",
      TRUE ~ NA_character_
    )
  }
  
  allowed_treatments <- c(
    "Control",
    "Lime_2",
    "Bolsdorfer_50",
    "Huhnerberg_50",
    "Eifelgold_50"
  )
  
  dat <- integrated %>%
    dplyr::filter(Treatment %in% allowed_treatments) %>%
    dplyr::mutate(
      sample = suppressWarnings(as.integer(as.character(sample))),
      Treatment_model = normalise_treatment(Treatment)
    ) %>%
    dplyr::filter(
      !is.na(sample),
      !is.na(Treatment_model)
    )
  
  message("Original model dataset reloaded: ", nrow(dat), " observations.")
}

# Confirm that the expected response is present.
if (!"d14C_2026_fPOM" %in% names(dat)) {
  stop("Response d14C_2026_fPOM is not present in dat.")
}

# -----------------------------------------------------------------------------
# 12B. LOAD POSTERIOR CORRELATION RESULTS
# -----------------------------------------------------------------------------

posterior_correlations <- readr::read_csv(
  posterior_correlation_file,
  show_col_types = FALSE
)

# Exclude diagonal correlations and retain pairs above the threshold.
posterior_pairs <- posterior_correlations %>%
  dplyr::filter(
    Predictor_1 != Predictor_2,
    Absolute_posterior_correlation >= posterior_threshold
  ) %>%
  dplyr::distinct(
    response,
    Predictor_1,
    Predictor_2,
    .keep_all = TRUE
  )

if (nrow(posterior_pairs) == 0L) {
  message("No posterior coefficient pairs exceed the selected threshold.")
} else {
  message(
    "\nPosterior coefficient pairs selected for investigation: ",
    nrow(posterior_pairs)
  )
  print(posterior_pairs)
}


# -----------------------------------------------------------------------------
# 12C. HELPER FUNCTION
# -----------------------------------------------------------------------------

calculate_raw_correlations <- function(data, response, predictor_1, predictor_2) {
  
  required_columns <- c(response, predictor_1, predictor_2)
  
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0L) {
    return(tibble::tibble(
      n_response = NA_integer_,
      n_complete = NA_integer_,
      missing_predictor_1 = NA_integer_,
      missing_predictor_2 = NA_integer_,
      spearman_rho = NA_real_,
      pearson_r = NA_real_,
      spearman_p = NA_real_,
      pearson_p = NA_real_
    ))
  }
  
  y <- suppressWarnings(as.numeric(data[[response]]))
  x1 <- suppressWarnings(as.numeric(data[[predictor_1]]))
  x2 <- suppressWarnings(as.numeric(data[[predictor_2]]))
  
  # Restrict to observations with a finite response.
  response_rows <- is.finite(y)
  
  # Among these observations, identify complete predictor pairs.
  complete_rows <- response_rows & is.finite(x1) & is.finite(x2)
  
  n_response <- sum(response_rows)
  n_complete <- sum(complete_rows)
  
  missing_1 <- sum(response_rows & !is.finite(x1))
  missing_2 <- sum(response_rows & !is.finite(x2))
  
  rho <- NA_real_
  pearson <- NA_real_
  rho_p <- NA_real_
  pearson_p <- NA_real_
  
  if (n_complete >= 3L) {
    
    x1_complete <- x1[complete_rows]
    x2_complete <- x2[complete_rows]
    
    if (stats::sd(x1_complete) > 0 &&
        stats::sd(x2_complete) > 0) {
      
      rho <- suppressWarnings(
        stats::cor(
          x1_complete,
          x2_complete,
          method = "spearman"
        )
      )
      
      pearson <- stats::cor(
        x1_complete,
        x2_complete,
        method = "pearson"
      )
      
      rho_p <- tryCatch(
        suppressWarnings(
          stats::cor.test(
            x1_complete,
            x2_complete,
            method = "spearman",
            exact = FALSE
          )$p.value
        ),
        error = function(e) NA_real_
      )
      
      pearson_p <- tryCatch(
        stats::cor.test(
          x1_complete,
          x2_complete,
          method = "pearson"
        )$p.value,
        error = function(e) NA_real_
      )
    }
  }
  
  tibble::tibble(
    n_response = as.integer(n_response),
    n_complete = as.integer(n_complete),
    missing_predictor_1 = as.integer(missing_1),
    missing_predictor_2 = as.integer(missing_2),
    spearman_rho = as.numeric(rho),
    pearson_r = as.numeric(pearson),
    spearman_p = as.numeric(rho_p),
    pearson_p = as.numeric(pearson_p)
  )
}


# -----------------------------------------------------------------------------
# 12D. CALCULATE RAW CORRELATIONS FOR EACH POSTERIOR PAIR
# -----------------------------------------------------------------------------

if (nrow(posterior_pairs) > 0L) {
  
  raw_correlation_results <- purrr::pmap_dfr(
    posterior_pairs %>%
      dplyr::select(
        response,
        Predictor_1,
        Predictor_2,
        Posterior_correlation,
        Absolute_posterior_correlation,
        Probability_opposite_sign
      ),
    
    function(
    response,
    Predictor_1,
    Predictor_2,
    Posterior_correlation,
    Absolute_posterior_correlation,
    Probability_opposite_sign
    ) {
      
      result <- calculate_raw_correlations(
        data = dat,
        response = response,
        predictor_1 = Predictor_1,
        predictor_2 = Predictor_2
      )
      
      dplyr::bind_cols(
        tibble::tibble(
          response = response,
          Predictor_1 = Predictor_1,
          Predictor_2 = Predictor_2,
          Posterior_correlation = Posterior_correlation,
          Absolute_posterior_correlation =
            Absolute_posterior_correlation,
          Probability_opposite_sign = Probability_opposite_sign
        ),
        result
      )
    }
  )
  
  
  # ---------------------------------------------------------------------------
  # 12E. CLASSIFY PAIRS FOR FURTHER INVESTIGATION
  # ---------------------------------------------------------------------------
  
  raw_correlation_results <- raw_correlation_results %>%
    dplyr::mutate(
      
      Absolute_raw_spearman =
        abs(spearman_rho),
      
      Absolute_raw_pearson =
        abs(pearson_r),
      
      recommendation = dplyr::case_when(
        
        is.na(spearman_rho) ~
          "Check data availability and predictor variation",
        
        Absolute_posterior_correlation >= 0.70 &
          Absolute_raw_spearman >= raw_correlation_threshold ~
          "Consider reduced-model sensitivity analysis",
        
        Absolute_posterior_correlation >= 0.70 &
          Absolute_raw_spearman < raw_correlation_threshold ~
          "Investigate posterior dependence; raw collinearity is not strong",
        
        Absolute_posterior_correlation >= posterior_threshold &
          Absolute_raw_spearman >= raw_correlation_threshold ~
          "Investigate predictor collinearity; consider sensitivity analysis",
        
        TRUE ~
          "No refit indicated by correlation evidence alone"
      )
    ) %>%
    dplyr::arrange(
      dplyr::desc(Absolute_posterior_correlation)
    )
  
  
  # ---------------------------------------------------------------------------
  # 12F. SAVE RESULTS
  # ---------------------------------------------------------------------------
  
  readr::write_csv(
    raw_correlation_results,
    file.path(
      best_table_dir,
      "posterior_vs_raw_predictor_correlations.csv"
    )
  )
  
  message("\nPosterior versus raw predictor correlations:")
  print(
    raw_correlation_results %>%
      dplyr::select(
        response,
        Predictor_1,
        Predictor_2,
        Posterior_correlation,
        spearman_rho,
        pearson_r,
        n_complete,
        Probability_opposite_sign,
        recommendation
      ),
    n = Inf
  )
}


# -----------------------------------------------------------------------------
# 12G. SPECIFIC CHECK: TEMPORAL MANGANESE CEC MEASUREMENTS
#
# Mn_Oct2023 and CEC_Mn represent measurements at different time points.
# Their correlation can indicate temporal persistence, but does not by itself
# establish that either measurement is redundant.
# -----------------------------------------------------------------------------

mn_temporal_check <- calculate_raw_correlations(
  data = dat,
  response = "d14C_2026_fPOM",
  predictor_1 = "Mn_Oct2023",
  predictor_2 = "CEC_Mn"
)

mn_temporal_check <- dplyr::bind_cols(
  tibble::tibble(
    comparison = "Mn_Oct2023 versus CEC_Mn",
    response = "d14C_2026_fPOM"
  ),
  mn_temporal_check
)

readr::write_csv(
  mn_temporal_check,
  file.path(
    best_table_dir,
    "temporal_Mn_CEC_correlation_check.csv"
  )
)

message("\nTemporal Mn CEC comparison:")
print(mn_temporal_check)


# -----------------------------------------------------------------------------
# 12H. FINAL SUMMARY
# -----------------------------------------------------------------------------

message("\n============================================================")
message("CORRELATION DIAGNOSTIC COMPLETE")
message("============================================================")

message(
  "\nResults saved to: ",
  best_table_dir
)

message(
  "\nNo Bayesian models have been refitted or modified."
)


# =============================================================================
# MALLE ERW: d14C_fPOM FIVE-MODEL SENSITIVITY ANALYSIS
#
# PURPOSE
#   Compare the existing selected d14C_2026_fPOM model against four
#   reduced versions, each excluding one potentially correlated predictor.
#
# MODELS
#   M1: Full selected best model
#   M2: Exclude Mn_Oct2023
#   M3: Exclude CEC_Mn
#   M4: Exclude FTIR_delta_aromatic
#   M5: Exclude FTIR_delta_Carbonate_bend
#
# IMPORTANT
#   - Uses the existing projpred-selected best-model predictor set.
#   - Uses identical observations and preprocessing for all five models.
#   - Uses the original Bayesian regression specification and priors.
#   - Fits five models using CmdStan.
#   - Calculates PSIS-LOO from Stan's generated posterior mean draws.
#   - Saves cumulative checkpoints after each completed model.
#   - Does not rerun projpred or any other response models.
#   - In-sample metrics are descriptive, not independent validation.
#
# =============================================================================


# =============================================================================
# 1. PACKAGES AND PATHS
# =============================================================================

required_packages <- c(
  "tidyverse",
  "here",
  "cmdstanr",
  "posterior",
  "loo"
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

if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not configured.")
}

project_dir <- here::here()

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

best_dir <- file.path(out_dir, "best_model")
best_model_dir <- file.path(best_dir, "models")
best_table_dir <- file.path(best_dir, "tables")

sensitivity_dir <- file.path(
  best_dir,
  "d14C_fPOM_predictor_sensitivity"
)

sensitivity_model_dir <- file.path(sensitivity_dir, "models")
sensitivity_table_dir <- file.path(sensitivity_dir, "tables")
sensitivity_diag_dir <- file.path(sensitivity_dir, "diagnostics")
sensitivity_plot_dir <- file.path(sensitivity_dir, "plots")

for (d in c(
  sensitivity_dir,
  sensitivity_model_dir,
  sensitivity_table_dir,
  sensitivity_diag_dir,
  sensitivity_plot_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# 2. SETTINGS
# =============================================================================

response <- "d14C_2026_fPOM"

# These match the original Bayesian model settings.
n_chains <- 4L
n_warmup <- 1500L
n_sampling <- 2000L
n_cores <- min(n_chains, parallel::detectCores(logical = FALSE))

adapt_delta <- 0.99
max_treedepth <- 12L
slope_prior_sd <- 0.5

set.seed(20261123)


# =============================================================================
# 3. LOAD THE ORIGINAL DATA
# =============================================================================

if (!file.exists(input_file)) {
  stop("Integrated dataset not found: ", input_file)
}

integrated <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  na = c("", "NA", "#VALUE!")
)

normalise_treatment <- function(x) {
  
  z <- stringr::str_to_lower(
    stringr::str_squish(as.character(x))
  )
  
  z <- stringr::str_replace_all(z, "ü", "u")
  z <- stringr::str_replace_all(z, "[^a-z0-9]+", "")
  
  dplyr::case_when(
    stringr::str_detect(z, "^control$") ~ "Control",
    stringr::str_detect(z, "^lime") ~ "Lime",
    stringr::str_detect(z, "^bolsdorfer") ~ "Bolsdorfer",
    stringr::str_detect(z, "^huhnerberg") ~ "Huhnerberg",
    stringr::str_detect(z, "^eifelgold") ~ "Eifelgold",
    TRUE ~ NA_character_
  )
}

allowed_treatments <- c(
  "Control",
  "Lime_2",
  "Bolsdorfer_50",
  "Huhnerberg_50",
  "Eifelgold_50"
)

dat <- integrated %>%
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(as.integer(as.character(sample))),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(
    !is.na(sample),
    !is.na(Treatment_model)
  )

if (!response %in% names(dat)) {
  stop("Response column not found: ", response)
}


# =============================================================================
# 4. LOAD THE EXISTING SELECTED BEST-MODEL PREDICTORS
# =============================================================================

selected_predictor_file <- file.path(
  best_table_dir,
  paste0(
    "selected_predictors_",
    response,
    "_bestmod.csv"
  )
)

if (!file.exists(selected_predictor_file)) {
  stop(
    "Selected predictor file not found: ",
    selected_predictor_file
  )
}

selected_tbl <- readr::read_csv(
  selected_predictor_file,
  show_col_types = FALSE
)

if (!"predictor" %in% names(selected_tbl)) {
  stop("The selected predictor file has no 'predictor' column.")
}

best_predictors <- as.character(selected_tbl$predictor)
best_predictors <- best_predictors[
  !is.na(best_predictors) & best_predictors != ""
]

best_predictors <- unique(best_predictors)

message("\nExisting selected best-model predictors:")
print(best_predictors)


# =============================================================================
# 5. DEFINE THE FIVE MODEL VERSIONS
# =============================================================================

sensitivity_variables <- c(
  "Mn_Oct2023",
  "CEC_Mn",
  "FTIR_delta_aromatic",
  "FTIR_delta_Carbonate_bend"
)

missing_from_best_model <- setdiff(
  sensitivity_variables,
  best_predictors
)

if (length(missing_from_best_model) > 0L) {
  
  stop(
    "\nThe following variables are not in the selected best model:\n",
    paste(missing_from_best_model, collapse = ", "),
    "\n\nNo sensitivity analysis was run. Check the selected predictor CSV."
  )
}

missing_from_data <- setdiff(
  best_predictors,
  names(dat)
)

if (length(missing_from_data) > 0L) {
  stop(
    "Selected predictors missing from integrated dataset: ",
    paste(missing_from_data, collapse = ", ")
  )
}

model_specs <- list(
  
  M1_full = list(
    label = "Full selected best model",
    exclude = character(0)
  ),
  
  M2_without_Mn_Oct2023 = list(
    label = "Without Mn_Oct2023",
    exclude = "Mn_Oct2023"
  ),
  
  M3_without_CEC_Mn = list(
    label = "Without CEC_Mn",
    exclude = "CEC_Mn"
  ),
  
  M4_without_FTIR_aromatic = list(
    label = "Without FTIR_delta_aromatic",
    exclude = "FTIR_delta_aromatic"
  ),
  
  M5_without_FTIR_carbonate = list(
    label = "Without FTIR_delta_Carbonate_bend",
    exclude = "FTIR_delta_Carbonate_bend"
  )
)

model_predictors <- lapply(
  model_specs,
  function(x) setdiff(best_predictors, x$exclude)
)

model_design_table <- purrr::imap_dfr(
  model_predictors,
  function(vars, model_id) {
    tibble(
      model_id = model_id,
      n_predictors = length(vars),
      predictors = paste(vars, collapse = " + "),
      excluded = paste(
        model_specs[[model_id]]$exclude,
        collapse = " + "
      )
    )
  }
)

readr::write_csv(
  model_design_table,
  file.path(
    sensitivity_table_dir,
    "five_model_predictor_sets.csv"
  )
)

print(model_design_table)


# =============================================================================
# 6. PREPARE ONE SHARED DATASET FOR ALL FIVE MODELS
# =============================================================================
#
# Important:
#   All models use exactly the same response observations.
#   All predictors are imputed and standardised once using this shared dataset.
#   Reduced models simply omit one column from the shared predictor matrix.
#
# This prevents differences in sample size, imputation or scaling from
# confounding the comparison.
# =============================================================================

y_raw_all <- suppressWarnings(
  as.numeric(dat[[response]])
)

keep_y <- is.finite(y_raw_all)

work <- dat[keep_y, , drop = FALSE]
y_raw <- y_raw_all[keep_y]

if (length(y_raw) < 8L) {
  stop("Fewer than eight finite response observations.")
}

# Check that all selected predictors have sufficient data and variation.
predictor_check <- purrr::map_dfr(
  best_predictors,
  function(v) {
    
    x <- suppressWarnings(as.numeric(work[[v]]))
    
    tibble(
      predictor = v,
      observed_fraction = mean(is.finite(x)),
      n_missing = sum(!is.finite(x)),
      sd = stats::sd(x, na.rm = TRUE),
      eligible =
        mean(is.finite(x)) >= 0.70 &&
        is.finite(stats::sd(x, na.rm = TRUE)) &&
        stats::sd(x, na.rm = TRUE) > 0
    )
  }
)

readr::write_csv(
  predictor_check,
  file.path(
    sensitivity_table_dir,
    "shared_predictor_data_check.csv"
  )
)

if (any(!predictor_check$eligible)) {
  stop(
    "At least one selected predictor fails the original missingness/variation criteria. ",
    "See shared_predictor_data_check.csv."
  )
}

# Median imputation and standardisation, performed once.
X_shared <- as.data.frame(
  lapply(best_predictors, function(v) {
    
    x <- suppressWarnings(
      as.numeric(work[[v]])
    )
    
    med <- stats::median(
      x[is.finite(x)],
      na.rm = TRUE
    )
    
    x[!is.finite(x)] <- med
    
    as.numeric(scale(x))
  })
)

names(X_shared) <- best_predictors

if (!all(vapply(
  X_shared,
  function(x) all(is.finite(x)),
  logical(1)
))) {
  stop("Non-finite values remain in the shared predictor matrix.")
}

y_mean <- mean(y_raw)
y_sd <- stats::sd(y_raw)

if (!is.finite(y_sd) || y_sd <= 0) {
  stop("Response has no finite variation.")
}

y_std <- (y_raw - y_mean) / y_sd

shared_data <- tibble(
  sample = work$sample,
  Treatment_model = as.character(work$Treatment_model),
  observed = y_raw,
  y_std = y_std
) %>%
  bind_cols(as_tibble(X_shared))

message(
  "\nShared observations: ", nrow(shared_data),
  "\nShared predictors: ", length(best_predictors)
)

# Save the exact shared preprocessing information.
preprocessing_summary <- predictor_check %>%
  mutate(
    median_imputed = n_missing > 0,
    standardisation = "z-score across shared response-complete observations"
  )

readr::write_csv(
  preprocessing_summary,
  file.path(
    sensitivity_table_dir,
    "shared_preprocessing_summary.csv"
  )
)


# =============================================================================
# 7. COMPILE THE ORIGINAL BAYESIAN REGRESSION MODEL
# =============================================================================

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, SLOPE_PRIOR_SD);
  sigma ~ exponential(1);

  y ~ normal(alpha + X * beta, sigma);
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;

  mu = alpha + X * beta;

  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
}
"

stan_code <- sub(
  "SLOPE_PRIOR_SD",
  format(slope_prior_sd, scientific = FALSE),
  stan_code,
  fixed = TRUE
)

stan_file <- file.path(
  sensitivity_model_dir,
  "d14C_fPOM_sensitivity_model.stan"
)

writeLines(stan_code, stan_file, useBytes = TRUE)

stan_model <- cmdstanr::cmdstan_model(
  stan_file,
  quiet = TRUE
)


# =============================================================================
# 8. FIT EACH MODEL
# =============================================================================

all_fits <- list()
all_draw_matrices <- list()
all_predictions <- list()
all_coefficients <- list()
all_metrics <- list()
all_loo <- list()
all_diagnostics <- list()

for (i in seq_along(model_specs)) {
  
  model_id <- names(model_specs)[i]
  model_label <- model_specs[[model_id]]$label
  predictors <- model_predictors[[model_id]]
  
  message("\n============================================================")
  message("FITTING: ", model_id)
  message("============================================================")
  message("Predictors: ", paste(predictors, collapse = ", "))
  
  X <- as.matrix(
    shared_data[, predictors, drop = FALSE]
  )
  
  storage.mode(X) <- "double"
  
  stan_data <- list(
    N = nrow(X),
    K = ncol(X),
    X = X,
    y = as.numeric(shared_data$y_std)
  )
  
  csv_dir <- file.path(
    sensitivity_model_dir,
    paste0("cmdstan_csv_", model_id)
  )
  
  dir.create(
    csv_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  fit <- stan_model$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_cores,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = 20261123 + i,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 200,
    output_dir = csv_dir
  )
  
  draws_array <- fit$draws(
    format = "draws_array"
  )
  
  draw_matrix <- posterior::as_draws_matrix(
    draws_array
  )
  
  all_fits[[model_id]] <- fit
  all_draw_matrices[[model_id]] <- draw_matrix
  
  saveRDS(
    draws_array,
    file.path(
      sensitivity_model_dir,
      paste0("draws_", model_id, ".rds")
    )
  )
  
  saveRDS(
    fit,
    file.path(
      sensitivity_model_dir,
      paste0("fit_", model_id, ".rds")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 8A. POSTERIOR PREDICTIONS
  # ---------------------------------------------------------------------------
  
  mu_cols <- paste0("mu[", seq_len(nrow(shared_data)), "]")
  rep_cols <- paste0("y_rep[", seq_len(nrow(shared_data)), "]")
  
  if (!all(mu_cols %in% colnames(draw_matrix))) {
    stop("Stan output is missing one or more mu columns in ", model_id, ".")
  }
  if (!all(rep_cols %in% colnames(draw_matrix))) {
    stop("Stan output is missing one or more y_rep columns in ", model_id, ".")
  }
  
  mu_draws <- draw_matrix[, mu_cols, drop = FALSE]
  rep_draws <- draw_matrix[, rep_cols, drop = FALSE]
  
  pred <- tibble(
    model_id = model_id,
    model = model_label,
    sample = shared_data$sample,
    Treatment_model = shared_data$Treatment_model,
    observed = shared_data$observed,
    
    predicted_mean =
      colMeans(mu_draws) * y_sd + y_mean,
    
    predicted_median =
      apply(mu_draws, 2, median) * y_sd + y_mean,
    
    fitted_lower_95 =
      apply(mu_draws, 2, quantile, probs = 0.025) * y_sd + y_mean,
    
    fitted_upper_95 =
      apply(mu_draws, 2, quantile, probs = 0.975) * y_sd + y_mean,
    
    predictive_lower_95 =
      apply(rep_draws, 2, quantile, probs = 0.025) * y_sd + y_mean,
    
    predictive_upper_95 =
      apply(rep_draws, 2, quantile, probs = 0.975) * y_sd + y_mean
  )
  
  all_predictions[[model_id]] <- pred
  
  
  # ---------------------------------------------------------------------------
  # 8B. POSTERIOR COEFFICIENT SUMMARIES
  # ---------------------------------------------------------------------------
  
  coef_names <- c(
    "alpha",
    paste0("beta[", seq_along(predictors), "]")
  )
  
  coef_labels <- c(
    "Intercept",
    predictors
  )
  
  coef_tbl <- purrr::map2_dfr(
    coef_names,
    coef_labels,
    function(coef_name, term_name) {
      
      vals <- draw_matrix[, coef_name]
      
      tibble(
        model_id = model_id,
        model = model_label,
        term = term_name,
        posterior_mean = mean(vals),
        posterior_sd = sd(vals),
        lower_95 = quantile(vals, 0.025),
        median = median(vals),
        upper_95 = quantile(vals, 0.975),
        probability_positive = mean(vals > 0),
        probability_negative = mean(vals < 0)
      )
    }
  ) %>%
    mutate(
      estimate_original_response_per_predictor_SD =
        posterior_mean * y_sd,
      lower_original_response_per_predictor_SD =
        lower_95 * y_sd,
      upper_original_response_per_predictor_SD =
        upper_95 * y_sd
    )
  
  all_coefficients[[model_id]] <- coef_tbl
  
  
  # ---------------------------------------------------------------------------
  # 8C. IN-SAMPLE DESCRIPTIVE METRICS
  # ---------------------------------------------------------------------------
  
  rmse <- sqrt(
    mean((pred$observed - pred$predicted_mean)^2)
  )
  
  mae <- mean(
    abs(pred$observed - pred$predicted_mean)
  )
  
  r2 <- if (
    stats::sd(pred$observed) > 0 &&
    stats::sd(pred$predicted_mean) > 0
  ) {
    cor(pred$observed, pred$predicted_mean)^2
  } else {
    NA_real_
  }
  
  all_metrics[[model_id]] <- tibble(
    model_id = model_id,
    model = model_label,
    n = nrow(pred),
    n_predictors = length(predictors),
    predictors = paste(predictors, collapse = " + "),
    RMSE_in_sample = rmse,
    MAE_in_sample = mae,
    R2_in_sample = r2,
    fitted_interval_coverage_95 = mean(
      pred$observed >= pred$fitted_lower_95 &
        pred$observed <= pred$fitted_upper_95
    ),
    predictive_interval_coverage_95 = mean(
      pred$observed >= pred$predictive_lower_95 &
        pred$observed <= pred$predictive_upper_95
    )
  )
  
  
  # Save predictions, coefficients and descriptive metrics before PSIS-LOO.
  # These remain available even if a later diagnostic step fails.
  readr::write_csv(
    bind_rows(all_predictions),
    file.path(sensitivity_table_dir, "predictions_checkpoint.csv")
  )
  readr::write_csv(
    bind_rows(all_coefficients),
    file.path(sensitivity_table_dir, "coefficients_checkpoint.csv")
  )
  readr::write_csv(
    bind_rows(all_metrics),
    file.path(sensitivity_table_dir, "in_sample_metrics_checkpoint.csv")
  )
  saveRDS(
    all_draw_matrices,
    file.path(sensitivity_model_dir, "draw_matrices_checkpoint.rds")
  )
  
  
  # ---------------------------------------------------------------------------
  # 8D. PSIS-LOO LOG-LIKELIHOOD
  # ---------------------------------------------------------------------------
  # Use Stan's generated-quantity mu draws directly. This avoids rebuilding
  # alpha + X beta and prevents draw/observation dimension mismatches.
  
  sigma_draws <- as.numeric(draw_matrix[, "sigma"])
  mu_std <- mu_draws
  
  if (nrow(mu_std) != length(sigma_draws) ||
      ncol(mu_std) != nrow(X)) {
    stop(
      "Unexpected dimensions for posterior means in ", model_id,
      ": mu = ", paste(dim(mu_std), collapse = " x "),
      "; sigma draws = ", length(sigma_draws),
      "; observations = ", nrow(X)
    )
  }
  
  log_lik <- matrix(
    NA_real_,
    nrow = nrow(mu_std),
    ncol = nrow(X)
  )
  
  for (j in seq_len(nrow(X))) {
    log_lik[, j] <- dnorm(
      as.numeric(shared_data$y_std[j]),
      mean = mu_std[, j],
      sd = sigma_draws,
      log = TRUE
    )
  }
  
  colnames(log_lik) <- paste0(
    "observation_",
    seq_len(ncol(log_lik))
  )
  
  loo_result <- loo::loo(log_lik)
  all_loo[[model_id]] <- loo_result
  
  saveRDS(
    all_loo,
    file.path(sensitivity_model_dir, "psis_loo_objects_checkpoint.rds")
  )
  
  # ---------------------------------------------------------------------------
  # 8E. CONVERGENCE DIAGNOSTICS
  # ---------------------------------------------------------------------------
  
  fit_summary <- fit$summary(
    variables = c("alpha", "beta", "sigma")
  )
  
  sampler_diag <- fit$sampler_diagnostics(
    format = "df"
  )
  
  divergences <- sum(
    sampler_diag$divergent__
  )
  
  max_treedepth_hits <- sum(
    sampler_diag$treedepth__ >= max_treedepth
  )
  
  all_diagnostics[[model_id]] <- tibble(
    model_id = model_id,
    model = model_label,
    max_Rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_bulk_ESS = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_tail_ESS = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = divergences,
    max_treedepth_hits = max_treedepth_hits,
    n_observations = nrow(X),
    n_predictors = length(predictors)
  )
  
  # Save cumulative outputs after EVERY completed model. If a later model
  # or downstream step fails, completed results are still available.
  readr::write_csv(
    bind_rows(all_predictions),
    file.path(sensitivity_table_dir, "predictions_checkpoint.csv")
  )
  readr::write_csv(
    bind_rows(all_coefficients),
    file.path(sensitivity_table_dir, "coefficients_checkpoint.csv")
  )
  readr::write_csv(
    bind_rows(all_metrics),
    file.path(sensitivity_table_dir, "in_sample_metrics_checkpoint.csv")
  )
  readr::write_csv(
    bind_rows(all_diagnostics),
    file.path(sensitivity_diag_dir, "sampling_diagnostics_checkpoint.csv")
  )
  saveRDS(
    all_draw_matrices,
    file.path(sensitivity_model_dir, "draw_matrices_checkpoint.rds")
  )
  
  message(
    "Completed ", model_id,
    " | max Rhat = ",
    round(max(fit_summary$rhat, na.rm = TRUE), 3),
    " | divergences = ", divergences
  )
}


# =============================================================================
# 9. COMBINE AND SAVE RESULTS
# =============================================================================

predictions_all <- bind_rows(all_predictions)
coefficients_all <- bind_rows(all_coefficients)
metrics_all <- bind_rows(all_metrics)
diagnostics_all <- bind_rows(all_diagnostics)

readr::write_csv(
  predictions_all,
  file.path(sensitivity_table_dir, "predictions_all_five_models.csv")
)

readr::write_csv(
  coefficients_all,
  file.path(sensitivity_table_dir, "coefficients_all_five_models.csv")
)

readr::write_csv(
  metrics_all,
  file.path(sensitivity_table_dir, "in_sample_metrics_all_five_models.csv")
)

readr::write_csv(
  diagnostics_all,
  file.path(sensitivity_diag_dir, "sampling_diagnostics_all_five_models.csv")
)

saveRDS(
  all_draw_matrices,
  file.path(sensitivity_model_dir, "all_five_model_draw_matrices.rds")
)


# =============================================================================
# 10. PSIS-LOO COMPARISON
# =============================================================================

saveRDS(
  all_loo,
  file.path(sensitivity_model_dir, "psis_loo_objects_all_five_models.rds")
)

# Validate and retain the model IDs explicitly. The names on all_loo are
# the canonical IDs used throughout this script; do not infer them from
# labels or rely on positional matching after model ranking.
if (is.null(names(all_loo)) || any(names(all_loo) == "")) {
  stop("The saved PSIS-LOO objects do not have model IDs as list names.")
}

if (!setequal(names(all_loo), names(model_specs))) {
  stop(
    "The PSIS-LOO model IDs do not match model_specs.\n",
    "LOO IDs: ", paste(names(all_loo), collapse = ", "), "\n",
    "Expected IDs: ", paste(names(model_specs), collapse = ", ")
  )
}

model_labels <- vapply(
  model_specs,
  function(spec) spec$label,
  character(1)
)

# loo_compare() returns rows ranked from best to worst, but depending on
# the loo package version, row names may be positional (1, 2, ...) rather
# than the names of the input list. Therefore, map rows using each model's
# elpd_loo ranking instead of relying on row names.
loo_comparison <- loo::loo_compare(all_loo)
loo_comparison_df <- as.data.frame(loo_comparison)

# Canonical model ranking: highest elpd_loo first, matching loo_compare().
model_elpd <- vapply(
  all_loo,
  function(x) unname(x$estimates["elpd_loo", "Estimate"]),
  numeric(1)
)
model_ids_ranked <- names(sort(model_elpd, decreasing = TRUE))

if (nrow(loo_comparison_df) != length(model_ids_ranked)) {
  stop(
    "The number of loo_compare() rows does not match the number of fitted models.\n",
    "Rows returned: ", nrow(loo_comparison_df), "\n",
    "Models available: ", length(model_ids_ranked)
  )
}

# Validate the row mapping against loo_compare()'s elpd_diff column.
# loo_compare() reports ELPD(model) - ELPD(best), so differences are
# zero for the best model and negative for models with lower ELPD.
if ("elpd_diff" %in% names(loo_comparison_df)) {
  # loo_compare() defines elpd_diff as the model's ELPD minus the
  # best model's ELPD. Thus the best model is 0 and poorer models are
  # negative; the sign here must not be reversed.
  expected_elpd_diff <- model_elpd[model_ids_ranked] - max(model_elpd)
  observed_elpd_diff <- loo_comparison_df$elpd_diff
  
  if (!isTRUE(all.equal(
    as.numeric(observed_elpd_diff),
    as.numeric(expected_elpd_diff),
    tolerance = 1e-6,
    check.attributes = FALSE
  ))) {
    stop(
      "Could not safely map loo_compare() rows to models using elpd_loo.\n",
      "Observed elpd_diff: ", paste(round(observed_elpd_diff, 6), collapse = ", "), "\n",
      "Expected elpd_diff: ", paste(round(expected_elpd_diff, 6), collapse = ", ")
    )
  }
}

loo_comparison_df$model_id <- model_ids_ranked

loo_comparison_tbl <- loo_comparison_df %>%
  mutate(
    model = unname(model_labels[model_id])
  ) %>%
  relocate(model_id, model)

readr::write_csv(
  loo_comparison_tbl,
  file.path(sensitivity_table_dir, "psis_loo_model_comparison.csv")
)

# Extract ELPD and Pareto-k diagnostics for each model.
loo_summary_tbl <- purrr::imap_dfr(
  all_loo,
  function(x, model_id) {
    
    if (!model_id %in% names(model_labels)) {
      stop("Unknown model ID in all_loo: ", model_id)
    }
    
    pointwise <- loo::pareto_k_values(x)
    
    tibble(
      model_id = model_id,
      model = unname(model_labels[model_id]),
      elpd_loo = x$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = x$estimates["elpd_loo", "SE"],
      p_loo = x$estimates["p_loo", "Estimate"],
      se_p_loo = x$estimates["p_loo", "SE"],
      looic = x$estimates["looic", "Estimate"],
      se_looic = x$estimates["looic", "SE"],
      max_pareto_k = max(pointwise),
      n_pareto_k_above_0.7 = sum(pointwise > 0.7),
      n_pareto_k_above_1 = sum(pointwise > 1)
    )
  }
)

readr::write_csv(
  loo_summary_tbl,
  file.path(sensitivity_table_dir, "psis_loo_diagnostics.csv")
)


# =============================================================================
# 11. PLOT OBSERVED VERSUS PREDICTED VALUES
# =============================================================================

p_prediction <- ggplot(
  predictions_all,
  aes(
    x = observed,
    y = predicted_mean,
    colour = model
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ model,
    ncol = 2
  ) +
  coord_equal() +
  labs(
    title = expression("Observed versus predicted " * delta^{14} * "C in fPOM"),
    x = expression("Observed " * delta^{14} * "C"),
    y = expression("Predicted " * delta^{14} * "C")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(sensitivity_plot_dir, "observed_vs_predicted_five_models.png"),
  p_prediction,
  width = 9,
  height = 7,
  dpi = 320
)


# =============================================================================
# 12. PLOT POSTERIOR COEFFICIENTS
# =============================================================================

coefficient_plot_data <- coefficients_all %>%
  filter(term != "Intercept") %>%
  mutate(
    term = factor(
      term,
      levels = rev(best_predictors)
    )
  )

p_coefficients <- ggplot(
  coefficient_plot_data,
  aes(
    x = posterior_mean,
    y = term,
    xmin = lower_95,
    xmax = upper_95,
    colour = model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.65)
  ) +
  facet_wrap(
    ~ model,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Posterior regression coefficients",
    subtitle = "Points = posterior means; lines = 95% credible intervals",
    x = "Standardised regression coefficient",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(sensitivity_plot_dir, "posterior_coefficients_five_models.png"),
  p_coefficients,
  width = 9,
  height = 10,
  dpi = 320
)


# =============================================================================
# 13. FINAL SUMMARY
# =============================================================================

final_summary <- metrics_all %>%
  left_join(
    loo_summary_tbl %>%
      select(
        model_id,
        elpd_loo,
        se_elpd_loo,
        max_pareto_k,
        n_pareto_k_above_0.7,
        n_pareto_k_above_1
      ),
    by = "model_id"
  ) %>%
  left_join(
    diagnostics_all %>%
      select(
        model_id,
        max_Rhat,
        min_bulk_ESS,
        min_tail_ESS,
        divergences,
        max_treedepth_hits
      ),
    by = "model_id"
  )

readr::write_csv(
  final_summary,
  file.path(sensitivity_table_dir, "final_five_model_comparison.csv")
)

message("\n============================================================")
message("FIVE-MODEL d14C_fPOM SENSITIVITY ANALYSIS COMPLETED")
message("============================================================")

message("\nFinal comparison:")
print(final_summary)

message("\nPSIS-LOO comparison:")
print(loo_comparison)

message("\nResults saved to:")
message(sensitivity_dir)
