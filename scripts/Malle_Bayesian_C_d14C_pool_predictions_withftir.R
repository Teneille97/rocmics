
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
