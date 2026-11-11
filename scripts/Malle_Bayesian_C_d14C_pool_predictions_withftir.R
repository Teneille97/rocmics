
# =============================================================================
# Malle ERW: Bayesian prediction of C and delta-14C for Bulk, MAOM, oPOM, fPOM
#
# Purpose:
#   Fit one regularised Bayesian multiple-regression model for each available
#   response: C_Bulk, C_MAOM, C_oPOM, C_fPOM and corresponding d14C_2026_*.
#
#   Save:
#     - Model objects and posterior draws
#     - Posterior predictions and coefficient summaries
#     - Model diagnostics and posterior predictive checks
#     - projpred predictor-selection results
#     - Observed-versus-predicted plots with 95% uncertainty intervals
#     - Relative predictor-importance bar plots underneath each prediction plot
#
# Important:
#   * Predictors are standardised; responses are also standardised for fitting.
#   * Regularising Normal priors reduce coefficient overfitting.
#   * Predictions are in-sample fitted predictions, not cross-validation.
#   * Relative importance = abs(standardised coefficient) / sum(abs(coefficients))
#   * Relative importance is not variance explained or causal contribution.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(posterior)
  library(bayesplot)
  library(patchwork)
})

set.seed(20261028)

# -----------------------------------------------------------------------------
# 0. Settings and paths
# -----------------------------------------------------------------------------

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

for (d in c(
  out_dir, model_dir, plot_dir, table_dir, diag_dir,
  selection_dir, selection_plot_dir, selection_table_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Package 'cmdstanr' is required. Install it and configure CmdStan first.")
}

if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop(
    "CmdStan is not configured. Check cmdstanr::cmdstan_path() ",
    "and the C++ toolchain."
  )
}

if (!requireNamespace("projpred", quietly = TRUE)) {
  stop(
    "Package 'projpred' is required. ",
    "Install it with install.packages('projpred') and rerun."
  )
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "Package 'patchwork' is required for the combined prediction ",
    "and relative-importance plots. Install it using ",
    "install.packages('patchwork')."
  )
}

message(
  "Using CmdStan ", as.character(cmdstanr::cmdstan_version()),
  " via cmdstanr ", as.character(utils::packageVersion("cmdstanr")),
  "; brms/rstan are not used."
)

message(
  "Using projpred ",
  as.character(utils::packageVersion("projpred")),
  " for projection predictive predictor selection."
)

# Sampling settings
n_chains <- 4L
n_warmup <- 1500L
n_sampling <- 2000L
n_cores <- min(n_chains, parallel::detectCores(logical = FALSE))
adapt_delta <- 0.99
max_treedepth <- 12L

# Predictor filtering
min_predictor_observed_fraction <- 0.70
max_abs_predictor_correlation <- 0.95

# Model priors on standardised variables
slope_prior_sd <- 0.5


# -----------------------------------------------------------------------------
# 1. Read data and validate key columns
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Cannot find input file: ", input_file,
    "\nExpected outputs/integrated_geochemistry_for_bayesian.csv."
  )
}

integrated <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  na = c("", "NA", "#VALUE!")
)

if (!"sample" %in% names(integrated)) {
  stop("The integrated CSV must contain a 'sample' plot identifier column.")
}

if (!"Treatment" %in% names(integrated)) {
  stop("The integrated CSV must contain a 'Treatment' column.")
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
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(as.integer(as.character(sample))),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(
    !is.na(sample),
    !is.na(Treatment_model)
  )

if (nrow(dat) == 0L) {
  stop(
    "No rows remain after filtering to the five expected treatments. ",
    "Check the Treatment labels in the integrated CSV."
  )
}

readr::write_csv(
  dat %>% count(Treatment, Treatment_model, name = "n_plots"),
  file.path(table_dir, "treatment_sample_counts.csv")
)


# -----------------------------------------------------------------------------
# 2. Define response variables
# -----------------------------------------------------------------------------

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

responses <- intersect(
  c(c_responses, d14c_responses),
  names(dat)
)

if (length(responses) == 0L) {
  stop(
    "None of the expected response columns were found. Expected C_Bulk, ",
    "C_MAOM, C_oPOM, C_fPOM and/or d14C_2026_Bulk/MAOM/oPOM/fPOM."
  )
}

missing_responses <- setdiff(
  c(c_responses, d14c_responses),
  names(dat)
)

if (length(missing_responses) > 0L) {
  message(
    "Response columns absent from the CSV (will be skipped): ",
    paste(missing_responses, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# 3. Define candidate predictors and prevent response leakage
# -----------------------------------------------------------------------------

explicit_geochem_predictors <- c(
  
  # Aggregate stability: H2O indices
  "GMD_h2o",
  "D_h2o",
  "MWD_h2o",
  
  # Aggregate stability: oxalate:H2O ratios
  "GMD_oxalate_over_h2o",
  "D_oxalate_over_h2o",
  "MWD_oxalate_over_h2o",
  
  # Oxalate-extractable mineral ratios
  "Fe_ox_over_DCB",
  "Al_ox_over_DCB",
  "Mn_ox_over_DCB",
  "Si_ox_over_DCB",
  
  # Secondary-mineral ratios
  "Fe_ox_over_Si_ox",
  "Al_ox_over_Si_ox",
  
  # Tessier exchangeable fractions
  "Tessier_Acetic_Ca_mg_g",
  "Tessier_Acetic_K_mg_g",
  "Tessier_Acetic_Na_mg_g",
  "Tessier_Acetic_Mg_mg_g",
  
  # Tessier reducible fractions
  "Tessier_Hydroxylamine_Fe_mg_g",
  "Tessier_Hydroxylamine_Al_mg_g",
  "Tessier_Hydroxylamine_Mn_mg_g",
  "Tessier_Hydroxylamine_Si_mg_g",
  
  # February 2026 soil chemistry
  "soil_pH",
  "soil_EC",
  
  # October 2023 soil chemistry and exchangeable cations
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


# Prevent carbon and radiocarbon outcomes from entering as predictors.
response_like <- names(dat)[
  stringr::str_detect(
    names(dat),
    "^C_(Bulk|MAOM|oPOM|fPOM)$"
  ) |
    stringr::str_detect(names(dat), "^d14C")
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
  stop(
    "No numeric candidate geochemical predictors were found. ",
    "Check column names in the integrated CSV."
  )
}


# Check availability of updated predictors.
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

updated_predictor_check <- tibble(
  predictor = updated_predictors,
  available_in_dataset = updated_predictors %in% names(dat),
  included_as_candidate = updated_predictors %in% candidate_predictors
)

readr::write_csv(
  updated_predictor_check,
  file.path(table_dir, "updated_predictor_availability_check.csv")
)

message("\nUpdated aggregate-stability and mineral-ratio predictors:")
print(updated_predictor_check)


# Report detected FTIR predictors.
ftir_predictors <- candidate_predictors[
  stringr::str_starts(candidate_predictors, "FTIR_delta_")
]

message("\nFTIR delta predictors detected: ", length(ftir_predictors))

if (length(ftir_predictors) > 0L) {
  print(ftir_predictors)
} else {
  warning(
    "No numeric FTIR_delta_* predictors were detected. ",
    "Check the integrated CSV and FTIR column names."
  )
}


readr::write_csv(
  tibble(candidate_predictor = candidate_predictors),
  file.path(
    table_dir,
    "candidate_predictors_before_response_filtering.csv"
  )
)

message(
  "\nTotal candidate geochemical predictors available: ",
  length(candidate_predictors)
)

print(candidate_predictors)


# -----------------------------------------------------------------------------
# 4. Helper functions
# -----------------------------------------------------------------------------

# Prepare response-specific model data:
# - Filter finite response observations.
# - Filter predictors based on missingness and variation.
# - Median-impute missing predictor values.
# - Standardise predictors and response.
# - Remove extreme predictor-predictor collinearity.

prepare_model_data <- function(data, response, candidates) {
  
  y_raw <- suppressWarnings(as.numeric(data[[response]]))
  keep_y <- is.finite(y_raw)
  
  work <- data[keep_y, , drop = FALSE]
  y_raw <- y_raw[keep_y]
  
  if (length(y_raw) < 8L) {
    stop(
      response, ": only ", length(y_raw),
      " finite response observations; at least 8 are required."
    )
  }
  
  eligible <- candidates[
    vapply(candidates, function(v) {
      
      x <- suppressWarnings(as.numeric(work[[v]]))
      observed_fraction <- mean(is.finite(x))
      sx <- stats::sd(x, na.rm = TRUE)
      
      is.finite(observed_fraction) &&
        observed_fraction >= min_predictor_observed_fraction &&
        is.finite(sx) &&
        sx > 0
      
    }, logical(1))
  ]
  
  if (length(eligible) == 0L) {
    stop(
      response,
      ": no eligible predictors after missingness/variance filtering."
    )
  }
  
  X <- as.data.frame(lapply(eligible, function(v) {
    
    x <- suppressWarnings(as.numeric(work[[v]]))
    med <- stats::median(x[is.finite(x)], na.rm = TRUE)
    x[!is.finite(x)] <- med
    
    as.numeric(scale(x))
    
  }))
  
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
    stop(response, ": no variable predictors remain after standardisation.")
  }
  
  # Thin extreme collinearity.
  dropped_correlated <- character(0)
  
  repeat {
    
    if (ncol(X) < 2L) break
    
    cm <- suppressWarnings(
      stats::cor(X, use = "pairwise.complete.obs")
    )
    
    diag(cm) <- 0
    
    max_cor <- max(abs(cm), na.rm = TRUE)
    
    if (
      !is.finite(max_cor) ||
      max_cor <= max_abs_predictor_correlation
    ) {
      break
    }
    
    idx <- which(abs(cm) == max_cor, arr.ind = TRUE)[1, ]
    drop_name <- colnames(cm)[idx[2]]
    
    dropped_correlated <- c(dropped_correlated, drop_name)
    
    X <- X[, setdiff(names(X), drop_name), drop = FALSE]
  }
  
  # Standardise response.
  y_mean <- mean(y_raw)
  y_sd <- stats::sd(y_raw)
  
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop(response, ": response has no finite variation.")
  }
  
  y_std <- (y_raw - y_mean) / y_sd
  
  model_data <- as_tibble(X)
  model_data$y_std <- y_std
  model_data$sample <- work$sample
  model_data$Treatment_model <- factor(work$Treatment_model)
  
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


safe_name <- function(x) {
  
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}


# Generate fitted-mean and posterior-predictive summaries.
make_predictions <- function(draw_matrix, model_data, prep, response) {
  
  mu_cols <- grep("^mu\\[", colnames(draw_matrix), value = TRUE)
  rep_cols <- grep("^y_rep\\[", colnames(draw_matrix), value = TRUE)
  
  if (
    length(mu_cols) != nrow(model_data) ||
    length(rep_cols) != nrow(model_data)
  ) {
    stop("Could not extract all mu/y_rep draws for ", response)
  }
  
  mu_draws <- draw_matrix[, mu_cols, drop = FALSE]
  rep_draws <- draw_matrix[, rep_cols, drop = FALSE]
  
  pred_mean <- apply(mu_draws, 2, mean)
  pred_median <- apply(mu_draws, 2, median)
  
  pred_lower <- apply(
    mu_draws, 2, stats::quantile, probs = 0.025
  )
  
  pred_upper <- apply(
    mu_draws, 2, stats::quantile, probs = 0.975
  )
  
  pred_lower_pi <- apply(
    rep_draws, 2, stats::quantile, probs = 0.025
  )
  
  pred_upper_pi <- apply(
    rep_draws, 2, stats::quantile, probs = 0.975
  )
  
  tibble(
    response = response,
    sample = model_data$sample,
    Treatment_model = as.character(model_data$Treatment_model),
    observed = prep$y_raw,
    predicted_mean = pred_mean * prep$y_sd + prep$y_mean,
    predicted_median = pred_median * prep$y_sd + prep$y_mean,
    fitted_lower_95 = pred_lower * prep$y_sd + prep$y_mean,
    fitted_upper_95 = pred_upper * prep$y_sd + prep$y_mean,
    predictive_lower_95 = pred_lower_pi * prep$y_sd + prep$y_mean,
    predictive_upper_95 = pred_upper_pi * prep$y_sd + prep$y_mean
  )
}


# -----------------------------------------------------------------------------
# 5. Relative predictor importance
# -----------------------------------------------------------------------------

# The coefficient CSVs use:
#   term    = predictor name
#   Estimate = posterior mean coefficient on standardised scales
#
# The intercept is excluded.
#
# Relative importance:
#   abs(Estimate) / sum(abs(Estimate)) * 100
#
# This gives percentages that sum to 100% within each model.

calculate_relative_importance <- function(coef_tbl, response) {
  
  predictor_col <- "term"
  coefficient_col <- "Estimate"
  
  if (!predictor_col %in% names(coef_tbl)) {
    stop(
      "Could not identify predictor column for ", response,
      ". Available columns: ",
      paste(names(coef_tbl), collapse = ", ")
    )
  }
  
  if (!coefficient_col %in% names(coef_tbl)) {
    stop(
      "Could not identify standardised coefficient column for ", response,
      ". Available columns: ",
      paste(names(coef_tbl), collapse = ", ")
    )
  }
  
  importance_tbl <- coef_tbl %>%
    transmute(
      Predictor = as.character(.data[[predictor_col]]),
      Standardised_coefficient = as.numeric(.data[[coefficient_col]])
    ) %>%
    filter(
      !is.na(Predictor),
      !is.na(Standardised_coefficient)
    ) %>%
    filter(
      !stringr::str_detect(
        tolower(Predictor),
        "intercept|^alpha$|^sigma$|residual|^lp__$"
      )
    )
  
  if (nrow(importance_tbl) == 0L) {
    stop(
      "No predictor coefficients available for relative importance: ",
      response
    )
  }
  
  total_absolute_coefficient <- sum(
    abs(importance_tbl$Standardised_coefficient)
  )
  
  if (
    !is.finite(total_absolute_coefficient) ||
    total_absolute_coefficient == 0
  ) {
    stop(
      "Cannot calculate relative importance for ", response,
      ": total absolute standardised coefficient is zero or invalid."
    )
  }
  
  importance_tbl %>%
    mutate(
      Relative_importance =
        abs(Standardised_coefficient) /
        total_absolute_coefficient * 100
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
# 6. Combined observed-vs-predicted and relative-importance plot
# -----------------------------------------------------------------------------

plot_observed_predicted <- function(
    pred,
    coef_tbl,
    response,
    output_path,
    model_label = "full reference model") {
  
  # ---------------------------------------------------------------------------
  # Upper panel: observed versus predicted
  # ---------------------------------------------------------------------------
  
  pdat <- pred %>%
    arrange(observed) %>%
    mutate(row_order = row_number())
  
  p_obs <- ggplot(
    pdat,
    aes(x = observed, y = predicted_mean)
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
  
  
  # ---------------------------------------------------------------------------
  # Lower panel: relative predictor importance
  # ---------------------------------------------------------------------------
  
  importance_tbl <- calculate_relative_importance(
    coef_tbl = coef_tbl,
    response = response
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
        label = sprintf("%.1f%%", Relative_importance)
      ),
      hjust = -0.15,
      size = 3.0
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
        "Absolute standardised coefficient as a percentage ",
        "of the total absolute coefficient magnitude"
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
  
  
  # ---------------------------------------------------------------------------
  # Combine panels
  # ---------------------------------------------------------------------------
  
  combined_plot <- p_obs / p_importance +
    patchwork::plot_layout(
      heights = c(
        3,
        max(1.1, 0.20 * n_predictors_plot)
      )
    ) +
    patchwork::plot_annotation(
      caption = paste0(
        "Relative importance is based on absolute posterior mean ",
        "standardised coefficients. Percentages sum to 100% within ",
        "each model; they do not represent variance explained."
      ),
      theme = theme(
        plot.caption = element_text(
          size = 8,
          colour = "grey35",
          hjust = 0
        )
      )
    )
  
  # Increase height for models with many predictors.
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
  
  return(combined_plot)
}


# -----------------------------------------------------------------------------
# 7. Define Stan model and compile
# -----------------------------------------------------------------------------

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
      "Could not compile the direct CmdStan model. ",
      "Check the C++ toolchain. Details: ",
      conditionMessage(e)
    )
  }
)


# -----------------------------------------------------------------------------
# 8. Fit models
# -----------------------------------------------------------------------------

all_predictions <- list()
all_metrics <- list()
all_coefficients <- list()
all_model_info <- list()
all_draws <- list()
all_projpred_info <- list()


for (response in responses) {
  
  message("\n============================================================")
  message("Fitting direct-cmdstanr Bayesian model for: ", response)
  message("============================================================")
  
  prep <- tryCatch(
    prepare_model_data(
      dat,
      response,
      candidate_predictors
    ),
    error = function(e) {
      message(
        "Skipping ", response, ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(prep)) next
  
  message(
    "n = ", prep$n,
    "; predictors = ", length(prep$predictors)
  )
  
  message(
    "Predictors: ",
    paste(prep$predictors, collapse = ", ")
  )
  
  if (length(prep$dropped_correlated) > 0L) {
    message(
      "Dropped for |r| > ",
      max_abs_predictor_correlation,
      ": ",
      paste(prep$dropped_correlated, collapse = ", ")
    )
  }
  
  model_id <- safe_name(response)
  
  X <- as.matrix(
    prep$data[, prep$predictors, drop = FALSE]
  )
  
  storage.mode(X) <- "double"
  
  stan_data <- list(
    N = nrow(X),
    K = ncol(X),
    X = X,
    y = as.numeric(prep$data$y_std)
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
  
  draws_file <- file.path(
    model_dir,
    paste0("draws_", model_id, ".rds")
  )
  
  fit <- NULL
  
  # Reuse saved full-model posterior draws when available.
  if (file.exists(draws_file)) {
    
    message(
      "Reusing saved full-model posterior draws for ",
      response,
      "; skipping the full-model MCMC refit."
    )
    
    draws_array <- readRDS(draws_file)
    
  } else {
    
    fit <- tryCatch(
      stan_model$sample(
        data = stan_data,
        chains = n_chains,
        parallel_chains = n_cores,
        iter_warmup = n_warmup,
        iter_sampling = n_sampling,
        seed = 20261028 + match(response, responses),
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = 200,
        output_dir = csv_dir
      ),
      error = function(e) {
        message(
          "Model failed for ", response, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(fit)) next
    
    draws_array <- fit$draws(
      format = "draws_array"
    )
    
    saveRDS(
      draws_array,
      draws_file
    )
  }
  
  all_draws[[response]] <- draws_array
  
  draw_matrix <- posterior::as_draws_matrix(
    draws_array
  )
  
  
  # ---------------------------------------------------------------------------
  # 8A. Projection predictive variable selection
  # ---------------------------------------------------------------------------
  
  tryCatch({
    
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
        newdata[, fit$predictors, drop = FALSE]
      )
      
      storage.mode(X_new) <- "double"
      
      beta_names <- paste0(
        "beta[",
        seq_along(fit$predictors),
        "]"
      )
      
      beta_draws <- as.matrix(
        fit$draw_matrix[, beta_names, drop = FALSE]
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
        paste(prep$predictors, collapse = " + ")
      )
    )
    
    ref_obj <- projpred::init_refmodel(
      object = ref_fit,
      data = as.data.frame(prep$data),
      formula = full_formula,
      family = stats::gaussian(),
      ref_predfun = ref_predfun,
      dis = as.numeric(draw_matrix[, "sigma"])
    )
    
    
    max_terms <- min(
      5L,
      length(prep$predictors)
    )
    
    vsel <- projpred::varsel(
      ref_obj,
      method = "forward",
      nterms_max = max_terms,
      nclusters = min(20L, nrow(draw_matrix)),
      nclusters_pred = min(20L, nrow(draw_matrix)),
      seed = 20261028 + match(response, responses),
      verbose = 1
    )
    
    saveRDS(
      vsel,
      file.path(
        selection_dir,
        paste0("projpred_", model_id, ".rds")
      )
    )
    
    capture.output(
      print(summary(vsel)),
      file = file.path(
        selection_table_dir,
        paste0("summary_", model_id, ".txt")
      )
    )
    
    capture.output(
      print(projpred::ranking(vsel)),
      file = file.path(
        selection_table_dir,
        paste0("ranking_", model_id, ".txt")
      )
    )
    
    p_vsel <- plot(vsel) +
      ggplot2::ggtitle(
        paste("projpred submodel performance:", response)
      ) +
      ggplot2::theme_bw(base_size = 12)
    
    ggplot2::ggsave(
      file.path(
        selection_plot_dir,
        paste0("projpred_performance_", model_id, ".png")
      ),
      p_vsel,
      width = 7.5,
      height = 5.5,
      dpi = 300
    )
    
    
    # -------------------------------------------------------------------------
    # 8B. Project and refit 3-, 4- and 5-predictor submodels
    # -------------------------------------------------------------------------
    
    for (k in intersect(
      c(3L, 4L, 5L),
      seq_len(max_terms)
    )) {
      
      prj_k <- projpred::project(
        vsel,
        nterms = k,
        ndraws = min(400L, nrow(draw_matrix)),
        seed = 20261028 +
          100L * match(response, responses) + k,
        verbose = 0
      )
      
      selected_k <- as.character(
        projpred::predictor_terms(prj_k)
      )
      
      selected_k <- intersect(
        selected_k,
        prep$predictors
      )
      
      if (length(selected_k) != k) {
        stop(
          "projpred returned ",
          length(selected_k),
          " unique usable predictors for the ",
          k,
          "-predictor model of ",
          response,
          ". Returned: ",
          paste(selected_k, collapse = ", ")
        )
      }
      
      all_projpred_info[[
        paste(response, k, sep = "__")
      ]] <- tibble::tibble(
        response = response,
        n_predictors = k,
        selection_order = seq_along(selected_k),
        predictor = selected_k
      )
      
      saveRDS(
        prj_k,
        file.path(
          selection_dir,
          paste0(
            "projection_",
            model_id,
            "_",
            k,
            "predictors.rds"
          )
        )
      )
      
      
      # Fit a separate CmdStan regression using selected predictors.
      sub_id <- paste0(
        model_id,
        "__projpred_",
        k,
        "predictors"
      )
      
      X_sub <- as.matrix(
        prep$data[, selected_k, drop = FALSE]
      )
      
      storage.mode(X_sub) <- "double"
      
      sub_stan_data <- list(
        N = nrow(X_sub),
        K = ncol(X_sub),
        X = X_sub,
        y = as.numeric(prep$data$y_std)
      )
      
      sub_csv_dir <- file.path(
        model_dir,
        paste0("cmdstan_csv_", sub_id)
      )
      
      dir.create(
        sub_csv_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      message(
        "Fitting selected ",
        k,
        "-predictor submodel for ",
        response,
        ": ",
        paste(selected_k, collapse = " + ")
      )
      
      sub_fit <- stan_model$sample(
        data = sub_stan_data,
        chains = n_chains,
        parallel_chains = n_cores,
        iter_warmup = n_warmup,
        iter_sampling = n_sampling,
        seed = 20261028 +
          1000L * match(response, responses) + k,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = 200,
        output_dir = sub_csv_dir
      )
      
      saveRDS(
        sub_fit,
        file.path(
          model_dir,
          paste0("fit_", sub_id, ".rds")
        )
      )
      
      sub_draws_array <- sub_fit$draws(
        format = "draws_array"
      )
      
      saveRDS(
        sub_draws_array,
        file.path(
          model_dir,
          paste0("draws_", sub_id, ".rds")
        )
      )
      
      sub_draw_matrix <- posterior::as_draws_matrix(
        sub_draws_array
      )
      
      all_draws[[sub_id]] <- sub_draws_array
      
      
      # -----------------------------------------------------------------------
      # Calculate predictions
      # -----------------------------------------------------------------------
      
      sub_prep <- prep
      sub_prep$predictors <- selected_k
      
      sub_data <- prep$data[
        ,
        c(
          selected_k,
          "y_std",
          "sample",
          "Treatment_model"
        ),
        drop = FALSE
      ]
      
      sub_pred <- make_predictions(
        sub_draw_matrix,
        sub_data,
        sub_prep,
        response
      ) %>%
        mutate(
          model = paste0(k, "-predictor refit"),
          n_predictors = k,
          predictors = paste(selected_k, collapse = " + "),
          .before = 1
        )
      
      all_predictions[[sub_id]] <- sub_pred
      
      readr::write_csv(
        sub_pred,
        file.path(
          table_dir,
          paste0("predictions_", sub_id, ".csv")
        )
      )
      
      
      # -----------------------------------------------------------------------
      # Calculate coefficient summaries BEFORE plotting
      # -----------------------------------------------------------------------
      
      sub_coef_names <- c(
        "alpha",
        paste0("beta[", seq_along(selected_k), "]")
      )
      
      sub_coef_labels <- c(
        "Intercept",
        selected_k
      )
      
      sub_coef <- purrr::map_dfr(
        seq_along(sub_coef_names),
        function(j) {
          
          vals <- sub_draw_matrix[
            ,
            sub_coef_names[j]
          ]
          
          tibble(
            term = sub_coef_labels[j],
            Estimate = mean(vals),
            Est.Error = stats::sd(vals),
            Q2.5 = unname(
              stats::quantile(vals, 0.025)
            ),
            Q50 = unname(
              stats::quantile(vals, 0.5)
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
          model = paste0(k, "-predictor refit"),
          n_predictors = k,
          predictors = paste(selected_k, collapse = " + "),
          estimate_original_response_per_predictor_SD =
            Estimate * prep$y_sd,
          lower_original_response_per_predictor_SD =
            Q2.5 * prep$y_sd,
          upper_original_response_per_predictor_SD =
            Q97.5 * prep$y_sd
        )
      
      all_coefficients[[sub_id]] <- sub_coef
      
      readr::write_csv(
        sub_coef,
        file.path(
          table_dir,
          paste0("coefficients_", sub_id, ".csv")
        )
      )
      
      
      # -----------------------------------------------------------------------
      # Combined observed-vs-predicted and relative-importance plot
      # -----------------------------------------------------------------------
      
      plot_observed_predicted(
        pred = sub_pred,
        coef_tbl = sub_coef,
        response = response,
        output_path = file.path(
          plot_dir,
          paste0("observed_vs_predicted_", sub_id, ".png")
        ),
        model_label = paste0(
          "projpred-selected ",
          k,
          "-predictor refit"
        )
      )
      
      
      # -----------------------------------------------------------------------
      # Diagnostics
      # -----------------------------------------------------------------------
      
      sub_summary <- sub_fit$summary(
        variables = c("alpha", "beta", "sigma")
      )
      
      sub_sampler <- as.data.frame(
        sub_fit$sampler_diagnostics(
          format = "draws_df"
        )
      )
      
      sub_diag <- tibble(
        response = response,
        model = paste0(k, "-predictor refit"),
        n = prep$n,
        n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        max_Rhat = if (
          all(is.na(sub_summary$rhat))
        ) {
          NA_real_
        } else {
          max(sub_summary$rhat, na.rm = TRUE)
        },
        min_bulk_ESS = if (
          all(is.na(sub_summary$ess_bulk))
        ) {
          NA_real_
        } else {
          min(sub_summary$ess_bulk, na.rm = TRUE)
        },
        min_tail_ESS = if (
          all(is.na(sub_summary$ess_tail))
        ) {
          NA_real_
        } else {
          min(sub_summary$ess_tail, na.rm = TRUE)
        },
        divergences = sum(
          sub_sampler$divergent__ == 1,
          na.rm = TRUE
        ),
        treedepth_hits = sum(
          sub_sampler$treedepth__ >= max_treedepth,
          na.rm = TRUE
        )
      )
      
      readr::write_csv(
        sub_diag,
        file.path(
          diag_dir,
          paste0("diagnostics_", sub_id, ".csv")
        )
      )
      
      capture.output(
        sub_fit$summary(),
        file = file.path(
          diag_dir,
          paste0("summary_", sub_id, ".txt")
        )
      )
      
      
      # Posterior predictive check
      sub_rep_cols <- grep(
        "^y_rep\\[",
        colnames(sub_draw_matrix),
        value = TRUE
      )
      
      sub_yrep <- sub_draw_matrix[
        seq_len(min(100L, nrow(sub_draw_matrix))),
        sub_rep_cols,
        drop = FALSE
      ]
      
      sub_pp <- bayesplot::ppc_dens_overlay(
        y = as.numeric(prep$data$y_std),
        yrep = sub_yrep
      ) +
        ggplot2::ggtitle(
          paste(
            "Posterior predictive check:",
            response,
            "-",
            k,
            "predictors"
          )
        ) +
        ggplot2::theme_bw(base_size = 12)
      
      ggplot2::ggsave(
        file.path(
          diag_dir,
          paste0("pp_check_", sub_id, ".png")
        ),
        sub_pp,
        width = 7,
        height = 5.5,
        dpi = 300
      )
      
      
      # In-sample performance metrics
      sub_rmse <- sqrt(
        mean(
          (sub_pred$observed - sub_pred$predicted_mean)^2,
          na.rm = TRUE
        )
      )
      
      sub_mae <- mean(
        abs(sub_pred$observed - sub_pred$predicted_mean),
        na.rm = TRUE
      )
      
      sub_r2 <- suppressWarnings(
        stats::cor(
          sub_pred$observed,
          sub_pred$predicted_mean,
          use = "complete.obs"
        )^2
      )
      
      sub_cov_mean <- mean(
        sub_pred$observed >= sub_pred$fitted_lower_95 &
          sub_pred$observed <= sub_pred$fitted_upper_95,
        na.rm = TRUE
      )
      
      sub_cov_pred <- mean(
        sub_pred$observed >= sub_pred$predictive_lower_95 &
          sub_pred$observed <= sub_pred$predictive_upper_95,
        na.rm = TRUE
      )
      
      all_metrics[[sub_id]] <- tibble(
        response = response,
        model = paste0(k, "-predictor refit"),
        n = prep$n,
        n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        RMSE_in_sample = sub_rmse,
        MAE_in_sample = sub_mae,
        squared_correlation_in_sample = sub_r2,
        observed_coverage_fitted_mean_95 = sub_cov_mean,
        observed_coverage_predictive_95 = sub_cov_pred
      )
      
      all_model_info[[sub_id]] <- tibble(
        response = response,
        model = paste0(k, "-predictor refit"),
        n = prep$n,
        n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        dropped_correlated = paste(
          prep$dropped_correlated,
          collapse = " + "
        )
      )
      
    }
    
    message(
      "projpred selection completed for ",
      response,
      "; inspect projpred_selection/tables and projpred_selection/plots."
    )
    
  }, error = function(e) {
    
    message(
      "projpred selection failed for ",
      response,
      ": ",
      conditionMessage(e),
      "\nThe full CmdStan model and its predictions will still be saved."
    )
    
  })
  
  
  # ---------------------------------------------------------------------------
  # 8C. Full reference model diagnostics
  # ---------------------------------------------------------------------------
  
  diag_summary <- if (!is.null(fit)) {
    fit$summary()
  } else {
    posterior::summarise_draws(
      posterior::as_draws_array(draws_array)
    )
  }
  
  capture.output(
    diag_summary,
    file = file.path(
      diag_dir,
      paste0("summary_", model_id, ".txt")
    )
  )
  
  readr::write_csv(
    diag_summary,
    file.path(
      diag_dir,
      paste0("posterior_draw_diagnostics_", model_id, ".csv")
    )
  )
  
  
  # Posterior predictive check
  rep_cols <- grep(
    "^y_rep\\[",
    colnames(draw_matrix),
    value = TRUE
  )
  
  n_pp <- min(100L, nrow(draw_matrix))
  
  yrep_check <- draw_matrix[
    seq_len(n_pp),
    rep_cols,
    drop = FALSE
  ]
  
  pp_plot <- bayesplot::ppc_dens_overlay(
    y = as.numeric(prep$data$y_std),
    yrep = yrep_check
  ) +
    ggtitle(
      paste("Posterior predictive check:", response)
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(
      diag_dir,
      paste0("pp_check_", model_id, ".png")
    ),
    pp_plot,
    width = 7,
    height = 5.5,
    dpi = 300
  )
  
  
  # ---------------------------------------------------------------------------
  # 8D. Full reference model coefficient summaries
  # ---------------------------------------------------------------------------
  
  coef_names <- c(
    "alpha",
    paste0("beta[", seq_along(prep$predictors), "]")
  )
  
  coef_labels <- c(
    "Intercept",
    prep$predictors
  )
  
  coef_summary <- lapply(
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
          stats::quantile(vals, 0.5)
        ),
        Q97.5 = unname(
          stats::quantile(vals, 0.975)
        )
      )
    }
  ) %>%
    bind_rows() %>%
    mutate(
      response = response,
      model = "full reference model",
      n_predictors = length(prep$predictors),
      predictors = paste(prep$predictors, collapse = " + "),
      estimate_original_response_per_predictor_SD =
        Estimate * prep$y_sd,
      lower_original_response_per_predictor_SD =
        Q2.5 * prep$y_sd,
      upper_original_response_per_predictor_SD =
        Q97.5 * prep$y_sd
    )
  
  all_coefficients[[response]] <- coef_summary
  
  readr::write_csv(
    coef_summary,
    file.path(
      table_dir,
      paste0("coefficients_", model_id, ".csv")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 8E. Full reference model predictions
  # ---------------------------------------------------------------------------
  
  pred <- make_predictions(
    draw_matrix,
    prep$data,
    prep,
    response
  ) %>%
    mutate(
      model = "full reference model",
      n_predictors = length(prep$predictors),
      predictors = paste(prep$predictors, collapse = " + "),
      .before = 1
    )
  
  all_predictions[[response]] <- pred
  
  readr::write_csv(
    pred,
    file.path(
      table_dir,
      paste0("predictions_", model_id, ".csv")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Full reference model combined plot
  # ---------------------------------------------------------------------------
  
  plot_observed_predicted(
    pred = pred,
    coef_tbl = coef_summary,
    response = response,
    output_path = file.path(
      plot_dir,
      paste0("observed_vs_predicted_", model_id, ".png")
    ),
    model_label = "full reference model"
  )
  
  
  # ---------------------------------------------------------------------------
  # Full reference model performance metrics
  # ---------------------------------------------------------------------------
  
  rmse <- sqrt(
    mean(
      (pred$observed - pred$predicted_mean)^2,
      na.rm = TRUE
    )
  )
  
  mae <- mean(
    abs(pred$observed - pred$predicted_mean),
    na.rm = TRUE
  )
  
  r2_corr <- suppressWarnings(
    cor(
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
  
  all_metrics[[response]] <- tibble(
    response = response,
    model = "full reference model",
    n = prep$n,
    n_predictors = length(prep$predictors),
    predictors = paste(prep$predictors, collapse = " + "),
    predictors_dropped_for_collinearity = paste(
      prep$dropped_correlated,
      collapse = " + "
    ),
    RMSE_in_sample = rmse,
    MAE_in_sample = mae,
    squared_correlation_in_sample = r2_corr,
    observed_coverage_fitted_mean_95 = coverage_mean,
    observed_coverage_predictive_95 = coverage_predictive
  )
  
  all_model_info[[response]] <- tibble(
    response = response,
    model = "full reference model",
    n = prep$n,
    predictors = paste(prep$predictors, collapse = " + "),
    dropped_correlated = paste(
      prep$dropped_correlated,
      collapse = " + "
    )
  )
  
}


# -----------------------------------------------------------------------------
# 9. Save combined outputs
# -----------------------------------------------------------------------------

if (length(all_predictions) == 0L) {
  stop(
    "No response models completed successfully. ",
    "Check messages and diagnostics."
  )
}

predictions_all <- bind_rows(all_predictions)
metrics_all <- bind_rows(all_metrics)
coefficients_all <- bind_rows(all_coefficients)
model_info_all <- bind_rows(all_model_info)

readr::write_csv(
  predictions_all,
  file.path(table_dir, "predictions_all_responses.csv")
)

readr::write_csv(
  metrics_all,
  file.path(table_dir, "model_performance_in_sample.csv")
)

readr::write_csv(
  coefficients_all,
  file.path(table_dir, "coefficients_all_responses.csv")
)

readr::write_csv(
  model_info_all,
  file.path(table_dir, "model_predictor_sets.csv")
)


# Save posterior draws.
saveRDS(
  all_draws,
  file.path(
    model_dir,
    "all_bayesian_C_d14C_posterior_draws.rds"
  )
)


# Save complete results bundle.
saveRDS(
  list(
    predictions = predictions_all,
    metrics = metrics_all,
    coefficients = coefficients_all,
    model_info = model_info_all,
    candidate_predictors = candidate_predictors,
    settings = list(
      min_predictor_observed_fraction =
        min_predictor_observed_fraction,
      max_abs_predictor_correlation =
        max_abs_predictor_correlation,
      slope_prior_sd = slope_prior_sd,
      backend = "cmdstanr_direct",
      n_chains = n_chains,
      n_warmup = n_warmup,
      n_sampling = n_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  ),
  file.path(
    out_dir,
    "bayesian_C_d14C_results_bundle.rds"
  )
)


# Save projpred selected predictors.
if (length(all_projpred_info) > 0L) {
  
  projpred_terms_all <- dplyr::bind_rows(
    all_projpred_info
  )
  
  readr::write_csv(
    projpred_terms_all,
    file.path(
      selection_table_dir,
      "selected_predictors_for_3_4_and_5_term_models.csv"
    )
  )
}


# -----------------------------------------------------------------------------
# 10. Completion messages
# -----------------------------------------------------------------------------

message("\n============================================================")
message("BAYESIAN C AND d14C MODELLING COMPLETED")
message("============================================================")

message("\nOutputs saved under: ", out_dir)
message("Model files and CmdStan CSVs: ", model_dir)
message("Observed-vs-predicted and relative-importance plots: ", plot_dir)
message("Tables: ", table_dir)
message("Diagnostics: ", diag_dir)
message("projpred outputs: ", selection_dir)

message(
  "\nEach observed-vs-predicted plot now includes a relative-importance ",
  "bar plot underneath it."
)

message(
  "\nRelative importance is calculated as the absolute standardised ",
  "posterior mean coefficient divided by the sum of the absolute ",
  "standardised predictor coefficients, multiplied by 100."
)

message(
  "\nReminder: relative importance is not variance explained, ",
  "and projpred selection/performance in this script is not cross-validation."
)

message(
  "\nTreat selected 3-, 4- or 5-predictor models as candidate ",
  "parsimonious models; validate before making strong predictive claims."
)
