# =============================================================================
# Malle ERW: Bayesian prediction of C and delta-14C for Bulk, MAOM, oPOM, fPOM
#
# Purpose:
#   Fit one regularised Bayesian multiple-regression model for each available
#   response: C_Bulk, C_MAOM, C_oPOM, C_fPOM and the corresponding d14C_2026_*.
#   Save model objects, summaries, posterior predictions, diagnostics, and
#   observed-vs-predicted plots with shaded 95% predictive uncertainty.
#
# Run from the RStudio project root. Expected input:
#   outputs/integrated_geochemistry_for_bayesian.csv
#
# Important:
#   * C and d14C responses are never used as predictors.
#   * Treatment labels, plot IDs, and metadata flags are never used as predictors.
#   * Predictors are standardised; responses are also standardised for fitting.
#   * Regularising Normal priors reduce coefficient overfitting. They do not
#     remove the limitations of a small dataset or strongly correlated predictors.
#   * Predictions shown are in-sample fitted predictions, not cross-validation.
#     Use the intervals to assess model uncertainty, not as external validation.
#   * Models are specified in Stan and fitted directly with cmdstanr.
#     Neither brms nor rstan is loaded or required.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(posterior)
  library(bayesplot)
})

set.seed(20261028)

# -----------------------------------------------------------------------------
# 0. Settings and paths
# -----------------------------------------------------------------------------

input_file <- here::here("outputs", "integrated_geochemistry_for_bayesian.csv")
out_dir <- here::here("outputs", "bayesian_C_d14C_predictions")
model_dir <- file.path(out_dir, "models")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
diag_dir <- file.path(out_dir, "diagnostics")
selection_dir <- file.path(out_dir, "projpred_selection")
selection_plot_dir <- file.path(selection_dir, "plots")
selection_table_dir <- file.path(selection_dir, "tables")

for (d in c(out_dir, model_dir, plot_dir, table_dir, diag_dir,
             selection_dir, selection_plot_dir, selection_table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Fit models directly with cmdstanr; do not load brms or rstan.
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Package 'cmdstanr' is required. Install it and configure CmdStan first.")
}
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not configured. Check cmdstanr::cmdstan_path() and the C++ toolchain.")
}
message("Using CmdStan ", as.character(cmdstanr::cmdstan_version()),
        " via cmdstanr ", as.character(utils::packageVersion("cmdstanr")),
        "; brms/rstan are not used.")
if (!requireNamespace("projpred", quietly = TRUE)) {
  stop("Package 'projpred' is required. Install it with install.packages('projpred') and rerun.")
}
message("Using projpred ", as.character(utils::packageVersion("projpred")),
        " for projection predictive predictor selection.")

# Sampling settings. Increase iter_sampling for final publication analyses.
n_chains <- 4L
n_warmup <- 1500L
n_sampling <- 2000L
n_cores <- min(n_chains, parallel::detectCores(logical = FALSE))
adapt_delta <- 0.99
max_treedepth <- 12L

# Predictor filtering.
# Predictors must have at least this fraction of observations available for
# the response-specific dataset. Missing predictor values are median-imputed.
min_predictor_observed_fraction <- 0.70
# If two predictors are correlated more strongly than this, keep only one of
# them to reduce severe collinearity. This is a preprocessing rule, not
# response-driven variable selection.
max_abs_predictor_correlation <- 0.95

# Model priors on standardised variables.
# Normal(0, 0.5) regularises slopes; it is intentionally more conservative
# than an unregularised multiple regression when n is small.
slope_prior_sd <- 0.5

# -----------------------------------------------------------------------------
# 1. Read data and validate key columns
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file,
       "\nExpected outputs/integrated_geochemistry_for_bayesian.csv.")
}

integrated <- readr::read_csv(input_file, show_col_types = FALSE,
                              na = c("", "NA", "#VALUE!"))

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

# Keep the five main treatment groups used in the Malle Bayesian analyses.
# If the CSV contains additional treatment rows, they are excluded here.
allowed_treatments <- c(
  "Control", "Lime_2", "Bolsdorfer_50", "Huhnerberg_50", "Eifelgold_50"
)

dat <- integrated %>%
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(as.integer(as.character(sample))),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(!is.na(sample), !is.na(Treatment_model))

if (nrow(dat) == 0L) {
  stop("No rows remain after filtering to the five expected treatments. ",
       "Check the Treatment labels in the integrated CSV.")
}

readr::write_csv(
  dat %>% count(Treatment, Treatment_model, name = "n_plots"),
  file.path(table_dir, "treatment_sample_counts.csv")
)

# Response columns: use the known C and d14C columns that actually exist.
c_responses <- c("C_Bulk", "C_MAOM", "C_oPOM", "C_fPOM")
d14c_responses <- c(
  "d14C_2026_Bulk", "d14C_2026_MAOM",
  "d14C_2026_oPOM", "d14C_2026_fPOM"
)
responses <- intersect(c(c_responses, d14c_responses), names(dat))

if (length(responses) == 0L) {
  stop("None of the expected response columns were found. Expected C_Bulk, ",
       "C_MAOM, C_oPOM, C_fPOM and/or d14C_2026_Bulk/MAOM/oPOM/fPOM.")
}
missing_responses <- setdiff(c(c_responses, d14c_responses), names(dat))
if (length(missing_responses)) {
  message("Response columns absent from the CSV (will be skipped): ",
          paste(missing_responses, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 2. Define candidate predictors and explicitly prevent response leakage
# -----------------------------------------------------------------------------

# Reuse the combined geochemical predictor family considered in the integrated
# correlation / MLR workflow. Add October 2023 variables only when present.
# The script intersects this list with the columns actually in the CSV.
explicit_geochem_predictors <- c(
  "Agg stab", "GMD_h2o", "GMD_oxalate", "MWD_h2o", "MWD_oxalate",
  "D_h2o", "D_oxalate", "Dh2o", "Doxalate",
  "Fe_ox_over_DCB", "Al_ox_over_DCB", "Mn_ox_over_DCB", "Si_ox_over_DCB",
  "Tessier_Acetic_Ca_mg_g", "Tessier_Acetic_K_mg_g",
  "Tessier_Acetic_Na_mg_g", "Tessier_Acetic_Mg_mg_g",
  "Tessier_Hydroxylamine_Fe_mg_g", "Tessier_Hydroxylamine_Al_mg_g",
  "Tessier_Hydroxylamine_Mn_mg_g", "Tessier_Hydroxylamine_Si_mg_g",
  "soil_pH", "soil_EC",
  "pH_Oct2023", "EC_Oct2023", "H_Oct2023", "Ca_Oct2023",
  "K_Oct2023", "Mg_Oct2023", "Na_Oct2023", "Al_Oct2023",
  "Fe_Oct2023", "Mn_Oct2023", "CEC_Oct2023", "Base_saturation_Oct2023"
)

# Include additional numeric geochemical columns in known families, but not
# response columns or metadata. This catches extra CEC/Tessier/FTIR columns.
numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
pattern_predictors <- numeric_cols[
  stringr::str_starts(numeric_cols, "Tessier_") |
    stringr::str_starts(numeric_cols, "CEC") |
    stringr::str_starts(numeric_cols, "FTIR_delta") |
    stringr::str_starts(numeric_cols, "Fe_ox") |
    stringr::str_starts(numeric_cols, "Al_ox") |
    stringr::str_starts(numeric_cols, "Mn_ox") |
    stringr::str_starts(numeric_cols, "Si_ox")
]

candidate_predictors <- unique(c(
  intersect(explicit_geochem_predictors, names(dat)),
  intersect(pattern_predictors, names(dat))
))

# Never use carbon / radiocarbon outcomes as predictors, including other pools.
response_like <- names(dat)[
  stringr::str_detect(names(dat), "^C_(Bulk|MAOM|oPOM|fPOM)$") |
    stringr::str_detect(names(dat), "^d14C")
]
metadata_cols <- c(
  "sample", "Treatment", "Treatment_model", "Basalt_type_exp",
  "Dose_resp_exp", "Date", "year", "depth", "fraction", "SOM_fraction"
)
candidate_predictors <- setdiff(
  candidate_predictors,
  c(response_like, metadata_cols)
)
# Exclude encoded metadata / treatment-design columns that may share the
# prefixes used for numeric geochemical variables.
candidate_predictors <- candidate_predictors[
  !stringr::str_detect(
    candidate_predictors,
    regex("dataset|dose_resp_exp|basalt_type_exp|treatment|^sample$", ignore_case = TRUE)
  )
]
candidate_predictors <- candidate_predictors[
  vapply(candidate_predictors, function(v) is.numeric(dat[[v]]), logical(1))
]

if (length(candidate_predictors) == 0L) {
  stop("No numeric candidate geochemical predictors were found. ",
       "Check column names in the integrated CSV and update explicit_geochem_predictors.")
}

readr::write_csv(
  tibble(candidate_predictor = candidate_predictors),
  file.path(table_dir, "candidate_predictors_before_response_filtering.csv")
)
message("Candidate geochemical predictors available: ", length(candidate_predictors))
print(candidate_predictors)

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

# Convert to numeric safely; remove constant / very sparse predictors, median
# impute remaining missing values, standardise. Highly correlated predictors
# are thinned using predictor-predictor correlations only (not response data).
prepare_model_data <- function(data, response, candidates) {
  y_raw <- suppressWarnings(as.numeric(data[[response]]))
  keep_y <- is.finite(y_raw)
  work <- data[keep_y, , drop = FALSE]
  y_raw <- y_raw[keep_y]

  if (length(y_raw) < 8L) {
    stop(response, ": only ", length(y_raw),
         " finite response observations; at least 8 are required.")
  }

  # Exclude predictors with too many missing values or no finite variation.
  eligible <- candidates[vapply(candidates, function(v) {
    x <- suppressWarnings(as.numeric(work[[v]]))
    observed_fraction <- mean(is.finite(x))
    sx <- stats::sd(x, na.rm = TRUE)
    is.finite(observed_fraction) &&
      observed_fraction >= min_predictor_observed_fraction &&
      is.finite(sx) && sx > 0
  }, logical(1))]

  if (length(eligible) == 0L) {
    stop(response, ": no eligible predictors after missingness/variance filtering.")
  }

  X <- as.data.frame(lapply(eligible, function(v) {
    x <- suppressWarnings(as.numeric(work[[v]]))
    med <- stats::median(x[is.finite(x)], na.rm = TRUE)
    x[!is.finite(x)] <- med
    as.numeric(scale(x))
  }))
  names(X) <- eligible

  # Remove any columns that became non-finite or constant after standardising.
  valid <- vapply(X, function(x) all(is.finite(x)) &&
                    is.finite(stats::sd(x)) && stats::sd(x) > 0, logical(1))
  X <- X[, valid, drop = FALSE]

  if (ncol(X) == 0L) {
    stop(response, ": no variable predictors remain after standardisation.")
  }

  # Thin extreme collinearity. The retained variable is the first in the
  # predeclared candidate list; this is not selection by association with y.
  dropped_correlated <- character(0)
  repeat {
    if (ncol(X) < 2L) break
    cm <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
    diag(cm) <- 0
    max_cor <- max(abs(cm), na.rm = TRUE)
    if (!is.finite(max_cor) || max_cor <= max_abs_predictor_correlation) break
    idx <- which(abs(cm) == max_cor, arr.ind = TRUE)[1, ]
    drop_name <- colnames(cm)[idx[2]]
    dropped_correlated <- c(dropped_correlated, drop_name)
    X <- X[, setdiff(names(X), drop_name), drop = FALSE]
  }

  # Standardise response for stable priors and sampling; save scale for inverse.
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
  x %>% stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# Calculate fitted-mean and posterior-predictive intervals from draws produced
# directly by CmdStan. `mu[i]` is the fitted mean and `y_rep[i]` a replicated
# observation, both on the standardised response scale.
make_predictions <- function(draw_matrix, model_data, prep, response) {
  mu_cols <- grep("^mu\\[", colnames(draw_matrix), value = TRUE)
  rep_cols <- grep("^y_rep\\[", colnames(draw_matrix), value = TRUE)
  if (length(mu_cols) != nrow(model_data) || length(rep_cols) != nrow(model_data)) {
    stop("Could not extract all mu/y_rep draws for ", response)
  }
  mu_draws <- draw_matrix[, mu_cols, drop = FALSE]
  rep_draws <- draw_matrix[, rep_cols, drop = FALSE]

  pred_mean <- apply(mu_draws, 2, mean)
  pred_median <- apply(mu_draws, 2, median)
  pred_lower <- apply(mu_draws, 2, stats::quantile, probs = 0.025)
  pred_upper <- apply(mu_draws, 2, stats::quantile, probs = 0.975)
  pred_lower_pi <- apply(rep_draws, 2, stats::quantile, probs = 0.025)
  pred_upper_pi <- apply(rep_draws, 2, stats::quantile, probs = 0.975)

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

plot_observed_predicted <- function(pred, response, output_path, model_label = "full reference model") {
  # Sort by observed values so uncertainty ribbons form a readable band.
  pdat <- pred %>% arrange(observed) %>%
    mutate(row_order = row_number())

  p <- ggplot(pdat, aes(x = observed, y = predicted_mean)) +
    geom_ribbon(
      aes(ymin = fitted_lower_95, ymax = fitted_upper_95),
      fill = "grey70", alpha = 0.55
    ) +
    geom_errorbar(
      aes(ymin = predictive_lower_95, ymax = predictive_upper_95),
      width = 0, alpha = 0.35, colour = "grey40"
    ) +
    geom_point(aes(shape = Treatment_model), size = 2.8, alpha = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                linewidth = 0.7, colour = "black") +
    coord_equal() +
    labs(
      title = paste("Observed vs Bayesian predicted:", response),
      subtitle = paste0(model_label, " | ",
        "Shaded band = 95% credible interval for fitted mean; whiskers = ",
        "95% posterior predictive interval"
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

  ggsave(output_path, p, width = 7.5, height = 6.2, dpi = 320)
  p
}

# -----------------------------------------------------------------------------
# 4. Define a Stan model and fit one regression per response with cmdstanr
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
stan_code <- sub("SLOPE_PRIOR_SD", format(slope_prior_sd, scientific = FALSE), stan_code, fixed = TRUE)
stan_file <- file.path(model_dir, "regularised_multiple_regression.stan")
writeLines(stan_code, stan_file, useBytes = TRUE)
stan_model <- tryCatch(
  cmdstanr::cmdstan_model(stan_file, quiet = TRUE),
  error = function(e) stop(
    "Could not compile the direct CmdStan model. Check the C++ toolchain. Details: ",
    conditionMessage(e)
  )
)

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
    prepare_model_data(dat, response, candidate_predictors),
    error = function(e) {
      message("Skipping ", response, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(prep)) next

  message("n = ", prep$n, "; predictors = ", length(prep$predictors))
  message("Predictors: ", paste(prep$predictors, collapse = ", "))
  if (length(prep$dropped_correlated)) {
    message("Dropped for |r| > ", max_abs_predictor_correlation, ": ",
            paste(prep$dropped_correlated, collapse = ", "))
  }

  model_id <- safe_name(response)
  X <- as.matrix(prep$data[, prep$predictors, drop = FALSE])
  storage.mode(X) <- "double"
  stan_data <- list(
    N = nrow(X), K = ncol(X), X = X,
    y = as.numeric(prep$data$y_std)
  )
  csv_dir <- file.path(model_dir, paste0("cmdstan_csv_", model_id))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  draws_file <- file.path(model_dir, paste0("draws_", model_id, ".rds"))
  fit <- NULL

  # Reuse saved full-model posterior draws when available. This is particularly
  # useful when the full model completed but projpred failed in a previous run.
  if (file.exists(draws_file)) {
    message("Reusing saved full-model posterior draws for ", response,
            "; skipping the full-model MCMC refit.")
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
        message("Model failed for ", response, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) next

    # Save draws to reuse them if a later selection/plotting step fails.
    draws_array <- fit$draws(format = "draws_array")
    saveRDS(draws_array, draws_file)
  }

  all_draws[[response]] <- draws_array
  draw_matrix <- posterior::as_draws_matrix(draws_array)

  # ---------------------------------------------------------------------------
  # Projection predictive variable selection (projpred)
  # ---------------------------------------------------------------------------
  # Custom reference-model adapter: use posterior draws from direct cmdstanr
  # fits, without loading brms or rstan.
  tryCatch({
    ref_fit <- list(
      draw_matrix = as.matrix(draw_matrix),
      predictors = prep$predictors,
      data = as.data.frame(prep$data)
    )

    ref_predfun <- function(fit, newdata = NULL) {
      if (is.null(newdata)) newdata <- fit$data
      X_new <- as.matrix(newdata[, fit$predictors, drop = FALSE])
      storage.mode(X_new) <- "double"
      beta_names <- paste0("beta[", seq_along(fit$predictors), "]")
      beta_draws <- as.matrix(fit$draw_matrix[, beta_names, drop = FALSE])
      alpha_draws <- as.numeric(fit$draw_matrix[, "alpha"])
      eta <- X_new %*% t(beta_draws)
      sweep(eta, MARGIN = 2, STATS = alpha_draws, FUN = "+")
    }

    full_formula <- stats::as.formula(
      paste("y_std ~", paste(prep$predictors, collapse = " + "))
    )
    ref_obj <- projpred::init_refmodel(
      object = ref_fit,
      data = as.data.frame(prep$data),
      formula = full_formula,
      family = stats::gaussian(),
      ref_predfun = ref_predfun,
      dis = as.numeric(draw_matrix[, "sigma"])
    )

    max_terms <- min(5L, length(prep$predictors))
    vsel <- projpred::varsel(
      ref_obj,
      method = "forward",
      nterms_max = max_terms,
      nclusters = min(20L, nrow(draw_matrix)),
      nclusters_pred = min(20L, nrow(draw_matrix)),
      seed = 20261028 + match(response, responses),
      verbose = 1
    )

    saveRDS(vsel, file.path(selection_dir, paste0("projpred_", model_id, ".rds")))
    capture.output(
      print(summary(vsel)),
      file = file.path(selection_table_dir, paste0("summary_", model_id, ".txt"))
    )
    capture.output(
      print(projpred::ranking(vsel)),
      file = file.path(selection_table_dir, paste0("ranking_", model_id, ".txt"))
    )

    p_vsel <- plot(vsel) +
      ggplot2::ggtitle(paste("projpred submodel performance:", response)) +
      ggplot2::theme_bw(base_size = 12)
    ggplot2::ggsave(
      file.path(selection_plot_dir, paste0("projpred_performance_", model_id, ".png")),
      p_vsel, width = 7.5, height = 5.5, dpi = 300
    )

    # Project onto 3- and 4-predictor submodels and save the exact terms.
    # This is exploratory selection, not cross-validated performance assessment.
    for (k in intersect(c(3L, 4L, 5L), seq_len(max_terms))) {
      prj_k <- projpred::project(
        vsel,
        nterms = k,
        ndraws = min(400L, nrow(draw_matrix)),
        seed = 20261028 + 100L * match(response, responses) + k,
        verbose = 0
      )
      selected_k <- as.character(projpred::predictor_terms(prj_k))
      selected_k <- intersect(selected_k, prep$predictors)
      if (length(selected_k) != k) {
        stop("projpred returned ", length(selected_k), " unique usable predictors for the ",
             k, "-predictor model of ", response, ". Returned: ",
             paste(selected_k, collapse = ", "))
      }

      all_projpred_info[[paste(response, k, sep = "__")]] <- tibble::tibble(
        response = response,
        n_predictors = k,
        selection_order = seq_along(selected_k),
        predictor = selected_k
      )
      saveRDS(prj_k, file.path(selection_dir,
                               paste0("projection_", model_id, "_", k, "predictors.rds")))

      # Fit a directly sampled CmdStan regression using exactly the predictors
      # selected by projpred. This gives the 3/4/5-predictor models their own
      # posterior draws, coefficient summaries, diagnostics and predictions.
      # These are refitted submodels (not the full model's fitted values and
      # not the projection object's draws).
      sub_id <- paste0(model_id, "__projpred_", k, "predictors")
      X_sub <- as.matrix(prep$data[, selected_k, drop = FALSE])
      storage.mode(X_sub) <- "double"
      sub_stan_data <- list(
        N = nrow(X_sub), K = ncol(X_sub), X = X_sub,
        y = as.numeric(prep$data$y_std)
      )
      sub_csv_dir <- file.path(model_dir, paste0("cmdstan_csv_", sub_id))
      dir.create(sub_csv_dir, recursive = TRUE, showWarnings = FALSE)
      message("Fitting selected ", k, "-predictor submodel for ", response,
              ": ", paste(selected_k, collapse = " + "))
      sub_fit <- stan_model$sample(
        data = sub_stan_data,
        chains = n_chains,
        parallel_chains = n_cores,
        iter_warmup = n_warmup,
        iter_sampling = n_sampling,
        seed = 20261028 + 1000L * match(response, responses) + k,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = 200,
        output_dir = sub_csv_dir
      )
      saveRDS(sub_fit, file.path(model_dir, paste0("fit_", sub_id, ".rds")))
      sub_draws_array <- sub_fit$draws(format = "draws_array")
      saveRDS(sub_draws_array, file.path(model_dir, paste0("draws_", sub_id, ".rds")))
      sub_draw_matrix <- posterior::as_draws_matrix(sub_draws_array)
      all_draws[[sub_id]] <- sub_draws_array

      sub_prep <- prep
      sub_prep$predictors <- selected_k
      sub_data <- prep$data[, c(selected_k, "y_std", "sample", "Treatment_model"), drop = FALSE]
      sub_pred <- make_predictions(sub_draw_matrix, sub_data, sub_prep, response) %>%
        mutate(model = paste0(k, "-predictor refit"), n_predictors = k,
               predictors = paste(selected_k, collapse = " + "), .before = 1)
      all_predictions[[sub_id]] <- sub_pred
      readr::write_csv(sub_pred,
        file.path(table_dir, paste0("predictions_", sub_id, ".csv")))
      plot_observed_predicted(
        sub_pred, response,
        file.path(plot_dir, paste0("observed_vs_predicted_", sub_id, ".png")),
        model_label = paste0("projpred-selected ", k, "-predictor refit")
      )

      # Coefficients are on standardised predictor and standardised response
      # scales; the additional columns convert the response scale back to the
      # original response units per 1 SD increase in each predictor.
      sub_coef_names <- c("alpha", paste0("beta[", seq_along(selected_k), "]"))
      sub_coef_labels <- c("Intercept", selected_k)
      sub_coef <- purrr::map_dfr(seq_along(sub_coef_names), function(j) {
        vals <- sub_draw_matrix[, sub_coef_names[j]]
        tibble(
          term = sub_coef_labels[j],
          Estimate = mean(vals),
          Est.Error = stats::sd(vals),
          Q2.5 = unname(stats::quantile(vals, 0.025)),
          Q50 = unname(stats::quantile(vals, 0.5)),
          Q97.5 = unname(stats::quantile(vals, 0.975)),
          prob_positive = mean(vals > 0),
          prob_negative = mean(vals < 0)
        )
      }) %>% mutate(
        response = response, model = paste0(k, "-predictor refit"),
        n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        estimate_original_response_per_predictor_SD = Estimate * prep$y_sd,
        lower_original_response_per_predictor_SD = Q2.5 * prep$y_sd,
        upper_original_response_per_predictor_SD = Q97.5 * prep$y_sd
      )
      all_coefficients[[sub_id]] <- sub_coef
      readr::write_csv(sub_coef,
        file.path(table_dir, paste0("coefficients_", sub_id, ".csv")))

      sub_summary <- sub_fit$summary(variables = c("alpha", "beta", "sigma"))
      sub_sampler <- as.data.frame(sub_fit$sampler_diagnostics(format = "draws_df"))
      sub_diag <- tibble(
        response = response, model = paste0(k, "-predictor refit"),
        n = prep$n, n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        max_Rhat = if (all(is.na(sub_summary$rhat))) NA_real_ else max(sub_summary$rhat, na.rm = TRUE),
        min_bulk_ESS = if (all(is.na(sub_summary$ess_bulk))) NA_real_ else min(sub_summary$ess_bulk, na.rm = TRUE),
        min_tail_ESS = if (all(is.na(sub_summary$ess_tail))) NA_real_ else min(sub_summary$ess_tail, na.rm = TRUE),
        divergences = sum(sub_sampler$divergent__ == 1, na.rm = TRUE),
        treedepth_hits = sum(sub_sampler$treedepth__ >= max_treedepth, na.rm = TRUE)
      )
      readr::write_csv(sub_diag,
        file.path(diag_dir, paste0("diagnostics_", sub_id, ".csv")))
      capture.output(sub_fit$summary(),
        file = file.path(diag_dir, paste0("summary_", sub_id, ".txt")))
      sub_rep_cols <- grep("^y_rep\\[", colnames(sub_draw_matrix), value = TRUE)
      sub_yrep <- sub_draw_matrix[seq_len(min(100L, nrow(sub_draw_matrix))), sub_rep_cols, drop = FALSE]
      sub_pp <- bayesplot::ppc_dens_overlay(y = as.numeric(prep$data$y_std), yrep = sub_yrep) +
        ggplot2::ggtitle(paste("Posterior predictive check:", response, "-", k, "predictors")) +
        ggplot2::theme_bw(base_size = 12)
      ggplot2::ggsave(file.path(diag_dir, paste0("pp_check_", sub_id, ".png")),
        sub_pp, width = 7, height = 5.5, dpi = 300)

      sub_rmse <- sqrt(mean((sub_pred$observed - sub_pred$predicted_mean)^2, na.rm = TRUE))
      sub_mae <- mean(abs(sub_pred$observed - sub_pred$predicted_mean), na.rm = TRUE)
      sub_r2 <- suppressWarnings(stats::cor(sub_pred$observed, sub_pred$predicted_mean,
                                            use = "complete.obs")^2)
      sub_cov_mean <- mean(sub_pred$observed >= sub_pred$fitted_lower_95 &
                           sub_pred$observed <= sub_pred$fitted_upper_95, na.rm = TRUE)
      sub_cov_pred <- mean(sub_pred$observed >= sub_pred$predictive_lower_95 &
                           sub_pred$observed <= sub_pred$predictive_upper_95, na.rm = TRUE)
      all_metrics[[sub_id]] <- tibble(
        response = response, model = paste0(k, "-predictor refit"),
        n = prep$n, n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        RMSE_in_sample = sub_rmse, MAE_in_sample = sub_mae,
        squared_correlation_in_sample = sub_r2,
        observed_coverage_fitted_mean_95 = sub_cov_mean,
        observed_coverage_predictive_95 = sub_cov_pred
      )
      all_model_info[[sub_id]] <- tibble(
        response = response, model = paste0(k, "-predictor refit"),
        n = prep$n, n_predictors = k,
        predictors = paste(selected_k, collapse = " + "),
        dropped_correlated = paste(prep$dropped_correlated, collapse = " + ")
      )
    }
    message("projpred selection completed for ", response,
            "; inspect projpred_selection/tables and projpred_selection/plots.")
  }, error = function(e) {
    message("projpred selection failed for ", response, ": ", conditionMessage(e),
            "\nThe full CmdStan model and its predictions will still be saved.")
  })
  # Summaries are computed from the posterior draws so this section also works
  # when the saved draws were reused and no live CmdStan fit object exists.
  diag_summary <- if (!is.null(fit)) {
    fit$summary()
  } else {
    posterior::summarise_draws(posterior::as_draws_array(draws_array))
  }
  capture.output(diag_summary,
                 file = file.path(diag_dir, paste0("summary_", model_id, ".txt")))
  readr::write_csv(diag_summary,
                   file.path(diag_dir, paste0("posterior_draw_diagnostics_", model_id, ".csv")))

  # Posterior predictive check on standardised response values.
  rep_cols <- grep("^y_rep\\[", colnames(draw_matrix), value = TRUE)
  n_pp <- min(100L, nrow(draw_matrix))
  yrep_check <- draw_matrix[seq_len(n_pp), rep_cols, drop = FALSE]
  pp_plot <- bayesplot::ppc_dens_overlay(
    y = as.numeric(prep$data$y_std), yrep = yrep_check
  ) +
    ggtitle(paste("Posterior predictive check:", response)) +
    theme_bw(base_size = 12)
  ggsave(file.path(diag_dir, paste0("pp_check_", model_id, ".png")),
         pp_plot, width = 7, height = 5.5, dpi = 300)

  # Summarise intercept and slopes on standardised and original response scales.
  coef_names <- c("alpha", paste0("beta[", seq_along(prep$predictors), "]"))
  coef_labels <- c("Intercept", prep$predictors)
  coef_summary <- lapply(seq_along(coef_names), function(j) {
    vals <- draw_matrix[, coef_names[j]]
    tibble(
      term = coef_labels[j],
      Estimate = mean(vals),
      Est.Error = stats::sd(vals),
      Q2.5 = unname(stats::quantile(vals, 0.025)),
      Q50 = unname(stats::quantile(vals, 0.5)),
      Q97.5 = unname(stats::quantile(vals, 0.975))
    )
  }) %>% bind_rows() %>%
    mutate(
      response = response,
      model = "full reference model",
      n_predictors = length(prep$predictors),
      predictors = paste(prep$predictors, collapse = " + "),
      estimate_original_response_per_predictor_SD = Estimate * prep$y_sd,
      lower_original_response_per_predictor_SD = Q2.5 * prep$y_sd,
      upper_original_response_per_predictor_SD = Q97.5 * prep$y_sd
    )
  all_coefficients[[response]] <- coef_summary
  readr::write_csv(coef_summary,
                   file.path(table_dir, paste0("coefficients_", model_id, ".csv")))

  pred <- make_predictions(draw_matrix, prep$data, prep, response) %>%
    mutate(model = "full reference model",
           n_predictors = length(prep$predictors),
           predictors = paste(prep$predictors, collapse = " + "), .before = 1)
  all_predictions[[response]] <- pred
  readr::write_csv(pred, file.path(table_dir, paste0("predictions_", model_id, ".csv")))
  plot_observed_predicted(
    pred, response,
    file.path(plot_dir, paste0("observed_vs_predicted_", model_id, ".png"))
  )

  rmse <- sqrt(mean((pred$observed - pred$predicted_mean)^2, na.rm = TRUE))
  mae <- mean(abs(pred$observed - pred$predicted_mean), na.rm = TRUE)
  r2_corr <- suppressWarnings(cor(pred$observed, pred$predicted_mean,
                                 use = "complete.obs")^2)
  coverage_mean <- mean(
    pred$observed >= pred$fitted_lower_95 &
      pred$observed <= pred$fitted_upper_95, na.rm = TRUE
  )
  coverage_predictive <- mean(
    pred$observed >= pred$predictive_lower_95 &
      pred$observed <= pred$predictive_upper_95, na.rm = TRUE
  )

  all_metrics[[response]] <- tibble(
    response = response,
    model = "full reference model",
    n = prep$n,
    n_predictors = length(prep$predictors),
    predictors = paste(prep$predictors, collapse = " + "),
    predictors_dropped_for_collinearity = paste(prep$dropped_correlated, collapse = " + "),
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
    dropped_correlated = paste(prep$dropped_correlated, collapse = " + ")
  )
}

# -----------------------------------------------------------------------------
# 5. Save combined outputs
# -----------------------------------------------------------------------------

if (length(all_predictions) == 0L) {
  stop("No response models completed successfully. Check messages and diagnostics.")
}

predictions_all <- bind_rows(all_predictions)
metrics_all <- bind_rows(all_metrics)
coefficients_all <- bind_rows(all_coefficients)
model_info_all <- bind_rows(all_model_info)

readr::write_csv(predictions_all, file.path(table_dir, "predictions_all_responses.csv"))
readr::write_csv(metrics_all, file.path(table_dir, "model_performance_in_sample.csv"))
readr::write_csv(coefficients_all, file.path(table_dir, "coefficients_all_responses.csv"))
readr::write_csv(model_info_all, file.path(table_dir, "model_predictor_sets.csv"))
# Save posterior draws (portable and sufficient for downstream summaries) and
# fit objects for the current R session. CmdStan CSV output is retained per model.
saveRDS(all_draws, file.path(model_dir, "all_bayesian_C_d14C_posterior_draws.rds"))
saveRDS(
  list(
    predictions = predictions_all,
    metrics = metrics_all,
    coefficients = coefficients_all,
    model_info = model_info_all,
    candidate_predictors = candidate_predictors,
    settings = list(
      min_predictor_observed_fraction = min_predictor_observed_fraction,
      max_abs_predictor_correlation = max_abs_predictor_correlation,
      slope_prior_sd = slope_prior_sd,
      backend = "cmdstanr_direct",
      n_chains = n_chains,
      n_warmup = n_warmup,
      n_sampling = n_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  ),
  file.path(out_dir, "bayesian_C_d14C_results_bundle.rds")
)

message("\nFinished. Outputs saved under: ", out_dir)
message("Model files and CmdStan CSVs: ", model_dir)
message("Backend: direct cmdstanr (no brms/rstan dependency).")
message("Existing full-model draws are reused when available; selected 3/4/5-predictor submodels are fitted separately.")
message("Observed-vs-predicted plots: ", plot_dir)
message("Tables: ", table_dir)
message("Diagnostics: ", diag_dir)
if (length(all_projpred_info)) {
  projpred_terms_all <- dplyr::bind_rows(all_projpred_info)
  readr::write_csv(
    projpred_terms_all,
    file.path(selection_table_dir, "selected_predictors_for_3_4_and_5_term_models.csv")
  )
}

message("\nprojpred outputs: ", selection_dir)
message("\nReminder: projpred selection/performance in this script is not cross-validation.")
message("Treat selected 3- or 4-predictor models as candidate parsimonious models; validate before strong predictive claims.")
