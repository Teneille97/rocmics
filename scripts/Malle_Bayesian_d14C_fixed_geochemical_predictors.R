# =============================================================================
# Malle ERW: Bayesian d14C models using pre-specified geochemical predictors
#
# Fixed model specifications:
#   d14C_2026_MAOM = f(soil_pH, Fe_ox_over_DCB, Si_ox_over_DCB)
#   d14C_2026_oPOM = f(soil_EC, CEC_Ca, GMD_h2o)
#   d14C_2026_fPOM = f(Mn_ox_over_DCB, Tessier_Acetic_Mg_mg_g)
#
# Run from the RStudio project root. Input:
#   outputs/integrated_geochemistry_for_bayesian.csv
#
# This script does NOT perform projpred or choose predictors from the data.
# It fits exactly the predictor sets listed above (unless a required variable
# is missing/constant, in which case it stops rather than silently changing the
# model). Predictors and responses are standardised for fitting; outputs and
# predictions are returned to the original d14C response scale.
#
# The observed-versus-predicted plots and performance metrics are in-sample;
# they are not cross-validation or independent validation.
# Models are fitted directly with cmdstanr (no brms/rstan dependency).
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
out_dir <- here::here("outputs", "bayesian_d14C_fixed_predictors")
model_dir <- file.path(out_dir, "models")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
diag_dir <- file.path(out_dir, "diagnostics")

for (d in c(out_dir, model_dir, plot_dir, table_dir, diag_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Package 'cmdstanr' is required. Install it and configure CmdStan first.")
}
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not configured. Check cmdstanr::cmdstan_path() and the C++ toolchain.")
}
message("Using CmdStan ", as.character(cmdstanr::cmdstan_version()),
        " via cmdstanr ", as.character(utils::packageVersion("cmdstanr")),
        "; brms/rstan are not used.")

# Sampling settings (consistent with the preceding direct-CmdStan workflow).
n_chains <- 4L
n_warmup <- 1500L
n_sampling <- 2000L
n_cores <- min(n_chains, parallel::detectCores(logical = FALSE))
adapt_delta <- 0.99
max_treedepth <- 12L
slope_prior_sd <- 0.5

# Explicit, fixed predictor sets. Do not reorder these unless you also update
# the coefficient labels in the summaries (the script uses this order).
model_specs <- list(
  d14C_2026_MAOM = c("soil_pH", "Fe_ox_over_DCB", "Si_ox_over_DCB"),
  d14C_2026_oPOM = c("soil_EC", "CEC_Ca", "GMD_h2o"),
  d14C_2026_fPOM = c("Mn_ox_over_DCB", "Tessier_Acetic_Mg_mg_g")
)

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file,
       "\nExpected outputs/integrated_geochemistry_for_bayesian.csv.")
}
integrated <- readr::read_csv(input_file, show_col_types = FALSE,
                              na = c("", "NA", "#VALUE!"))

required_base <- c("sample", "Treatment")
missing_base <- setdiff(required_base, names(integrated))
if (length(missing_base)) {
  stop("Required column(s) missing from input CSV: ", paste(missing_base, collapse = ", "))
}

required_cols <- unique(c(names(model_specs), unlist(model_specs, use.names = FALSE)))
missing_cols <- setdiff(required_cols, names(integrated))
if (length(missing_cols)) {
  stop("Required response/predictor column(s) missing from input CSV: ",
       paste(missing_cols, collapse = ", "),
       "\nThe script will not substitute other predictors automatically.")
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

allowed_treatments <- c("Control", "Lime_2", "Bolsdorfer_50",
                        "Huhnerberg_50", "Eifelgold_50")
dat <- integrated %>%
  filter(Treatment %in% allowed_treatments) %>%
  mutate(
    sample = suppressWarnings(as.integer(as.character(sample))),
    Treatment_model = normalise_treatment(Treatment)
  ) %>%
  filter(!is.na(sample), !is.na(Treatment_model))

if (nrow(dat) == 0L) {
  stop("No rows remain after filtering to the five expected treatments. Check Treatment labels.")
}
readr::write_csv(
  dat %>% count(Treatment, Treatment_model, name = "n_plots"),
  file.path(table_dir, "treatment_sample_counts.csv")
)

# -----------------------------------------------------------------------------
# 1. Prepare each response with only its fixed predictors
# -----------------------------------------------------------------------------
prepare_fixed_model_data <- function(data, response, predictors) {
  y_raw_all <- suppressWarnings(as.numeric(data[[response]]))
  keep_y <- is.finite(y_raw_all)
  work <- data[keep_y, , drop = FALSE]
  y_raw <- y_raw_all[keep_y]

  if (length(y_raw) < 8L) {
    stop(response, ": only ", length(y_raw),
         " finite response observations; at least 8 are required.")
  }

  X <- list()
  imputation <- list()
  predictor_scaling <- list()
  for (v in predictors) {
    x <- suppressWarnings(as.numeric(work[[v]]))
    finite <- is.finite(x)
    if (!any(finite)) stop(response, ": predictor ", v, " has no finite values.")
    med <- stats::median(x[finite])
    x[!finite] <- med
    x_sd <- stats::sd(x)
    if (!is.finite(x_sd) || x_sd <= 0) {
      stop(response, ": predictor ", v, " is constant or has no variation; fixed model cannot be fitted.")
    }
    x_mean <- mean(x)
    X[[v]] <- (x - x_mean) / x_sd
    imputation[[v]] <- tibble(
      response = response, predictor = v,
      n_missing_imputed = sum(!finite), imputation_median = med,
      predictor_mean = x_mean, predictor_sd = x_sd
    )
    predictor_scaling[[v]] <- c(mean = x_mean, sd = x_sd, median = med)
  }

  y_mean <- mean(y_raw)
  y_sd <- stats::sd(y_raw)
  if (!is.finite(y_sd) || y_sd <= 0) stop(response, ": response has no finite variation.")
  model_data <- as_tibble(X)
  model_data$y_std <- (y_raw - y_mean) / y_sd
  model_data$sample <- work$sample
  model_data$Treatment_model <- factor(work$Treatment_model)

  list(
    data = model_data, y_raw = y_raw, y_mean = y_mean, y_sd = y_sd,
    predictors = predictors, n = nrow(model_data),
    imputation = bind_rows(imputation), predictor_scaling = predictor_scaling
  )
}

safe_name <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

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
    response = response, sample = model_data$sample,
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

plot_observed_predicted <- function(pred, response, predictors, output_path) {
  p <- ggplot(pred %>% arrange(observed), aes(x = observed, y = predicted_mean)) +
    geom_ribbon(aes(ymin = fitted_lower_95, ymax = fitted_upper_95),
                fill = "grey70", alpha = 0.55) +
    geom_errorbar(aes(ymin = predictive_lower_95, ymax = predictive_upper_95),
                  width = 0, alpha = 0.35, colour = "grey40") +
    geom_point(aes(shape = Treatment_model), size = 2.8, alpha = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7) +
    coord_equal() +
    labs(
      title = paste("Observed vs Bayesian predicted:", response),
      subtitle = paste0("Fixed predictors: ", paste(predictors, collapse = " + "),
                        " | shaded = 95% fitted-mean credible interval; whiskers = 95% predictive interval"),
      x = paste("Observed", response), y = paste("Predicted", response),
      shape = "Treatment"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(output_path, p, width = 7.8, height = 6.2, dpi = 320)
  p
}

# -----------------------------------------------------------------------------
# 2. Define and compile the regularised Gaussian regression in Stan
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
  for (n in 1:N) y_rep[n] = normal_rng(mu[n], sigma);
}
"
stan_code <- sub("SLOPE_PRIOR_SD", format(slope_prior_sd, scientific = FALSE),
                  stan_code, fixed = TRUE)
stan_file <- file.path(model_dir, "regularised_fixed_predictor_regression.stan")
writeLines(stan_code, stan_file, useBytes = TRUE)
stan_model <- tryCatch(
  cmdstanr::cmdstan_model(stan_file, quiet = TRUE),
  error = function(e) stop("Could not compile CmdStan model. Check C++ toolchain. Details: ",
                           conditionMessage(e))
)

all_predictions <- list()
all_metrics <- list()
all_coefficients <- list()
all_model_info <- list()
all_draws <- list()
all_diagnostics <- list()
all_imputation <- list()

# -----------------------------------------------------------------------------
# 3. Fit the three pre-specified models
# -----------------------------------------------------------------------------
for (response in names(model_specs)) {
  predictors <- model_specs[[response]]
  message("\n============================================================")
  message("Fitting fixed-predictor model: ", response)
  message("Formula: ", response, " ~ ", paste(predictors, collapse = " + "))
  message("============================================================")

  prep <- prepare_fixed_model_data(dat, response, predictors)
  all_imputation[[response]] <- prep$imputation
  message("n = ", prep$n, "; fixed predictors = ", paste(prep$predictors, collapse = ", "))

  X <- as.matrix(prep$data[, predictors, drop = FALSE])
  storage.mode(X) <- "double"
  stan_data <- list(N = nrow(X), K = ncol(X), X = X,
                    y = as.numeric(prep$data$y_std))
  model_id <- safe_name(response)
  csv_dir <- file.path(model_dir, paste0("cmdstan_csv_", model_id))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- stan_model$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_cores,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    seed = 20261028 + match(response, names(model_specs)),
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 200,
    output_dir = csv_dir
  )

  saveRDS(fit, file.path(model_dir, paste0("fit_", model_id, ".rds")))
  draws_array <- fit$draws(format = "draws_array")
  saveRDS(draws_array, file.path(model_dir, paste0("draws_", model_id, ".rds")))
  all_draws[[response]] <- draws_array
  draw_matrix <- posterior::as_draws_matrix(draws_array)

  # Posterior coefficient summaries. Slopes are for a 1-SD increase in each
  # predictor; the response-scale columns convert coefficients to d14C units.
  coef_names <- c("alpha", paste0("beta[", seq_along(predictors), "]"))
  coef_labels <- c("Intercept", predictors)
  coef_summary <- purrr::map_dfr(seq_along(coef_names), function(j) {
    vals <- draw_matrix[, coef_names[j]]
    tibble(
      term = coef_labels[j], Estimate = mean(vals), Est.Error = stats::sd(vals),
      Q2.5 = unname(stats::quantile(vals, 0.025)),
      Q50 = unname(stats::quantile(vals, 0.5)),
      Q97.5 = unname(stats::quantile(vals, 0.975)),
      prob_positive = mean(vals > 0), prob_negative = mean(vals < 0)
    )
  }) %>% mutate(
    response = response, model = "fixed predictor set",
    n_predictors = length(predictors), predictors = paste(predictors, collapse = " + "),
    estimate_original_response_per_predictor_SD = Estimate * prep$y_sd,
    lower_original_response_per_predictor_SD = Q2.5 * prep$y_sd,
    upper_original_response_per_predictor_SD = Q97.5 * prep$y_sd
  )
  all_coefficients[[response]] <- coef_summary
  readr::write_csv(coef_summary, file.path(table_dir, paste0("coefficients_", model_id, ".csv")))

  pred <- make_predictions(draw_matrix, prep$data, prep, response) %>%
    mutate(model = "fixed predictor set", n_predictors = length(predictors),
           predictors = paste(predictors, collapse = " + "), .before = 1)
  all_predictions[[response]] <- pred
  readr::write_csv(pred, file.path(table_dir, paste0("predictions_", model_id, ".csv")))
  plot_observed_predicted(
    pred, response, predictors,
    file.path(plot_dir, paste0("observed_vs_predicted_", model_id, "__fixed_predictors.png"))
  )

  # Posterior predictive check on standardised response values.
  rep_cols <- grep("^y_rep\\[", colnames(draw_matrix), value = TRUE)
  yrep_check <- draw_matrix[seq_len(min(100L, nrow(draw_matrix))), rep_cols, drop = FALSE]
  pp_plot <- bayesplot::ppc_dens_overlay(
    y = as.numeric(prep$data$y_std), yrep = yrep_check
  ) + ggplot2::ggtitle(paste("Posterior predictive check:", response)) +
    ggplot2::theme_bw(base_size = 12)
  ggsave(file.path(diag_dir, paste0("pp_check_", model_id, ".png")),
         pp_plot, width = 7, height = 5.5, dpi = 300)

  diag_summary <- fit$summary(variables = c("alpha", "beta", "sigma"))
  readr::write_csv(diag_summary,
                   file.path(diag_dir, paste0("posterior_draw_diagnostics_", model_id, ".csv")))
  capture.output(fit$summary(), file = file.path(diag_dir, paste0("summary_", model_id, ".txt")))
  sampler <- as.data.frame(fit$sampler_diagnostics(format = "draws_df"))
  diag_row <- tibble(
    response = response, n = prep$n, n_predictors = length(predictors),
    predictors = paste(predictors, collapse = " + "),
    max_Rhat = if (all(is.na(diag_summary$rhat))) NA_real_ else max(diag_summary$rhat, na.rm = TRUE),
    min_bulk_ESS = if (all(is.na(diag_summary$ess_bulk))) NA_real_ else min(diag_summary$ess_bulk, na.rm = TRUE),
    min_tail_ESS = if (all(is.na(diag_summary$ess_tail))) NA_real_ else min(diag_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__ == 1, na.rm = TRUE),
    treedepth_hits = sum(sampler$treedepth__ >= max_treedepth, na.rm = TRUE)
  )
  all_diagnostics[[response]] <- diag_row
  readr::write_csv(diag_row, file.path(diag_dir, paste0("diagnostics_", model_id, ".csv")))

  rmse <- sqrt(mean((pred$observed - pred$predicted_mean)^2, na.rm = TRUE))
  mae <- mean(abs(pred$observed - pred$predicted_mean), na.rm = TRUE)
  r2_corr <- suppressWarnings(stats::cor(pred$observed, pred$predicted_mean,
                                         use = "complete.obs")^2)
  coverage_mean <- mean(pred$observed >= pred$fitted_lower_95 &
                        pred$observed <= pred$fitted_upper_95, na.rm = TRUE)
  coverage_predictive <- mean(pred$observed >= pred$predictive_lower_95 &
                              pred$observed <= pred$predictive_upper_95, na.rm = TRUE)
  all_metrics[[response]] <- tibble(
    response = response, model = "fixed predictor set", n = prep$n,
    n_predictors = length(predictors), predictors = paste(predictors, collapse = " + "),
    RMSE_in_sample = rmse, MAE_in_sample = mae,
    squared_correlation_in_sample = r2_corr,
    observed_coverage_fitted_mean_95 = coverage_mean,
    observed_coverage_predictive_95 = coverage_predictive
  )
  all_model_info[[response]] <- tibble(
    response = response, n = prep$n, n_predictors = length(predictors),
    predictors = paste(predictors, collapse = " + "),
    missing_predictor_values_imputed = sum(prep$imputation$n_missing_imputed)
  )
}

# -----------------------------------------------------------------------------
# 4. Save combined outputs and a reproducibility bundle
# -----------------------------------------------------------------------------
predictions_all <- bind_rows(all_predictions)
metrics_all <- bind_rows(all_metrics)
coefficients_all <- bind_rows(all_coefficients)
model_info_all <- bind_rows(all_model_info)
diagnostics_all <- bind_rows(all_diagnostics)
imputation_all <- bind_rows(all_imputation)

readr::write_csv(predictions_all, file.path(table_dir, "predictions_all_fixed_models.csv"))
readr::write_csv(metrics_all, file.path(table_dir, "model_performance_in_sample.csv"))
readr::write_csv(coefficients_all, file.path(table_dir, "coefficients_all_fixed_models.csv"))
readr::write_csv(model_info_all, file.path(table_dir, "model_predictor_sets.csv"))
readr::write_csv(diagnostics_all, file.path(diag_dir, "diagnostics_all_models.csv"))
readr::write_csv(imputation_all, file.path(table_dir, "predictor_imputation_and_scaling.csv"))
saveRDS(all_draws, file.path(model_dir, "all_fixed_predictor_d14C_posterior_draws.rds"))
saveRDS(
  list(
    model_specs = model_specs,
    predictions = predictions_all,
    metrics = metrics_all,
    coefficients = coefficients_all,
    model_info = model_info_all,
    diagnostics = diagnostics_all,
    predictor_imputation_and_scaling = imputation_all,
    settings = list(
      backend = "cmdstanr_direct", n_chains = n_chains,
      n_warmup = n_warmup, n_sampling = n_sampling,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      slope_prior_sd = slope_prior_sd
    )
  ),
  file.path(out_dir, "bayesian_d14C_fixed_predictors_results_bundle.rds")
)

message("\nFinished. Outputs saved under: ", out_dir)
message("Observed-vs-predicted plots: ", plot_dir)
message("Tables: ", table_dir)
message("Diagnostics and posterior predictive checks: ", diag_dir)
message("No projpred selection was performed; each model used only its pre-specified predictor set.")
