# =============================================================================
# Malle ERW: Bayesian multivariable models for SoilR decomposition/transfer
# parameters using treatment-level geochemistry
#
# Run from the RStudio project root.
#
# ONE multivariable model per SoilR parameter:
#   kf        ~ Mn_ox_over_DCB + soil_pH
#   ki        ~ GMD_h2o + CEC_Ca + soil_EC
#   ks        ~ Fe_ox_over_DCB + soil_pH + Si_ox_over_DCB
#   alpha_fi  ~ GMD_h2o + CEC_Ca
#   alpha_is  ~ Fe_ox_over_DCB + CEC_Ca + Si_ox_over_DCB
#
# Important limitations:
# - There are only five independent treatment-level SoilR fits. Each model
#   therefore has five response observations, regardless of how many posterior
#   draws are available from each SoilR fit.
# - Geochemistry is summarised across the three plots within each treatment.
# - These multivariable regressions are exploratory. With n = 5 and up to three
#   predictors, estimates can be weakly identified and sensitive to priors.
# - The script loads existing SoilR MCMC files; it does not refit SoilR.
# - No predictor-vs-response plots are produced. The only plots compare the
#   observed and fitted SoilR parameter estimates across the five treatments.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(posterior)
})

set.seed(20261027)

# ---------------------------------------------------------------------------
# 0. Paths and settings
# ---------------------------------------------------------------------------

integrated_file <- here::here(
  "outputs", "integrated_geochemistry_for_bayesian.csv"
)
mcmc_dir <- here::here("mod_runs")
out_dir <- here::here("outputs", "bayesian_rate_transfer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

burnin <- 20000L
n_posterior_draws <- 1000L
n_chains <- 4L
n_warmup <- 2000L
n_sampling <- 2000L

mcmc_files <- c(
  Control = "MCMC_3p_Control.Rdata",
  Lime = "MCMC_3p_Lime.Rdata",
  Bolsdorfer = "MCMC_3p_Bolsdorfer.Rdata",
  Huhnerberg = "MCMC_3p_Huhnerberg.Rdata",
  Eifelgold = "MCMC_3p_Eifelgold.Rdata"
)
model_treatments <- c("Control", "Lime", "Bolsdorfer", "Huhnerberg", "Eifelgold")

# Each response has ONE model containing ALL predictors listed for it.
model_specs <- list(
  kf = c("Mn_ox_over_DCB", "soil_pH"),
  ki = c("GMD_h2o", "CEC_Ca", "soil_EC"),
  ks = c("Fe_ox_over_DCB", "soil_pH", "Si_ox_over_DCB"),
  alpha_fi = c("GMD_h2o", "CEC_Ca"),
  alpha_is = c("Fe_ox_over_DCB", "CEC_Ca", "Si_ox_over_DCB")
)

response_transform <- c(
  kf = "log",
  ki = "log",
  ks = "log",
  alpha_fi = "logit",
  alpha_is = "logit"
)

# ---------------------------------------------------------------------------
# 1. Read and validate integrated geochemistry
# ---------------------------------------------------------------------------

if (!file.exists(integrated_file)) {
  stop("Cannot find ", integrated_file,
       ". Expected outputs/integrated_geochemistry_for_bayesian.csv.")
}

integrated <- readr::read_csv(integrated_file, show_col_types = FALSE)
if (!all(c("sample", "Treatment") %in% names(integrated))) {
  stop("The integrated CSV must contain columns named sample and Treatment.")
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

integrated <- integrated %>%
  filter(Treatment %in% c(
    "Control", "Lime_2", "Bolsdorfer_50",
    "Huhnerberg_50", "Eifelgold_50"
  )) %>%
  mutate(Treatment_model = normalise_treatment(Treatment))

if (anyNA(integrated$Treatment_model)) {
  stop("Unmapped treatment labels: ",
       paste(unique(integrated$Treatment[is.na(integrated$Treatment_model)]),
             collapse = ", "))
}

treatment_check <- integrated %>%
  count(Treatment, Treatment_model, name = "n_plots")
print(treatment_check, n = Inf)

if (nrow(treatment_check) != 5L || any(treatment_check$n_plots != 3L) ||
    !setequal(unique(integrated$Treatment_model), model_treatments)) {
  stop("Expected exactly five modeled treatments with three plots each.")
}

# ---------------------------------------------------------------------------
# 2. Validate and summarise geochemical predictors by treatment
# ---------------------------------------------------------------------------

all_predictors <- unique(unlist(model_specs, use.names = FALSE))
missing_predictors <- setdiff(all_predictors, names(integrated))
if (length(missing_predictors) > 0L) {
  stop(
    "Missing primary predictor columns in integrated CSV: ",
    paste(missing_predictors, collapse = ", ")
  )
}

integrated <- integrated %>%
  mutate(across(all_of(all_predictors), ~ suppressWarnings(as.numeric(.x))))

geochem_summary <- integrated %>%
  group_by(Treatment_model) %>%
  summarise(
    across(
      all_of(all_predictors),
      list(
        mean = ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE),
        sd = ~ if (sum(is.finite(.x)) < 2L) NA_real_ else sd(.x, na.rm = TRUE),
        n = ~ sum(is.finite(.x))
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

readr::write_csv(geochem_summary,
                 file.path(out_dir, "geochemistry_treatment_summary.csv"))

geochem_long <- geochem_summary %>%
  pivot_longer(
    cols = -Treatment_model,
    names_to = c("predictor", "stat"),
    names_pattern = "^(.*)__(mean|sd|n)$"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(se = if_else(n > 0, sd / sqrt(n), NA_real_))

geochem_wide_means <- geochem_long %>%
  select(Treatment_model, predictor, mean) %>%
  pivot_wider(names_from = predictor, values_from = mean)

geochem_wide_ses <- geochem_long %>%
  select(Treatment_model, predictor, se) %>%
  pivot_wider(names_from = predictor, values_from = se) %>%
  rename_with(~ paste0(.x, "_se"), -Treatment_model)

geochem_prepped <- left_join(geochem_wide_means, geochem_wide_ses,
                             by = "Treatment_model")

# ---------------------------------------------------------------------------
# 3. Extract posterior draws of the five SoilR parameters from saved MCMC
# ---------------------------------------------------------------------------

get_par_draws <- function(path) {
  if (!file.exists(path)) stop("MCMC file not found: ", path)
  e <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = e)
  if (!"MCMC" %in% obj_names || is.null(e$MCMC$pars)) {
    stop(path, " must contain an object named MCMC with MCMC$pars.")
  }

  pars <- as.data.frame(e$MCMC$pars)
  wanted <- c("kf", "ki", "ks", "alpha_fi", "alpha_is")
  if (all(wanted %in% names(pars))) {
    pars <- pars[, wanted, drop = FALSE]
  } else {
    if (ncol(pars) < 5L) stop("Fewer than five parameter columns in ", path)
    warning("Parameter names not all found in ", basename(path),
            "; using first five columns in assumed order: ",
            paste(names(pars)[1:5], collapse = ", "))
    pars <- pars[, 1:5, drop = FALSE]
    names(pars) <- wanted
  }

  if (nrow(pars) <= burnin) {
    stop("burnin = ", burnin, " leaves no draws in ", basename(path),
         ". Check chain length and adjust burnin.")
  }
  pars <- pars[(burnin + 1L):nrow(pars), , drop = FALSE]
  pars <- pars[complete.cases(pars) &
                 apply(pars, 1, function(x) all(is.finite(x))), , drop = FALSE]
  if (nrow(pars) == 0L) stop("No finite posterior draws in ", path)

  pars <- pars[sample(seq_len(nrow(pars)), min(n_posterior_draws, nrow(pars))),
               , drop = FALSE]
  pars
}

parameter_draws <- purrr::map_dfr(names(mcmc_files), function(tr) {
  path <- file.path(mcmc_dir, mcmc_files[[tr]])
  message("Reading saved SoilR MCMC parameters for ", tr, ": ", path)
  get_par_draws(path) %>% mutate(Treatment_model = tr, draw = row_number(), .before = 1)
})

readr::write_csv(parameter_draws,
                 file.path(out_dir, "soilR_parameter_posterior_draws_sampled.csv"))

parameter_summary <- parameter_draws %>%
  pivot_longer(cols = all_of(names(response_transform)),
               names_to = "outcome", values_to = "value") %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    # The plotted observation is the posterior mean parameter estimate.
    # Vertical observed error bars show +/- one posterior SD (a concise
    # summary of uncertainty in the saved SoilR parameter distribution).
    y_mean = mean(value, na.rm = TRUE),
    y_sd = sd(value, na.rm = TRUE),
    observed_lower = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE),
    observed_upper = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE),
    q2.5 = quantile(value, 0.025, na.rm = TRUE),
    q97.5 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(parameter_summary,
                 file.path(out_dir, "soilR_parameter_summary.csv"))

# ---------------------------------------------------------------------------
# 4. Transform response posterior draws and prepare backup data
# ---------------------------------------------------------------------------

transform_response <- function(x, method) {
  if (method == "log") {
    if (any(x <= 0)) stop("A decomposition-rate posterior draw was <= 0.")
    log(x)
  } else if (method == "logit") {
    if (any(x <= 0 | x >= 1)) {
      stop("Transfer-fraction draws must be strictly between 0 and 1 for logit.")
    }
    qlogis(x)
  } else {
    stop("Unknown response transform: ", method)
  }
}

transformed_parameter_summary <- parameter_draws %>%
  pivot_longer(cols = all_of(names(response_transform)),
               names_to = "outcome", values_to = "value") %>%
  mutate(
    transform = unname(response_transform[outcome]),
    transformed_value = map2_dbl(value, transform, transform_response)
  ) %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    y_mean_trans = mean(transformed_value, na.rm = TRUE),
    y_sd_trans = sd(transformed_value, na.rm = TRUE),
    .groups = "drop"
  )

if (any(!is.finite(transformed_parameter_summary$y_mean_trans)) ||
    any(!is.finite(transformed_parameter_summary$y_sd_trans))) {
  stop("Non-finite transformed parameter summaries; check parameter draws.")
}

model_data_long <- parameter_summary %>%
  left_join(geochem_prepped, by = "Treatment_model") %>%
  arrange(match(Treatment_model, model_treatments), outcome)
readr::write_csv(model_data_long,
                 file.path(out_dir, "treatment_level_model_data_long.csv"))

dir.create(file.path(out_dir, "session_backup"), recursive = TRUE,
            showWarnings = FALSE)
saveRDS(model_data_long, file.path(out_dir, "session_backup", "model_data_long.rds"))
saveRDS(parameter_draws, file.path(out_dir, "session_backup", "parameter_draws.rds"))
saveRDS(geochem_summary, file.path(out_dir, "session_backup", "geochem_summary.rds"))

# ---------------------------------------------------------------------------
# 5. Configure CmdStan multivariable measurement-error model
# ---------------------------------------------------------------------------

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr is required. Install/configure cmdstanr and CmdStan first.")
}
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not configured. Check cmdstanr::cmdstan_path() and toolchain.")
}

stan_file <- file.path(out_dir, "bayesian_multivariable_measurement_error.stan")
stan_code <- "
data {
  int<lower=2> N;
  int<lower=1> P;
  vector[N] y;
  vector<lower=0>[N] y_se;
  matrix[N, P] x_obs;
  matrix<lower=0>[N, P] x_se;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> sigma;
  matrix[N, P] z_x;
}
transformed parameters {
  matrix[N, P] x_true;
  for (n in 1:N)
    for (p in 1:P)
      x_true[n, p] = x_obs[n, p] + x_se[n, p] * z_x[n, p];
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  sigma ~ exponential(1);
  to_vector(z_x) ~ std_normal();
  for (n in 1:N)
    y[n] ~ normal(alpha + x_true[n] * beta,
                  sqrt(square(y_se[n]) + square(sigma)));
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  for (n in 1:N) {
    mu[n] = alpha + x_true[n] * beta;
    y_rep[n] = normal_rng(mu[n], sqrt(square(y_se[n]) + square(sigma)));
  }
}
"
writeLines(stan_code, stan_file)
stan_model <- cmdstanr::cmdstan_model(stan_file, quiet = FALSE)

summarise_draws_simple <- function(x) {
  tibble(
    mean = mean(x), median = median(x), sd = sd(x),
    q2.5 = unname(quantile(x, 0.025)),
    q97.5 = unname(quantile(x, 0.975)),
    Pr_positive = mean(x > 0),
    Pr_negative = mean(x < 0)
  )
}

inverse_response_transform <- function(x, method) {
  if (method == "log") exp(x) else if (method == "logit") plogis(x) else
    stop("Unknown response transform: ", method)
}

# ---------------------------------------------------------------------------
# 6. Fit exactly ONE multivariable model per outcome
# ---------------------------------------------------------------------------

fit_one_multivariable <- function(outcome, predictors) {
  dat <- transformed_parameter_summary %>%
    filter(.data$outcome == .env$outcome) %>%
    left_join(geochem_prepped, by = "Treatment_model") %>%
    arrange(match(Treatment_model, model_treatments))

  if (nrow(dat) != 5L || dplyr::n_distinct(dat$Treatment_model) != 5L) {
    stop(outcome, ": expected one row for each of the five treatments.")
  }

  x_se_names <- paste0(predictors, "_se")
  required <- c("y_mean_trans", "y_sd_trans", predictors, x_se_names)
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(outcome, ": missing model columns: ", paste(missing, collapse = ", "))
  }

  x_raw <- as.matrix(dat[, predictors, drop = FALSE])
  x_se_raw <- as.matrix(dat[, x_se_names, drop = FALSE])
  storage.mode(x_raw) <- "double"
  storage.mode(x_se_raw) <- "double"

  if (any(!is.finite(x_raw)) || any(!is.finite(x_se_raw)) ||
      any(!is.finite(dat$y_mean_trans)) || any(!is.finite(dat$y_sd_trans))) {
    stop(outcome, ": non-finite predictor or response values in model data.")
  }

  x_centres <- colMeans(x_raw)
  x_scales <- apply(x_raw, 2, sd)
  if (any(!is.finite(x_scales)) || any(x_scales <= 0)) {
    stop(outcome, ": at least one predictor has no between-treatment variation.")
  }

  x_obs <- sweep(sweep(x_raw, 2, x_centres, "-"), 2, x_scales, "/")
  x_se <- sweep(x_se_raw, 2, x_scales, "/")

  y_center <- mean(dat$y_mean_trans)
  y_scale <- sd(dat$y_mean_trans)
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop(outcome, ": transformed response has no between-treatment variation.")
  }
  y <- (dat$y_mean_trans - y_center) / y_scale
  y_se <- pmax(dat$y_sd_trans / y_scale, 1e-6)

  stan_data <- list(
    N = nrow(dat), P = length(predictors),
    y = as.numeric(y), y_se = as.numeric(y_se),
    x_obs = x_obs, x_se = x_se
  )

  model_id <- paste0(outcome, "__multivariable")
  sample_dir <- file.path(out_dir, paste0("cmdstan_csv_", model_id))
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  message("\nFitting ONE multivariable model for ", outcome, ": ",
          paste(predictors, collapse = " + "),
          " (five treatment-level observations).")

  fit <- stan_model$sample(
    data = stan_data,
    seed = 20261027,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    adapt_delta = 0.99,
    max_treedepth = 12,
    refresh = 500,
    output_dir = sample_dir
  )

  saveRDS(fit, file.path(out_dir, paste0("cmdstan_", model_id, ".rds")))
  draws <- posterior::as_draws_df(fit$draws(inc_warmup = FALSE))
  saveRDS(draws, file.path(out_dir, paste0("posterior_draws_", model_id, ".rds")))

  # Coefficients back-transformed to the transformed response units per original
  # predictor unit. Intercept is at all original predictor values equal to zero.
  beta_cols <- paste0("beta[", seq_along(predictors), "]")
  beta_std <- as.matrix(draws[, beta_cols, drop = FALSE])
  beta_orig <- sweep(beta_std, 2, y_scale / x_scales, "*")
  colnames(beta_orig) <- predictors

  coefficient_rows <- purrr::map_dfr(seq_along(predictors), function(j) {
    summarise_draws_simple(beta_std[, j]) %>%
      mutate(outcome = outcome, predictor = predictors[j],
             coefficient_scale = "standardised", .before = 1) %>%
      bind_rows(
        summarise_draws_simple(beta_orig[, j]) %>%
          mutate(outcome = outcome, predictor = predictors[j],
                 coefficient_scale = "transformed response units per original predictor unit",
                 .before = 1)
      )
  })
  readr::write_csv(coefficient_rows,
                   file.path(out_dir, paste0("coefficients_", model_id, ".csv")))

  smry <- fit$summary(variables = c("alpha", "beta", "sigma"))
  sampler <- as.data.frame(fit$sampler_diagnostics(format = "draws_df"))
  diagnostics <- tibble(
    model = model_id,
    outcome = outcome,
    predictors = paste(predictors, collapse = " + "),
    n_treatments = nrow(dat),
    n_predictors = length(predictors),
    residual_df_nominal = nrow(dat) - length(predictors) - 1L,
    max_Rhat = if (all(is.na(smry$rhat))) NA_real_ else max(smry$rhat, na.rm = TRUE),
    min_bulk_ESS = if (all(is.na(smry$ess_bulk))) NA_real_ else min(smry$ess_bulk, na.rm = TRUE),
    min_tail_ESS = if (all(is.na(smry$ess_tail))) NA_real_ else min(smry$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__ == 1, na.rm = TRUE),
    treedepth_hits = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
  )
  readr::write_csv(diagnostics,
                   file.path(out_dir, paste0("diagnostics_", model_id, ".csv")))

  # Posterior fitted means on the transformed scale, then back-transform every
  # draw to preserve asymmetry in the rate/fraction uncertainty intervals.
  mu_cols <- paste0("mu[", seq_len(nrow(dat)), "]")
  mu_std <- as.matrix(draws[, mu_cols, drop = FALSE])
  mu_trans <- mu_std * y_scale + y_center
  mu_original <- apply(mu_trans, 2, inverse_response_transform,
                       method = response_transform[[outcome]])
  observed_original <- parameter_draws %>%
    select(Treatment_model, all_of(outcome)) %>%
    group_by(Treatment_model) %>%
    summarise(
      observed_mean = mean(.data[[outcome]], na.rm = TRUE),
      observed_sd = sd(.data[[outcome]], na.rm = TRUE),
      observed_lower = mean(.data[[outcome]], na.rm = TRUE) - sd(.data[[outcome]], na.rm = TRUE),
      observed_upper = mean(.data[[outcome]], na.rm = TRUE) + sd(.data[[outcome]], na.rm = TRUE),
      .groups = "drop"
    )

  pred <- purrr::map_dfr(seq_len(nrow(dat)), function(i) {
    tibble(
      outcome = outcome,
      predictors = paste(predictors, collapse = " + "),
      Treatment_model = dat$Treatment_model[i],
      observed_mean = observed_original$observed_mean[
        match(dat$Treatment_model[i], observed_original$Treatment_model)],
      observed_sd = observed_original$observed_sd[
        match(dat$Treatment_model[i], observed_original$Treatment_model)],
      observed_lower = observed_original$observed_lower[
        match(dat$Treatment_model[i], observed_original$Treatment_model)],
      observed_upper = observed_original$observed_upper[
        match(dat$Treatment_model[i], observed_original$Treatment_model)],
      predicted_mean = mean(mu_original[, i]),
      predicted_median = median(mu_original[, i]),
      predicted_lower = unname(quantile(mu_original[, i], 0.025)),
      predicted_upper = unname(quantile(mu_original[, i], 0.975))
    )
  })
  readr::write_csv(pred,
                   file.path(out_dir, paste0("predictions_", model_id, ".csv")))

  list(
    fit = fit, draws = draws, data = dat,
    coefficients = coefficient_rows,
    diagnostics = diagnostics,
    predictions = pred
  )
}

# ---------------------------------------------------------------------------
# 7. Fit the five multivariable models
# ---------------------------------------------------------------------------

fits <- list()
for (outcome in names(model_specs)) {
  fits[[outcome]] <- fit_one_multivariable(outcome, model_specs[[outcome]])
}

saveRDS(fits, file.path(out_dir, "all_bayesian_rate_transfer_multivariable_fits.rds"))

all_coefficients <- purrr::map_dfr(fits, "coefficients")
all_diagnostics <- purrr::map_dfr(fits, "diagnostics")
all_predictions <- purrr::map_dfr(fits, "predictions") %>%
  mutate(Treatment_model = factor(Treatment_model, levels = model_treatments),
         outcome = factor(outcome, levels = names(model_specs))) %>%
  arrange(outcome, Treatment_model)

readr::write_csv(all_coefficients,
                 file.path(out_dir, "all_primary_multivariable_coefficients.csv"))
readr::write_csv(all_diagnostics,
                 file.path(out_dir, "all_primary_multivariable_diagnostics.csv"))
readr::write_csv(all_predictions,
                 file.path(out_dir, "all_primary_multivariable_predictions.csv"))

# ---------------------------------------------------------------------------
# 8. Observed vs fitted plots: one plot for EACH of the five parameters
#
# Plot design:
# - x-axis is the five treatment-level observations, ordered from low to high
#   by the observed posterior mean for that parameter.
# - coloured circles are observed posterior means; error bars are +/- 1 SD
#   across the saved SoilR posterior draws.
# - black open diamonds are model-fitted means; their error bars are 95%
#   credible intervals for the fitted mean.
# - observed and fitted points are nudged slightly apart horizontally.
# - no ribbon, connecting line, or predictor-vs-response plot is drawn.
# ---------------------------------------------------------------------------

plot_data <- all_predictions %>%
  mutate(
    outcome = as.character(outcome),
    outcome_label = case_when(
      outcome %in% c("kf", "ki", "ks") ~ paste0(outcome, " (decomposition rate)"),
      TRUE ~ paste0(outcome, " (transfer fraction)")
    ),
    Treatment_model = factor(Treatment_model, levels = model_treatments)
  )

treatment_palette <- c(
  Control = "#333333",
  Lime = "#D8A900",
  Bolsdorfer = "#0072B2",
  Huhnerberg = "#009E73",
  Eifelgold = "#CC79A7"
)

plot_observed_fitted <- function(dat) {
  this_outcome <- unique(dat$outcome)[1]
  this_label <- unique(dat$outcome_label)[1]
  predictor_text <- unique(dat$predictors)[1]

  dat <- dat %>%
    arrange(observed_mean) %>%
    mutate(
      x_order = row_number(),
      x_observed = x_order - 0.10,
      x_predicted = x_order + 0.10,
      Treatment_model = factor(as.character(Treatment_model),
                               levels = names(treatment_palette))
    )

  p <- ggplot() +
    geom_errorbar(
      data = dat,
      aes(x = x_observed, ymin = observed_lower, ymax = observed_upper,
          colour = Treatment_model),
      width = 0.08, linewidth = 0.7, show.legend = FALSE
    ) +
    geom_point(
      data = dat,
      aes(x = x_observed, y = observed_mean, colour = Treatment_model,
          shape = "Observed mean"),
      size = 3.2, stroke = 0.9
    ) +
    geom_errorbar(
      data = dat,
      aes(x = x_predicted, ymin = predicted_lower, ymax = predicted_upper),
      width = 0.08, linewidth = 0.7, colour = "black"
    ) +
    geom_point(
      data = dat,
      aes(x = x_predicted, y = predicted_mean, shape = "Model-fitted mean"),
      size = 3.0, colour = "black", fill = "white", stroke = 1.0
    ) +
    scale_colour_manual(name = "Treatment", values = treatment_palette,
                        breaks = names(treatment_palette), drop = FALSE) +
    scale_shape_manual(
      name = "Point type",
      values = c("Observed mean" = 16, "Model-fitted mean" = 23)
    ) +
    scale_x_continuous(
      breaks = dat$x_order,
      labels = as.character(dat$Treatment_model),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      title = paste("Observed vs model-fitted", this_label),
      subtitle = paste0("Predictors in joint model: ", predictor_text),
      x = "Treatment observation (ordered by observed mean, low to high)",
      y = this_outcome,
      caption = paste(
        "Observed: coloured circle with +/- 1 posterior SD.",
        "Model-fitted: open black diamond with 95% credible interval.",
        "Only five independent treatment-level observations are available."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  p
}

plot_list <- setNames(vector("list", length(model_specs)), names(model_specs))
for (outcome in names(model_specs)) {
  dat <- plot_data %>% filter(.data$outcome == outcome)
  if (nrow(dat) != 5L) {
    warning("Plot for ", outcome, " has ", nrow(dat),
            " rows; expected five treatment observations.")
    next
  }
  p <- plot_observed_fitted(dat)
  plot_list[[outcome]] <- p
  ggsave(file.path(out_dir, paste0("observed_vs_fitted_", outcome, ".png")),
         p, width = 9, height = 5.8, dpi = 300)
  ggsave(file.path(out_dir, paste0("observed_vs_fitted_", outcome, ".pdf")),
         p, width = 9, height = 5.8)
}

valid_plots <- plot_list[!vapply(plot_list, is.null, logical(1))]
if (length(valid_plots) > 0L) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(valid_plots, ncol = 2)
  } else {
    # Save individual parameter plots regardless of whether patchwork is installed.
    combined_plot <- NULL
  }
  if (!is.null(combined_plot)) {
    ggsave(file.path(out_dir, "observed_vs_fitted_all_parameters.png"),
           combined_plot, width = 16, height = 13, dpi = 300)
    ggsave(file.path(out_dir, "observed_vs_fitted_all_parameters.pdf"),
           combined_plot, width = 16, height = 13)
  }
}

saveRDS(plot_list, file.path(out_dir, "observed_vs_fitted_plot_objects.rds"))
# Save all five observed-versus-fitted plots as PNG
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Save all five observed-versus-fitted plots as PNG
# =============================================================================

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

outcomes_to_plot <- c(
  "kf", "ki", "ks", "alpha_fi", "alpha_is"
)

for (outcome_name in outcomes_to_plot) {
  
  dat <- plot_data %>%
    dplyr::filter(.data[["outcome"]] == outcome_name)
  
  message(
    "Processing ", outcome_name,
    ": ", nrow(dat), " rows"
  )
  
  if (nrow(dat) != 5L) {
    warning(
      "Skipping ", outcome_name,
      ": expected 5 rows, found ", nrow(dat)
    )
    next
  }
  
  p <- plot_observed_fitted(dat)
  
  filename <- file.path(
    out_dir,
    paste0("observed_vs_fitted_", outcome_name, ".png")
  )
  
  ggplot2::ggsave(
    filename = filename,
    plot = p,
    width = 9,
    height = 5.8,
    units = "in",
    dpi = 300
  )
  
  message("Saved: ", filename)
}

# Verify that the PNGs exist
list.files(
  out_dir,
  pattern = "^observed_vs_fitted_.*\\.png$",
  full.names = TRUE
)
message(
  "\nFinished fitting the five primary multivariable models.",
  "\nOne model per response; all listed predictors enter simultaneously.",
  "\nNo predictor-vs-response plots were created.",
  "\nOutput directory: ", out_dir,
  "\nKey outputs:",
  "\n- all_primary_multivariable_coefficients.csv",
  "\n- all_primary_multivariable_diagnostics.csv",
  "\n- all_primary_multivariable_predictions.csv",
  "\n- observed_vs_fitted_kf.png/.pdf",
  "\n- observed_vs_fitted_ki.png/.pdf",
  "\n- observed_vs_fitted_ks.png/.pdf",
  "\n- observed_vs_fitted_alpha_fi.png/.pdf",
  "\n- observed_vs_fitted_alpha_is.png/.pdf",
  "\n- observed_vs_fitted_all_parameters.png/.pdf (if patchwork is installed)",
  "\n- all_bayesian_rate_transfer_multivariable_fits.rds",
  "\n\nInterpret cautiously: n = 5 treatment-level observations for each model,",
  "\nwith two or three predictors. Priors regularise the coefficients but cannot",
  "\nreplace independent treatment replication. Treat all results as exploratory."
)
