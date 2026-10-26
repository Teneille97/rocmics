# =============================================================================
# Malle ERW: Bayesian SoilR persistence vs treatment-level geochemistry
#
# Run from the RStudio project root.
#
# IMPORTANT:
# - Only five independent SoilR treatment outcomes exist.
# - Plot-level geochemistry is summarised within treatment.
# - Never replicate one SoilR outcome across three plots and count those plots
#   as independent response observations.
# - All regressions are exploratory; four models each use only five rows.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(SoilR)
  library(brms)
  library(posterior)
  library(dplyr)
})

set.seed(20261026)

integrated_file <- here::here(
  "outputs", "integrated_geochemistry_for_bayesian.csv"
)
mcmc_dir <- here::here("mod_runs")
out_dir <- here::here("outputs", "bayesian_persistence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

burnin <- 20000L
n_posterior_draws <- 1000L
age_grid <- seq(0, 500, by = 1)
mean_C_inputs <- 100

# Edit filenames if your saved MCMC files differ.
mcmc_files <- c(
  Control = "MCMC_3p_Control.Rdata",
  Lime = "MCMC_3p_Lime.Rdata",
  Bolsdorfer = "MCMC_3p_Bolsdorfer.Rdata",
  Huhnerberg = "MCMC_3p_Huhnerberg.Rdata",
  Eifelgold = "MCMC_3p_Eifelgold.Rdata"
)

model_treatments <- c(
  "Control", "Lime", "Bolsdorfer", "Huhnerberg", "Eifelgold"
)

# ---------------------------------------------------------------------------
# 1. Treatment labels and integrated geochemistry
# ---------------------------------------------------------------------------

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

if (!file.exists(integrated_file)) {
  stop(
    "Cannot find ", integrated_file,
    ". Export the integrated dataset to ",
    "outputs/integrated_geochemistry_for_bayesian.csv."
  )
}

integrated <- readr::read_csv(
  integrated_file,
  show_col_types = FALSE
)

if (!all(c("sample", "Treatment") %in% names(integrated))) {
  stop("The integrated CSV must contain columns named sample and Treatment.")
}

# Map labels once, then retain only the five treatments represented by the
# five SoilR MCMC files. This excludes other Huhnerberg application rates.
integrated <- integrated %>%
  mutate(Treatment_model = normalise_treatment(Treatment)) %>%
  filter(
    Treatment %in% c(
      "Control",
      "Lime_2",
      "Bolsdorfer_50",
      "Huhnerberg_50",
      "Eifelgold_50"
    )
  )

if (anyNA(integrated$Treatment_model)) {
  stop(
    "Unmapped treatment labels: ",
    paste(
      unique(integrated$Treatment[is.na(integrated$Treatment_model)]),
      collapse = ", "
    )
  )
}

treatment_check <- integrated %>%
  count(Treatment, Treatment_model, name = "n_plots")
print(treatment_check, n = Inf)

if (nrow(treatment_check) != 5L || any(treatment_check$n_plots != 3L)) {
  stop(
    "Expected exactly five modeled treatments with three plots each. ",
    "Check the treatment labels and integrated CSV."
  )
}

# ---------------------------------------------------------------------------
# 2. Predictor specifications
# ---------------------------------------------------------------------------

model_specs <- list(
  MAOM_age = c("soil_pH", "Fe_ox_over_DCB", "Si_ox_over_DCB"),
  oPOM_age = c("soil_EC", "CEC_Ca", "GMD_h2o"),
  fPOM_age = c(
    "Mn_ox_over_DCB",
    "Tessier_Acetic_Mg_mg_g"
  ),
  system_transit_time = c("soil_pH", "Fe_ox_over_DCB", "CEC_Ca")
)

all_predictors <- unique(unlist(model_specs))

missing_predictors <- setdiff(all_predictors, names(integrated))
if (length(missing_predictors) > 0L) {
  stop(
    "Missing predictor columns in integrated CSV: ",
    paste(missing_predictors, collapse = ", "),
    ". Check exact names and edit model_specs."
  )
}

for (v in all_predictors) {
  integrated[[v]] <- suppressWarnings(as.numeric(integrated[[v]]))
}

# Treatment means, between-plot SDs, sample sizes and standard errors.
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

if (nrow(geochem_summary) != 5L ||
    !setequal(geochem_summary$Treatment_model, model_treatments)) {
  stop("The geochemistry summary does not contain the expected five treatments.")
}

readr::write_csv(
  geochem_summary,
  file.path(out_dir, "geochemistry_treatment_summary.csv")
)

# ---------------------------------------------------------------------------
# 3. Reconstruct SoilR persistence metrics from posterior parameter draws
# ---------------------------------------------------------------------------

# Parameter order assumed to match the supplied three-pool SoilR script:
# kf, ki, ks, alpha_fi, alpha_is.
build_A_u <- function(pars) {
  pars <- as.numeric(pars)
  if (length(pars) < 5L) {
    stop("Expected five SoilR parameters: kf, ki, ks, alpha_fi, alpha_is.")
  }
  
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[5]
  
  A <- diag(-c(kf, ki, ks))
  A[2, 1] <- kf * alpha_fi
  A[3, 2] <- ki * alpha_is
  
  list(
    A = A,
    u = matrix(c(mean_C_inputs, 0, 0), ncol = 1)
  )
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

mean_age <- function(age, density) {
  d <- as.numeric(density)
  d[!is.finite(d) | d < 0] <- 0
  
  area <- trapz(age, d)
  if (!is.finite(area) || area <= 0) return(NA_real_)
  
  trapz(age, age * d) / area
}

one_persistence_draw <- function(pars) {
  au <- build_A_u(pars)
  
  sa <- SoilR::systemAge(
    A = au$A,
    u = au$u,
    a = age_grid
  )
  tt <- SoilR::transitTime(
    A = au$A,
    u = au$u,
    a = age_grid
  )
  
  if (is.null(sa$poolAgeDensity) ||
      ncol(sa$poolAgeDensity) < 3L ||
      is.null(tt$transitTimeDensity)) {
    stop(
      "Unexpected output from SoilR::systemAge() or transitTime(). ",
      "Check the installed SoilR version and the returned object structure."
    )
  }
  
  tibble(
    MAOM_age = mean_age(age_grid, sa$poolAgeDensity[, 3]),
    oPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 2]),
    fPOM_age = mean_age(age_grid, sa$poolAgeDensity[, 1]),
    system_transit_time = mean_age(age_grid, tt$transitTimeDensity)
  )
}

get_par_draws <- function(path) {
  if (!file.exists(path)) {
    stop("MCMC file not found: ", path)
  }
  
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
    if (ncol(pars) < 5L) {
      stop("Fewer than five parameter columns in ", path)
    }
    warning(
      "Using the first five MCMC columns in assumed parameter order for ",
      basename(path), ": ",
      paste(names(pars)[1:5], collapse = ", ")
    )
    pars <- pars[, 1:5, drop = FALSE]
  }
  
  if (nrow(pars) <= burnin) {
    stop(
      "burnin = ", burnin, " leaves no draws in ", basename(path),
      ". Check chain length and adjust burnin."
    )
  }
  
  pars <- pars[(burnin + 1L):nrow(pars), , drop = FALSE]
  pars <- pars[
    complete.cases(pars) &
      apply(pars, 1, function(x) all(is.finite(x))),
    ,
    drop = FALSE
  ]
  
  if (nrow(pars) == 0L) {
    stop("No finite posterior draws in ", path)
  }
  
  pars[
    sample(seq_len(nrow(pars)), min(n_posterior_draws, nrow(pars))),
    ,
    drop = FALSE
  ]
}

# Run full SoilR persistence calculation across posterior draws
persistence_draws_list <- list()

for (tr in names(mcmc_files)) {
  mcmc_path <- file.path(mcmc_dir, mcmc_files[[tr]])
  message("Processing MCMC file for treatment: ", tr, " (", mcmc_path, ")")
  
  pars <- get_par_draws(mcmc_path)
  
  draw_res <- purrr::map_dfr(seq_len(nrow(pars)), function(i) {
    one_persistence_draw(pars[i, ])
  }) %>%
    mutate(Treatment_model = tr, draw = row_number())
  
  persistence_draws_list[[tr]] <- draw_res
}

persistence_draws <- bind_rows(persistence_draws_list)

# ---------------------------------------------------------------------------
# 4. Assemble analysis dataset
# ---------------------------------------------------------------------------

soilr_summary <- persistence_draws %>%
  pivot_longer(
    cols = c(MAOM_age, oPOM_age, fPOM_age, system_transit_time),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(Treatment_model, outcome) %>%
  summarise(
    y_mean = mean(value, na.rm = TRUE),
    y_se   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  soilr_summary,
  file.path(out_dir, "soilR_persistence_summary.csv")
)

# Reshape geochem summary for merging with outcomes
geochem_long <- geochem_summary %>%
  pivot_longer(
    cols = -Treatment_model,
    names_to = c("predictor", "stat"),
    names_pattern = "^(.*)__(mean|sd|n)$"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(se = sd / sqrt(n))

geochem_wide_means <- geochem_long %>%
  select(Treatment_model, predictor, mean) %>%
  pivot_wider(names_from = predictor, values_from = mean)

geochem_wide_ses <- geochem_long %>%
  select(Treatment_model, predictor, se) %>%
  pivot_wider(names_from = predictor, values_from = se, names_prefix = "") %>%
  rename_with(~ paste0(.x, "_se"), -Treatment_model)

geochem_prepped <- left_join(geochem_wide_means, geochem_wide_ses, by = "Treatment_model")

model_data_long <- soilr_summary %>%
  left_join(geochem_prepped, by = "Treatment_model")

readr::write_csv(
  model_data_long,
  file.path(out_dir, "treatment_level_model_data_long.csv")
)

# Optional session backup
dir.create("outputs/bayesian_persistence/session_backup",
           recursive = TRUE, showWarnings = FALSE)

saveRDS(
  model_data_long,
  "outputs/bayesian_persistence/session_backup/model_data_long.rds"
)

saveRDS(
  geochem_summary,
  "outputs/bayesian_persistence/session_backup/geochem_summary.rds"
)

# Commented out reading saved data for clean first-time run
# model_data_long <- readRDS(
#   here::here(
#     "outputs", "bayesian_persistence",
#     "session_backup", "model_data_long.rds"
#   )
# )

# ---------------------------------------------------------------------------
# 5. Bayesian measurement-error regressions using cmdstanr directly
# ---------------------------------------------------------------------------
#
# ANALYSIS:
# - Fits all four outcomes specified in model_specs.
# - Uses exactly five independent treatment-level observations per outcome.
# - Treats geochemistry treatment means as uncertain measurements.
# - Incorporates the posterior SD of each SoilR-derived response.
# - Standardises predictors AND responses for numerically sensible priors.
# - Reports coefficients on standardised and original scales.
#
# Model:
#   x_true[n,k] ~ Normal(x_observed[n,k], x_SE[n,k])
#   y[n] ~ Normal(alpha + sum(beta[k] * x_true[n,k]),
#                 sqrt(y_SE[n]^2 + sigma^2))
#
# Priors on standardised scales:
#   alpha ~ Normal(0, 1)
#   beta  ~ Normal(0, 0.5)
#   sigma ~ Exponential(1)
# ---------------------------------------------------------------------------

# Confirm the real analysis inputs are available.
stopifnot(
  exists("model_data_long"),
  exists("model_specs"),
  exists("out_dir")
)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr is required. Install it from the Stan R-universe if needed.")
}

if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop(
    "CmdStan is not configured. Check cmdstanr::cmdstan_path() ",
    "and cmdstanr::check_cmdstan_toolchain()."
  )
}

# Save the actual Stan model to disk so the model's purpose and specification
# remain clear when you return to this project later.
stan_file <- file.path(out_dir, "bayesian_measurement_error.stan")


stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  vector<lower=0>[N] y_se;
  matrix[N, K] x_obs;
  matrix<lower=0>[N, K] x_se;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
  matrix[N, K] z_x;
}
transformed parameters {
  matrix[N, K] x_true;

  for (n in 1:N) {
    for (k in 1:K) {
      x_true[n, k] = x_obs[n, k] + x_se[n, k] * z_x[n, k];
    }
  }
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  sigma ~ exponential(1);

  to_vector(z_x) ~ std_normal();

  for (n in 1:N) {
    y[n] ~ normal(
      alpha + dot_product(beta, to_vector(x_true[n])),
      sqrt(square(y_se[n]) + square(sigma))
    );
  }
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;

  for (n in 1:N) {
    mu[n] = alpha + dot_product(beta, to_vector(x_true[n]));
    y_rep[n] = normal_rng(
      mu[n],
      sqrt(square(y_se[n]) + square(sigma))
    );
  }
}
"
writeLines(stan_code, stan_file)

stan_model <- cmdstanr::cmdstan_model(
  stan_file = stan_file,
  quiet = FALSE
)
  
  # Compile the actual model. The compiled executable is retained with the
  # project output, rather than using a disposable diagnostic model.
  stan_model <- cmdstanr::cmdstan_model(
    stan_file = stan_file,
    quiet = FALSE
  )
  
  # Summarise a vector of posterior draws.
  summarise_parameter <- function(x) {
    tibble(
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      q2.5 = unname(stats::quantile(x, 0.025)),
      q97.5 = unname(stats::quantile(x, 0.975)),
      Pr_positive = mean(x > 0),
      Pr_negative = mean(x < 0)
    )
  }
  
  fit_one <- function(outcome, predictors) {
    
    # Select the five treatment-level responses for this outcome.
    dat <- model_data_long %>%
      filter(.data$outcome == .env$outcome) %>%
      arrange(match(Treatment_model, model_treatments))
    
    message(
      "\nFitting real model: ", outcome,
      "\nTreatment rows: ", nrow(dat),
      "\nPredictors: ", paste(predictors, collapse = ", "),
      "\nTreatments: ", paste(dat$Treatment_model, collapse = ", ")
    )
    
    if (nrow(dat) != 5L ||
        dplyr::n_distinct(dat$Treatment_model) != 5L ||
        !setequal(dat$Treatment_model, model_treatments)) {
      stop(
        outcome,
        ": expected exactly one row for each of the five model treatments."
      )
    }
    
    required <- c(
      "y_mean", "y_se",
      predictors,
      paste0(predictors, "_se")
    )
    
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
      stop(outcome, ": missing columns: ", paste(missing, collapse = ", "))
    }
    
    if (!all(vapply(dat[required], is.numeric, logical(1)))) {
      stop(outcome, ": all response, predictor and SE columns must be numeric.")
    }
    
    if (!all(vapply(dat[required], function(x) all(is.finite(x)), logical(1)))) {
      stop(outcome, ": non-finite values detected in the model data.")
    }
    
    if (any(dat$y_se <= 0)) {
      stop(
        outcome,
        ": response SDs must be positive. Check the SoilR posterior draws."
      )
    }
    
    predictor_se <- as.matrix(dat[paste0(predictors, "_se")])
    
    if (any(predictor_se < 0)) {
      stop(outcome, ": predictor SEs cannot be negative.")
    }
    
    # Standardise predictors and their SEs using the same scale.
    x_means <- vapply(
      predictors,
      function(p) mean(dat[[p]]),
      numeric(1)
    )
    
    x_sds <- vapply(
      predictors,
      function(p) sd(dat[[p]]),
      numeric(1)
    )
    
    if (any(!is.finite(x_sds)) || any(x_sds <= 0)) {
      stop(outcome, ": at least one predictor has no between-treatment variation.")
    }
    
    x_obs <- sweep(
      as.matrix(dat[predictors]),
      MARGIN = 2,
      STATS = x_means,
      FUN = "-"
    )
    
    x_obs <- sweep(x_obs, 2, x_sds, "/")
    
    x_se <- sweep(
      as.matrix(dat[paste0(predictors, "_se")]),
      MARGIN = 2,
      STATS = x_sds,
      FUN = "/"
    )
    
    # Standardise the response too, so the priors are on a sensible scale.
    y_center <- mean(dat$y_mean)
    y_scale <- sd(dat$y_mean)
    
    if (!is.finite(y_scale) || y_scale <= 0) {
      stop(outcome, ": response means have no variation.")
    }
    
    y_std <- (dat$y_mean - y_center) / y_scale
    y_se_std <- dat$y_se / y_scale
    
    # Stan data: five treatment rows and the specified predictors.
    stan_data <- list(
      N = nrow(dat),
      K = length(predictors),
      y = as.numeric(y_std),
      y_se = as.numeric(y_se_std),
      x_obs = unname(x_obs),
      x_se = unname(x_se)
    )
    
    # Keep all sampling CSVs in a permanent output folder.
    sample_dir <- file.path(out_dir, paste0("cmdstan_csv_", outcome))
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Four-chain posterior sampling for the actual outcome.
    fit <- stan_model$sample(
      data = stan_data,
      seed = 20261026,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.99,
      max_treedepth = 12,
      refresh = 500,
      output_dir = sample_dir
    )
    
    # Save the fitted object and posterior draws immediately after sampling.
    saveRDS(
      fit,
      file.path(out_dir, paste0("cmdstan_", outcome, ".rds"))
    )
    
    draws <- posterior::as_draws_df(
      fit$draws(inc_warmup = FALSE)
    )
    
    saveRDS(
      draws,
      file.path(out_dir, paste0("posterior_draws_", outcome, ".rds"))
    )
    
    # Save a standard CmdStan summary for alpha, beta and sigma.
    parameter_summary <- fit$summary(
      variables = c("alpha", "beta", "sigma")
    )
    
    readr::write_csv(
      parameter_summary,
      file.path(out_dir, paste0("cmdstan_summary_", outcome, ".csv"))
    )
    
    # Convert standardised coefficients to original response/predictor units.
    # beta_original = beta_standardised * SD(response) / SD(predictor)
    beta_original_draws <- lapply(seq_along(predictors), function(k) {
      draws[[paste0("beta[", k, "]")]] * y_scale / x_sds[k]
    })
    names(beta_original_draws) <- predictors
    
    # The original-unit intercept depends on the centring of all predictors.
    alpha_original_draws <- y_center + y_scale * draws$alpha
    
    for (p in predictors) {
      alpha_original_draws <- alpha_original_draws -
        beta_original_draws[[p]] * x_means[p]
    }
    
    coefficient_rows <- list()
    
    coefficient_rows[["Intercept_standardised"]] <- summarise_parameter(
      draws$alpha
    ) %>%
      mutate(
        outcome = outcome,
        parameter = "Intercept",
        scale = "standardised"
      )
    
    coefficient_rows[["Intercept_original"]] <- summarise_parameter(
      alpha_original_draws
    ) %>%
      mutate(
        outcome = outcome,
        parameter = "Intercept",
        scale = "original response units"
      )
    
    for (p in predictors) {
      k <- match(p, predictors)
      
      coefficient_rows[[paste0(p, "_standardised")]] <-
        summarise_parameter(draws[[paste0("beta[", k, "]")]]) %>%
        mutate(
          outcome = outcome,
          parameter = p,
          scale = "standardised"
        )
      
      coefficient_rows[[paste0(p, "_original")]] <-
        summarise_parameter(beta_original_draws[[p]]) %>%
        mutate(
          outcome = outcome,
          parameter = p,
          scale = "original response units per predictor unit"
        )
    }
    
    coefficients <- bind_rows(coefficient_rows) %>%
      select(
        outcome, parameter, scale,
        mean, median, sd, q2.5, q97.5,
        Pr_positive, Pr_negative
      )
    
    readr::write_csv(
      coefficients,
      file.path(out_dir, paste0("coefficients_", outcome, ".csv"))
    )
    
    # Posterior predictive check:
    # compare observed treatment-level estimates with replicated responses
    # drawn from the fitted model. This is not a placeholder or test dataset.
    ppc_rows <- lapply(seq_len(nrow(dat)), function(i) {
      yrep_name <- paste0("y_rep[", i, "]")
      yrep_original <- y_center + y_scale * draws[[yrep_name]]
      
      tibble(
        Treatment_model = dat$Treatment_model[i],
        observed = dat$y_mean[i],
        observed_se = dat$y_se[i],
        predictive_median = median(yrep_original),
        predictive_lower = unname(quantile(yrep_original, 0.025)),
        predictive_upper = unname(quantile(yrep_original, 0.975))
      )
    })
    
    ppc_data <- bind_rows(ppc_rows)
    
    ppc_plot <- ggplot(
      ppc_data,
      aes(y = reorder(Treatment_model, observed))
    ) +
      geom_linerange(
        aes(xmin = predictive_lower, xmax = predictive_upper),
        linewidth = 0.8
      ) +
      geom_point(
        aes(x = predictive_median, shape = "Posterior predictive median"),
        size = 2.5
      ) +
      geom_errorbarh(
        aes(
          xmin = observed - observed_se,
          xmax = observed + observed_se
        ),
        height = 0.15
      ) +
      geom_point(aes(x = observed), size = 2.5) +
      labs(
        title = paste("Posterior predictive check:", outcome),
        x = "Persistence metric (original units)",
        y = "Treatment",
        shape = NULL
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    ggsave(
      file.path(out_dir, paste0("pp_check_", outcome, ".png")),
      ppc_plot,
      width = 8,
      height = 5,
      dpi = 300
    )
    
    # Trace plots are optional for visualisation, not for model fitting.
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      trace_draws <- fit$draws(
        variables = c("alpha", "beta", "sigma"),
        format = "draws_array"
      )
      
      trace_plot <- bayesplot::mcmc_trace(trace_draws)
      
      ggsave(
        file.path(out_dir, paste0("trace_", outcome, ".png")),
        trace_plot,
        width = 9,
        height = 7,
        dpi = 300
      )
    }
    
    # Keep enough information for downstream diagnostics and reproducibility.
    list(
      fit = fit,
      draws = draws,
      data = dat,
      coefficients = coefficients,
      ppc = ppc_data,
      x_means = x_means,
      x_sds = x_sds,
      y_center = y_center,
      y_scale = y_scale
    )
  }
  
  # Fit all four real models. Each model uses its own actual predictor set.
  fits <- list()
  
  for (nm in names(model_specs)) {
    fits[[nm]] <- fit_one(nm, model_specs[[nm]])
  }
  
  saveRDS(
    fits,
    file.path(out_dir, "all_bayesian_persistence_fits.rds")
  )
  
  # ---------------------------------------------------------------------------
  # 6. Sampling diagnostics for all four models
  # ---------------------------------------------------------------------------
  
  diagnostics <- purrr::map_dfr(names(fits), function(nm) {
    
    fit <- fits[[nm]]$fit
    
    s <- fit$summary(
      variables = c("alpha", "beta", "sigma")
    )
    
    sampler <- as.data.frame(
      fit$sampler_diagnostics(format = "draws_df")
    )
    
    tibble(
      model = nm,
      n_treatments = nrow(fits[[nm]]$data),
      max_Rhat = if (all(is.na(s$rhat))) {
        NA_real_
      } else {
        max(s$rhat, na.rm = TRUE)
      },
      min_bulk_ESS = if (all(is.na(s$ess_bulk))) {
        NA_real_
      } else {
        min(s$ess_bulk, na.rm = TRUE)
      },
      min_tail_ESS = if (all(is.na(s$ess_tail))) {
        NA_real_
      } else {
        min(s$ess_tail, na.rm = TRUE)
      },
      divergences = sum(sampler$divergent__ == 1, na.rm = TRUE),
      treedepth_hits = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
    )
  })
  
  readr::write_csv(
    diagnostics,
    file.path(out_dir, "model_diagnostics.csv")
  )
  
  print(diagnostics, n = Inf)
  
  # ---------------------------------------------------------------------------
  # 7. Completion message
  # ---------------------------------------------------------------------------
  
  message(
    "\nFinished the four real Bayesian measurement-error regressions.",
    "\nOutputs saved to: ", out_dir,
    "\nKey files:",
    "\n- bayesian_measurement_error.stan",
    "\n- coefficients_*.csv",
    "\n- cmdstan_summary_*.csv",
    "\n- posterior_draws_*.rds",
    "\n- cmdstan_*.rds",
    "\n- model_diagnostics.csv",
    "\n- pp_check_*.png",
    "\n- trace_*.png (if bayesplot is installed)",
    "\n- all_bayesian_persistence_fits.rds",
    "\n\nInterpret cautiously: each regression has five independent treatment-level",
    "\nobservations and three predictors, so posterior uncertainty will be substantial."
  )
  
  
  # =============================================================================
  # 8. OBSERVED VS BAYESIAN-PREDICTED PERSISTENCE PLOTS
  # =============================================================================
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # Output folder
  plot_dir <- file.path(out_dir, "observed_vs_predicted")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Model labels for plot titles and axes
  model_labels <- c(
    MAOM_age = "MAOM mean age",
    oPOM_age = "oPOM mean age",
    fPOM_age = "fPOM mean age",
    system_transit_time = "System transit time"
  )
  
  # Colours
  observed_colour <- "black"
  predicted_colour <- "#2878B5"
  ribbon_colour <- "#2878B5"
  
  # ---------------------------------------------------------------------------
  # Function to extract posterior predictions and create a plot
  # ---------------------------------------------------------------------------
  
  
  plot_observed_predicted <- function(model_name, fit_obj) {
    
    dat <- fit_obj$data
    
    # Accept either treatment column name
    if ("Treatment_model" %in% names(dat)) {
      dat$Treatment <- as.character(dat$Treatment_model)
    } else if ("Treatment" %in% names(dat)) {
      dat$Treatment <- as.character(dat$Treatment)
    } else {
      stop(
        model_name,
        ": no Treatment_model or Treatment column found. Available columns: ",
        paste(names(dat), collapse = ", ")
      )
    }
    
    # Check response columns
    required_cols <- c("Treatment", "y_mean", "y_se")
    
    if (!all(required_cols %in% names(dat))) {
      stop(
        model_name,
        ": missing required columns. Available columns: ",
        paste(names(dat), collapse = ", ")
      )
    }
    
    # Extract posterior draws
    posterior_draws <- as.data.frame(
      fit_obj$fit$draws(format = "draws_df")
    )
    
    mu_cols <- grep("^mu\\[[0-9]+\\]$", names(posterior_draws), value = TRUE)
    
    if (length(mu_cols) != nrow(dat)) {
      stop(
        model_name, ": found ", length(mu_cols),
        " mu columns but ", nrow(dat),
        " treatment rows. Check model-data ordering."
      )
    }
    
    mu_matrix <- as.matrix(posterior_draws[, mu_cols, drop = FALSE])
    
    # Convert predictions back to original response units
    mu_original <- mu_matrix * fit_obj$y_scale + fit_obj$y_center
    
    # Summarise predictions for each treatment
    pred_summary <- dplyr::bind_rows(
      lapply(seq_len(nrow(dat)), function(i) {
        tibble::tibble(
          Treatment = dat$Treatment[i],
          observed = dat$y_mean[i],
          observed_se = dat$y_se[i],
          predicted_median = median(mu_original[, i], na.rm = TRUE),
          predicted_lower = unname(quantile(
            mu_original[, i], 0.025, na.rm = TRUE
          )),
          predicted_upper = unname(quantile(
            mu_original[, i], 0.975, na.rm = TRUE
          ))
        )
      })
    )
    
    pred_summary$Treatment <- factor(
      pred_summary$Treatment,
      levels = dat$Treatment
    )
    
    # Plot observed estimates and posterior predictions
    p <- ggplot2::ggplot(pred_summary, ggplot2::aes(x = Treatment)) +
      
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = predicted_lower,
          ymax = predicted_upper,
          group = 1
        ),
        fill = "#2878B5",
        alpha = 0.20
      ) +
      
      ggplot2::geom_line(
        ggplot2::aes(y = predicted_median, group = 1),
        colour = "#2878B5",
        linewidth = 1
      ) +
      
      ggplot2::geom_point(
        ggplot2::aes(y = predicted_median),
        colour = "#2878B5",
        size = 3
      ) +
      
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = observed - observed_se,
          ymax = observed + observed_se
        ),
        width = 0.12,
        colour = "black",
        linewidth = 0.7,
        na.rm = TRUE
      ) +
      
      ggplot2::geom_point(
        ggplot2::aes(y = observed),
        colour = "black",
        size = 3
      ) +
      
      ggplot2::labs(
        title = paste(
          "Observed vs Bayesian-predicted",
          model_labels[[model_name]]
        ),
        subtitle = paste(
          "Black = observed SoilR estimate ± 1 SD;",
          "blue = posterior median prediction;",
          "shading = 95% credible interval"
        ),
        x = "Treatment",
        y = model_labels[[model_name]]
      ) +
      
      ggplot2::theme_bw(base_size = 12) +
      
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(
          angle = 35, hjust = 1, vjust = 1
        ),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    # Save plots
    ggplot2::ggsave(
      file.path(plot_dir, paste0(model_name, "_observed_vs_predicted.png")),
      plot = p, width = 9, height = 6, dpi = 300
    )
    
    ggplot2::ggsave(
      file.path(plot_dir, paste0(model_name, "_observed_vs_predicted.pdf")),
      plot = p, width = 9, height = 6
    )
    
    list(plot = p, predictions = pred_summary)
  }
  
  
  # ---------------------------------------------------------------------------
  # Create plots for all fitted models
  # ---------------------------------------------------------------------------
  
  observed_predicted_plots <- list()
  
  for (nm in names(model_specs)) {
    
    if (!nm %in% names(fits)) {
      warning("Skipping ", nm, ": no fitted model found.")
      next
    }
    
    message("Plotting observed vs predicted: ", nm)
    
    observed_predicted_plots[[nm]] <- plot_observed_predicted(
      model_name = nm,
      fit_obj = fits[[nm]]
    )
    
    print(observed_predicted_plots[[nm]]$plot)
  }
  
  # Save prediction summaries for all models
  prediction_summary <- bind_rows(
    lapply(names(observed_predicted_plots), function(nm) {
      
      observed_predicted_plots[[nm]]$predictions %>%
        mutate(model = nm, .before = 1)
      
    })
  )
  
  readr::write_csv(
    prediction_summary,
    file.path(plot_dir, "observed_vs_predicted_summary.csv")
  )
  
  message("Observed vs predicted plots saved to: ", plot_dir)
  