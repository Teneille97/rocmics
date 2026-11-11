# =============================================================================
# MALLE ERW: SoilR posterior treatment comparisons
#
# PURPOSE:
#   Compare SoilR decomposition rates, transfer coefficients, pool mean ages
#   and system mean transit time between treatments and Control.
#
# IMPORTANT STATISTICAL NOTE:
#   The existing SoilR MCMC models were fitted separately at treatment level.
#   Posterior draws are NOT independent experimental replicates.
#
#   Therefore, this script calculates posterior treatment contrasts directly.
#   It does not treat posterior draws as additional plot-level observations.
#
#   A joint hierarchical SoilR model would require plot-level SoilR fits or
#   refitting the underlying SoilR model using plot-level data.
#
# TREATMENTS:
#   Control
#   Lime
#   Huhnerberg (HB 50)
#   Eifelgold (EG 50)
#   Bolsdorfer (BD 50)
#
# OUTCOMES:
#   Decomposition rates: kf, ki, ks
#   Transfer coefficients: alpha_fi, alpha_is
#   Pool mean ages: fPOM_age, oPOM_age, MAOM_age
#   System mean transit time: system_transit_time
#
# PRIMARY COMPARISONS:
#   Lime - Control
#   HB 50 - Control
#   EG 50 - Control
#   BD 50 - Control
#
# OUTPUTS:
#   Posterior treatment means
#   Primary treatment-versus-Control contrasts
#   Posterior contrast draws
#   Grouped forest plots
#   Grouped posterior distribution plots
#   Model/persistence draw summaries
#   RDS backup
#
# No SoilR model refitting.
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(here)
  library(SoilR)
  library(posterior)
  
})


set.seed(20261110)


# =============================================================================
# 2. PATHS AND SETTINGS
# =============================================================================

mcmc_dir <- here::here("mod_runs")

out_dir <- here::here(
  "outputs",
  "Bayesian_SoilR_Treatment_Comparisons"
)

models_dir <- file.path(out_dir, "posterior_draws")
plots_dir <- file.path(out_dir, "plots")
csv_dir <- file.path(out_dir, "csv_files")

dir.create(
  models_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  csv_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# SoilR settings matching the existing persistence script

burnin <- 20000L

n_posterior_draws <- 1000L

age_grid <- seq(0, 500, by = 1)

mean_C_inputs <- 100


# =============================================================================
# 3. TREATMENTS AND OUTCOMES
# =============================================================================

mcmc_files <- c(
  
  Control = "MCMC_3p_Control.Rdata",
  
  Lime = "MCMC_3p_Lime.Rdata",
  
  Bolsdorfer = "MCMC_3p_Bolsdorfer.Rdata",
  
  Huhnerberg = "MCMC_3p_Huhnerberg.Rdata",
  
  Eifelgold = "MCMC_3p_Eifelgold.Rdata"
  
)


model_treatments <- names(mcmc_files)


treatment_labels <- c(
  
  Control = "Control",
  
  Lime = "Lime",
  
  Bolsdorfer = "BD 50",
  
  Huhnerberg = "HB 50",
  
  Eifelgold = "EG 50"
  
)


plot_treatment_order <- c(
  "Control",
  "Lime",
  "HB 50",
  "EG 50",
  "BD 50"
)


# SoilR parameters

soilR_parameters <- c(
  
  "kf",
  "ki",
  "ks",
  "alpha_fi",
  "alpha_is"
  
)


# Derived persistence outcomes

age_outcomes <- c(
  
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
  
)


all_outcomes <- c(
  soilR_parameters,
  age_outcomes
)


# Group outcomes for plotting

plot_groups <- list(
  
  Decomposition_Rates = c(
    "kf",
    "ki",
    "ks"
  ),
  
  Transfer_Coefficients = c(
    "alpha_fi",
    "alpha_is"
  ),
  
  Pool_Mean_Ages = c(
    "fPOM_age",
    "oPOM_age",
    "MAOM_age"
  ),
  
  System_Transit_Time = c(
    "system_transit_time"
  )
  
)


outcome_labels <- outcome_display_labels


# Primary comparisons ONLY

primary_comparisons <- tibble(
  
  Treatment = c(
    "Lime",
    "Huhnerberg",
    "Eifelgold",
    "Bolsdorfer"
  ),
  
  Reference = "Control"
  
) %>%
  
  mutate(
    
    Comparison = paste(
      unname(treatment_labels[Treatment]),
      "minus Control"
    )
    
  )


print(primary_comparisons)



# =============================================================================
# PLOT DISPLAY LABELS, UNITS AND TREATMENT COLOURS
# =============================================================================

treatment_colours <- c(
  "Control" = "#777777",
  "Lime" = "#78A641",
  "HB 50" = "#C65A9E",
  "EG 50" = "#D99A16",
  "BD 50" = "#2878B5"
)

treatment_display <- c(
  "Control" = "Control",
  "Lime" = "Lime",
  "HB 50" = "Hühnerberg",
  "EG 50" = "Eifelgold",
  "BD 50" = "Bolsdorfer"
)

group_display <- c(
  Decomposition_Rates = "Decomposition rates",
  Transfer_Coefficients = "Transfer coefficients",
  Pool_Mean_Ages = "Pool mean ages",
  System_Transit_Time = "System mean transit time"
)

outcome_display_labels <- c(
  kf = "fPOM decomposition rate, kf (yr⁻¹)",
  ki = "oPOM decomposition rate, ki (yr⁻¹)",
  ks = "MAOM decomposition rate, ks (yr⁻¹)",
  alpha_fi = "fPOM → oPOM transfer coefficient, αfi (no units)",
  alpha_is = "oPOM → MAOM transfer coefficient, αis (no units)",
  fPOM_age = "fPOM mean age (years)",
  oPOM_age = "oPOM mean age (years)",
  MAOM_age = "MAOM mean age (years)",
  system_transit_time = "System mean transit time (years)"
)

display_group <- function(x) {
  out <- unname(group_display[x])
  ifelse(is.na(out), gsub("_", " ", x), out)
}


# =============================================================================
# 4. CHECK MCMC FILES
# =============================================================================

missing_files <- file.path(
  mcmc_dir,
  mcmc_files
) %>%
  setdiff(
    file.path(mcmc_dir, mcmc_files)[
      file.exists(file.path(mcmc_dir, mcmc_files))
    ]
  )


if (length(missing_files) > 0) {
  
  stop(
    "Missing MCMC files:\n",
    paste(missing_files, collapse = "\n")
  )
  
}


cat("\nAll five SoilR MCMC files found.\n")


# =============================================================================
# 5. SOILR PERSISTENCE FUNCTIONS
# =============================================================================

trapz <- function(x, y) {
  
  sum(
    diff(x) *
      (head(y, -1) + tail(y, -1)) / 2,
    na.rm = TRUE
  )
  
}


mean_age <- function(age, density) {
  
  density <- as.numeric(density)
  
  density[
    !is.finite(density) | density < 0
  ] <- 0
  
  area <- trapz(age, density)
  
  if (!is.finite(area) || area <= 0) {
    return(NA_real_)
  }
  
  trapz(
    age,
    age * density
  ) / area
  
}


build_A_u <- function(pars) {
  
  pars <- as.numeric(pars)
  
  if (length(pars) < 5L) {
    stop("Five SoilR parameters are required.")
  }
  
  kf <- pars[1]
  
  ki <- pars[2]
  
  ks <- pars[3]
  
  alpha_fi <- pars[4]
  
  alpha_is <- pars[5]
  
  
  A <- diag(
    -c(kf, ki, ks)
  )
  
  
  A[2, 1] <- kf * alpha_fi
  
  A[3, 2] <- ki * alpha_is
  
  
  list(
    
    A = A,
    
    u = matrix(
      c(mean_C_inputs, 0, 0),
      ncol = 1
    )
    
  )
  
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
  
  
  if (
    is.null(sa$poolAgeDensity) ||
    ncol(sa$poolAgeDensity) < 3L ||
    is.null(tt$transitTimeDensity)
  ) {
    
    stop(
      "Unexpected SoilR systemAge/transitTime output structure."
    )
    
  }
  
  
  tibble(
    
    fPOM_age = mean_age(
      age_grid,
      sa$poolAgeDensity[, 1]
    ),
    
    oPOM_age = mean_age(
      age_grid,
      sa$poolAgeDensity[, 2]
    ),
    
    MAOM_age = mean_age(
      age_grid,
      sa$poolAgeDensity[, 3]
    ),
    
    system_transit_time = mean_age(
      age_grid,
      tt$transitTimeDensity
    )
    
  )
  
}


# =============================================================================
# 6. EXTRACT POSTERIOR PARAMETER DRAWS
# =============================================================================

get_par_draws <- function(
    path,
    n_draws = n_posterior_draws
) {
  
  if (!file.exists(path)) {
    
    stop("MCMC file not found: ", path)
    
  }
  
  
  e <- new.env(parent = emptyenv())
  
  obj_names <- load(
    path,
    envir = e
  )
  
  
  if (
    !"MCMC" %in% obj_names ||
    is.null(e$MCMC$pars)
  ) {
    
    stop(
      basename(path),
      " must contain MCMC$pars."
    )
    
  }
  
  
  pars <- as.data.frame(
    e$MCMC$pars
  )
  
  
  if (
    all(soilR_parameters %in% names(pars))
  ) {
    
    pars <- pars[
      ,
      soilR_parameters,
      drop = FALSE
    ]
    
  } else {
    
    if (ncol(pars) < 5L) {
      
      stop(
        "Fewer than five parameter columns in ",
        basename(path)
      )
      
    }
    
    warning(
      "Parameter names not all identified in ",
      basename(path),
      ". Using first five columns in assumed order."
    )
    
    pars <- pars[
      ,
      seq_along(soilR_parameters),
      drop = FALSE
    ]
    
    names(pars) <- soilR_parameters
    
  }
  
  
  if (nrow(pars) <= burnin) {
    
    stop(
      "Burn-in exceeds available draws in ",
      basename(path)
    )
    
  }
  
  
  pars <- pars[
    (burnin + 1L):nrow(pars),
    ,
    drop = FALSE
  ]
  
  
  pars <- pars[
    complete.cases(pars) &
      apply(
        pars,
        1,
        function(x) all(is.finite(x))
      ),
    ,
    drop = FALSE
  ]
  
  
  if (nrow(pars) == 0L) {
    
    stop(
      "No valid posterior draws in ",
      basename(path)
    )
    
  }
  
  
  idx <- sample.int(
    
    nrow(pars),
    
    size = n_draws,
    
    replace = n_draws > nrow(pars)
    
  )
  
  
  pars[
    idx,
    ,
    drop = FALSE
  ]
  
}


# =============================================================================
# 7. GENERATE PARAMETER AND PERSISTENCE POSTERIOR DRAWS
# =============================================================================

persistence_draws <- purrr::map_dfr(
  
  model_treatments,
  
  function(tr) {
    
    cat(
      "\nProcessing treatment:",
      tr,
      "\n"
    )
    
    
    path <- file.path(
      mcmc_dir,
      mcmc_files[[tr]]
    )
    
    
    pars <- get_par_draws(path)
    
    
    # Derive ages and transit times from each sampled parameter combination
    
    ages <- purrr::map_dfr(
      
      seq_len(nrow(pars)),
      
      function(i) {
        
        one_persistence_draw(
          pars[i, ]
        )
        
      }
      
    )
    
    
    bind_cols(
      
      tibble(
        
        Treatment_model = tr,
        
        Treatment = unname(
          treatment_labels[tr]
        ),
        
        draw = seq_len(nrow(pars))
        
      ),
      
      pars %>%
        as_tibble(),
      
      ages
      
    )
    
  }
  
) %>%
  
  mutate(
    
    Treatment = factor(
      Treatment,
      levels = plot_treatment_order
    )
    
  ) %>%
  
  arrange(
    Treatment,
    draw
  )


# Save full posterior draws

write_csv(
  
  persistence_draws,
  
  file.path(
    models_dir,
    "SoilR_parameter_and_persistence_posterior_draws.csv"
  )
  
)


saveRDS(
  
  persistence_draws,
  
  file.path(
    models_dir,
    "SoilR_parameter_and_persistence_posterior_draws.rds"
  )
  
)


cat(
  "\nPosterior draws generated successfully.\n"
)


# =============================================================================
# 8. SUMMARISE POSTERIOR TREATMENT MEANS
# =============================================================================

posterior_treatment_means <- persistence_draws %>%
  
  pivot_longer(
    
    cols = all_of(all_outcomes),
    
    names_to = "Response",
    
    values_to = "Value"
    
  ) %>%
  
  filter(
    is.finite(Value)
  ) %>%
  
  group_by(
    Treatment,
    Response
  ) %>%
  
  summarise(
    
    n_draws = n(),
    
    Posterior_Mean = mean(Value),
    
    Posterior_Median = median(Value),
    
    Posterior_SD = sd(Value),
    
    Lower_95 = quantile(
      Value,
      0.025
    ),
    
    Upper_95 = quantile(
      Value,
      0.975
    ),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    Response = factor(
      Response,
      levels = all_outcomes
    )
    
  )


write_csv(
  
  posterior_treatment_means,
  
  file.path(
    csv_dir,
    "SoilR_Bayesian_posterior_treatment_means.csv"
  )
  
)


# =============================================================================
# 9. GENERATE PRIMARY TREATMENT-VERSUS-CONTROL CONTRASTS
# =============================================================================

# For each posterior iteration:
#
#   Difference = Treatment posterior draw - Control posterior draw
#
# Independent draws are sampled from the treatment-specific posterior
# distributions. This estimates uncertainty in the difference between
# treatment-level SoilR quantities.
#
# These are posterior contrasts, not hierarchical treatment effects.


n_contrast_draws <- 10000L


get_outcome_draws <- function(
    treatment_name,
    outcome
) {
  
  persistence_draws %>%
    
    filter(
      Treatment_model == treatment_name
    ) %>%
    
    pull(
      all_of(outcome)
    ) %>%
    
    .[is.finite(.)]
  
}


contrast_draws_list <- list()

idx <- 1L


for (outcome in all_outcomes) {
  
  cat(
    "\nCalculating contrasts for:",
    outcome,
    "\n"
  )
  
  
  control_draws <- get_outcome_draws(
    "Control",
    outcome
  )
  
  
  if (length(control_draws) < 100L) {
    
    stop(
      "Insufficient Control posterior draws for ",
      outcome
    )
    
  }
  
  
  for (tr in primary_comparisons$Treatment) {
    
    
    treatment_draws <- get_outcome_draws(
      tr,
      outcome
    )
    
    
    if (length(treatment_draws) < 100L) {
      
      warning(
        "Insufficient draws for ",
        tr,
        " and ",
        outcome
      )
      
      next
      
    }
    
    
    # Independently sample from the two posterior distributions
    
    control_sample <- sample(
      
      control_draws,
      
      size = n_contrast_draws,
      
      replace = TRUE
      
    )
    
    
    treatment_sample <- sample(
      
      treatment_draws,
      
      size = n_contrast_draws,
      
      replace = TRUE
      
    )
    
    
    difference <- treatment_sample - control_sample
    
    
    contrast_draws_list[[idx]] <- tibble(
      
      Response = outcome,
      
      Treatment_model = tr,
      
      Treatment = unname(
        treatment_labels[tr]
      ),
      
      Comparison = paste(
        unname(treatment_labels[tr]),
        "minus Control"
      ),
      
      draw = seq_len(n_contrast_draws),
      
      Treatment_draw = treatment_sample,
      
      Control_draw = control_sample,
      
      Difference = difference
      
    )
    
    
    idx <- idx + 1L
    
  }
  
}


contrast_draws <- bind_rows(
  contrast_draws_list
)


write_csv(
  
  contrast_draws,
  
  file.path(
    models_dir,
    "SoilR_primary_treatment_contrast_draws.csv"
  )
  
)


# =============================================================================
# 10. SUMMARISE PRIMARY CONTRASTS
# =============================================================================

primary_contrasts <- contrast_draws %>%
  
  group_by(
    Response,
    Treatment_model,
    Treatment,
    Comparison
  ) %>%
  
  summarise(
    
    n_draws = n(),
    
    Difference_mean = mean(Difference),
    
    Difference_median = median(Difference),
    
    Difference_SD = sd(Difference),
    
    Lower_95 = quantile(
      Difference,
      0.025
    ),
    
    Upper_95 = quantile(
      Difference,
      0.975
    ),
    
    Probability_Treatment_higher = mean(
      Difference > 0
    ),
    
    Probability_Treatment_lower = mean(
      Difference < 0
    ),
    
    Probability_difference_positive = mean(
      Difference > 0
    ),
    
    Probability_difference_negative = mean(
      Difference < 0
    ),
    
    # Probability that the difference exceeds 0.5 absolute units
    # is deliberately NOT calculated here because units differ
    # substantially among outcomes.
    
    CrI_excludes_zero = (
      Lower_95 > 0 |
        Upper_95 < 0
    ),
    
    Evidence_Flag = if_else(
      CrI_excludes_zero,
      "*",
      ""
    ),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    Response = factor(
      Response,
      levels = all_outcomes
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Lime",
        "HB 50",
        "EG 50",
        "BD 50"
      )
    )
    
  ) %>%
  
  arrange(
    Response,
    Treatment
  )


write_csv(
  
  primary_contrasts,
  
  file.path(
    csv_dir,
    "SoilR_primary_treatment_vs_Control_contrasts.csv"
  )
  
)


# =============================================================================
# 11. EXPORT A COMPACT RESULTS TABLE
# =============================================================================

primary_results_compact <- primary_contrasts %>%
  
  select(
    
    Response,
    
    Treatment,
    
    Comparison,
    
    Difference_mean,
    
    Difference_median,
    
    Lower_95,
    
    Upper_95,
    
    Probability_Treatment_higher,
    
    Probability_Treatment_lower,
    
    Evidence_Flag
    
  ) %>%
  
  mutate(
    
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
    
  )


write_csv(
  
  primary_results_compact,
  
  file.path(
    csv_dir,
    "SoilR_primary_contrasts_compact.csv"
  )
  
)


print(
  primary_results_compact,
  n = Inf
)


# =============================================================================
# 12. PLOT THE POSTERIOR TREATMENT DISTRIBUTIONS
# =============================================================================

# Each violin represents the posterior uncertainty of the treatment-level
# SoilR estimate. It does NOT represent variation among individual plots.


posterior_plot_data <- persistence_draws %>%
  
  pivot_longer(
    
    cols = all_of(all_outcomes),
    
    names_to = "Response",
    
    values_to = "Value"
    
  ) %>%
  
  filter(
    is.finite(Value)
  ) %>%
  
  mutate(
    
    Treatment = factor(
      Treatment,
      levels = plot_treatment_order
    ),
    
    Response = factor(
      Response,
      levels = all_outcomes,
      labels = unname(
        outcome_labels[all_outcomes]
      )
    )
    
  )


make_posterior_violin <- function(
    data,
    responses,
    group_name
) {
  
  selected_labels <- unname(
    outcome_labels[responses]
  )
  
  
  plot_data <- data %>%
    
    filter(
      as.character(Response) %in% selected_labels
    ) %>%
    
    mutate(
      
      Response = factor(
        Response,
        levels = selected_labels
      )
      
    )
  
  
  if (nrow(plot_data) == 0L) {
    
    return(NULL)
    
  }
  
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = Treatment,
      y = Value,
      fill = Treatment
    )
    
  ) +
    
    geom_violin(
      
      alpha = 0.65,
      
      trim = FALSE,
      
      scale = "width",
      
      linewidth = 0.35
      
    ) +
    
    geom_boxplot(
      
      width = 0.12,
      
      colour = "black",
      
      alpha = 0.8,
      
      outlier.shape = NA,
      
      linewidth = 0.35
      
    ) +
    
    facet_wrap(
      
      ~ Response,
      
      scales = "free_y",
      
      ncol = 2
      
    ) +
    
    scale_x_discrete(labels = treatment_display) +
    
    scale_fill_manual(
      values = treatment_colours,
      breaks = plot_treatment_order,
      labels = treatment_display[plot_treatment_order],
      drop = FALSE
    ) +
    
    
    labs(
      
      title = paste("SoilR posterior distributions:", display_group(group_name)),
      
      subtitle = paste(
        "Treatment-specific posterior distributions;",
        "not independent plot-level observations"
      ),
      
      x = NULL,
      
      y = "Posterior value (units specified in each panel heading)"
      
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      
      legend.position = "none",
      
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      ),
      
      plot.title = element_text(
        face = "bold"
      ),
      
      panel.grid.minor = element_blank()
      
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    
    plots_dir,
    
    paste0(
      "Posterior_Violin_",
      group_name,
      ".png"
    )
    
  )
  
  
  ggsave(
    
    filename,
    
    p,
    
    width = 13,
    
    height = height,
    
    dpi = 300,
    
    bg = "white"
    
  )
  
  
  ggsave(
    
    sub(
      "\\.png$",
      ".pdf",
      filename
    ),
    
    p,
    
    width = 13,
    
    height = height,
    
    bg = "white"
    
  )
  
  
  print(p)
  
  p
  
}


posterior_violin_plots <- list()


for (group_name in names(plot_groups)) {
  
  posterior_violin_plots[[group_name]] <-
    
    make_posterior_violin(
      
      data = posterior_plot_data,
      
      responses = plot_groups[[group_name]],
      
      group_name = group_name
      
    )
  
}


# =============================================================================
# 13. GROUPED PRIMARY TREATMENT-VERSUS-CONTROL FOREST PLOTS
# =============================================================================

make_primary_forest_plot <- function(
    data,
    responses,
    group_name
) {
  
  
  plot_data <- data %>%
    
    filter(
      as.character(Response) %in% responses
    ) %>%
    
    mutate(
      
      Response = factor(
        as.character(Response),
        levels = responses,
        labels = unname(
          outcome_labels[responses]
        )
      ),
      
      Treatment = factor(
        as.character(Treatment),
        levels = rev(c("Lime", "HB 50", "EG 50", "BD 50")),
        labels = rev(c("Lime - Control", "HB 50 - Control", "EG 50 - Control", "BD 50 - Control"))
      )
      
    )
  
  
  if (nrow(plot_data) == 0L) {
    
    message(
      "Skipping forest plot: ",
      group_name
    )
    
    return(NULL)
    
  }
  
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = Difference_median,
      y = Treatment,
      colour = Treatment
    )
    
  ) +
    
    geom_vline(
      
      xintercept = 0,
      
      linetype = "dashed",
      
      colour = "grey40",
      
      linewidth = 0.5
      
    ) +
    
    geom_errorbar(
      
      aes(
        xmin = Lower_95,
        xmax = Upper_95
      ),
      
      orientation = "y",
      
      width = 0.18,
      
      linewidth = 0.65
      
    ) +
    
    geom_point(
      
      size = 2.5
      
    ) +
    
    scale_colour_manual(
      values = c(
        "Lime - Control" = "#78A641",
        "HB 50 - Control" = "#C65A9E",
        "EG 50 - Control" = "#D99A16",
        "BD 50 - Control" = "#2878B5"
      ),
      breaks = c("Lime - Control", "HB 50 - Control", "EG 50 - Control", "BD 50 - Control"),
      name = "Treatment"
    ) +
    
    facet_wrap(
      
      ~ Response,
      
      scales = "free_x",
      
      ncol = 2
      
    ) +
    
    labs(
      
      title = paste("SoilR primary treatment contrasts:", display_group(group_name)),
      
      subtitle = paste(
        "Treatment minus Control;",
        "posterior median differences and 95% credible intervals"
      ),
      
      x = "Treatment difference (units specified in each panel heading)",
      
      y = NULL,
      
      caption = paste(
        "* 95% credible interval excludes zero.",
        "Contrasts propagate uncertainty from the saved treatment-level",
        "SoilR posterior distributions."
      )
      
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      
      axis.text.y = element_text(
        size = 9
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      
      plot.subtitle = element_text(
        size = 9
      ),
      
      plot.caption = element_text(
        size = 8
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.spacing = grid::unit(
        1,
        "lines"
      )
      
    )
  
  
  height <- max(
    
    4.5,
    
    ceiling(length(responses) / 2) * 3
    
  )
  
  
  filename <- file.path(
    
    plots_dir,
    
    paste0(
      "Primary_Forest_",
      group_name,
      ".png"
    )
    
  )
  
  
  ggsave(
    
    filename,
    
    p,
    
    width = 13,
    
    height = height,
    
    dpi = 300,
    
    bg = "white"
    
  )
  
  
  ggsave(
    
    sub(
      "\\.png$",
      ".pdf",
      filename
    ),
    
    p,
    
    width = 13,
    
    height = height,
    
    bg = "white"
    
  )
  
  
  print(p)
  
  p
  
}


primary_forest_plots <- list()


for (group_name in names(plot_groups)) {
  
  primary_forest_plots[[group_name]] <-
    
    make_primary_forest_plot(
      
      data = primary_contrasts,
      
      responses = plot_groups[[group_name]],
      
      group_name = group_name
      
    )
  
}


# =============================================================================
# 14. PLOT POSTERIOR CONTRAST DISTRIBUTIONS
# =============================================================================

# These plots show the full posterior distribution of treatment differences,
# rather than only the median and credible interval.


contrast_plot_data <- contrast_draws %>%
  
  mutate(
    
    Treatment = factor(
      
      Treatment,
      
      levels = c(
        "Lime",
        "HB 50",
        "EG 50",
        "BD 50"
      )
      
    ),
    
    Response = factor(
      
      Response,
      
      levels = all_outcomes,
      
      labels = unname(
        outcome_labels[all_outcomes]
      )
      
    )
    
  )


make_contrast_distribution_plot <- function(
    data,
    responses,
    group_name
) {
  
  
  selected_labels <- unname(
    outcome_labels[responses]
  )
  
  
  plot_data <- data %>%
    
    filter(
      as.character(Response) %in% selected_labels
    ) %>%
    
    mutate(
      
      Response = factor(
        Response,
        levels = selected_labels
      )
      
    )
  
  
  if (nrow(plot_data) == 0L) {
    
    return(NULL)
    
  }
  
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = Difference,
      fill = Treatment
    )
    
  ) +
    
    geom_density(
      alpha = 0.35,
      linewidth = 0.55
    ) +
    
    geom_vline(
      
      xintercept = 0,
      
      linetype = "dashed",
      
      colour = "black"
      
    ) +
    
    facet_wrap(
      
      ~ Response,
      
      scales = "free",
      
      ncol = 2
      
    ) +
    
    scale_fill_manual(
      values = treatment_colours[c("Lime", "HB 50", "EG 50", "BD 50")],
      breaks = c("Lime", "HB 50", "EG 50", "BD 50"),
      labels = treatment_display[c("Lime", "HB 50", "EG 50", "BD 50")]
    ) +
    
    
    labs(
      
      title = paste("Posterior treatment differences:", display_group(group_name)),
      
      subtitle = "Each distribution represents treatment minus Control",
      
      x = "Posterior treatment difference (units specified in each panel heading)",
      
      y = "Density",
      
      fill = "Treatment"
      
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      
      strip.background = element_rect(
        fill = "grey92",
        colour = "grey40"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      plot.title = element_text(
        face = "bold"
      ),
      
      panel.grid.minor = element_blank()
      
    )
  
  
  height <- max(
    5,
    ceiling(length(responses) / 2) * 3
  )
  
  
  filename <- file.path(
    
    plots_dir,
    
    paste0(
      "Posterior_Contrast_Distributions_",
      group_name,
      ".png"
    )
    
  )
  
  
  ggsave(
    
    filename,
    
    p,
    
    width = 13,
    
    height = height,
    
    dpi = 300,
    
    bg = "white"
    
  )
  
  
  ggsave(
    
    sub(
      "\\.png$",
      ".pdf",
      filename
    ),
    
    p,
    
    width = 13,
    
    height = height,
    
    bg = "white"
    
  )
  
  
  print(p)
  
  p
  
}


contrast_distribution_plots <- list()


for (group_name in names(plot_groups)) {
  
  contrast_distribution_plots[[group_name]] <-
    
    make_contrast_distribution_plot(
      
      data = contrast_plot_data,
      
      responses = plot_groups[[group_name]],
      
      group_name = group_name
      
    )
  
}


# =============================================================================
# 15. STRONGEST TREATMENT DIFFERENCES
# =============================================================================

# Descriptive ranking by probability that a treatment differs in either
# direction from Control.
#
# This is NOT a significance test or a multiple-comparison correction.


contrast_strength <- primary_contrasts %>%
  
  mutate(
    
    Maximum_direction_probability = pmax(
      
      Probability_Treatment_higher,
      
      Probability_Treatment_lower
      
    ),
    
    Absolute_median_difference = abs(
      Difference_median
    )
    
  ) %>%
  
  arrange(
    Response,
    desc(Maximum_direction_probability)
  )


write_csv(
  
  contrast_strength,
  
  file.path(
    csv_dir,
    "SoilR_treatment_contrast_probability_ranking.csv"
  )
  
)


# =============================================================================
# 16. COMPARISON AND OUTPUT CHECKS
# =============================================================================

cat("\n\n============================================================\n")

cat("SOILR TREATMENT COMPARISON ANALYSIS COMPLETE\n")

cat("============================================================\n\n")


cat(
  "Treatments:",
  length(model_treatments),
  "\n"
)


cat(
  "Outcomes:",
  length(all_outcomes),
  "\n"
)


cat(
  "Posterior draws per treatment:",
  n_posterior_draws,
  "\n"
)


cat(
  "Posterior contrast draws per comparison:",
  n_contrast_draws,
  "\n"
)


cat(
  "Primary contrasts:",
  nrow(primary_contrasts),
  "\n"
)


cat(
  "\nPosterior treatment means:\n"
)

print(
  posterior_treatment_means,
  n = Inf
)


cat(
  "\nPrimary treatment-versus-Control contrasts:\n"
)

print(
  primary_results_compact,
  n = Inf
)


cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)


# =============================================================================
# 17. SAVE COMPLETE ANALYSIS
# =============================================================================

soilR_treatment_analysis <- list(
  
  settings = list(
    
    burnin = burnin,
    
    n_posterior_draws = n_posterior_draws,
    
    n_contrast_draws = n_contrast_draws,
    
    age_grid = age_grid,
    
    mean_C_inputs = mean_C_inputs
    
  ),
  
  persistence_draws = persistence_draws,
  
  posterior_treatment_means = posterior_treatment_means,
  
  contrast_draws = contrast_draws,
  
  primary_contrasts = primary_contrasts,
  
  contrast_strength = contrast_strength,
  
  outcomes = all_outcomes,
  
  plot_groups = plot_groups,
  
  treatments = model_treatments
  
)


saveRDS(
  
  soilR_treatment_analysis,
  
  file.path(
    out_dir,
    "SoilR_Bayesian_Treatment_Comparison_Complete.rds"
  )
  
)


cat("\nAll CSV, PNG, PDF and RDS outputs saved.\n")

cat("\n============================================================\n")
cat("IMPORTANT INTERPRETATION\n")
cat("============================================================\n")

cat(
  "\nThese results compare treatment-level SoilR posterior distributions.",
  "\nThe posterior draws quantify uncertainty in the fitted SoilR parameters.",
  "\nThey do not represent additional independent plot replicates.",
  "\nA joint hierarchical SoilR model requires plot-level model fitting",
  "\nor a joint model of the underlying plot-level observations.",
  "\n"
)