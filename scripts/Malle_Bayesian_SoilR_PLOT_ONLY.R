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

# =============================================================================
# PLOT-ONLY MODE: RESTORE THE SAVED COMPLETE ANALYSIS
# =============================================================================
# This script does NOT load the original MCMC_3p_*.Rdata files, call SoilR
# persistence functions, generate posterior draws, or calculate new contrasts.
# It restores the saved analysis object created by the full script.

saved_analysis_file <- file.path(
  out_dir,
  "SoilR_Bayesian_Treatment_Comparison_Complete.rds"
)

if (!file.exists(saved_analysis_file)) {
  stop(
    "Cannot find the saved SoilR analysis object:\n",
    saved_analysis_file,
    "\n\nRun the full SoilR treatment-comparison script once first."
  )
}

soilR_treatment_analysis <- readRDS(saved_analysis_file)

required_objects <- c(
  "persistence_draws",
  "posterior_treatment_means",
  "contrast_draws",
  "primary_contrasts",
  "contrast_strength",
  "plot_groups",
  "treatments",
  "settings"
)
missing_objects <- setdiff(required_objects, names(soilR_treatment_analysis))
if (length(missing_objects) > 0) {
  stop(
    "The saved RDS is missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    ". Please regenerate it using the updated full SoilR script."
  )
}

persistence_draws <- soilR_treatment_analysis$persistence_draws
posterior_treatment_means <- soilR_treatment_analysis$posterior_treatment_means
contrast_draws <- soilR_treatment_analysis$contrast_draws
primary_contrasts <- soilR_treatment_analysis$primary_contrasts
contrast_strength <- soilR_treatment_analysis$contrast_strength
plot_groups <- soilR_treatment_analysis$plot_groups
model_treatments <- soilR_treatment_analysis$treatments

burnin <- soilR_treatment_analysis$settings$burnin
n_posterior_draws <- soilR_treatment_analysis$settings$n_posterior_draws
n_contrast_draws <- soilR_treatment_analysis$settings$n_contrast_draws
age_grid <- soilR_treatment_analysis$settings$age_grid
mean_C_inputs <- soilR_treatment_analysis$settings$mean_C_inputs

all_outcomes <- soilR_treatment_analysis$outcomes
outcome_labels <- outcome_display_labels

cat("\nLoaded saved SoilR posterior draws and contrasts. No models will be fitted.\n")
cat("Treatments restored: ", paste(model_treatments, collapse = ", "), "\n", sep = "")


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



cat("\nSOILR PLOT-ONLY RUN COMPLETE\n")
cat("Plots saved in: ", plots_dir, "\n", sep = "")