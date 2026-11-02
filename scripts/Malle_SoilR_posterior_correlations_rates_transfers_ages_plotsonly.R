
# =============================================================================
# MALLE ERW: PLOT-ONLY SCRIPT
# Load existing posterior correlation analysis and regenerate plots
#
# NO model fitting
# NO SoilR age/transit-time calculations
# NO correlation recalculation
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# -----------------------------------------------------------------------------
# 0. Load existing statistical results
# -----------------------------------------------------------------------------

out_dir <- here::here(
  "outputs",
  "soilR_posterior_correlations_age_transfer"
)

analysis_file <- file.path(
  out_dir,
  "soilR_rates_transfer_ages_correlation_analysis.rds"
)

analysis <- readRDS(analysis_file)

# Restore previously calculated objects
correlation_summary <- analysis$correlation_summary
correlation_draws <- analysis$correlation_draws
persistence_draws <- analysis$persistence_draws
geochem_summary <- analysis$geochemistry_summary

all_outcomes <- analysis$outcomes
valid_predictors <- analysis$predictors

# Existing definitions from the original analysis
methods <- c("spearman", "kendall", "pearson")

age_outcomes <- c(
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

message("Existing statistical results successfully loaded.")

# -----------------------------------------------------------------------------
# 1. Outcome labels
# -----------------------------------------------------------------------------

outcome_labels <- c(
  kf = "kf: fPOM decomposition rate",
  ki = "ki: oPOM decomposition rate",
  ks = "ks: MAOM decomposition rate",
  alpha_fi = "alpha_fi: fPOM -> oPOM transfer",
  alpha_is = "alpha_is: oPOM -> MAOM transfer",
  fPOM_age = "fPOM mean age",
  oPOM_age = "oPOM mean age",
  MAOM_age = "MAOM mean age",
  system_transit_time = "System mean transit time"
)

# -----------------------------------------------------------------------------
# 2. Heatmap function
# -----------------------------------------------------------------------------

plot_heatmap <- function(
    method_name,
    outcome_subset,
    filename_stub,
    title_text
) {
  
  plot_dat <- correlation_summary %>%
    filter(
      method == method_name,
      outcome %in% outcome_subset
    ) %>%
    mutate(
      outcome = factor(
        outcome,
        levels = rev(outcome_subset),
        labels = rev(unname(outcome_labels[outcome_subset]))
      ),
      predictor = factor(
        predictor,
        levels = valid_predictors
      )
    )
  
  p <- ggplot(
    plot_dat,
    aes(
      x = predictor,
      y = outcome,
      fill = correlation_median
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.25
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = ifelse(
        method_name == "spearman",
        "Median rho",
        ifelse(
          method_name == "kendall",
          "Median tau",
          "Median r"
        )
      )
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Five treatment-level observations; median ",
        method_name,
        " correlation across sampled SoilR posterior draws"
      ),
      x = "Geochemical predictor",
      y = NULL,
      caption = paste(
        "Exploratory only: posterior draws propagate SoilR uncertainty,",
        "but do not increase the number of independent treatments."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        colour = "black"
      ),
      axis.text.x = element_text(
        colour = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.grid = element_blank()
    )
  
  width <- max(
    12,
    min(24, 0.28 * length(valid_predictors) + 5)
  )
  
  height <- max(
    5.5,
    0.48 * length(outcome_subset) + 2.5
  )
  
  ggsave(
    file.path(out_dir, paste0(filename_stub, ".png")),
    p,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    file.path(out_dir, paste0(filename_stub, ".pdf")),
    p,
    width = width,
    height = height,
    limitsize = FALSE,
    bg = "white"
  )
  
  return(p)
}

# -----------------------------------------------------------------------------
# 3. Generate heatmaps
# -----------------------------------------------------------------------------

for (method_name in methods) {
  
  method_label <- stringr::str_to_title(method_name)
  
  # Combined rates, transfers, ages and transit time
  p_all <- plot_heatmap(
    method_name,
    all_outcomes,
    paste0(
      "heatmap_", method_name,
      "_all_rates_transfers_ages"
    ),
    paste0(
      "SoilR rates, transfer coefficients and ages vs geochemistry: ",
      method_label
    )
  )
  
  # Ages and transit time only
  p_ages <- plot_heatmap(
    method_name,
    age_outcomes,
    paste0(
      "heatmap_", method_name,
      "_ages_transit_only"
    ),
    paste0(
      "SoilR pool ages and system transit time vs geochemistry: ",
      method_label
    )
  )
  
  print(p_all)
  print(p_ages)
}

# -----------------------------------------------------------------------------
# 4. System transit time: weak Spearman correlations only
# -----------------------------------------------------------------------------

system_transit_plot_data <- correlation_summary %>%
  filter(
    outcome == "system_transit_time",
    method == "spearman",
    abs(correlation_median) <= 0.5
  ) %>%
  mutate(
    predictor = forcats::fct_reorder(
      predictor,
      correlation_median
    )
  )

if (nrow(system_transit_plot_data) > 0) {
  
  p_system_transit_weak <- ggplot(
    system_transit_plot_data,
    aes(
      x = predictor,
      y = correlation_median
    )
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = correlation_q025,
        ymax = correlation_q975
      ),
      width = 0.2
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_hline(
      yintercept = c(-0.5, 0.5),
      linetype = "dotted",
      colour = "grey50"
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25)
    ) +
    labs(
      title = "System transit time vs geochemistry",
      subtitle = paste(
        "Spearman correlations with |median rho| <= 0.5;",
        "posterior median and 95% interval"
      ),
      x = "Geochemical predictor",
      y = "Spearman correlation (rho)",
      caption = paste(
        "Five independent treatment-level observations.",
        "Intervals propagate SoilR posterior uncertainty;",
        "they are not significance intervals."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.background = element_rect(
        fill = "white", colour = NA
      ),
      panel.background = element_rect(
        fill = "white", colour = NA
      ),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(
        face = "bold", colour = "black"
      ),
      plot.subtitle = element_text(colour = "black"),
      plot.caption = element_text(colour = "black")
    )
  
  ggsave(
    file.path(
      out_dir,
      "system_transit_time_spearman_absrho_le_0.5.png"
    ),
    p_system_transit_weak,
    width = 10,
    height = max(
      5,
      0.35 * nrow(system_transit_plot_data) + 2
    ),
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "system_transit_time_spearman_absrho_le_0.5.pdf"
    ),
    p_system_transit_weak,
    width = 10,
    height = max(
      5,
      0.35 * nrow(system_transit_plot_data) + 2
    ),
    bg = "white"
  )
  
  print(p_system_transit_weak)
  
} else {
  
  message(
    "No system transit time predictors with |median Spearman rho| <= 0.5."
  )
}

message("\nPlot-only workflow completed.")

# -----------------------------------------------------------------------------
# System transit time: correlation plot for weak correlations only
# -----------------------------------------------------------------------------

system_transit_plot_data <- correlation_summary %>%
  filter(
    outcome == "system_transit_time",
    method == "spearman",
    abs(correlation_median) <= 0.5
  ) %>%
  mutate(
    predictor = forcats::fct_reorder(
      predictor,
      correlation_median
    )
  )

if (nrow(system_transit_plot_data) > 0) {
  
  p_system_transit_weak <- ggplot(
    system_transit_plot_data,
    aes(
      x = predictor,
      y = correlation_median
    )
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = correlation_q025,
        ymax = correlation_q975
      ),
      width = 0.2
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_hline(
      yintercept = c(-0.5, 0.5),
      linetype = "dotted",
      colour = "grey50"
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25)
    ) +
    labs(
      title = "System transit time vs geochemistry",
      subtitle = paste(
        "Spearman correlations with |median rho| <= 0.5;",
        "posterior median and 95% interval"
      ),
      x = "Geochemical predictor",
      y = "Spearman correlation (rho)",
      caption = paste(
        "Five independent treatment-level observations.",
        "Intervals propagate SoilR posterior uncertainty;",
        "they are not significance intervals."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      plot.subtitle = element_text(colour = "black"),
      plot.caption = element_text(colour = "black")
    )
  
  ggsave(
    file.path(
      out_dir,
      "system_transit_time_spearman_absrho_le_0.5.png"
    ),
    p_system_transit_weak,
    width = 10,
    height = max(5, 0.35 * nrow(system_transit_plot_data) + 2),
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "system_transit_time_spearman_absrho_le_0.5.pdf"
    ),
    p_system_transit_weak,
    width = 10,
    height = max(5, 0.35 * nrow(system_transit_plot_data) + 2),
    bg = "white"
  )
  
  print(p_system_transit_weak)
  
} else {
  message(
    "No system transit time predictors with |median Spearman rho| <= 0.5."
  )
}


# -----------------------------------------------------------------------------
# SoilR ages and system transit time:
# Combined heatmap retaining predictors with strong correlations
# (|median Spearman rho| >= 0.5) for at least one outcome
# -----------------------------------------------------------------------------

# Outcomes to include
persistence_outcomes <- c(
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

# -----------------------------------------------------------------------------
# 1. Identify predictors with at least one strong correlation
# -----------------------------------------------------------------------------

strong_persistence_predictors <- correlation_summary %>%
  filter(
    outcome %in% persistence_outcomes,
    method == "spearman"
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

# -----------------------------------------------------------------------------
# 2. Retain ALL correlations for the selected predictors
# -----------------------------------------------------------------------------

persistence_heatmap_data <- correlation_summary %>%
  filter(
    outcome %in% persistence_outcomes,
    method == "spearman",
    predictor %in% strong_persistence_predictors
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = rev(persistence_outcomes),
      labels = rev(c(
        "fPOM mean age",
        "oPOM mean age",
        "MAOM mean age",
        "System mean transit time"
      ))
    ),
    predictor = factor(
      predictor,
      levels = valid_predictors
    )
  )

# -----------------------------------------------------------------------------
# 3. Plot combined heatmap:
#    Decomposition rates, transfer coefficients, pool ages and transit time
# -----------------------------------------------------------------------------

# All SoilR outcomes to include
soilR_outcomes <- c(
  "kf",
  "ki",
  "ks",
  "alpha_fi",
  "alpha_is",
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

# Human-readable outcome labels
soilR_outcome_labels <- c(
  "kf"                  = "fPOM decomposition rate (kf)",
  "ki"                  = "oPOM decomposition rate (ki)",
  "ks"                  = "MAOM decomposition rate (ks)",
  "alpha_fi"             = "fPOM -> oPOM transfer (alpha_fi)",
  "alpha_is"             = "oPOM -> MAOM transfer (alpha_is)",
  "fPOM_age"             = "fPOM mean age",
  "oPOM_age"             = "oPOM mean age",
  "MAOM_age"             = "MAOM mean age",
  "system_transit_time"  = "System mean transit time"
)

# -----------------------------------------------------------------------------
# 3.1 Identify geochemical predictors with at least one strong correlation
# -----------------------------------------------------------------------------

strong_soilR_predictors <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

# -----------------------------------------------------------------------------
# 3.2 Prepare combined heatmap data
#     Show ALL correlations for retained predictors
# -----------------------------------------------------------------------------

soilR_heatmap_data <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    predictor %in% strong_soilR_predictors
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = rev(soilR_outcomes),
      labels = rev(
        unname(soilR_outcome_labels[soilR_outcomes])
      )
    ),
    predictor = factor(
      predictor,
      levels = strong_soilR_predictors
    )
  )

# -----------------------------------------------------------------------------
# 3. Plot combined heatmap:
#    Decomposition rates, transfer coefficients, pool ages and transit time
# -----------------------------------------------------------------------------

# All SoilR outcomes to include
soilR_outcomes <- c(
  "kf",
  "ki",
  "ks",
  "alpha_fi",
  "alpha_is",
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

# Human-readable outcome labels
soilR_outcome_labels <- c(
  "kf"                  = "fPOM decomposition rate (kf)",
  "ki"                  = "oPOM decomposition rate (ki)",
  "ks"                  = "MAOM decomposition rate (ks)",
  "alpha_fi"             = "fPOM -> oPOM transfer (alpha_fi)",
  "alpha_is"             = "oPOM -> MAOM transfer (alpha_is)",
  "fPOM_age"             = "fPOM mean age",
  "oPOM_age"             = "oPOM mean age",
  "MAOM_age"             = "MAOM mean age",
  "system_transit_time"  = "System mean transit time"
)

# -----------------------------------------------------------------------------
# 3.1 Identify geochemical predictors with at least one strong correlation
# -----------------------------------------------------------------------------

strong_soilR_predictors <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

# -----------------------------------------------------------------------------
# 3.2 Prepare combined heatmap data
#     Show ALL correlations for retained predictors
# -----------------------------------------------------------------------------

soilR_heatmap_data <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    predictor %in% strong_soilR_predictors
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = rev(soilR_outcomes),
      labels = rev(
        unname(soilR_outcome_labels[soilR_outcomes])
      )
    ),
    predictor = factor(
      predictor,
      levels = strong_soilR_predictors
    )
  )
# -----------------------------------------------------------------------------
# 3.1 Identify strong geochemical predictors
#     EXCLUDING all FTIR-derived variables
# -----------------------------------------------------------------------------

# Identify FTIR-derived predictor names
ftir_predictors <- correlation_summary %>%
  distinct(predictor) %>%
  filter(
    grepl(
      "FTIR|ftir|DRIFT|drift|aromatic|carbonate|carbonyl|carboxyl|aliphatic|spectral|peak",
      predictor
    )
  ) %>%
  pull(predictor)

# Select predictors with at least one strong correlation
# across the nine SoilR outcomes, excluding FTIR variables
strong_soilR_predictors <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    !predictor %in% ftir_predictors
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

message(
  "\nExcluded ", length(ftir_predictors),
  " FTIR-derived predictors:"
)

print(ftir_predictors)


# -----------------------------------------------------------------------------
# 3. Plot combined heatmap:
#    Decomposition rates, transfer coefficients, pool ages and transit time
# -----------------------------------------------------------------------------

# All SoilR outcomes to include
soilR_outcomes <- c(
  "kf",
  "ki",
  "ks",
  "alpha_fi",
  "alpha_is",
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

# Human-readable outcome labels
soilR_outcome_labels <- c(
  "kf"                  = "fPOM decomposition rate (kf)",
  "ki"                  = "oPOM decomposition rate (ki)",
  "ks"                  = "MAOM decomposition rate (ks)",
  "alpha_fi"             = "fPOM -> oPOM transfer (alpha_fi)",
  "alpha_is"             = "oPOM -> MAOM transfer (alpha_is)",
  "fPOM_age"             = "fPOM mean age",
  "oPOM_age"             = "oPOM mean age",
  "MAOM_age"             = "MAOM mean age",
  "system_transit_time"  = "System mean transit time"
)

# -----------------------------------------------------------------------------
# 3.1 Identify strong geochemical predictors
#     EXCLUDING FTIR-derived variables
# -----------------------------------------------------------------------------

# Identify FTIR-derived predictor names
ftir_predictors <- correlation_summary %>%
  distinct(predictor) %>%
  filter(
    grepl(
      "FTIR|DRIFT|aromatic|carbonate|carbonyl|carboxyl|aliphatic|spectral|peak",
      predictor,
      ignore.case = TRUE
    )
  ) %>%
  pull(predictor)

# Identify retained predictors with at least one strong correlation
strong_soilR_predictors <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    !predictor %in% ftir_predictors
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

# Report excluded predictors
message(
  "\nExcluded FTIR-derived predictors:"
)

print(ftir_predictors)

message(
  "\nRetained non-FTIR predictors:"
)

print(strong_soilR_predictors)

# -----------------------------------------------------------------------------
# 3.2 Prepare combined heatmap data
#     Show all correlations for retained non-FTIR predictors
# -----------------------------------------------------------------------------

soilR_heatmap_data <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    predictor %in% strong_soilR_predictors,
    !predictor %in% ftir_predictors
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = rev(soilR_outcomes),
      labels = rev(
        unname(soilR_outcome_labels[soilR_outcomes])
      )
    ),
    predictor = factor(
      predictor,
      levels = strong_soilR_predictors
    )
  )

# Explicit check that no FTIR predictors remain
stopifnot(
  !any(
    as.character(soilR_heatmap_data$predictor) %in%
      ftir_predictors
  )
)

message(
  "\nNumber of non-FTIR predictors in heatmap: ",
  n_distinct(soilR_heatmap_data$predictor)
)


# =============================================================================
# 4. STRONG CORRELATIONS ONLY — EXCLUDING FTIR
# =============================================================================

# -----------------------------------------------------------------------------
# 4.1 Identify FTIR-derived predictors
# -----------------------------------------------------------------------------

ftir_predictors <- correlation_summary %>%
  distinct(predictor) %>%
  filter(
    grepl(
      "FTIR|DRIFT|aromatic|carbonate|carbonyl|carboxyl|aliphatic|spectral|peak",
      predictor,
      ignore.case = TRUE
    )
  ) %>%
  pull(predictor)

# -----------------------------------------------------------------------------
# 4.2 Identify non-FTIR predictors with at least one strong correlation
# -----------------------------------------------------------------------------

strong_soilR_predictors <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    !predictor %in% ftir_predictors
  ) %>%
  group_by(predictor) %>%
  summarise(
    max_abs_rho = max(
      abs(correlation_median),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    max_abs_rho >= 0.5
  ) %>%
  pull(predictor)

# -----------------------------------------------------------------------------
# 4.3 Retain ONLY strong correlations
# -----------------------------------------------------------------------------

significant_heatmap_data <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    predictor %in% strong_soilR_predictors,
    !predictor %in% ftir_predictors
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = rev(soilR_outcomes),
      labels = rev(
        unname(soilR_outcome_labels[soilR_outcomes])
      )
    ),
    predictor = factor(
      predictor,
      levels = strong_soilR_predictors
    )
  )

# -----------------------------------------------------------------------------
# 4.4 Plot strong correlations only
# -----------------------------------------------------------------------------

if (nrow(significant_heatmap_data) > 0) {
  
  p_significant_only <- ggplot(
    significant_heatmap_data,
    aes(
      x = predictor,
      y = outcome,
      fill = correlation_median
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.25
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = "Median rho"
    ) +
    labs(
      title = "Strong persistence parameter correlations with geochemistry",
      subtitle = paste(
        "Only correlations with |median Spearman rho| >= 0.5",
        "are displayed"
      ),
      x = "Geochemical predictor",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.key = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.grid = element_blank(),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      axis.ticks = element_line(
        colour = "black"
      ),
      plot.title = element_text(
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        colour = "black"
      ),
      plot.caption = element_text(
        colour = "black"
      ),
      legend.text = element_text(
        colour = "black"
      ),
      legend.title = element_text(
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = 65,
        hjust = 1,
        vjust = 1
      )
    )
  
  width <- max(
    12,
    min(
      24,
      0.28 * length(strong_soilR_predictors) + 5
    )
  )
  
  height <- 9
  
  # Save PNG
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_significant_only_no_FTIR.png"
    ),
    p_significant_only,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  # Save PDF
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_significant_only_no_FTIR.pdf"
    ),
    p_significant_only,
    width = width,
    height = height,
    limitsize = FALSE,
    bg = "white"
  )
  
  print(p_significant_only)
  
  message(
    "\nStrong non-FTIR predictors retained: ",
    length(strong_soilR_predictors)
  )
  
  message(
    "Strong correlation cells displayed: ",
    nrow(significant_heatmap_data)
  )
  
  print(strong_soilR_predictors)
  
} else {
  
  message(
    "No non-FTIR correlations met the |median Spearman rho| >= 0.5 threshold."
  )
}