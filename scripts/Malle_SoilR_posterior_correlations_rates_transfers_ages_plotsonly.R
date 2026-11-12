# =============================================================================
# MALLE ERW: PLOT-ONLY SCRIPT
# Load existing posterior correlation analysis and regenerate plots
#
# NO model fitting
# NO SoilR age/transit-time calculations
# NO correlation recalculation
#
# Heatmaps:
#   - All correlation cells retain their coloured tiles.
#   - Only cells with |median correlation| >= 0.5 display numerical labels.
#   - Labels show the signed correlation coefficient to two decimal places.
#
# Additional plots:
#   - System transit time: weak Spearman correlations only.
#   - Combined SoilR heatmap: all nine outcomes, retaining predictors
#     with at least one strong Spearman correlation.
#   - Combined SoilR heatmap excluding FTIR-derived predictors.
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

correlation_summary <- analysis$correlation_summary
correlation_draws <- analysis$correlation_draws
persistence_draws <- analysis$persistence_draws
geochem_summary <- analysis$geochemistry_summary

all_outcomes <- analysis$outcomes
valid_predictors <- analysis$predictors

methods <- c("spearman", "kendall", "pearson")

age_outcomes <- c(
  "fPOM_age",
  "oPOM_age",
  "MAOM_age",
  "system_transit_time"
)

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

soilR_outcome_labels <- c(
  kf = "fPOM decomposition rate (kf)",
  ki = "oPOM decomposition rate (ki)",
  ks = "MAOM decomposition rate (ks)",
  alpha_fi = "fPOM -> oPOM transfer (alpha_fi)",
  alpha_is = "oPOM -> MAOM transfer (alpha_is)",
  fPOM_age = "fPOM mean age",
  oPOM_age = "oPOM mean age",
  MAOM_age = "MAOM mean age",
  system_transit_time = "System mean transit time"
)

# -----------------------------------------------------------------------------
# 2. General heatmap function
#
# All tiles are retained.
# Numerical labels are displayed only when |correlation| >= 0.5.
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
      ),
      rho_label = ifelse(
        !is.na(correlation_median) &
          abs(correlation_median) >= 0.5,
        sprintf("%.2f", correlation_median),
        ""
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
    geom_text(
      aes(label = rho_label),
      colour = "black",
      size = 2.7,
      fontface = "bold"
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
        " correlation across sampled SoilR posterior draws. ",
        "Numerical labels shown for |coefficient| >= 0.5."
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
# 3. Generate general heatmaps
# -----------------------------------------------------------------------------

for (method_name in methods) {
  
  method_label <- stringr::str_to_title(method_name)
  
  # All rates, transfers, ages and transit time
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
        fill = "white",
        colour = NA
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.grid.minor = element_blank(),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
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
      )
    )
  
  transit_height <- max(
    5,
    0.35 * nrow(system_transit_plot_data) + 2
  )
  
  ggsave(
    file.path(
      out_dir,
      "system_transit_time_spearman_absrho_le_0.5.png"
    ),
    p_system_transit_weak,
    width = 10,
    height = transit_height,
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
    height = transit_height,
    bg = "white"
  )
  
  print(p_system_transit_weak)
  
} else {
  
  message(
    "No system transit time predictors with |median Spearman rho| <= 0.5."
  )
}

# -----------------------------------------------------------------------------
# 5. Combined SoilR heatmap: all nine outcomes
#
# Retain predictors with at least one strong Spearman correlation.
# Display ALL their correlation tiles, including weak correlations.
# Only strong correlations receive numerical labels.
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
    ),
    rho_label = ifelse(
      !is.na(correlation_median) &
        abs(correlation_median) >= 0.5,
      sprintf("%.2f", correlation_median),
      ""
    )
  )

if (nrow(soilR_heatmap_data) > 0) {
  
  p_soilR_combined <- ggplot(
    soilR_heatmap_data,
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
    geom_text(
      aes(label = rho_label),
      colour = "black",
      size = 2.8,
      fontface = "bold"
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
      title = "SoilR parameters and persistence vs geochemistry",
      subtitle = paste(
        "All nine SoilR outcomes; predictors retained if at least one",
        "outcome has |median Spearman rho| >= 0.5.",
        "All correlation tiles displayed; strong correlations labelled."
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
      panel.grid = element_blank(),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = 65,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        colour = "black"
      ),
      legend.text = element_text(
        colour = "black"
      ),
      legend.title = element_text(
        colour = "black"
      )
    )
  
  n_predictors <- length(
    unique(soilR_heatmap_data$predictor)
  )
  
  n_outcomes <- length(
    unique(soilR_heatmap_data$outcome)
  )
  
  width <- max(
    12,
    min(30, 0.35 * n_predictors + 5)
  )
  
  height <- max(
    9,
    min(18, 0.35 * n_outcomes + 4)
  )
  
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_all_soilR_strong_predictors.png"
    ),
    p_soilR_combined,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_all_soilR_strong_predictors.pdf"
    ),
    p_soilR_combined,
    width = width,
    height = height,
    limitsize = FALSE,
    bg = "white"
  )
  
  print(p_soilR_combined)
  
  message(
    "\nCombined SoilR heatmap: ",
    length(strong_soilR_predictors),
    " predictors and ",
    nrow(soilR_heatmap_data),
    " correlation cells."
  )
  
} else {
  
  message(
    "No strong predictors found for the combined SoilR heatmap."
  )
}

# -----------------------------------------------------------------------------
# 6. Combined SoilR heatmap excluding FTIR-derived predictors
#
# Retain predictors with at least one strong correlation across the
# nine SoilR outcomes, excluding FTIR-derived variables.
#
# All correlation tiles remain visible.
# Only cells with |median rho| >= 0.5 are labelled.
# -----------------------------------------------------------------------------

ftir_predictors <- correlation_summary %>%
  distinct(predictor) %>%
  filter(
    grepl(
      "FTIR|ftir|DRIFT|drift|aromatic|carbonate|carbonyl|carboxyl|aliphatic|spectral|peak",
      predictor
    )
  ) %>%
  pull(predictor)

message(
  "\nExcluded FTIR-derived predictors: ",
  length(ftir_predictors)
)

print(ftir_predictors)

strong_soilR_predictors_no_ftir <- correlation_summary %>%
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

soilR_heatmap_data_no_ftir <- correlation_summary %>%
  filter(
    method == "spearman",
    outcome %in% soilR_outcomes,
    predictor %in% strong_soilR_predictors_no_ftir
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
      levels = strong_soilR_predictors_no_ftir
    ),
    rho_label = ifelse(
      !is.na(correlation_median) &
        abs(correlation_median) >= 0.5,
      sprintf("%.2f", correlation_median),
      ""
    )
  )

if (nrow(soilR_heatmap_data_no_ftir) > 0) {
  
  p_soilR_no_ftir <- ggplot(
    soilR_heatmap_data_no_ftir,
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
    geom_text(
      aes(label = rho_label),
      colour = "black",
      size = 2.8,
      fontface = "bold"
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
      title = "SoilR parameters and persistence vs geochemistry",
      subtitle = paste(
        "FTIR-derived predictors excluded.",
        "All correlation tiles displayed; numerical labels shown",
        "only for |median Spearman rho| >= 0.5."
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
      panel.grid = element_blank(),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = 65,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        colour = "black"
      ),
      legend.text = element_text(
        colour = "black"
      ),
      legend.title = element_text(
        colour = "black"
      )
    )
  
  n_predictors <- length(
    unique(soilR_heatmap_data_no_ftir$predictor)
  )
  
  n_outcomes <- length(
    unique(soilR_heatmap_data_no_ftir$outcome)
  )
  
  width <- max(
    12,
    min(30, 0.35 * n_predictors + 5)
  )
  
  height <- max(
    9,
    min(18, 0.35 * n_outcomes + 4)
  )
  
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_all_soilR_strong_predictors_no_FTIR.png"
    ),
    p_soilR_no_ftir,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "heatmap_spearman_all_soilR_strong_predictors_no_FTIR.pdf"
    ),
    p_soilR_no_ftir,
    width = width,
    height = height,
    limitsize = FALSE,
    bg = "white"
  )
  
  print(p_soilR_no_ftir)
  
  message(
    "\nNon-FTIR combined SoilR heatmap: ",
    length(strong_soilR_predictors_no_ftir),
    " predictors and ",
    nrow(soilR_heatmap_data_no_ftir),
    " correlation cells."
  )
  
} else {
  
  message(
    "No strong non-FTIR predictors found for the combined SoilR heatmap."
  )
}

# -----------------------------------------------------------------------------
# 7. Completion
# -----------------------------------------------------------------------------

message(
  "\nPlot-only workflow completed successfully."
)

message(
  "All heatmaps retain weak and strong coloured correlation tiles."
)

message(
  "Numerical correlation labels are displayed only when |correlation| >= 0.5."
)