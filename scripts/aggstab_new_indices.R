#load packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)
library(ggthemr)
#set theme
ggthemr('fresh')
library(rstatix)
library(emmeans)
library(multcomp)
library(multcompView)
library(here)
library(ARTool)
library(tidyr)

# Read treatment names once
treatment_names <- read.csv(here("csv_files", "treatment_names.csv"), header = TRUE)
tmt_names_vec <- rep(treatment_names$Treatment, each = 2)

#--------------------------------------------------
# Fractal dimension calculation
#--------------------------------------------------

calc_fractal_D <- function(f_gt2,
                           f_gt025,
                           f_gt0063,
                           f_lt0063){
  
  # cumulative mass finer than sieve size
  cum_mass <- c(
    f_lt0063,
    f_lt0063 + f_gt0063,
    f_lt0063 + f_gt0063 + f_gt025
  )
  
  R <- c(0.063, 0.25, 2)
  
  fit <- lm(
    log(cum_mass) ~ log(R/max(R))
  )
  
  slope <- coef(fit)[2]
  
  D <- 3 - slope
  
  return(as.numeric(D))
}

# data cleaning function
clean_aggstab_extract <- function(filename, extract_label) {
  
  # Import data
  df <- read.csv(here("csv_files", filename), header = TRUE)
  
  # Add treatment names
  df$treatment <- tmt_names_vec
  
  # Convert to factors
  df[c("extract", "rep", "treatment")] <- 
    lapply(df[c("extract", "rep", "treatment")], as.factor)
  
  # Select columns
  df <- df %>%
    dplyr::select(
      sample, extract, rep, recovery....,
      X...2mm.fraction,
      X...0.25mm.fraction,
      X...0.063mm.fraction,
      X...0.063mm.fraction.1,
      treatment
    )
  
  # Rename columns
  df <- df %>%
    rename(
      recovery = recovery....,
      `>2mm.fraction`   = X...2mm.fraction,
      `>0.25mm.fraction` = X...0.25mm.fraction,
      `>0.063mm.fraction` = X...0.063mm.fraction,
      `<0.063mm.fraction` = X...0.063mm.fraction.1
    )
  
  # Remove missing rows
  df <- df %>%
    filter(
      !is.na(`>2mm.fraction`),
      !is.na(`>0.25mm.fraction`),
      !is.na(`>0.063mm.fraction`),
      !is.na(`<0.063mm.fraction`)
    )
  
  #--------------------------------------------------
  # Aggregate metrics (fractions are PERCENTAGES → convert to proportions)
  #--------------------------------------------------
  
  df <- df %>%
    mutate(
      
      # Mean Weight Diameter
      MWD =
        (`>0.25mm.fraction` / 100) * 1.125 +
        (`>0.063mm.fraction` / 100) * 0.1565 +
        (`<0.063mm.fraction` / 100) * 0.0315,
      
      # Geometric Mean Diameter
      GMD = exp(
        (`>0.25mm.fraction` / 100)  * log(1.125) +
          (`>0.063mm.fraction` / 100) * log(0.1565) +
          (`<0.063mm.fraction` / 100) * log(0.0315)
      )
    )
  
  # Fractal dimension (unchanged — slope invariant to scaling)
  df$D <- mapply(
    calc_fractal_D,
    df$`>2mm.fraction`,
    df$`>0.25mm.fraction`,
    df$`>0.063mm.fraction`,
    df$`<0.063mm.fraction`
  )
  
  # Pivot to long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(
        `>2mm.fraction`,
        `>0.25mm.fraction`,
        `>0.063mm.fraction`,
        `<0.063mm.fraction`
      ),
      names_to = "fraction_type",
      values_to = "fraction"
    ) %>%
    mutate(
      fraction_type = case_when(
        fraction_type == ">2mm.fraction"    ~ ">2",
        fraction_type == ">0.25mm.fraction"  ~ ">0.25",
        fraction_type == ">0.063mm.fraction" ~ ">0.063",
        fraction_type == "<0.063mm.fraction" ~ "<0.063"
      ),
      fraction_type = factor(
        fraction_type,
        levels = c(">2", ">0.25", ">0.063", "<0.063")
      ),
      extract_type = factor(
        extract_label,
        levels = c("h2o", "oxalate", "h2o2", "naoh")
      )
    )
  
  return(as.data.frame(df_long))
}

#--------------------------------------------------
# Generic analysis function
# Compare treatments within extracts
#--------------------------------------------------
#--------------------------------------------------
# Generic analysis function
# Compare treatments within extracts (EMMEANS ONLY)
#--------------------------------------------------

analyze_metric_by_extract <- function(df,
                                      metric = "MWD",
                                      ylab = metric){
  
  results_list <- list()
  plot_list <- list()
  
  for(xtract in unique(df$extract_type)) {
    
    metric_df <- df %>%
      filter(extract_type == xtract)
    
    metric_df$response <- metric_df[[metric]]
    
    # ---- MODEL (no ANOVA decision tree anymore) ----
    fit <- lm(response ~ treatment, data = metric_df)
    
    emm <- emmeans(fit, ~ treatment)
    
    letters_df <- multcomp::cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    ) %>%
      as.data.frame() %>%
      dplyr::select(treatment, .group)
    
    # ---- SUMMARY ----
    summary_df <- metric_df %>%
      group_by(treatment) %>%
      summarise(
        mean_value = mean(response, na.rm = TRUE),
        se_value = sd(response, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      ) %>%
      left_join(letters_df, by = "treatment")
    
    # ---- PLOT ----
    p <- ggplot(
      summary_df,
      aes(treatment, mean_value, fill = treatment)
    ) +
      geom_col(width = 0.7) +
      geom_errorbar(
        aes(
          ymin = mean_value - se_value,
          ymax = mean_value + se_value
        ),
        width = 0.2
      ) +
      geom_text(
        aes(label = .group,
            y = mean_value + se_value),
        fontface = "bold"
      ) +
      scale_fill_viridis_d() +
      labs(
        y = ylab,
        x = NULL,
        title = paste("Extract:", xtract)
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    results_list[[xtract]] <- summary_df
    plot_list[[xtract]] <- p
  }
  
  list(results = results_list, plots = plot_list)
}

#--------------------------------------------------
# Compare extracts within treatments (EMMEANS ONLY)
#--------------------------------------------------

analyze_metric_by_treatment <- function(df,
                                        metric = "MWD") {
  
  results_list <- list()
  
  for(tmt in unique(df$treatment)) {
    
    metric_df <- df %>%
      filter(treatment == tmt)
    
    metric_df$response <- metric_df[[metric]]
    
    # ---- MODEL ----
    fit <- lm(response ~ extract_type, data = metric_df)
    
    emm <- emmeans(fit, ~ extract_type)
    
    letters_df <- multcomp::cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    ) %>%
      as.data.frame() %>%
      dplyr::select(extract_type, .group)
    
    # ---- SUMMARY ----
    summary_df <- metric_df %>%
      group_by(extract_type) %>%
      summarise(
        mean_value = mean(response, na.rm = TRUE),
        se_value = sd(response, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      ) %>%
      left_join(letters_df, by = "extract_type") %>%
      mutate(
        treatment = tmt,
        metric = metric
      )
    
    results_list[[tmt]] <- summary_df
  }
  
  results_list
}
aggstab_h2o <- clean_aggstab_extract("aggstab_2026_h2o.csv", "h2o")

aggstab_h2o <- aggstab_h2o %>% #remove faulty sample
  filter(!(sample == 4 &
             rep == "b" &
             treatment == "Bolsdorfer_50"))


aggstab_oxalate <- clean_aggstab_extract("aggstab_2026_oxalate.csv", "oxalate")
aggstab_h2o2    <- clean_aggstab_extract("aggstab_2026_h2o2.csv", "h2o2")
aggstab_naoh    <- clean_aggstab_extract("aggstab_2026_naoh.csv", "naoh")

aggstab_all <- bind_rows(
  aggstab_h2o,
  aggstab_oxalate,
  aggstab_h2o2,
  aggstab_naoh
)
aggstab_metrics <- aggstab_all %>%
  distinct(
    sample,
    rep,
    treatment,
    extract_type,
    MWD,
    GMD,
    D
  )
MWD_extract <- analyze_metric_by_extract(
  aggstab_metrics,
  metric = "MWD",
  ylab = "Mean Weight Diameter (mm)"
)

GMD_extract <- analyze_metric_by_extract(
  aggstab_metrics,
  metric = "GMD",
  ylab = "Geometric Mean Diameter (mm)"
)

D_extract <- analyze_metric_by_extract(
  aggstab_metrics,
  metric = "D",
  ylab = "Fractal Dimension (D)"
)
MWD_treatment <- analyze_metric_by_treatment(
  aggstab_metrics,
  metric = "MWD"
)

GMD_treatment <- analyze_metric_by_treatment(
  aggstab_metrics,
  metric = "GMD"
)

D_treatment <- analyze_metric_by_treatment(
  aggstab_metrics,
  metric = "D"
)


#--------------------------------------------------
# H2O-only comparison: MWD, GMD, D
#--------------------------------------------------

# Letters from treatment comparisons WITHIN H2O

MWD_letters <- MWD_extract$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "MWD")

GMD_letters <- GMD_extract$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "GMD")

D_letters <- D_extract$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "D")

letters_all <- bind_rows(
  MWD_letters,
  GMD_letters,
  D_letters
)

h2o_metrics <- aggstab_metrics %>%
  filter(extract_type == "h2o") %>%
  pivot_longer(
    cols = c(MWD, GMD, D),
    names_to = "metric",
    values_to = "value"
  )
h2o_summary <- h2o_metrics %>%
  group_by(treatment, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

h2o_summary <- h2o_summary %>%
  left_join(
    letters_all,
    by = c("treatment", "metric")
  )

ggplot(
  h2o_summary,
  aes(
    x = treatment,
    y = mean_value,
    fill = treatment
  )
) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = .group,
      y = mean_value + se_value + 0.02
    ),
    fontface = "bold",
    vjust = 0
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

### eliminate dose-response

keep_treatments <- c(
  "Control",
  "Eifelgold_50",
  "Bolsdorfer_50",
  "Huhnerberg_50"
)

aggstab_metrics_filtered <- aggstab_metrics %>%
  filter(treatment %in% keep_treatments)

aggstab_metrics_filtered %>%
  distinct(treatment) %>%
  arrange(treatment)

h2o_metrics_filt <- aggstab_metrics_filtered %>%
  filter(extract_type == "h2o") %>%
  pivot_longer(
    cols = c(MWD, GMD, D),
    names_to = "metric",
    values_to = "value"
  )

h2o_summary_filt <- h2o_metrics_filt %>%
  group_by(treatment, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  h2o_summary_filt,
  aes(
    x = treatment,
    y = mean_value,
    fill = treatment
  )
) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = .group,
      y = mean_value + se_value + 0.02
    ),
    fontface = "bold",
    vjust = 0
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )


MWD_extract_filt <- analyze_metric_by_extract(
  aggstab_metrics_filtered,
  metric = "MWD",
  ylab = "MWD"
)

GMD_extract_filt <- analyze_metric_by_extract(
  aggstab_metrics_filtered,
  metric = "GMD",
  ylab = "GMD"
)

D_extract_filt <- analyze_metric_by_extract(
  aggstab_metrics_filtered,
  metric = "D",
  ylab = "D"
)

MWD_letters_filt <- MWD_extract_filt$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "MWD")

GMD_letters_filt <- GMD_extract_filt$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "GMD")

D_letters_filt <- D_extract_filt$results$h2o %>%
  dplyr::select(treatment, .group) %>%
  mutate(metric = "D")

letters_all_filt <- bind_rows(
  MWD_letters_filt,
  GMD_letters_filt,
  D_letters_filt
)

h2o_summary_filt <- h2o_summary_filt %>%
  left_join(letters_all_filt, by = c("treatment", "metric"))

#--------------------------------------------------
# MWD faceted plot
#--------------------------------------------------

MWD_all_summary <- bind_rows(MWD_treatment)

MWD_all_summary$treatment <- factor(
  MWD_all_summary$treatment,
  levels = unique(MWD_all_summary$treatment)
)

ggplot(
  MWD_all_summary,
  aes(
    x = extract_type,
    y = mean_value,
    fill = extract_type
  )
) +
  geom_col(
    width = 0.7,
    colour = "grey20"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = .group,
      y = mean_value + se_value + 0.01
    ),
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ treatment,
    scales = "free_y"
  ) +
  scale_fill_viridis_d(
    option = "D",
    end = 0.9
  ) +
  labs(
    x = NULL,
    y = "Mean weight diameter (mm)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )


#--------------------------------------------------
# GMD faceted plot
#--------------------------------------------------

GMD_all_summary <- bind_rows(GMD_treatment)

GMD_all_summary$treatment <- factor(
  GMD_all_summary$treatment,
  levels = unique(GMD_all_summary$treatment)
)

ggplot(
  GMD_all_summary,
  aes(
    x = extract_type,
    y = mean_value,
    fill = extract_type
  )
) +
  geom_col(
    width = 0.7,
    colour = "grey20"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = .group,
      y = mean_value + se_value + 0.01
    ),
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ treatment,
    scales = "free_y"
  ) +
  scale_fill_viridis_d(
    option = "D",
    end = 0.9
  ) +
  labs(
    x = NULL,
    y = "Geometric mean diameter (mm)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

#--------------------------------------------------
# Fractal dimension faceted plot
#--------------------------------------------------

D_all_summary <- bind_rows(D_treatment)

D_all_summary$treatment <- factor(
  D_all_summary$treatment,
  levels = unique(D_all_summary$treatment)
)

ggplot(
  D_all_summary,
  aes(
    x = extract_type,
    y = mean_value,
    fill = extract_type
  )
) +
  geom_col(
    width = 0.7,
    colour = "grey20"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = .group,
      y = mean_value + se_value + 0.05
    ),
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ treatment,
    scales = "free_y"
  ) +
  scale_fill_viridis_d(
    option = "D",
    end = 0.9
  ) +
  labs(
    x = NULL,
    y = "Fractal dimension (D)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

#normalize all plots
normalize_to_control <- function(df, control_name) {
  
  control_vals <- df %>%
    filter(treatment == control_name) %>%
    dplyr::select(metric, mean_value) %>%
    rename(control_mean = mean_value)
  
  df %>%
    left_join(control_vals, by = "metric") %>%
    mutate(
      mean_norm = mean_value / control_mean,
      se_norm = se_value / control_mean
    )
}

h2o_summary_norm <- normalize_to_control(h2o_summary, control_name = "Control")

ggplot(h2o_summary_norm,
       aes(x = treatment, y = mean_norm, fill = treatment)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(ymin = mean_norm - se_norm,
        ymax = mean_norm + se_norm),
    width = 0.2
  ) +
  geom_text(
    aes(label = .group,
        y = mean_norm + se_norm + 0.02),
    fontface = "bold"
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(y = "Relative to control (Control = 1)", x = NULL)

MWD_all_summary <- bind_rows(MWD_treatment) %>%
  normalize_to_control(control_name = "Control")
MWD_all_summary <- MWD_all_summary %>%
  distinct(treatment, extract_type, .keep_all = TRUE)

ggplot(MWD_all_summary,
       aes(x = extract_type, y = mean_norm, fill = extract_type)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(ymin = mean_norm - se_norm,
        ymax = mean_norm + se_norm),
    width = 0.2
  ) +
  geom_text(
    aes(label = .group,
        y = mean_norm + se_norm + 0.01),
    fontface = "bold"
  ) +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(y = "MWD (relative to control)", x = NULL) +
  theme_minimal()

GMD_all_summary <- bind_rows(GMD_treatment) %>%
  normalize_to_control("Control")
GMD_all_summary <- GMD_all_summary %>%
  distinct(treatment, extract_type, .keep_all = TRUE)

ggplot(GMD_all_summary,
       aes(x = extract_type, y = mean_norm, fill = extract_type)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(ymin = mean_norm - se_norm,
        ymax = mean_norm + se_norm),
    width = 0.2
  ) +
  geom_text(
    aes(label = .group,
        y = mean_norm + se_norm + 0.01),
    fontface = "bold"
  ) +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(y = "GWD (relative to control)", x = NULL) +
  theme_minimal()

D_all_summary <- bind_rows(D_treatment) %>%
  normalize_to_control("Control")
D_all_summary <- D_all_summary %>%
  distinct(treatment, extract_type, .keep_all = TRUE)


ggplot(D_all_summary,
       aes(x = extract_type, y = mean_norm, fill = extract_type)) +
  geom_col(width = 0.7, colour = "grey20") +
  geom_errorbar(
    aes(ymin = mean_norm - se_norm,
        ymax = mean_norm + se_norm),
    width = 0.2
  ) +
  geom_text(
    aes(label = .group,
        y = mean_norm + se_norm + 0.01),
    fontface = "bold"
  ) +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(y = "D (relative to control)", x = NULL) +
  theme_minimal()