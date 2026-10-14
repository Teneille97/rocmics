# SOM decomposition models for Malle soil

# ============================================================
# LOAD LIBRARIES
# ============================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(emmeans)
library(stringr)


# ============================================================
# DATA CLEANING
# ============================================================

radiocarbon_samples <- read.csv(
  here("outputs", "radiocarbon_samples.csv")
)

aggstab_h2o <- read.csv(
  here("outputs", "aggstab_h2o.csv")
)

SOM_fraction_data <- radiocarbon_samples[
  radiocarbon_samples$Sampling_year == "2026",
]


# ============================================================
# LOAD PIPETTE DATA
# ============================================================

pipette_rawdata <- read.csv(
  here("csv_files", "Pipette_rawdata.csv")
)


# ============================================================
# LOAD METADATA
# ============================================================

metadata <- read.csv(
  here("csv_files", "treatment_names.csv")
) %>%
  mutate(
    analysis_group = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Eifelgold_50" ~ "Eifelgold_50",
      Treatment == "Bolsdorfer_50" ~ "Bolsdorfer_50",
      Treatment == "Lime_2" ~ "Lime_2",
      Treatment == "Huhnerberg_50" ~ "Huhnerberg_50",
      TRUE ~ NA_character_
    )
  )


# ============================================================
# JOIN METADATA TO PIPETTE DATA
# ============================================================

pipette_with_metadata <- merge(
  pipette_rawdata,
  metadata,
  by.x = "Sample",
  by.y = "sample",
  all.x = TRUE
)

# Remove clay outlier: Sample 16, Rep B, Energy III
pipette_with_metadata <- pipette_with_metadata %>%
  filter(
    !(Sample == 16 & Rep == "b" & Energy == "iii")
  ) %>%
  filter(
    !(X.clay_corrected ==- 0)
  )


# ============================================================
# LOAD SOM FRACTION DATA
# ============================================================

SOM_wide <- SOM_fraction_data %>%
  mutate(
    sample_id = as.integer(
      str_extract(sample, "^[0-9]+")
    )
  ) %>%
  pivot_wider(
    id_cols = sample_id,
    names_from = SOM_fraction,
    values_from = fraction_prop
  )


# ============================================================
# JOIN SOM FRACTIONS
# ============================================================

pipette_with_metadata <- pipette_with_metadata %>%
  left_join(
    SOM_wide,
    by = c("Sample" = "sample_id")
  )


# ============================================================
# PREPARE AGGREGATE STABILITY DATA
# ============================================================

aggstab_fractions <- aggstab_h2o %>%
  filter(
    fraction_type %in% c(
      ">0.25",
      ">0.063",
      "<0.063"
    )
  ) %>%
  group_by(
    sample,
    fraction_type
  ) %>%
  summarise(
    fraction_mean = mean(
      fraction,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    fraction_type = case_when(
      fraction_type == ">0.25" ~ "fraction_more_0.25",
      fraction_type == ">0.063" ~ "fraction_more_0.063",
      fraction_type == "<0.063" ~ "fraction_less_0.063"
    )
  ) %>%
  pivot_wider(
    names_from = fraction_type,
    values_from = fraction_mean
  ) %>%
  rename(
    sample_id = sample
  )


# ============================================================
# JOIN AGGREGATE STABILITY DATA
# ============================================================

pipette_with_metadata <- pipette_with_metadata %>%
  left_join(
    aggstab_fractions,
    by = c("Sample" = "sample_id")
  )

pipette_with_metadata$Treatment<-as.factor(pipette_with_metadata$Treatment)

# ============================================================
# SUMMARISE
# ============================================================

summary_data <- pipette_with_metadata %>%
  group_by(
    Treatment,
    Energy
  ) %>%
  summarise(
    
    mean_X.clay_corrected = mean(
      X.clay_corrected,
      na.rm = TRUE
    ),
    
    sd_X.clay_corrected = sd(
      X.clay_corrected,
      na.rm = TRUE
    ),
    
    n = sum(
      !is.na(X.clay_corrected)
    ),
    
    Bulk_prop = mean(
      Bulk,
      na.rm = TRUE
    ),
    
    fPOM_prop = mean(
      fPOM,
      na.rm = TRUE
    ),
    
    oPOM_prop = mean(
      oPOM,
      na.rm = TRUE
    ),
    
    MAOM_prop = mean(
      MAOM,
      na.rm = TRUE
    ),
    
    fraction_more_0.25 = mean(
      fraction_more_0.25,
      na.rm = TRUE
    ),
    
    fraction_more_0.063 = mean(
      fraction_more_0.063,
      na.rm = TRUE
    ),
    
    fraction_less_0.063 = mean(
      fraction_less_0.063,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

summary_data$Treatment<-as.factor(summary_data$Treatment)

# ============================================================
# FILTER INTO TREATMENT AND DOSE-RESPONSE DATAFRAMES
# ============================================================

pipette_app50_df <- summary_data %>%
  filter(
    !str_detect(Treatment, "Huhnerberg") |
      Treatment == "Huhnerberg_50"
  )
pipette_app50_df$Treatment<-as.factor(pipette_app50_df$Treatment)

dose_response_df <- summary_data %>%
  filter(
    Treatment == "Control" |
      str_detect(Treatment, "Huhnerberg")
  )
dose_response_df$Treatment<-as.factor(dose_response_df$Treatment)


# ============================================================
# PLOT CLAY
# ============================================================

ggplot(
  pipette_app50_df,
  aes(
    x = Energy,
    y = mean_X.clay_corrected,
    fill = Energy
  )
) +
  geom_col(
    color = "black",
    alpha = 0.8,
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(
        0,
        mean_X.clay_corrected - sd_X.clay_corrected
      ),
      ymax = mean_X.clay_corrected + sd_X.clay_corrected
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(~ Treatment) +
  labs(
    x = "Energy Level",
    y = "Mean % Clay (Corrected)",
    title = "Mean % Clay ± SD by Energy Level Across Treatments"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


# ============================================================
# CLAY NORMALISED TO oPOM FRACTION
# ============================================================

pipette_app50_df <- pipette_app50_df %>%
  mutate(
    clay_to_oPOM_ratio =
      mean_X.clay_corrected / oPOM_prop,
    
    sd_clay_to_oPOM_ratio =
      sd_X.clay_corrected / oPOM_prop
  )


ggplot(
  pipette_app50_df,
  aes(
    x = Energy,
    y = clay_to_oPOM_ratio,
    fill = Energy
  )
) +
  geom_col(
    color = "black",
    alpha = 0.8,
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(
        0,
        clay_to_oPOM_ratio - sd_clay_to_oPOM_ratio
      ),
      ymax =
        clay_to_oPOM_ratio + sd_clay_to_oPOM_ratio
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(~ Treatment) +
  labs(
    x = "Energy Level",
    y = "Mean % Clay / oPOM",
    title = "Clay Normalised to oPOM Fraction"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


# corr heatmaps

make_correlation_heatmap <- function(data, title_text) {
  
  correlation_data <- data %>%
    dplyr::select(
      X.clay_corrected,     # Raw sample clay %
      oPOM,                 # Raw sample oPOM fraction
      MAOM,                 # Raw sample MAOM fraction
      fraction_more_0.25,
      fraction_more_0.063,
      fraction_less_0.063
    )
  
  cor_matrix <- cor(
    correlation_data,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  cor_long <- as.data.frame(cor_matrix) %>%
    mutate(variable1 = rownames(.)) %>%
    pivot_longer(
      cols = -variable1,
      names_to = "variable2",
      values_to = "correlation"
    )
  
  variable_order <- c(
    "X.clay_corrected",
    "oPOM",
    "MAOM",
    "fraction_more_0.25",
    "fraction_more_0.063",
    "fraction_less_0.063"
  )
  
  variable_labels <- c(
    X.clay_corrected = "Clay (%)",
    oPOM = "oPOM",
    MAOM = "MAOM",
    fraction_more_0.25 = ">0.25 mm",
    fraction_more_0.063 = ">0.063 mm",
    fraction_less_0.063 = "<0.063 mm"
  )
  
  cor_long <- cor_long %>%
    mutate(
      variable1 = factor(variable1, levels = variable_order),
      variable2 = factor(variable2, levels = variable_order)
    )
  
  ggplot(cor_long, aes(x = variable1, y = variable2, fill = correlation)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.5) +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Spearman\nρ") +
    scale_x_discrete(labels = variable_labels) +
    scale_y_discrete(labels = variable_labels) +
    coord_equal() +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

# correlation heatmaps per energy level

pipette_with_metadata_50<-pipette_with_metadata%>%
  filter(
    !str_detect(Treatment, "Huhnerberg") |
      Treatment == "Huhnerberg_50"
  )

correlation_heatmap_Energy_I <- make_correlation_heatmap(
  pipette_with_metadata_50 %>% filter(Energy == "i"),
  "Correlations at Energy I (App50 Subset - Replicates)"
)

# view
correlation_heatmap_Energy_I

correlation_heatmap_Energy_II <- make_correlation_heatmap(
  pipette_with_metadata_50 %>% filter(Energy == "ii"),
  "Correlations at Energy II (App50 Subset - Replicates)"
)

# view
correlation_heatmap_Energy_II

correlation_heatmap_Energy_III <- make_correlation_heatmap(
  pipette_with_metadata_50 %>% filter(Energy == "iii"),
  "Correlations at Energy III (App50 Subset - Replicates)"
)

# view
correlation_heatmap_Energy_III

# plot clay energy vs opom
plot_clay_vs_oPOM_Energy_III <- pipette_with_metadata_50 %>%
  filter(Energy == "iii") %>%
  ggplot(aes(x = X.clay_corrected, y = oPOM)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") + # overall trend across replicates
  labs(
    x = "Clay (%)",
    y = "oPOM (%)",
    title = "Clay vs oPOM — Energy III (All Replicates)"
  ) +
  theme_minimal()

plot_clay_vs_oPOM_Energy_III

# plot clay energy vs maom
plot_clay_vs_MAOM_Energy_III <- pipette_with_metadata_50 %>%
  filter(Energy == "iii") %>%
  ggplot(aes(x = X.clay_corrected, y = MAOM)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") + # overall trend across replicates
  labs(
    x = "Clay (%)",
    y = "MAOM (%)",
    title = "Clay vs MAOM — Energy III (All Replicates)"
  ) +
  theme_minimal()

plot_clay_vs_MAOM_Energy_III

