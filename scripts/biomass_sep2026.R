
# =========================================================
# BIOMASS — SEPTEMBER 2026
# Conversion to kg m-2
# ANOVA + TUKEY HSD + MEAN ± SE PLOTS
# Plot area = 0.36 m2
# =========================================================


# ---------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(dplyr)
library(here)
library(multcompView)
library(ggthemr)

options(bitmapType = "cairo")

# ---------------------------------------------------------
# 2. Theme
# ---------------------------------------------------------

ggthemr("flat dark")

theme_update(
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8)
)


# ---------------------------------------------------------
# 3. Import raw biomass data
# ---------------------------------------------------------

biomass_Sep_2026 <- read.csv(
  here("csv_files", "biomass_Sep_2026.csv"),
  header = TRUE
)

# Ensure mass column is numeric
biomass_Sep_2026$Mass.Difference..g. <-
  as.numeric(biomass_Sep_2026$Mass.Difference..g.)


# ---------------------------------------------------------
# 4. Summarise biomass by plot and plant component
# ---------------------------------------------------------

# Conversion:
# g per plot / 360 = kg m-2
#
# 0.36 m2 * 1000 g kg-1 = 360 g m2 per kg

biomass_Sep_2026_clean <- biomass_Sep_2026 %>%
  
  group_by(Plot, Part) %>%
  
  summarise(
    Biomass_g = sum(Total.Mass..g., na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Convert long data to wide format
  pivot_wider(
    names_from = Part,
    values_from = Biomass_g,
    values_fill = 0
  ) %>%
  
  # Calculate total biomass in grams FIRST
  mutate(
    Total_Biomass_g = rowSums(
      across(-Plot),
      na.rm = TRUE
    )
  ) %>%
  
  # Convert all biomass components from g per plot
  # to kg m-2
  mutate(
    across(
      -Plot,
      ~ .x / 360
    )
  ) %>%
  
  # Rename total biomass column to reflect units
  dplyr::rename(
    Total_Biomass_kg_m2 = Total_Biomass_g
  ) %>%
  
  arrange(Plot)


# ---------------------------------------------------------
# 5. Import treatment metadata
# ---------------------------------------------------------

treatment_names <- read.csv(
  here("csv_files", "treatment_names.csv")
)


# ---------------------------------------------------------
# 6. Join treatment metadata
# ---------------------------------------------------------

biomass_Sep_2026_clean <- biomass_Sep_2026_clean %>%
  
  left_join(
    treatment_names %>%
      select(
        sample,
        Treatment,
        Basalt_type_exp,
        Dose_resp_exp
      ),
    by = c("Plot" = "sample")
  ) %>%
  
  relocate(
    Treatment,
    Basalt_type_exp,
    Dose_resp_exp,
    .after = Plot
  ) %>%
  
  arrange(Plot)


# ---------------------------------------------------------
# 7. Check joined data
# ---------------------------------------------------------

print(biomass_Sep_2026_clean)

# Check that all plots have a treatment
if (any(is.na(biomass_Sep_2026_clean$Treatment))) {
  warning("Some plots have no matching treatment metadata.")
}


# =========================================================
# ANOVA + TUKEY HSD + MEAN ± SE PLOTS
# =========================================================


# ---------------------------------------------------------
# 8. Prepare data for statistical analyses
# ---------------------------------------------------------

biomass_data <- biomass_Sep_2026_clean %>%
  
  mutate(
    Treatment = as.character(Treatment)
  )


# ---------------------------------------------------------
# 9. Define biomass properties
# ---------------------------------------------------------


# =========================================================
# BIOMASS — SEPTEMBER 2026
# TWO SEPARATE PLOTS
# 1. 50 t/ha + Lime comparison (ANOVA + Tukey HSD)
# 2. Huhnerberg dose-response (mean ± SE, no Tukey)
# Biomass units: kg m-2
# =========================================================


# ---------------------------------------------------------
# 1. Define biomass properties and facet labels
# ---------------------------------------------------------

biomass_properties <- c(
  "Corn",
  "Rem",
  "Stem",
  "TL",
  "Total_Biomass_kg_m2"
)

biomass_labels <- c(
  "Corn"                = "Corn",
  "Rem"                 = "Rem",
  "Stem"                = "Stem",
  "TL"                  = "TL",
  "Total_Biomass_kg_m2" = "Total~biomass~(kg~m^{-2})"
)


# =========================================================
# PLOT 1: 50 t/ha + LIME COMPARISON
# =========================================================


# ---------------------------------------------------------
# 2. Filter 50 t/ha treatments + lime
# ---------------------------------------------------------

biomass_50 <- biomass_data %>%
  filter(
    Treatment %in% c(
      "Control",
      "Bolsdorfer_50",
      "Eifelgold_50",
      "Huhnerberg_50",
      "Lime_2"
    )
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Bolsdorfer_50",
        "Eifelgold_50",
        "Huhnerberg_50",
        "Lime_2"
      )
    )
  )


# Check treatments included
print(table(biomass_50$Treatment, useNA = "ifany"))


# ---------------------------------------------------------
# 3. ANOVA + Tukey HSD for 50 t/ha comparison
# ---------------------------------------------------------

anova_50_results <- list()
tukey_50_results <- list()
letters_50_results <- list()

for (prop in biomass_properties) {
  
  formula <- as.formula(
    paste0("`", prop, "` ~ Treatment")
  )
  
  model <- aov(
    formula,
    data = biomass_50
  )
  
  # Store ANOVA
  anova_50_results[[prop]] <- summary(model)
  
  # Overall ANOVA p-value
  p_value <- summary(model)[[1]][["Pr(>F)"]][1]
  
  # Tukey HSD only if ANOVA is significant
  if (!is.na(p_value) && p_value < 0.05) {
    
    tukey <- TukeyHSD(model)
    
    tukey_50_results[[prop]] <- tukey
    
    letters <- multcompLetters4(
      model,
      tukey
    )
    
    letters_50_results[[prop]] <- data.frame(
      Treatment = names(letters$Treatment$Letters),
      Letters = as.character(letters$Treatment$Letters),
      Property = prop,
      stringsAsFactors = FALSE
    )
    
  } else {
    
    tukey_50_results[[prop]] <- NULL
    
    # No letters if ANOVA is not significant
    letters_50_results[[prop]] <- data.frame(
      Treatment = levels(biomass_50$Treatment),
      Letters = "",
      Property = prop,
      stringsAsFactors = FALSE
    )
  }
}


# ---------------------------------------------------------
# 4. Print ANOVA p-values — 50 t/ha comparison
# ---------------------------------------------------------

anova_50_pvalues <- data.frame(
  Property = biomass_properties,
  p_value = sapply(
    anova_50_results,
    function(x) x[[1]][["Pr(>F)"]][1]
  )
)

print(anova_50_pvalues)


# ---------------------------------------------------------
# 5. Calculate mean ± SE — 50 t/ha comparison
# ---------------------------------------------------------

biomass_50_summary <- biomass_50 %>%
  group_by(Treatment) %>%
  summarise(
    across(
      all_of(biomass_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )


biomass_50_long <- biomass_50_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  ) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property)
  )


# ---------------------------------------------------------
# 6. Prepare Tukey letters and positions
# ---------------------------------------------------------

letters_50_df <- bind_rows(letters_50_results) %>%
  mutate(
    Treatment = as.character(Treatment),
    Property = as.character(Property),
    Letters = as.character(Letters)
  )


letter_positions_50 <- letters_50_df %>%
  filter(
    !is.na(Letters),
    Letters != ""
  ) %>%
  left_join(
    biomass_50_long,
    by = c("Treatment", "Property")
  ) %>%
  group_by(Property) %>%
  mutate(
    y_max = max(mean + se, na.rm = TRUE),
    y_min = min(mean - se, na.rm = TRUE),
    y_range = y_max - y_min,
    y_range = if_else(
      !is.finite(y_range) | y_range == 0,
      1,
      y_range
    ),
    y_position = mean + se + 0.08 * y_range
  ) %>%
  ungroup() %>%
  filter(
    !is.na(mean),
    !is.na(se),
    is.finite(y_position)
  )


# ---------------------------------------------------------
# 7. Plot 1 — 50 t/ha + lime
#    Point symbols + SE + Tukey letters
# ---------------------------------------------------------

biomass_plot_50 <- ggplot(
  biomass_50_long,
  aes(
    x = Treatment,
    y = mean
  )
) +
  
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  geom_text(
    data = letter_positions_50,
    aes(
      x = Treatment,
      y = y_position,
      label = Letters
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    vjust = 0
  ) +
  
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      biomass_labels,
      default = label_parsed
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.25)
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Mean biomass " %+-% " SE (kg " * m^{-2} * ")"
    ),
    title = "Biomass — 50 t ha-1 Silicate Treatments + Lime"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(biomass_plot_50)


# =========================================================
# PLOT 2: HUHNERBERG DOSE-RESPONSE
# =========================================================


# ---------------------------------------------------------
# 8. Filter Control + Huhnerberg dose-response treatments
# ---------------------------------------------------------

biomass_dose <- biomass_data %>%
  filter(
    Treatment == "Control" |
      Treatment %in% c(
        "Huhnerberg_2",
        "Huhnerberg_4",
        "Huhnerberg_8",
        "Huhnerberg_12",
        "Huhnerberg_20",
        "Huhnerberg_30",
        "Huhnerberg_50"
      )
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Huhnerberg_2",
        "Huhnerberg_4",
        "Huhnerberg_8",
        "Huhnerberg_12",
        "Huhnerberg_20",
        "Huhnerberg_30",
        "Huhnerberg_50"
      ),
      labels = c(
        "Control",
        "2 t ha-1",
        "4 t ha-1",
        "8 t ha-1",
        "12 t ha-1",
        "20 t ha-1",
        "30 t ha-1",
        "50 t ha-1"
      )
    )
  )


# Check treatments included
print(table(biomass_dose$Treatment, useNA = "ifany"))


# ---------------------------------------------------------
# 9. Calculate mean ± SE — dose-response
# ---------------------------------------------------------

biomass_dose_summary <- biomass_dose %>%
  group_by(Treatment) %>%
  summarise(
    across(
      all_of(biomass_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )


biomass_dose_long <- biomass_dose_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "2 t ha-1",
        "4 t ha-1",
        "8 t ha-1",
        "12 t ha-1",
        "20 t ha-1",
        "30 t ha-1",
        "50 t ha-1"
      )
    ),
    Property = as.character(Property)
  )


# ---------------------------------------------------------
# 10. Plot 2 — Huhnerberg dose-response
#     Point symbols + SE, no Tukey letters
# ---------------------------------------------------------

biomass_plot_dose <- ggplot(
  biomass_dose_long,
  aes(
    x = Treatment,
    y = mean,
    group = 1
  )
) +
  
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      biomass_labels,
      default = label_parsed
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Mean biomass " %+-% " SE (kg " * m^{-2} * ")"
    ),
    title = "Biomass — Hühnerberg Dose-Response"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(biomass_plot_dose)

