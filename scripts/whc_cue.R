library(stringr)
library(reshape)
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------
metadata <- read.csv(here("csv_files", "treatment_names.csv"))
whc <- read.csv(here("csv_files", "WHC_Malle_Feb2026.csv"))

metadata <- metadata %>%
  mutate(
    sample = paste0("soil", sample),
    
    # Extract Hühnerberg application rate from Treatment
    Dose = case_when(
      str_detect(Treatment, "^Huhnerberg_") ~ 
        as.numeric(str_extract(Treatment, "(?<=_)\\d+")),
      TRUE ~ NA_real_
    )
  )

metadata <- metadata %>%
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

metadata_analysis <- metadata %>%
  filter(!is.na(analysis_group))

metadata_join <- metadata %>%
  select(
    sample,
    Treatment,
    Basalt_type_exp,
    Dose_resp_exp,
    analysis_group,
    Dose
  )

library(dplyr)

# Add a numeric sample ID to metadata to match whc$Sample
metadata_join2 <- metadata_join %>%
  mutate(Sample = as.numeric(gsub("soil", "", sample)))

# Join metadata with WHC data
whc_joined <- metadata_join2 %>%
  left_join(whc, by = "Sample")

# Calculate the average WHC values for each Treatment
# If there is only one sample for a Treatment, its original value is retained
whc_treatment_mean <- whc_joined %>%
  group_by(Treatment) %>%
  summarise(
    Water.content.at.60..WHC = mean(
      Water.content.at.60..WHC,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

whc_treatment_mean

whc_joined_mean <- whc_joined %>%
  group_by(Treatment) %>%
  mutate(
    Water.content.at.60..WHC_mean = mean(
      Water.content.at.60..WHC,
      na.rm = TRUE
    ),
    n_Treatment = n()
  ) %>%
  ungroup()

whc_joined_mean_df<-as.data.frame(whc_joined_mean)


# =========================================================
# WHC — FEBRUARY 2026
# TWO PLOTS:
# 1. 50 t/ha treatments + Lime
# 2. Huhnerberg dose-response
# Mean ± SE, point symbols + error bars
# =========================================================


# ---------------------------------------------------------
# 1. Prepare WHC data
# ---------------------------------------------------------

# Ensure WHC is numeric and treatment is character
whc_plot_data <- whc_joined %>%
  mutate(
    Treatment = as.character(Treatment),
    Water.content.at.60..WHC = as.numeric(
      as.character(Water.content.at.60..WHC)
    )
  )


# ---------------------------------------------------------
# 2. Helper function: calculate mean ± SE
# ---------------------------------------------------------

calculate_whc_summary <- function(data) {
  
  data %>%
    group_by(Treatment) %>%
    summarise(
      mean_WHC = mean(
        Water.content.at.60..WHC,
        na.rm = TRUE
      ),
      
      n_WHC = sum(
        !is.na(Water.content.at.60..WHC)
      ),
      
      se_WHC = ifelse(
        n_WHC > 1,
        sd(Water.content.at.60..WHC, na.rm = TRUE) /
          sqrt(n_WHC),
        NA_real_
      ),
      
      .groups = "drop"
    )
}


# =========================================================
# PLOT 1: 50 t/ha TREATMENTS + LIME
# =========================================================


# ---------------------------------------------------------
# 3. Filter comparison treatments
# ---------------------------------------------------------

whc_50 <- whc_plot_data %>%
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


# Check sample sizes
print(
  whc_50 %>%
    group_by(Treatment) %>%
    summarise(
      n_samples = sum(!is.na(Water.content.at.60..WHC)),
      .groups = "drop"
    )
)


# ---------------------------------------------------------
# 4. Calculate mean ± SE
# ---------------------------------------------------------

whc_50_summary <- calculate_whc_summary(whc_50) %>%
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


print(whc_50_summary)


# ---------------------------------------------------------
# 5. Plot 1 — WHC treatment comparison
# ---------------------------------------------------------

whc_plot_50 <- ggplot(
  whc_50_summary,
  aes(
    x = Treatment,
    y = mean_WHC
  )
) +
  
  # Mean point symbols
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  # Standard error bars
  geom_errorbar(
    aes(
      ymin = mean_WHC - se_WHC,
      ymax = mean_WHC + se_WHC
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  
  labs(
    x = "",
    y = "Water content at 60% WHC",
    title = "WHC — 50 t ha-1 Treatments + Lime"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(whc_plot_50)


# =========================================================
# PLOT 2: HUHNERBERG DOSE-RESPONSE
# =========================================================


# ---------------------------------------------------------
# 6. Filter Control + Huhnerberg dose-response treatments
# ---------------------------------------------------------

whc_dose <- whc_plot_data %>%
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


# Check sample sizes
print(
  whc_dose %>%
    group_by(Treatment) %>%
    summarise(
      n_samples = sum(!is.na(Water.content.at.60..WHC)),
      .groups = "drop"
    )
)


# ---------------------------------------------------------
# 7. Calculate mean ± SE
# ---------------------------------------------------------

whc_dose_summary <- calculate_whc_summary(whc_dose) %>%
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
    )
  )


print(whc_dose_summary)


# ---------------------------------------------------------
# 8. Plot 2 — Huhnerberg dose-response
# ---------------------------------------------------------

whc_plot_dose <- ggplot(
  whc_dose_summary,
  aes(
    x = Treatment,
    y = mean_WHC,
    group = 1
  )
) +
  
  # Mean point symbols
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  
  # Standard error bars
  geom_errorbar(
    aes(
      ymin = mean_WHC - se_WHC,
      ymax = mean_WHC + se_WHC
    ),
    width = 0.2,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  
  labs(
    x = "",
    y = "Water content at 60% WHC",
    title = "WHC — Hühnerberg Dose-Response"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )


print(whc_plot_dose)


# ---------------------------------------------------------
# 9. Optional: Save plots
# ---------------------------------------------------------

# ggsave(
#   here("WHC_50_t_ha_lime_Feb2026.png"),
#   whc_plot_50,
#   width = 9,
#   height = 6,
#   dpi = 300
# )

# ggsave(
#   here("WHC_Huhnerberg_dose_response_Feb2026.png"),
#   whc_plot_dose,
#   width = 10,
#   height = 6,
#   dpi = 300
# )