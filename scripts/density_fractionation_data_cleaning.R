library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)
library(ggrepel)


# import
density_fractionation <- read.csv(here("csv_files", "density_fractionation_data_final.csv")) #replace with final csv version after replacing originals with redo's
moisture_content <- read.csv(here("csv_files", "moisture_content_density_fractionation.csv"))
rbf_masses <- read.csv(here("csv_files", "rbf_masses_density_fractionation.csv"))
treatment_names <- read.csv(here("csv_files", "treatment_names.csv"))

# clean + join + compute
density_fractionation <- density_fractionation %>%
  
  # --- clean lookup tables ---
  left_join(
    moisture_content %>%
      rename(moisture_content = Water.content..m.m.1.),
    by = c("Plot" = "Sample",
           "Sampling_year" = "Sample.year")
  ) %>%
  
  left_join(
    rbf_masses %>%
      rename(rbf = RBF., mass_rbf = mass..g.),
    by = c("rbf.." = "rbf")
  ) %>%
  
  left_join(
    treatment_names %>%
      rename(Plot = sample),
    by = "Plot"
  ) %>%
  
  # --- drop junk columns safely ---
  select(-any_of(c("mass.rbf", "moisture.content", "recovery....",
                   "som.fraction....", "X", "X.1", "X.2", "Sample.date"))) %>%
  
  # --- fix treatment (2023 = Control) ---
  mutate(
    Treatment = if_else(Sampling_year == 2023, "Control", Treatment)
  ) %>%
  
  # --- split treatment ---
  separate(Treatment, into = c("Tmt", "App_rate"), sep = "_", fill = "right") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Tmt = factor(Tmt),
    App_rate = factor(App_rate)
  ) %>%
  
  # --- numeric fixes + calculations ---
  mutate(
    mass.som.on.filter.corr = as.numeric(mass.som.on.filter.corr),
    mass.som.on.filter.corr = replace_na(mass.som.on.filter.corr, 0),
    
    mass.soil.corr = mass.soil * (1 - moisture_content),
    
    mass.som.fraction =
      mass.rbf...som.fraction.freeze.dried +
      mass.som.on.filter.corr -
      mass_rbf,
    
    SOM_fraction = factor(SOM_fraction)
  )

# calculate total mass som fractions per sample
density_fractionation <- density_fractionation %>%
  # 1. Sort by Sampling_year, then Plot
  arrange(Sampling_year, Plot) %>%
  
  # 2. Group by the unique identifier for a physical soil sample
  group_by(Sampling_year, Plot) %>%
  
  # 3. Calculate the sum of mass.som.fraction and assign it to the column
  mutate(total.mass.som.fractions = sum(mass.som.fraction, na.rm = TRUE)) %>%
  
  # 4. Ungroup 
  ungroup()


#calculate recovery
density_fractionation$recovery_percent<-100*density_fractionation$total.mass.som.fractions/density_fractionation$mass.soil.corr

#plot outliers

# 1. Ensure recovery_percent is calculated
density_fractionation$recovery_percent <- 100 * (density_fractionation$total.mass.som.fractions / density_fractionation$mass.soil.corr)

# 2. Filter to one row per sample and identify outliers
# We define outliers as points outside 1.5 * IQR
plot_data <- density_fractionation %>%
  distinct(Sampling_year, Plot, recovery_percent) %>%
  mutate(
    # Logical check for outliers
    is_outlier = recovery_percent < (quantile(recovery_percent, 0.25, na.rm=T) - 1.5 * IQR(recovery_percent, na.rm=T)) | 
      recovery_percent > (quantile(recovery_percent, 0.75, na.rm=T) + 1.5 * IQR(recovery_percent, na.rm=T)),
    # Create the label only for outliers
    outlier_label = ifelse(is_outlier, paste0("Plot ", Plot, " (", Sampling_year, ")"), NA)
  )

# 3. Generate the plot
ggplot(plot_data, aes(x = factor(Sampling_year), y = recovery_percent)) +
  geom_boxplot(outlier.shape = NA) + # Hide default outliers to avoid double-plotting
  geom_jitter(width = 0.1, alpha = 0.4, color = "blue") + # Show individual data points
  geom_text_repel(aes(label = outlier_label), 
                  na.rm = TRUE, 
                  box.padding = 0.5, 
                  segment.color = 'grey50') +
  labs(title = "SOM Fraction Recovery Percentage",
       subtitle = "Outliers labeled by Plot number and Sampling Year",
       x = "Sampling Year",
       y = "Mass Recovery (%)") +
  theme_minimal()

# Calculate the proportion of each fraction within its own sample
df_proportions <- density_fractionation %>%
  mutate(
    fraction_prop = (mass.som.fraction / total.mass.som.fractions) * 100
  )

# Quick check: Do the replicates look similar?
# We'll group by the factors that define a "replicate"
replicate_summary <- df_proportions %>%
  group_by(Tmt, App_rate, Sampling_year, SOM_fraction) %>%
  summarise(
    mean_prop = mean(fraction_prop, na.rm = TRUE),
    sd_prop = sd(fraction_prop, na.rm = TRUE),
    cv = (sd_prop / mean_prop) * 100, # Coefficient of Variation
    .groups = "drop"
  )

print(replicate_summary, n = 63)

#visualize
# 1. Prepare data (filter for replicates and ensure proportions are calculated)
plot_data_faceted <- density_fractionation %>%
  mutate(fraction_prop = 100 * (mass.som.fraction / total.mass.som.fractions)) %>%
  group_by(Sampling_year, Tmt, App_rate) %>%
  filter(n_distinct(Plot) > 1) %>%
  ungroup()

# 2. Create the plot
ggplot(plot_data_faceted, aes(x = factor(Tmt), y = fraction_prop, fill = Tmt)) +
  # Using boxplots or points here is better than bars for comparing replicates
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.7) +
  # Facet by Fraction (Rows) and Year (Columns)
  # 'free_y' is the magic ingredient here
  facet_grid(SOM_fraction ~ Sampling_year, scales = "free_y") +
  labs(
    title = "Comparison of SOM Fraction Proportions",
    subtitle = "Note: Y-axes differ by fraction to highlight POM variability",
    x = "Treatment",
    y = "Percentage of Total SOM (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # Hide legend since X-axis labels cover it
  )