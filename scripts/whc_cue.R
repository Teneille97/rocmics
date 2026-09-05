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

