library(stringr)
library(reshape)
library(data.table)
library(nc)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(emmeans)
library(multcomp)
library(multcompView)
# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------

metadata <- read.csv(here("csv_files", "treatment_names.csv"))

metadata <- metadata %>%
  mutate(
    sample = sample
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

# import oxalate data

oxalate_data <- read.csv(here("csv_files", "oxalate_sesq_feb2026.csv"))
alldata_long <- cbind(oxalate_data, metadata)
alldata_long<-alldata_long%>%
  
  # --- split treatment ---
  separate(Treatment, into = c("Tmt", "App_rate"), sep = "_", fill = "right") %>%
  mutate(
    App_rate = replace_na(App_rate, "0"),
    Tmt = factor(Tmt),
    App_rate = factor(App_rate)
  ) 

# plot all data

library(tidyverse)

# Convert to long format
elements_long <- alldata_long %>%
  dplyr::select(Tmt, App_rate,
         Fe..mg.g.1.soil.,
         Al..mg.g.1.soil.,
         Mn..mg.g.1.soil.,
         Si..mg.g.1.soil.) %>%
  pivot_longer(
    cols = starts_with(c("Fe", "Al", "Mn", "Si")),
    names_to = "Element",
    values_to = "Concentration"
  ) %>%
  mutate(
    Element = recode(
      Element,
      "Fe..mg.g.1.soil." = "Fe",
      "Al..mg.g.1.soil." = "Al",
      "Mn..mg.g.1.soil." = "Mn",
      "Si..mg.g.1.soil." = "Si"
    ),
    App_rate = factor(App_rate,
                      levels = c("0", "2", "4", "8", "12", "20", "30", "50"))
  )

# Plot
ggplot(elements_long,
       aes(x = App_rate, y = Concentration)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_grid(Element ~ Tmt, scales = "free_y") +
  labs(
    x = "Application rate (t ha⁻¹)",
    y = expression(paste("Concentration (mg g"^-1, " soil)"))
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# app = 50

elements_long_50 <- alldata_long %>%
  dplyr::select(
    analysis_group, Tmt, App_rate,
    Fe..mg.g.1.soil.,
    Al..mg.g.1.soil.,
    Mn..mg.g.1.soil.,
    Si..mg.g.1.soil.
  ) %>%
  pivot_longer(
    cols = starts_with(c("Fe", "Al", "Mn", "Si")),
    names_to = "Element",
    values_to = "Concentration"
  ) %>%
  mutate(
    Element = recode(
      Element,
      "Fe..mg.g.1.soil." = "Fe",
      "Al..mg.g.1.soil." = "Al",
      "Mn..mg.g.1.soil." = "Mn",
      "Si..mg.g.1.soil." = "Si"
    )
  )

plot50 <- elements_long_50 %>%
  filter(
    analysis_group %in% c(
      "Control",
      "Bolsdorfer_50",
      "Eifelgold_50",
      "Huhnerberg_50",
      "Lime_2"
    )
  ) %>%
  mutate(
    analysis_group = factor(
      analysis_group,
      levels = c(
        "Control",
        "Bolsdorfer_50",
        "Eifelgold_50",
        "Huhnerberg_50",
        "Lime_2"
      )
    )
  )


# stats

letters_df <- plot50 %>%
  group_by(Element) %>%
  group_modify(~{
    
    mod <- aov(Concentration ~ analysis_group, data = .x)
    
    emm <- emmeans(mod, ~ analysis_group)
    
    cld_results <- cld(
      emm,
      Letters = letters,
      adjust = "tukey"
    )
    
    cld_results %>%
      as.data.frame() %>%
      dplyr::select(analysis_group, .group) %>%
      mutate(
        .group = gsub(" ", "", .group),
        y = max(.x$Concentration, na.rm = TRUE) * 1.08
      )
    
  }) %>%
  ungroup()

# add letters to plot
ggplot(plot50,
       aes(x = analysis_group,
           y = Concentration)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 2) +
  geom_text(data = letters_df,
            aes(y = y, label = .group),
            inherit.aes = FALSE,
            x = letters_df$analysis_group,
            size = 5) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression(paste("Concentration (mg g"^-1," soil)"))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# dose-response
dose_resp <- elements_long %>%
  filter(
    Tmt == "Control" |
      (Tmt == "Huhnerberg")
  ) %>%
  mutate(
    Treatment = ifelse(
      Tmt == "Control",
      "Control",
      paste0(App_rate, " t ha^-1")
    ),
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "2 t ha^-1",
        "4 t ha^-1",
        "8 t ha^-1",
        "12 t ha^-1",
        "20 t ha^-1",
        "30 t ha^-1",
        "50 t ha^-1"
      )
    )
  )

ggplot(dose_resp,
       aes(x = Treatment,
           y = Concentration)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 2) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression(paste("Concentration (mg g"^-1," soil)"))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(dose_resp,
       aes(Treatment, Concentration, group = 1)) +
  geom_point(size = 3) +
  facet_grid(Element ~ ., scales = "free_y") +
  labs(
    x = "",
    y = expression(paste("Concentration (mg g"^-1," soil)"))
  ) +
  theme_bw()

