# load packages
library(ggplot2)
library(gganimate)
library(ragg)
options(bitmapType = "cairo")
library(gifski)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(car)
library(zoo)
library(ggthemr)
theme_update(
  legend.title = element_text(size = 10),
  legend.text  = element_text(size = 8)
)
library(gridExtra)
library(here)

# set theme
ggthemr('flat dark')

CEC_Feb_26 <- read.csv(
  here("csv_files", "CEC_Feb_26.csv"),
  header = TRUE
)[1:252, ]

treatment_names <- read.csv(
  here("csv_files", "treatment_names.csv")
)

# clean and add metadata
CEC_Feb_26 <- CEC_Feb_26 %>%
  left_join(
    treatment_names %>%
      rename(Sample = sample),
    by = "Sample"
  ) %>%
  na.omit()

# Summarise
CEC_Feb_26_summary <- CEC_Feb_26 %>%
  group_by(Treatment) %>%
  summarise(
    across(
      c(
        H...cmolg.kg.1.,
        CALCIUM..cmolg.kg.1.,
        KALIUM..cmolg.kg.1.,
        MAGNESIUM..cmolg.kg.1.,
        NATRIUM..cmolg.kg.1.,
        Aluminium..cmolg.kg.1.,
        Ijzer..cmolg.kg.1.,
        Mangaan..cmolg.kg.1.,
        CEC..cmolg.kg.1.,
        Basesaturation....
      ),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se   = ~ sd(.x, na.rm = TRUE) / sqrt(n())
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )

# Plot data
CEC_long <- CEC_Feb_26_summary %>%
  pivot_longer(
    cols = -Treatment,
    names_to = c(".value", "Property"),
    names_pattern = "^(mean|se)_(.*)$"
  )
# ---------------------------------------------------------
# Publication-quality facet labels
# ---------------------------------------------------------

property_labels <- c(
  "H...cmolg.kg.1."        = "H^'+'~(cmol[c]~kg^{-1})",
  "CALCIUM..cmolg.kg.1."   = "Calcium~(cmol[c]~kg^{-1})",
  "KALIUM..cmolg.kg.1."    = "Potassium~(cmol[c]~kg^{-1})",
  "MAGNESIUM..cmolg.kg.1." = "Magnesium~(cmol[c]~kg^{-1})",
  "NATRIUM..cmolg.kg.1."   = "Sodium~(cmol[c]~kg^{-1})",
  "Aluminium..cmolg.kg.1." = "Aluminium~(cmol[c]~kg^{-1})",
  "Ijzer..cmolg.kg.1."     = "Iron~(cmol[c]~kg^{-1})",
  "Mangaan..cmolg.kg.1."   = "Manganese~(cmol[c]~kg^{-1})",
  "CEC..cmolg.kg.1."       = "CEC~(cmol[c]~kg^{-1})",
  "Basesaturation...."      = "Base~saturation~plain('%')"
)

# Plot
ggplot(CEC_long, aes(x = Treatment, y = mean)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ Property,
    scales = "free_y",
    labeller = as_labeller(
      property_labels,
      default = label_parsed
    )
  ) +
  labs(
    x = "Treatment",
    y = "Mean ± SE"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )