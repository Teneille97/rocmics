#load packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)
library(car)
library(lme4)
library(ggthemr)
#set theme
ggthemr('fresh')
library(rstatix)
library(emmeans)
library(multcompView)
library(here)

#import data
CN_total_data <- read.csv(here("csv_files", "total_CN_MAOM_bulk.csv"), header = TRUE)
C_inorganic_data <- read.csv(here("csv_files", "inorgCN_bulk_MAOM.csv"), header = TRUE) #negligble


CN_clean <- CN_total_data %>%
  select(Sample.name, SN2, SN3, Sample.year, Parameter, Result) %>%
  rename(
    treatment = SN2,
    dose_resp_exp = SN3
  ) %>%
  pivot_wider(
    names_from = Parameter,
    values_from = Result
  ) %>%
  rename(
    C_total = Ctotal,
    N_total = Ntotal
  ) %>%
  
  # --- FIXED parsing ---
  mutate(
    sample_num  = as.numeric(str_extract(Sample.name, "^[0-9]+")),
    sample_year = as.numeric(str_extract(Sample.name, "(?<=_)[0-9]{4}")),
    sample_type = str_extract(Sample.name, "(?<=_[0-9]{4}_).*")
  ) %>%
  
  # --- enforce your custom year order ---
  mutate(
    sample_year = factor(sample_year, levels = c(2023, 2026, 2025))
  ) %>%
  
  arrange(
    sample_year,
    sample_num
  ) %>%
  
  select(-sample_num, -sample_year, -sample_type)

CN_clean
write.csv(
  CN_clean,
  here("csv_files", "CN_2026_clean.csv"),
  row.names = FALSE
)


# also for 2025 just for submission of 14c samples

CN_nov2025 <- read.csv(here("csv_files", "CN_nov2025.csv"), header = TRUE)

library(dplyr)
library(tidyr)
library(stringr)

library(dplyr)

CN_2025_clean <- CN_nov2025 %>%
  
  # group by plot
  group_by(sample) %>%
  
  # average the 3 replicates
  summarise(
    C_total = mean(X.C, na.rm = TRUE),
    N_total = mean(X.N, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # rename for consistency with your other dataset
  rename(plot = sample) %>%
  
  # sort nicely
  arrange(plot)

write.csv(
  CN_2025_clean,
  here("csv_files", "CN_2025_clean.csv"),
  row.names = FALSE
)
