#load packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(tidyverse)
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
aggstab_data_h2o <- read.csv(here("csv_files", "aggstab_2026_h2o.csv"), header = TRUE)
treatment_names<- read.csv(here("csv_files", "treatment_names.csv"), header = TRUE)
tmt_names_vec <- rep(treatment_names$Treatment, each = 2) #create vector to add tmt names
aggstab_data_h2o$treatment<-tmt_names_vec #add tmt names to df


#data cleaning
aggstab_data_h2o[c("extract", "rep", "treatment")] <- lapply(aggstab_data[c("extract", "rep", "treatment")], as.factor)
aggstab_data_h2o <- aggstab_data_h2o %>%
  select(sample, extract, rep, recovery...., X...2mm.fraction, X...0.25mm.fraction, X...0.063mm.fraction, X...0.063mm.fraction.1, treatment)
aggstab_data_h2o <- aggstab_data_h2o %>% 
  rename(
    recovery = recovery....,
    `>2mm.fraction` = X...2mm.fraction,
    `>0.25mm.fraction` = X...0.25mm.fraction,
    `>0.063mm.fraction` = X...0.063mm.fraction,
    `<0.063mm.fraction` = X...0.063mm.fraction.1
  )

aggstab_data_h2o_long <- aggstab_data_h2o %>%
  pivot_longer(
    cols = c(`>2mm.fraction`,
             `>0.25mm.fraction`,
             `>0.063mm.fraction`,
             `<0.063mm.fraction`),
    names_to = "fraction_type",
    values_to = "fraction"
  ) %>%
  mutate(
    fraction_type = case_when(
      fraction_type == ">2mm.fraction" ~ ">2",
      fraction_type == ">0.25mm.fraction" ~ ">0.25",
      fraction_type == ">0.063mm.fraction" ~ ">0.063",
      fraction_type == "<0.063mm.fraction" ~ "<0.063"
    )
  )

aggstab_df_h2o_long<-as.data.frame(aggstab_data_h2o_long)