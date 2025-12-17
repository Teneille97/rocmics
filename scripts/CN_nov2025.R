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
CN_data <- read.csv( here("csv_files", "CN_nov2025.csv"), header = TRUE)

#remove NA
CN_data <- CN_data[1:63,-8]

# Helper vector for filtering treatments
valid_treatments <- tribble(
  ~Tmt,          ~App_rate,
  "Control",     "0",
  "Lime",        "02",
  "Bolsdorfer",  "50",
  "Eifelgold",   "50",
  "Huhnerberg",  "50"
)

# Data cleaning and transformation function 
clean_data <- function(df) {
  df %>%
    mutate(across(c(Treatment, Dose_response_exp.), as.factor)) %>%
    separate(Treatment, into = c("Tmt", "App_rate"), sep = "_") %>%
    mutate(
      App_rate = replace_na(App_rate, "0"),
      Tmt = as.factor(Tmt)
    )%>%
    rename(
      N = X.N,
      C = X.C
    )
}

CN_data <- clean_data(CN_data)
CN_data$CN<-CN_data$C/CN_data$N

# Summarise CN data
CN_data_summary <- CN_data %>%
  group_by(Tmt, App_rate) %>%
  summarise(across(
    c(N, C, CN),
    list(mean = ~ mean(.x, na.rm = TRUE),
         se   = ~ sd(.x, na.rm = TRUE) / sqrt(n())),
    .names = "{.fn}_{.col}"),
    .groups = "drop")


#rate to numeric
CN_data_summary <- CN_data_summary %>%
  mutate(
    App_rate_num = as.numeric(App_rate)
  )

ggplot(CN_data_summary,
       aes(x = Tmt, y = mean_C,
           colour = Tmt,
           alpha = App_rate_num)) +
  
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  
  geom_errorbar(
    aes(ymin = mean_C - se_C,
        ymax = mean_C + se_C),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  
  scale_alpha_continuous(
    range = c(0.4, 1),
    name = "Application rate"
  ) +
  
  labs(
    x = "Treatment",
    y = "Mean C (%)"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# subset for selected treatments and rates
CN_sel <- CN_data %>%
  filter(
    (Tmt == "Control" & App_rate == "0") |
      (Tmt == "Lime" & App_rate == "02") |
      (Tmt %in% c("Eifelgold", "Huhnerberg", "Bolsdorfer") & App_rate == "50")
  )

# ensure ordering
CN_sel$Tmt <- factor(
  CN_sel$Tmt,
  levels = c("Control", "Lime", "Eifelgold", "Huhnerberg", "Bolsdorfer")
)

# ANOVA
aov_C <- aov(C ~ Tmt, data = CN_sel)
summary(aov_C)

# Tukey HSD
tuk <- TukeyHSD(aov_C)

# compact letter display
letters_C <- multcompLetters4(aov_C, tuk)

# convert to dataframe
cld_df <- as.data.frame.list(letters_C$Tmt)
cld_df$Tmt <- rownames(cld_df)
colnames(cld_df)[1] <- "letters"

CN_sel_summary <- CN_sel %>%
  group_by(Tmt) %>%
  summarise(
    mean_C = mean(C, na.rm = TRUE),
    se_C   = sd(C, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(cld_df, by = "Tmt")

ggplot(CN_sel_summary,
       aes(x = Tmt, y = mean_C, colour = Tmt)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = mean_C - se_C,
        ymax = mean_C + se_C),
    width = 0.2
  ) +
  geom_text(
    aes(label = letters, y = mean_C + se_C + 0.15),
    size = 5,
    show.legend = FALSE
  ) +
  labs(
    x = "Treatment",
    y = "Mean C (%)",
    title = "Total C for control, lime, and basalt treatments (50)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

CN_dose <- CN_data %>%
  filter(Dose_response_exp. == "Y") %>%
  group_by(Tmt, App_rate) %>%
  summarise(
    mean_C = mean(C, na.rm = TRUE),
    se_C   = sd(C, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(App_rate_num = as.numeric(App_rate))

ggplot(CN_dose,
       aes(x = App_rate_num, y = mean_C,
           colour = Tmt, group = Tmt)) +
  
  geom_point(size = 3) +
  geom_line() +
  
  geom_errorbar(
    aes(ymin = mean_C - se_C,
        ymax = mean_C + se_C),
    width = 1
  ) +
  labs(
    x = "Application rate",
    y = "Mean C (%)",
    title = "Dose–response treatments only"
  ) +
  theme_classic()
