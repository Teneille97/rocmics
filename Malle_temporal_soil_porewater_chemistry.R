#load packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(zoo)
library(ggthemr)

#set theme
ggthemr('fresh')


#import files
Alldata_Rhizon<-read.csv("csv_files/Alldata_Rhizon.csv", header = T, sep = ";") #sep seems to change from my pc to the other
Alldata_Soil_phEC<-read.csv("csv_files/Alldata_Soil_phEC.csv", header = T)


#data cleaning
Alldata_Rhizon <- Alldata_Rhizon %>%
  mutate(across(c(Treatment, Dose_response_exp.), as.factor)) #convert chr to factor

Alldata_Rhizon$Sampled_on<-as.Date(Alldata_Rhizon$Sampled_on,format = "%Y/%m/%d") #date format

Alldata_Soil_phEC <- Alldata_Soil_phEC %>%
  mutate(across(c(Pot, Plot, Treatment, Dose_response_exp.), as.factor)) #convert chr to factor

Alldata_Soil_phEC$Date<-as.Date(Alldata_Soil_phEC$Date, format = "%Y/%m/%d")

Alldata_Soil_phEC <- Alldata_Soil_phEC %>%
  separate(
    col = Treatment,       # The column to split
    into = c("Tmt", "App_rate"), # The names of the new columns
    sep = "_"              # The separator character in your data
  ) %>%
  mutate(
    App_rate = replace_na(App_rate, "0") 
  )

Alldata_Soil_phEC$Tmt<-as.factor(Alldata_Soil_phEC$Tmt)

Alldata_Soil_phEC_summary <- Alldata_Soil_phEC %>%
  group_by(Date, Tmt, App_rate) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    se_pH   = sd(pH, na.rm = TRUE) / sqrt(n()),
    
    mean_EC = mean(EC, na.rm = TRUE),
    se_EC   = sd(EC, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  )

Alldata_Soil_phEC_summary$Date<-as.Date(Alldata_Soil_phEC_summary$Date, format = "%d/%m/%Y")

#plotting
soilpH_plot <- Alldata_Soil_phEC_summary %>%
  dplyr::filter(App_rate %in% c("0", "2", "50")) %>%
  ggplot(aes(x = Date, y = mean_pH, colour = Tmt, shape = App_rate)) +
  #geom_line(aes(group = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_pH - se_pH,
                    ymax = mean_pH + se_pH),
                width = 0.1)+
  scale_x_date(
    breaks = "3 months",  # Display a tick mark every month
    date_labels = "%b %Y"
  )+
  theme(
    # Specify the element for the x-axis text
    axis.text.x = element_text(
      angle = 45,           # Rotate labels by 45 degrees
      hjust = 1,            # align the right edge of the label with the tick mark
      vjust = 1             # move the label down slightly 
    )
  ) + labs(
    y = "Soil pH", # <- Change the y-axis label
    title = "Soil pH over time",
    shape = "Application rate (t/ha)",
    colour = "Treatment"
  )