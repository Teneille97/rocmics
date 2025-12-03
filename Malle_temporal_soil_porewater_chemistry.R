#load packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(zoo)
library(ggthemr)

#set theme
ggthemr('pale')


#import files
Alldata_Rhizon<-read.csv("csv_files/Alldata_Rhizon.csv", header = T, sep = ";")
Alldata_Soil_phEC<-read.csv("csv_files/Alldata_Soil_phEC.csv", header = T, sep = ";")


#data cleaning
Alldata_Rhizon <- Alldata_Rhizon %>%
  mutate(across(c(Treatment, Dose_response_exp.), as.factor)) #convert chr to factor

Alldata_Rhizon$Sampled_on<-as.Date(Alldata_Rhizon$Sampled_on,format = "%d/%m/%Y") #date format

Alldata_Soil_phEC <- Alldata_Soil_phEC %>%
  mutate(across(c(Pot, Plot, Treatment, Dose_response_exp.), as.factor)) #convert chr to factor

Alldata_Soil_phEC$Date<-as.Date(Alldata_Soil_phEC$Date, format = "%d/%m/%Y")

Alldata_Soil_phEC_summary <- Alldata_Soil_phEC %>%
  group_by(Date, Treatment, Dose_response_exp.) %>%
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
  dplyr::filter(Dose_response_exp. == "N") %>%
  ggplot(aes(x = Date, y = mean_pH, colour = Treatment)) +
  geom_line(aes(group = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_pH - se_pH,
                    ymax = mean_pH + se_pH),
                width = 0.1)+
  scale_x_date(
    breaks = "3 months"  # Display a tick mark every month
  )