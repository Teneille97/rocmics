# ============================================================
# Atmospheric Delta 14C Plot (2023–2026)
# ============================================================

# Load required libraries
library(here)
library(SoilR)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Load and Combine Atmospheric 14C Data
# ------------------------------------------------------------

# Extract historical atmospheric data from SoilR (Hua et al. 2021, NH Zone 1)
Atm14C <- Hua2021$NHZone1[, 1:2]

# Load the forecast CSV file
fAtm14C <- read.csv(here("csv_files", "NHZ1forecast.csv"))

# Append forecast data and convert F14C to Delta 14C
Atm14C <- rbind(
  Atm14C,
  data.frame(
    Year = fAtm14C$time,
    mean.Delta14C = Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)
  )
)

# Filter for the 2023 to 2026 timeframe
Atm14C_subset <- Atm14C %>%
  filter(Year >= 2023 & Year <= 2026)

# ------------------------------------------------------------
# 2. Plot Atmospheric Delta 14C
# ------------------------------------------------------------

ggplot(Atm14C_subset, aes(x = Year, y = mean.Delta14C)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "navy", size = 2) +
  labs(
    title = "Atmospheric Delta 14C (2023–2026)",
    x = "Year",
    y = expression(paste("Atmospheric ", Delta^{14}, "C (‰)"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 3. Plot Long-Term Atmospheric Delta 14C (1950–2026)
# ------------------------------------------------------------

# Filter for the bomb spike era up to present (1950 to 2026)
Atm14C_longterm <- Atm14C %>%
  filter(Year >= 1950 & Year <= 2026)

# Long-term plot showing the bomb peak and decline
ggplot(Atm14C_longterm, aes(x = Year, y = mean.Delta14C)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  labs(
    title = "Atmospheric Delta 14C Bomb Peak and Decline (1950–2026)",
    subtitle = "Hua et al. (2021) NH Zone 1 data combined with post-2019 forecast",
    x = "Year",
    y = expression(paste("Atmospheric ", Delta^{14}, "C (‰)"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 4. Plot Long-Term Atmospheric Delta 14C (~500 Years to 2026)
# ------------------------------------------------------------

# Filter from approximately 500 years ago (1526) to 2026
# (Using min(Atm14C$Year) ensures it grabs the earliest available historical record)
start_year <- max(1526, min(Atm14C$Year, na.rm = TRUE))

Atm14C_500yr <- Atm14C %>%
  filter(Year >= start_year & Year <= 2026)

# Plot long-term historical trend including the pre-bomb baseline
ggplot(Atm14C_500yr, aes(x = Year, y = mean.Delta14C)) +
  geom_line(color = "darkgreen", linewidth = 0.6) +
  labs(
    title = "Atmospheric Delta 14C Over the Last ~500 Years (to 2026)",
    subtitle = "Showing pre-bomb baseline variability and the 20th-century bomb peak",
    x = "Year",
    y = expression(paste("Atmospheric ", Delta^{14}, "C (‰)"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(face = "bold")
  )