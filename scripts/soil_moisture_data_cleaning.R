process_plot <- function(folder_name) {
  
  folder_path <- file.path(base_dir, folder_name)
  
  file_path <- list.files(
    folder_path,
    pattern = "^data.*\\.csv$",
    full.names = TRUE
  )[1]
  
  df <- read_delim(
    file_path,
    delim = ";",
    col_names = FALSE,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Handle column mismatch
  if (ncol(df) == 10) {
    # Drop last column (usually empty or junk)
    df <- df[, 1:9]
  } else if (ncol(df) != 9) {
    stop(paste("Unexpected number of columns in", file_path, ":", ncol(df)))
  }
  
  # Assign names
  names(df) <- c(
    "index", "date_time", "time_zone",
    "T1", "T2", "T3",
    "moisture_count", "shake", "flag"
  )
  
  df <- df %>%
    mutate(
      moisture_count = as.numeric(moisture_count),
      soil_moisture_content = -0.0000000134 * moisture_count^2 +
        0.000249622 * moisture_count -
        0.157889
    ) %>%
    separate(date_time, into = c("date", "time"), sep = " ", remove = FALSE) %>%
    mutate(
      date = str_replace_all(date, "\\.", "/")
    )
  
  return(df)
}
base_dir <- "C:/Program Files (x86)/Lolly/data"
plot_folders <- sprintf("tms_%02d", 1:21)
plot_data_list <- map(plot_folders, process_plot)
names(plot_data_list) <- sprintf("plot_%02d", 1:21)

# OPTIONAL: assign each to global environment (if you want individual dataframes)
#list2env(plot_data_list, envir = .GlobalEnv)

# OPTIONAL: combine all into one dataframe with plot ID
combined_data <- bind_rows(plot_data_list, .id = "plot_id")

combined_data <- combined_data %>%
  mutate(datetime = as.POSIXct(paste(date, time), format = "%Y/%m/%d %H:%M")) %>%
  filter(datetime >= as.POSIXct("2023-01-01"))

# Read metadata
metadata <- read_csv("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics/csv_files/treatment_names.csv")
# Check structure
print(metadata)

# Add plot number to combined data
combined_data <- combined_data %>%
  mutate(
    plot_num = as.numeric(str_extract(plot_id, "\\d+"))
  )

# Join metadata
combined_data <- combined_data %>%
  left_join(metadata, by = c("plot_num" = "sample"))

#plot by treatment
plot_by_treatment<-ggplot(combined_data, aes(x = as.POSIXct(paste(date, time)),
                          y = soil_moisture_content,
                          group = plot_id)) +
  geom_line(alpha = 0.6, color = "steelblue") +
  facet_wrap(~ Treatment, scales = "free_y") +
  labs(
    x = "Date",
    y = "Soil moisture content"
  ) +
  theme_minimal()


# Define a realistic minimum threshold 
# Also filter out the extreme "shake" flags if your sensor records movement
combined_data_clean <- combined_data %>%
  filter(soil_moisture_content > 0.1) %>%  # Removes zeros and negatives
  filter(!is.na(soil_moisture_content))     # Removes NAs

#view summary
treatment_summary <- combined_data_clean %>%
  group_by(Treatment, datetime) %>%
  summarise(
    mean_moisture = mean(soil_moisture_content, na.rm = TRUE),
    sd_moisture = sd(soil_moisture_content, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(treatment_summary, aes(x = datetime, y = mean_moisture, color = Treatment)) +
  geom_line(size = 1) +
  labs(title = "Average Soil Moisture by Treatment",
       y = "Volumetric Water Content",
       x = "Time") +
  theme_minimal()

#view boxplots

# Calculate one mean value per plot to avoid temporal autocorrelation issues in simple tests
plot_means <- combined_data_clean %>%
  group_by(plot_id, Treatment) %>%
  summarise(avg_smc = mean(soil_moisture_content), .groups = "drop")

ggplot(plot_means, aes(x = Treatment, y = avg_smc, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) + # Shows the individual 21 plots
  labs(title = "Consistency Check: Mean Moisture per Plot")
