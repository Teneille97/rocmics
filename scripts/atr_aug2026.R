library(stringr)
library(reshape)
library(data.table)
library(nc)
library(baseline)
library(ggplot2)
library(tidyr)
library(pracma)
library(dplyr)
library(prospectr)
library(broom)
library(here)
library(FSA)
library(car)
library(multcompView)
library(FSA)
# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------
metadata <- read.csv(here("csv_files", "treatment_names.csv"))

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
# ---------------------------------------------------------
# PROCESS SPECTRA
# ---------------------------------------------------------

process_spectra <- function(path){
  
  # ---- import ----
  long_table <- data.table(
    filename = Sys.glob(path)
  )[by = filename,
    j = fread(filename, col.names = c("wavenumber", "absorbance"))
  ]
  
  # ---- metadata ----
  meta <- nc::capture_first_vec(
    long_table$filename,
    sample_no = "soil\\d+",
    rep_no = "[abc]"
  )
  
  meta_long <- cbind(long_table, meta)
  
  meta_long <- meta_long %>%
    mutate(
      wavenumber = round(wavenumber, 0),
      sample = gsub("\\.", "", sample_no),
      sample_rep = paste(sample, rep_no, sep = "_")
    ) %>%
    rename(
      wave = wavenumber,
      abs = absorbance,
      filename = filename
    )
  
  # ---------------------------------------------------------
  # BUILD WIDE MATRIX (raw)
  # ---------------------------------------------------------
  spec_wide <- dcast(meta_long, sample_rep ~ wave, value.var = "abs")
  
  ids <- spec_wide$sample_rep
  spec_mat <- as.matrix(spec_wide[, -1, with = FALSE])
  waves <- as.numeric(colnames(spec_wide)[-1])
  
  # ---------------------------------------------------------
  # BASELINE CORRECTION
  # ---------------------------------------------------------
  bl <- baseline::baseline(spec_mat, method = "als")
  spec_corr <- getCorrected(bl)
  
  corr_mat <- as.data.frame(spec_corr)
  corr_mat$sample_rep <- ids
  
  corr_long <- data.table::melt(
    setDT(corr_mat),
    id.vars = "sample_rep",
    variable.name = "wave",
    value.name = "abs_corr"
  )
  
  corr_long[, wave := as.numeric(as.character(wave))]
  
  # ---------------------------------------------------------
  # DEFINE COMMON GRID HERE
  # ---------------------------------------------------------
  common_grid <- sort(unique(corr_long$wave))
  
  # force consistent grid ordering
  corr_long <- corr_long %>%
    filter(wave %in% common_grid)
  
  
  list(
    raw = meta_long,
    corrected = corr_long,
    peaks = peaks,
    grid = common_grid
  )
}

# ---------------------------------------------------------
# RUN
# ---------------------------------------------------------

before <- process_spectra("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics/dpt_files/aug2026_dpt_before/*.dpt")
after  <- process_spectra("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics/dpt_files/aug2026_dpt_after/*.dpt")

# ---------------------------------------------------------
# FILTER REGION
# ---------------------------------------------------------

before_filt <- before$corrected %>%
  filter(wave >= 800, wave <= 3000)

after_filt <- after$corrected %>%
  filter(wave >= 800, wave <= 3000)

# ---------------------------------------------------------
# NORMALIZATION
# ---------------------------------------------------------

normalize_spectra <- function(df){
  
  df %>%
    group_by(sample_rep) %>%
    mutate(
      abs_norm = abs_corr / sqrt(sum(abs_corr^2, na.rm = TRUE))
    ) %>%
    ungroup()
}

before_filt <- normalize_spectra(before_filt)
after_filt  <- normalize_spectra(after_filt)

# ---------------------------------------------------------
# KEEP INDIVIDUAL SCANS FOR INTEGRATION
# ---------------------------------------------------------

before_individual <- before_filt %>%
  mutate(sample = sub("_[a-z]+$", "", sample_rep))

after_individual <- after_filt %>%
  mutate(sample = sub("_[a-z]+$", "", sample_rep))

# ---------------------------------------------------------
# SECOND DERIVATIVE
# ---------------------------------------------------------

compute_second_derivative <- function(df){
  
  scans <- unique(df$sample_rep)
  grid <- sort(unique(df$wave))
  n <- length(grid)
  
  out <- list()
  
  for(scan in scans){
    
    tmp <- df %>%
      filter(sample_rep == scan) %>%
      arrange(wave)
    
    y <- tmp$abs_corr
    
    y2 <- savitzkyGolay(
      matrix(y, nrow = 1),
      p = 3,
      w = 11,
      m = 2
    )[1, ]
    
    
    if(length(y2) != n){
      
      pad <- rep(NA, n)
      
      offset <- floor((n - length(y2)) / 2)
      
      pad[(offset + 1):(offset + length(y2))] <- y2
      
      y2 <- pad
    }
    
    
    out[[scan]] <- data.frame(
      sample_rep = scan,
      sample = unique(tmp$sample),
      wave = grid,
      deriv2 = as.numeric(y2)
    )
  }
  
  bind_rows(out)
}


# plot overlaid sample spectra

samples <- unique(before_individual$sample)
overlay_all <- bind_rows(
  before_individual %>% mutate(dataset = "Before"),
  after_individual %>% mutate(dataset = "After")
) %>%
  mutate(
    sample = factor(
      sample,
      levels = paste0("soil", 1:21)
    )
  )

ggplot(
  overlay_all,
  aes(
    x = wave,
    y = abs_corr,
    colour = dataset,
    group = interaction(dataset, sample_rep)
  )
) +
  geom_line(alpha = 0.6, linewidth = 0.6) +
  scale_x_reverse() +
  coord_cartesian(xlim = c(3000, 800)) +
  facet_wrap(~ sample, ncol = 4) +
  theme_bw() +
  labs(
    x = "Wavenumber (cm⁻¹)",
    y = "Corrected absorbance",
    colour = "Dataset"
  ) +
  theme(
    strip.background = element_rect(colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "bottom"
  )


#### overlay by treatment group 
# Join treatment metadata to individual corrected spectra
overlay_treatment <- overlay_all %>%
  left_join(metadata_join, by = "sample") %>%
  filter(!is.na(analysis_group))

# Option 1: Individual scans faceted by dataset and colored by treatment
ggplot(
  overlay_treatment,
  aes(
    x = wave,
    y = abs_corr,
    color = analysis_group,
    group = interaction(dataset, sample_rep)
  )
) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  scale_x_reverse() +
  facet_wrap(~ dataset, ncol = 2) +
  theme_bw() +
  labs(
    x = expression("Wavenumber (cm"^-1*")"),
    y = "Corrected Absorbance",
    color = "Treatment Group",
    title = "Spectral Overlay by Treatment Group"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

# Option 2: Average spectrum per treatment group
overlay_treatment_mean <- overlay_treatment %>%
  mutate(dataset = factor(dataset, levels = c("Before", "After"))) %>%
  group_by(dataset, analysis_group, wave) %>%
  summarise(
    mean_abs = mean(abs_corr, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  overlay_treatment_mean,
  aes(x = wave, y = mean_abs, color = analysis_group)
) +
  geom_line(linewidth = 0.8) +
  scale_x_reverse() +
  facet_wrap(~ dataset, ncol = 2) +
  theme_bw() +
  labs(
    x = expression("Wavenumber (cm"^-1*")"),
    y = "Mean Baseline-Corrected Absorbance",
    color = "Treatment Group",
    title = "Average Spectra by Treatment Group"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )


# ---------------------------------------------------------
# Function to integrate spectral regions per scan
# ---------------------------------------------------------

integrate_bands <- function(df){
  
  bands <- list(
    COO_sym  = c(1380,1425),
    COO_asym = c(1580,1650),
    Silicate = c(1000,1120),
    poorly_crystalline = c(950,1050),
    Carbonate_bend = c(870,880),
    carbonate_stretch = c(1415,1480),
    CH_assym = c(2898,2976),
    CH_sym = c(2839,2870),
    aromatic = c(1500,1550),
    carbonyl = c(1570,1710),
    polysaccharide = c(1030,1080)
  )
  
  
  results <- lapply(unique(df$sample_rep), function(scan){
    
    sub <- df %>%
      filter(sample_rep == scan)
    
    out <- data.frame(
      sample_rep = scan,
      sample = unique(sub$sample)
    )
    
    for(b in names(bands)){
      
      rng <- bands[[b]]
      
      tmp <- sub %>%
        filter(wave >= rng[1],
               wave <= rng[2]) %>%
        arrange(wave)
      
      # avoid errors if region missing
      if(nrow(tmp) > 1){
        area <- trapz(tmp$wave, tmp$abs_norm)
      } else {
        area <- NA
      }
      
      out[[b]] <- area
    }
    
    out
  })
  
  bind_rows(results)
}

# ---------------------------------------------------------
# Calculate integrated areas per rep
# ---------------------------------------------------------
before_scan_areas <- integrate_bands(before_individual)

after_scan_areas <- integrate_bands(after_individual)

# ---------------------------------------------------------
# Calculate combined functional groups per scan
# ---------------------------------------------------------

before_scan_areas <- before_scan_areas %>%
  mutate(
    Aliphatic_total = CH_assym + CH_sym,
    COO_total = COO_sym + COO_asym
  )

after_scan_areas <- after_scan_areas %>%
  mutate(
    Aliphatic_total = CH_assym + CH_sym,
    COO_total = COO_sym + COO_asym
  )

before_scan_areas <- before_scan_areas %>%
  mutate(
    total_area = Silicate + poorly_crystalline +
      carbonate_stretch +
      carbonyl +
      Aliphatic_total
  )

after_scan_areas <- after_scan_areas %>%
  mutate(
    total_area = Silicate + poorly_crystalline +
      carbonate_stretch +
      carbonyl +
      Aliphatic_total
  )
# ---------------------------------------------------------
# Average integrated peak areas across replicate scans
# ---------------------------------------------------------

before_areas <- before_scan_areas %>%
  group_by(sample) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Aliphatic_total = CH_assym + CH_sym,
    COO_total = COO_sym + COO_asym,
    dataset = "before"
  ) %>%
  select(
    -CH_assym,
    -CH_sym,
    -COO_sym,
    -COO_asym
  )


after_areas <- after_scan_areas %>%
  group_by(sample) %>% 
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Aliphatic_total = CH_assym + CH_sym,
    COO_total = COO_sym + COO_asym,
    dataset = "after"
  ) %>%
  select(
    -CH_assym,
    -CH_sym,
    -COO_sym,
    -COO_asym
  )

before_areas <- before_areas %>%
  left_join(metadata_join, by = "sample")

after_areas <- after_areas %>%
  left_join(metadata_join, by = "sample")

# ---------------------------------------------------------
# Combine and calculate extraction effect
# ---------------------------------------------------------

band_compare <- before_areas %>%
  
  # Rename spectral variables from BEFORE
  rename_with(
    ~ paste0(.x, "_before"),
    -c(sample, Treatment, analysis_group, Dose)
  ) %>%
  
  # Join AFTER spectral data
  left_join(
    after_areas %>%
      select(
        sample,
        Treatment,
        analysis_group,
        everything()
      ) %>%
      select(-Dose) %>%
      rename_with(
        ~ paste0(.x, "_after"),
        -c(sample, Treatment, analysis_group)
      ),
    by = c("sample", "Treatment", "analysis_group")
  ) %>%
  
  # Calculate extraction-induced changes
  mutate(
    COO_change      = COO_total_after - COO_total_before,
    carbonyl_change = carbonyl_after - carbonyl_before,
    aromatic_change = aromatic_after - aromatic_before,
    CH_change       = Aliphatic_total_after - Aliphatic_total_before
  )


# inspect
band_compare_long <- band_compare %>%
  # Filter out NA treatments if needed
  filter(!is.na(Treatment), Treatment != "NA") %>%
  pivot_longer(
    cols = ends_with("_change"),
    names_to = "functional_group",
    values_to = "change"
  ) %>%
  # Clean up column names (removes "_change" for cleaner facet titles)
  mutate(functional_group = gsub("_change$", "", functional_group))

band_compare_long <- band_compare %>%
  filter(!is.na(analysis_group)) %>%
  pivot_longer(
    cols = ends_with("_change"),
    names_to = "functional_group",
    values_to = "change"
  ) %>%
  mutate(functional_group = gsub("_change$", "", functional_group))

ggplot(band_compare_long, aes(x = analysis_group, y = change)) +
  #geom_point(size = 2.5, alpha = 0.6, colour = "blue") +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se, 
    geom = "errorbar", 
    width = 0.2, 
    color = "black"
  ) +
  # scales = "free_y" lets each functional group have its own y-axis scale
  facet_wrap(~ functional_group, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Rock amendment",
    y = "Change in band area after extraction",
    title = "SOM Functional Group Changes by Treatment"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )



##### treatment summary

treatment_summary <- band_compare %>%
  # Filter out NA rows and isolated dose-response treatments
  filter(!is.na(Treatment), Treatment != "NA") %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    
    # Dynamically calculate Mean and SE for all *_change columns
    across(
      ends_with("_change"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se   = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )

treatment_summary %>% 
  print(n = Inf, width = Inf)


# ---------------------------------------------------------
# HÜHNERBERG DOSE-RESPONSE ANALYSIS
# ---------------------------------------------------------

huhnerberg_summary <- band_compare %>%
  filter(
    str_detect(Treatment, "^Huhnerberg_")
  ) %>%
  group_by(Dose) %>%
  summarise(
    n = n(),
    
    across(
      ends_with("_change"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se   = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )

huhnerberg_long <- band_compare %>%
  filter(
    str_detect(Treatment, "^Huhnerberg_")
  ) %>%
  pivot_longer(
    cols = ends_with("_change"),
    names_to = "functional_group",
    values_to = "change"
  ) %>%
  mutate(
    functional_group = gsub("_change$", "", functional_group)
  )

ggplot(
  huhnerberg_long,
  aes(x = Dose, y = change)
) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 1.5
  ) +
  facet_wrap(
    ~ functional_group,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    x = expression("Hühnerberg application rate (t ha"^-1*")"),
    y = "Change in band area after extraction",
    title = "Spectral changes across the Hühnerberg dose-response gradient"
  )


huhnerberg_summary %>% 
  print(n = Inf, width = Inf)

########### stats

band_compare_filtered <- band_compare_long %>%
  filter(!is.na(analysis_group))

# Generate compact letter display per functional group
get_cld <- function(df) {
  dt <- FSA::dunnTest(change ~ factor(analysis_group), data = df, method = "bh")$res
  p_vals <- dt$P.adj
  names(p_vals) <- gsub(" ", "", dt$Comparison)
  
  cld <- multcompLetters(p_vals)$Letters
  data.frame(
    analysis_group = names(cld),
    letter = as.character(cld),
    stringsAsFactors = FALSE
  )
}

cld_df <- band_compare_filtered %>%
  group_by(functional_group) %>%
  group_modify(~ get_cld(.x)) %>%
  ungroup()

# Position letters slightly above mean + SE
pos_df <- band_compare_filtered %>%
  group_by(functional_group, analysis_group) %>%
  summarise(
    y_pos = mean(change, na.rm = TRUE) + (sd(change, na.rm = TRUE) / sqrt(n())) + 0.00005,
    .groups = "drop"
  )

cld_df <- left_join(cld_df, pos_df, by = c("functional_group", "analysis_group"))

# Plot with statistical letters
ggplot(band_compare_filtered, aes(x = analysis_group, y = change)) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  geom_text(
    data = cld_df,
    aes(x = analysis_group, y = y_pos, label = letter),
    vjust = -0.5,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~ functional_group, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Rock amendment",
    y = "Change in band area after extraction",
    title = "SOM Functional Group Changes by Treatment (Dunn test, p < 0.05)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )


#checks 

# Levene's test per functional group
levene_results <- band_compare_long %>%
  filter(!is.na(analysis_group)) %>%
  group_by(functional_group) %>%
  group_modify(~ {
    res <- car::leveneTest(change ~ factor(analysis_group), data = .x)
    data.frame(
      F_stat  = res$`F value`[1],
      df1     = res$Df[1],
      df2     = res$Df[2],
      p_value = res$`Pr(>F)`[1]
    )
  }) %>%
  ungroup()

print(levene_results)


# =========================================================
# SECOND VERSION: PERCENTAGE CHANGE IN FUNCTIONAL GROUPS
# =========================================================
#
# Percentage change is calculated separately for each
# functional group:
#
#   % change = ((After - Before) / Before) * 100
#
# Thus, each functional group is expressed relative to its
# own BEFORE integrated area.
# =========================================================


# ---------------------------------------------------------
# Calculate percentage change for each functional group
# ---------------------------------------------------------

band_compare_pct <- before_areas %>%
  
  # Rename spectral variables from BEFORE
  rename_with(
    ~ paste0(.x, "_before"),
    -c(sample, Treatment, analysis_group, Dose)
  ) %>%
  
  # Join AFTER spectral data
  left_join(
    after_areas %>%
      select(
        sample,
        Treatment,
        analysis_group,
        everything()
      ) %>%
      select(-Dose) %>%
      rename_with(
        ~ paste0(.x, "_after"),
        -c(sample, Treatment, analysis_group)
      ),
    by = c("sample", "Treatment", "analysis_group")
  ) %>%
  
  # Calculate percentage change relative to BEFORE
  mutate(
    
    COO_change_pct =
      (COO_total_after - COO_total_before) /
      COO_total_before * 100,
    
    carbonyl_change_pct =
      (carbonyl_after - carbonyl_before) /
      carbonyl_before * 100,
    
    aromatic_change_pct =
      (aromatic_after - aromatic_before) /
      aromatic_before * 100,
    
    CH_change_pct =
      (Aliphatic_total_after - Aliphatic_total_before) /
      Aliphatic_total_before * 100
  )


# ---------------------------------------------------------
# Convert to long format for plotting
# ---------------------------------------------------------

band_compare_pct_long <- band_compare_pct %>%
  
  filter(
    !is.na(Treatment),
    Treatment != "NA",
    !is.na(analysis_group)
  ) %>%
  
  pivot_longer(
    cols = ends_with("_change_pct"),
    names_to = "functional_group",
    values_to = "change_pct"
  ) %>%
  
  mutate(
    functional_group =
      gsub("_change_pct$", "", functional_group)
  )


# ---------------------------------------------------------
# Treatment summary
# ---------------------------------------------------------

treatment_summary_pct <- band_compare_pct %>%
  
  filter(
    !is.na(Treatment),
    Treatment != "NA"
  ) %>%
  
  group_by(Treatment) %>%
  
  summarise(
    
    n = n(),
    
    across(
      ends_with("_change_pct"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )


# Print treatment summary
treatment_summary_pct %>%
  print(n = Inf, width = Inf)


# ---------------------------------------------------------
# Plot percentage change by treatment
# ---------------------------------------------------------

ggplot(
  band_compare_pct_long,
  aes(
    x = analysis_group,
    y = change_pct
  )
) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    color = "black"
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  
  facet_wrap(
    ~ functional_group,
    scales = "free_y"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  theme_bw() +
  
  labs(
    x = "Rock amendment",
    y = "Change in band area (%)",
    title = "SOM Functional Group Changes by Treatment"
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_rect(
      fill = "gray90"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )


# ---------------------------------------------------------
# HÜHNERBERG DOSE-RESPONSE ANALYSIS
# ---------------------------------------------------------

huhnerberg_summary_pct <- band_compare_pct %>%
  
  filter(
    str_detect(Treatment, "^Huhnerberg_")
  ) %>%
  
  group_by(Dose) %>%
  
  summarise(
    
    n = n(),
    
    across(
      ends_with("_change_pct"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )


huhnerberg_long_pct <- band_compare_pct %>%
  
  filter(
    str_detect(Treatment, "^Huhnerberg_")
  ) %>%
  
  pivot_longer(
    cols = ends_with("_change_pct"),
    names_to = "functional_group",
    values_to = "change_pct"
  ) %>%
  
  mutate(
    functional_group =
      gsub("_change_pct$", "", functional_group)
  )


# ---------------------------------------------------------
# Hühnerberg dose-response plot
# ---------------------------------------------------------

ggplot(
  huhnerberg_long_pct,
  aes(
    x = Dose,
    y = change_pct
  )
) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 1.5
  ) +
  
  facet_wrap(
    ~ functional_group,
    scales = "free_y"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  theme_bw() +
  
  labs(
    x = expression(
      "Hühnerberg application rate (t ha"^-1*")"
    ),
    y = "Change in band area (%)",
    title = "Spectral Changes Across the Hühnerberg Dose-Response Gradient"
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_rect(
      fill = "gray90"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )


# =========================================================
# STATISTICS FOR PERCENTAGE CHANGE
# =========================================================

# ---------------------------------------------------------
# Filter data for statistical analysis
# ---------------------------------------------------------

band_compare_pct_filtered <- band_compare_pct_long %>%
  filter(
    !is.na(analysis_group),
    !is.na(change_pct)
  )


# ---------------------------------------------------------
# Dunn test + compact letter display
# ---------------------------------------------------------

get_cld_pct <- function(df) {
  
  dt <- FSA::dunnTest(
    change_pct ~ factor(analysis_group),
    data = df,
    method = "bh"
  )$res
  
  p_vals <- dt$P.adj
  
  names(p_vals) <- gsub(
    " ",
    "",
    dt$Comparison
  )
  
  cld <- multcompLetters(p_vals)$Letters
  
  data.frame(
    analysis_group = names(cld),
    letter = as.character(cld),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------
# Generate significance letters separately for each
# functional group
# ---------------------------------------------------------

cld_df_pct <- band_compare_pct_filtered %>%
  
  group_by(functional_group) %>%
  
  group_modify(
    ~ get_cld_pct(.x)
  ) %>%
  
  ungroup()


# ---------------------------------------------------------
# Calculate positions for significance letters
# ---------------------------------------------------------
#
# Letters are positioned above the mean + SE.
# The offset is scaled to the range of each functional
# group so that it works appropriately with percentage data.
# ---------------------------------------------------------

pos_df_pct <- band_compare_pct_filtered %>%
  
  group_by(
    functional_group,
    analysis_group
  ) %>%
  
  summarise(
    
    mean_change = mean(
      change_pct,
      na.rm = TRUE
    ),
    
    se = sd(
      change_pct,
      na.rm = TRUE
    ) / sqrt(
      sum(!is.na(change_pct))
    ),
    
    .groups = "drop"
  ) %>%
  
  group_by(functional_group) %>%
  
  mutate(
    
    range_change = diff(
      range(
        band_compare_pct_filtered$change_pct[
          band_compare_pct_filtered$functional_group ==
            first(functional_group)
        ],
        na.rm = TRUE
      )
    ),
    
    y_pos = mean_change +
      se +
      0.05 * range_change
  ) %>%
  
  ungroup()


# ---------------------------------------------------------
# Add positions to significance letters
# ---------------------------------------------------------

cld_df_pct <- cld_df_pct %>%
  
  left_join(
    pos_df_pct %>%
      select(
        functional_group,
        analysis_group,
        y_pos
      ),
    by = c(
      "functional_group",
      "analysis_group"
    )
  )


# ---------------------------------------------------------
# Final percentage-change plot with significance letters
# ---------------------------------------------------------

ggplot(
  band_compare_pct_filtered,
  aes(
    x = analysis_group,
    y = change_pct
  )
) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3,
    color = "black"
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_text(
    data = cld_df_pct,
    aes(
      x = analysis_group,
      y = y_pos,
      label = letter
    ),
    vjust = -0.5,
    size = 3.5,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ functional_group,
    scales = "free_y"
  ) +
  
  theme_bw() +
  
  labs(
    x = "Rock amendment",
    y = "Change in band area (%)",
    title = "SOM Functional Group Changes by Treatment",
    subtitle = "Percentage change relative to the pre-extraction band area"
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_rect(
      fill = "gray90"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )


# ---------------------------------------------------------
# Levene's test for percentage changes
# ---------------------------------------------------------

levene_results_pct <- band_compare_pct_long %>%
  
  filter(
    !is.na(analysis_group),
    !is.na(change_pct)
  ) %>%
  
  group_by(functional_group) %>%
  
  group_modify(~ {
    
    res <- car::leveneTest(
      change_pct ~ factor(analysis_group),
      data = .x
    )
    
    data.frame(
      F_stat = res$`F value`[1],
      df1 = res$Df[1],
      df2 = res$Df[2],
      p_value = res$`Pr(>F)`[1]
    )
  }) %>%
  
  ungroup()


# Print Levene's test results
print(levene_results_pct)


# =========================================================
# OPTIONAL: INSPECT THE PERCENTAGE-CHANGE DATA
# =========================================================

band_compare_pct_long %>%
  print(n = Inf, width = Inf)
