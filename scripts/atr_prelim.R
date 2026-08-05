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
    sample_no = "\\d{1,2}",
    rep_no = "[a|b|c]{1,2}",
    nomatch.error = FALSE
  )
  
  meta_long <- cbind(long_table, meta)
  
  meta_long <- meta_long %>%
    mutate(
      wavenumber = round(wavenumber, 0),
      sample_rep = paste(sample_no, rep_no, sep = "_")
    ) %>%
    rename(
      wave = wavenumber,
      abs = absorbance,
      sample = filename
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
  
  # ---------------------------------------------------------
  # PEAK DETECTION (raw derivative space)
  # ---------------------------------------------------------
  
  spec_deriv <- prospectr::savitzkyGolay(
    X = spec_mat,
    p = 2,
    w = 11,
    m = 1
  )
  
  detect_peaks <- function(x, y, threshold = 0.0005){
    
    dy <- diff(y)
    peak_idx <- which(diff(sign(dy)) == -2) + 1
    peak_idx <- peak_idx[abs(y[peak_idx]) > threshold]
    
    data.frame(
      wave = x[peak_idx],
      intensity = y[peak_idx],
      index = peak_idx
    )
  }
  
  x <- waves
  
  peak_list <- lapply(seq_len(nrow(spec_deriv)), function(i){
    
    res <- detect_peaks(x, spec_deriv[i, ])
    
    if(is.null(res) || nrow(res) == 0) return(NULL)
    
    data.frame(
      sample_rep = ids[i],
      res
    )
  })
  
  peaks <- bind_rows(peak_list)
  
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

before <- process_spectra("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics/dpt_files/prelim_before/*.dpt")
after  <- process_spectra("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Documents/rocmics/dpt_files/prelim_after/*.dpt")

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

before_deriv2 <- compute_second_derivative(before_individual)
after_deriv2  <- compute_second_derivative(after_individual)

ggplot(before_deriv2, aes(wave, deriv2, colour = sample)) +
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  labs(title = "Second derivative spectra - BEFORE")

ggplot(after_deriv2, aes(wave, deriv2, colour = sample)) +
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  labs(title = "Second derivative spectra - AFTER")

# peak picking 
detect_peaks <- function(df,
                         value_col = "abs_corr",
                         smooth_window = 11,
                         threshold = 0.001,
                         use_second_derivative = FALSE){
  
  out <- list()
  
  samples <- unique(df$sample_rep)
  
  for(s in samples){
    
    tmp <- df %>%
      filter(sample_rep == s) %>%
      arrange(wave)
    
    x <- tmp$wave
    y <- tmp[[value_col]]
    
    # ----------------------------
    # smoothing (always applied)
    # ----------------------------
    y_smooth <- savitzkyGolay(
      matrix(y, nrow = 1),
      p = 3,
      w = smooth_window,
      m = 0
    )[1,]
    
    # ----------------------------
    # OPTION 1: normal peak picking
    # ----------------------------
    if(!use_second_derivative){
      
      dy <- diff(y_smooth)
      peak_idx <- which(diff(sign(dy)) == -2) + 1
      
      peak_idx <- peak_idx[
        y_smooth[peak_idx] > threshold
      ]
      
    } else {
      
      # ----------------------------
      # OPTION 2: second derivative minima
      # (better for overlapping bands)
      # ----------------------------
      
      y2 <- savitzkyGolay(
        matrix(y, nrow = 1),
        p = 3,
        w = smooth_window,
        m = 2
      )[1,]
      
      peak_idx <- which(diff(sign(diff(y2))) == -2) + 1
      
      # second derivative peaks are minima -> invert criterion
      peak_idx <- peak_idx[
        abs(y2[peak_idx]) > threshold
      ]
    }
    
    if(length(peak_idx) == 0) next
    
    out[[s]] <- data.frame(
      sample_rep = s,
      sample = unique(tmp$sample),
      wave = x[peak_idx],
      intensity = y_smooth[peak_idx]
    )
  }
  
  bind_rows(out)
}

before_peaks <- detect_peaks(before_individual,
                             value_col = "abs_corr",
                             use_second_derivative = FALSE)

after_peaks  <- detect_peaks(after_individual,
                             value_col = "abs_corr",
                             use_second_derivative = FALSE)

before_peaks_d2 <- detect_peaks(before_individual,
                                value_col = "abs_corr",
                                use_second_derivative = TRUE)

after_peaks_d2  <- detect_peaks(after_individual,
                                value_col = "abs_corr",
                                use_second_derivative = TRUE)

# checkpoint peak detection - individual scans

ggplot(before_individual,
       aes(wave, abs_corr, group = sample_rep)) +
  geom_line(alpha = 0.5) +
  geom_point(data = before_peaks,
             aes(wave, intensity),
             size = 2) +
  facet_wrap(~sample, scales = "free_y") +
  theme_bw() +
  coord_cartesian(xlim = c(800, 3000)) +
  labs(title = "Before extraction - individual spectra")


ggplot(after_individual,
       aes(wave, abs_corr, group = sample_rep)) +
  geom_line(alpha = 0.5) +
  geom_point(data = after_peaks,
             aes(wave, intensity),
             size = 2) +
  facet_wrap(~sample, scales = "free_y") +
  theme_bw() +
  coord_cartesian(xlim = c(800, 3000)) +
  labs(title = "After extraction - individual spectra")

# ---------------------------------------------------------
# Match before and after peaks
# ---------------------------------------------------------

match_peak_shifts <- function(before_peaks,
                              after_peaks,
                              tolerance = 30){
  
  results <- list()
  
  samples <- intersect(
    unique(before_peaks$sample),
    unique(after_peaks$sample)
  )
  
  for(s in samples){
    
    b <- before_peaks %>%
      filter(sample == s)
    
    a <- after_peaks %>%
      filter(sample == s)
    
    if(nrow(b) == 0 | nrow(a) == 0) next
    
    sample_matches <- lapply(seq_len(nrow(b)), function(i){
      
      d <- abs(a$wave - b$wave[i])
      
      if(min(d) > tolerance) return(NULL)
      
      j <- which.min(d)
      
      data.frame(
        sample = s,
        before_wave = b$wave[i],
        after_wave = a$wave[j],
        shift_cm1 = a$wave[j] - b$wave[i],
        before_intensity = b$intensity[i],
        after_intensity = a$intensity[j]
      )
    })
    
    results[[s]] <- bind_rows(sample_matches)
  }
  
  bind_rows(results)
}

peak_shifts_d2 <- match_peak_shifts(
  before_peaks_d2,
  after_peaks_d2,
  tolerance = 20
)

peak_shifts <- match_peak_shifts(
  before_peaks,
  after_peaks,
  tolerance = 30
)

interesting_shifts <- peak_shifts %>%
  filter(
    before_wave >= 1200,
    before_wave <= 1800
  ) %>%
  arrange(sample, before_wave)

interesting_shifts


#visualize shifts
ggplot(peak_shifts,
       aes(before_wave,
           shift_cm1,
           colour = sample)) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = "Before peak position (cm-1)",
    y = "Peak shift after extraction (cm-1)",
    title = "Peak position changes after cation removal"
  )


# plot overlaid sample spectra

samples <- unique(before_individual$sample)

for(s in samples){
  
  overlay <- bind_rows(
    before_individual %>%
      filter(sample == s) %>%
      mutate(dataset = "Before"),
    
    after_individual %>%
      filter(sample == s) %>%
      mutate(dataset = "After")
  )
  
  p <- ggplot(overlay,
              aes(wave,
                  abs_corr,
                  colour = dataset)) +
    geom_line(size = 1) +
    scale_x_reverse() +
    theme_bw() +
    coord_cartesian(xlim = c(3000, 800)) +
    labs(
      title = paste("Sample", s),
      x = "Wavenumber (cm-1)",
      y = "Corrected absorbance"
    )
  
  print(p)
  
}

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


# ---------------------------------------------------------
# Combine and calculate extraction effect
# ---------------------------------------------------------

band_compare <- before_areas %>%
  select(-dataset) %>%
  rename_with(~paste0(.x,"_before"), -sample) %>%
  left_join(
    after_areas %>%
      select(-dataset) %>%
      rename_with(~paste0(.x,"_after"), -sample),
    by = "sample"
  ) %>%
  mutate(
    COO_total_change =
      100*(COO_total_after - COO_total_before)/COO_total_before,
    
    Silicate_change =
      100*(Silicate_after - Silicate_before)/Silicate_before,
    
    poorly_crystalline_change =
      100*(poorly_crystalline_after - poorly_crystalline_before)/poorly_crystalline_before,
    
    Carbonate_bend_change =
      100*(Carbonate_bend_after - Carbonate_bend_before)/Carbonate_bend_before,
    
    carbonate_stretch_change =
      100*(carbonate_stretch_after - carbonate_stretch_before)/carbonate_stretch_before,
    
    aromatic_change =
      100*(aromatic_after - aromatic_before)/aromatic_before,
    
    carbonyl_change =
      100*(carbonyl_after - carbonyl_before)/carbonyl_before,
    
    Aliphatic_total_change =
      100*(Aliphatic_total_after - Aliphatic_total_before)/Aliphatic_total_before,
    
    polysaccharide_change =
      100*(polysaccharide_after - polysaccharide_before)/polysaccharide_before  
  )

# inspect
band_compare

library(broom)

# --------------------------------------------------
# Statistical testing of peak areas
# paired before vs after
# --------------------------------------------------

stat_test <- band_compare %>%
  select(
    sample,
    ends_with("_before"),
    ends_with("_after")
  ) %>%
  pivot_longer(
    cols = -sample,
    names_to = c("band", "dataset"),
    names_pattern = "(.*)_(before|after)$",
    values_to = "area"
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = area
  )


# paired t-test + Wilcoxon
results <- stat_test %>%
  group_by(band) %>%
  summarise(
    
    n = sum(!is.na(before) & !is.na(after)),
    
    mean_before = mean(before, na.rm = TRUE),
    mean_after = mean(after, na.rm = TRUE),
    
    percent_change =
      100 * (mean_after - mean_before) / mean_before,
    
    t_test_p = if(n >= 2) {
      t.test(after, before, paired = TRUE)$p.value
    } else {
      NA_real_
    },
    
    wilcox_p = if(n >= 2) {
      wilcox.test(after, before, paired = TRUE)$p.value
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  arrange(t_test_p)

# adjusted FDR 
results <- results %>%
  mutate(
    t_test_FDR = p.adjust(t_test_p, method = "BH"),
    wilcox_FDR = p.adjust(wilcox_p, method = "BH")
  )

#cohen's d
effect_sizes <- stat_test %>%
  group_by(band) %>%
  summarise(
    cohens_d =
      mean(after - before, na.rm = TRUE) /
      sd(after - before, na.rm = TRUE)
  )

results <- left_join(results, effect_sizes, by="band")

results #silicate has large change

# --------------------------------------------------
# Average technical scans within biological sample
# for difference spectra only
# --------------------------------------------------

before_avg_diff <- before_individual %>%
  group_by(sample, wave) %>%
  summarise(
    before = mean(abs_corr, na.rm = TRUE),
    .groups = "drop"
  )


after_avg_diff <- after_individual %>%
  group_by(sample, wave) %>%
  summarise(
    after = mean(abs_corr, na.rm = TRUE),
    .groups = "drop"
  )
# --------------------------------------------------
# Difference spectra
# after - before
# --------------------------------------------------

diff_spec <- before_avg_diff %>%
  inner_join(
    after_avg_diff,
    by = c("sample", "wave")
  ) %>%
  mutate(
    diff = after - before
  )

# --------------------------------------------------
# Plot all samples
# --------------------------------------------------

ggplot(diff_spec,
       aes(wave, diff, colour = factor(sample))) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(size = 0.8) +
  
  scale_x_reverse() +
  
  theme_bw() +
  
  labs(
    title = "Difference spectra (After − Before)",
    x = expression(Wavenumber~(cm^{-1})),
    y = "Difference absorbance",
    colour = "Sample"
  ) +
  
  coord_cartesian(
    xlim = c(1800, 800)
  )

ggplot(diff_spec,
       aes(wave, diff)) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(size = 0.8) +
  
  scale_x_reverse() +
  
  facet_wrap(~sample,
             ncol = 2,
             scales = "free_y") +
  
  theme_bw() +
  
  labs(
    title = "Difference spectra (After − Before)",
    x = expression(Wavenumber~(cm^{-1})),
    y = "Difference absorbance"
  ) +
  
  coord_cartesian(
    xlim = c(1800, 800)
  )

# pca
diff_wide <- diff_spec %>%
  select(sample, wave, diff) %>%
  pivot_wider(names_from = wave,
              values_from = diff)

diff_mat <- as.data.frame(diff_wide)
rownames(diff_mat) <- diff_mat$sample
diff_mat$sample <- NULL

diff_mat <- as.matrix(diff_mat)

# replace NA (important)
diff_mat[is.na(diff_mat)] <- 0

pca_diff <- prcomp(diff_mat, center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca_diff$x)
scores$sample <- rownames(scores)

ggplot(scores, aes(PC1, PC2, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5) +
  theme_bw() +
  labs(
    title = "PCA of difference spectra (After − Before)",
    x = "PC1",
    y = "PC2"
  )

loadings <- as.data.frame(pca_diff$rotation)
loadings$wave <- as.numeric(rownames(loadings))

ggplot(loadings, aes(wave, PC1)) +
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  labs(
    title = "PC1 loadings (spectral regions driving treatment effects)",
    x = "Wavenumber (cm-1)",
    y = "Loading"
  )
