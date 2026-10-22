# =============================================================================
# Integrated Malle soil correlations: MAOM/oPOM C, aggregate stability, FTIR
# peak-area changes after cation extraction, and DCB/oxalate-extractable Fe, Al,
# Mn and Si.
#
# Expected project layout:
#   csv_files/  : treatment_names.csv, total_CN_MAOM_bulk.csv,
#                 oxalate_sesq_feb2026.csv, DCB_Feb2026.csv, and aggregate CSVs
#   outputs/    : ftir_peak_changes_feb2026.csv
#
# Results are written to outputs/integrated_correlation_results/.
# Required packages: here, tidyverse, corrplot
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(corrplot)
})

# ---- Paths ------------------------------------------------------------------
csv_dir <- here::here("csv_files")
outputs_dir <- here::here("outputs")
out_dir <- file.path(outputs_dir, "integrated_correlation_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_csv_required <- function(filename) {
  path <- file.path(csv_dir, filename)
  if (!file.exists(path)) stop("Required input file not found: ", path)
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA", "#VALUE!"))
}
read_csv_optional <- function(filename) {
  path <- file.path(csv_dir, filename)
  if (!file.exists(path)) {
    message("Optional input not found; skipping: ", path)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA", "#VALUE!"))
}

# ---- Plot/treatment metadata -------------------------------------------------
metadata <- read_csv_required("treatment_names.csv") %>%
  transmute(
    sample = as.integer(sample),
    Treatment = as.character(Treatment),
    Basalt_type_exp = toupper(as.character(Basalt_type_exp)) == "Y",
    Dose_resp_exp = toupper(as.character(Dose_resp_exp)) == "Y"
  ) %>%
  distinct(sample, .keep_all = TRUE)

if (nrow(metadata) != 21 || !setequal(metadata$sample, 1:21)) {
  warning("Treatment metadata does not contain exactly plots 1-21.")
}

# ---- MAOM-C and oPOM-C (% C in fraction) ------------------------------------
cn_raw <- read_csv_required("total_CN_MAOM_bulk.csv")
cn_wide <- cn_raw %>%
  dplyr::rename(sample_name = `Sample name`, result = `Result (%)`) %>%
  mutate(
    sample_name = as.character(sample_name),
    result = readr::parse_number(as.character(result)),
    Parameter = str_trim(as.character(Parameter))
  ) %>%
  filter(Parameter == "Ctotal") %>%
  extract(
    sample_name,
    into = c("sample", "year", "fraction"),
    regex = "^(\\d+)_((?:20)?\\d{2})_(MAOM|oPOM|fPOM|Bulk)$",
    remove = FALSE,
    convert = TRUE
  ) %>%
  filter(!is.na(sample), year == 2026,
         fraction %in% c("MAOM", "oPOM", "fPOM", "Bulk")) %>%
  group_by(sample, fraction) %>%
  summarise(C_percent = mean(result, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = fraction, values_from = C_percent, names_prefix = "C_")

# These are C concentrations (%) within each fraction, not stocks per unit area.
if (nrow(cn_wide) == 0) {
  warning("No 2026 Ctotal fraction rows matched the expected sample-name pattern.")
}

# ---- Aggregate stability: MWD, GMD and fractal D -----------------------------
calc_D <- function(gt2, gt025, gt0063, lt0063) {
  vals <- c(gt2, gt025, gt0063, lt0063)
  if (any(!is.finite(vals)) || any(vals <= 0)) return(NA_real_)
  cum_finer <- c(lt0063, lt0063 + gt0063,
                 lt0063 + gt0063 + gt025)
  R <- c(0.063, 0.25, 2)
  fit <- try(lm(log(cum_finer) ~ log(R / max(R))), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  as.numeric(3 - coef(fit)[2])
}

read_aggstab <- function(filename, extract_label) {
  dat <- read_csv_optional(filename)
  if (is.null(dat)) return(NULL)
  
  required_cols <- c("sample", "% >0.25mm fraction",
                     "% >0.063mm fraction", "% <0.063mm fraction",
                     "% >2mm fraction")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) {
    warning(filename, " missing expected columns: ",
            paste(missing_cols, collapse = ", "), "; skipping.")
    return(NULL)
  }
  
  dat %>%
    mutate(
      sample = as.integer(sample),
      extract_method = extract_label,
      MWD = (`% >0.25mm fraction` / 100) * 1.125 +
        (`% >0.063mm fraction` / 100) * 0.1565 +
        (`% <0.063mm fraction` / 100) * 0.0315,
      GMD = exp((`% >0.25mm fraction` / 100) * log(1.125) +
                  (`% >0.063mm fraction` / 100) * log(0.1565) +
                  (`% <0.063mm fraction` / 100) * log(0.0315)),
      D = pmap_dbl(
        list(`% >2mm fraction`, `% >0.25mm fraction`,
             `% >0.063mm fraction`, `% <0.063mm fraction`), calc_D
      )
    ) %>%
    group_by(sample, extract_method) %>%
    summarise(
      across(c(MWD, GMD, D), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
      n_agg_reps = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = extract_method,
      values_from = c(MWD, GMD, D),
      names_glue = "{.value}_{extract_method}"
    )
}

agg_files <- c(
  h2o = "aggstab_2026_h2o.csv",
  oxalate = "aggstab_2026_oxalate.csv",
  h2o2 = "aggstab_2026_h2o2.csv",
  naoh = "aggstab_2026_naoh.csv"
)
agg_tables <- imap(agg_files, ~read_aggstab(.x, .y)) %>% compact()
agg_data <- if (length(agg_tables)) reduce(agg_tables, full_join, by = "sample") else
  tibble(sample = integer())

# ---- DCB and oxalate-extractable elements -----------------------------------
ox_raw <- read_csv_required("oxalate_sesq_feb2026.csv")
ox <- ox_raw %>%
  transmute(
    sample = as.integer(Sample),
    ox_Fe_mg_g = `Fe (mg g-1 soil)`,
    ox_Al_mg_g = `Al (mg g-1 soil)`,
    ox_Mn_mg_g = `Mn (mg g-1 soil)`,
    ox_Si_mg_g = `Si (mg g-1 soil)`
  )

dcb_raw <- read_csv_required("DCB_Feb2026.csv")
dcb <- dcb_raw %>%
  dplyr::rename(sample = Plot) %>%
  mutate(sample = as.integer(sample)) %>%
  transmute(
    sample,
    dcb_Fe_mg_g = `%Fe` * 10,
    dcb_Al_mg_g = `%Al` * 10,
    dcb_Mn_mg_g = `%Mn` * 10,
    dcb_Si_mg_g = `%Si` * 10
  )

mineral_data <- full_join(ox, dcb, by = "sample")
for (el in c("Fe", "Al", "Mn", "Si")) {
  ox_col <- paste0("ox_", el, "_mg_g")
  dcb_col <- paste0("dcb_", el, "_mg_g")
  mineral_data[[paste0(el, "_DCB_minus_ox_mg_g")]] <-
    mineral_data[[dcb_col]] - mineral_data[[ox_col]]
  mineral_data[[paste0(el, "_ox_over_DCB")]] <-
    mineral_data[[ox_col]] / mineral_data[[dcb_col]]
}


# ---- Tessier sequential extraction ------------------------------------------
# Expected files in csv_files/:
# tessier_feb2026_nh4oac.csv
# tessier_feb2026_acetic.csv
# tessier_feb2026_hydroxylamine.csv
# tessier_feb2026_h2o2.csv

read_tessier <- function(filename, extract_label) {
  dat <- read_csv_required(filename)
  
  # Accommodate common variations in imported CSV column names.
  names(dat) <- make.names(names(dat), unique = TRUE)
  nm <- names(dat)
  
  sample_col <- intersect(c("sample", "Sample", "Plot", "plot"), nm)[1]
  if (is.na(sample_col)) {
    stop(filename, ": cannot identify the sample/plot column.")
  }
  
  element_map <- c(
    Fe = "Fe..mg.L.",
    Al = "Al..mg.L.",
    Mn = "Mn..mg.L.",
    Si = "Si..mg.L.",
    Ca = "Ca..mg.L.",
    K  = "K..mg.L.",
    Na = "Na..mg.L.",
    Mg = "Mg..mg.L."
  )
  
  # Also support alternative column names after make.names().
  find_col <- function(candidates) {
    found <- intersect(candidates, names(dat))
    if (length(found)) found[1] else NA_character_
  }
  
  volume_col <- find_col(c(
    "Corrected.extract.volume..mL.",
    "Corrected.extract.volume.mL.",
    "Corrected.extract.volume..mL",
    "Extract.volume..mL.",
    "Extract.volume.mL."
  ))
  
  mass_col <- find_col(c(
    "Mass.soil..g.",
    "Mass.soil.g.",
    "Soil.mass..g.",
    "Soil.mass.g."
  ))
  
  if (is.na(volume_col) || is.na(mass_col)) {
    stop(filename,
         ": cannot identify corrected extract volume and soil mass columns. ",
         "Check names(dat) and update volume_col/mass_col candidates.")
  }
  
  # Match element concentration columns flexibly.
  find_element_col <- function(element) {
    candidates <- c(
      paste0(element, "..mg.L."),
      paste0(element, ".mg.L."),
      paste0(element, "..mg.L"),
      paste0(element, ".mg.L")
    )
    find_col(candidates)
  }
  
  out <- tibble(
    sample = as.integer(dat[[sample_col]]),
    extract_volume_mL = readr::parse_number(
      as.character(dat[[volume_col]])
    ),
    soil_mass_g = readr::parse_number(
      as.character(dat[[mass_col]])
    )
  )
  
  for (el in names(element_map)) {
    el_col <- find_element_col(el)
    
    if (!is.na(el_col)) {
      conc_mg_L <- readr::parse_number(as.character(dat[[el_col]]))
      
      # mg/L × mL / 1000 = mg extracted; divide by soil mass (g).
      out[[paste0(el, "_mg_g")]] <- ifelse(
        is.finite(conc_mg_L) &
          is.finite(out$extract_volume_mL) &
          out$extract_volume_mL > 0 &
          is.finite(out$soil_mass_g) &
          out$soil_mass_g > 0,
        conc_mg_L * out$extract_volume_mL / 1000 / out$soil_mass_g,
        NA_real_
      )
    }
  }
  
  value_cols <- setdiff(
    names(out), c("sample", "extract_volume_mL", "soil_mass_g")
  )
  if (!length(value_cols)) {
    stop(filename, ": no recognised element concentration columns found.")
  }
  
  # Average replicate measurements if present; output one row per plot.
  out %>%
    group_by(sample) %>%
    summarise(
      across(all_of(value_cols),
             ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
      n_tessier_reps = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(value_cols),
      names_to = "element",
      values_to = "concentration_mg_g"
    ) %>%
    mutate(
      element = str_remove(element, "_mg_g$"),
      variable = paste0("Tessier_", extract_label, "_", element, "_mg_g")
    ) %>%
    select(sample, variable, concentration_mg_g) %>%
    pivot_wider(names_from = variable, values_from = concentration_mg_g)
}

tessier_files <- c(
  NH4OAc = "tessier_feb2026_nh4oac.csv",
  Acetic = "tessier_feb2026_acetic.csv",
  Hydroxylamine = "tessier_feb2026_hydroxylamine.csv",
  H2O2 = "tessier_feb2026_h2o2.csv"
)

tessier_tables <- purrr::imap(
  tessier_files,
  ~ read_tessier(.x, .y)
)

tessier_data <- purrr::reduce(
  tessier_tables, full_join, by = "sample"
)


# ---- CEC and exchangeable cations -------------------------------------------

cec_candidates <- c(
  "CEC_Feb_26.csv",
  "CEC_Feb2026.csv",
  "cec_feb2026.csv",
  "CEC_Feb_2026.csv"
)

cec_filename <- cec_candidates[
  file.exists(file.path(csv_dir, cec_candidates))
][1]

cec_data <- NULL

if (length(cec_filename) > 0 && !is.na(cec_filename)) {
  
  cec_raw <- read_csv_required(cec_filename)
  names(cec_raw) <- make.names(names(cec_raw), unique = TRUE)
  
  # Identify sample/plot column
  sample_col <- intersect(
    c("Sample", "sample", "Plot", "plot"),
    names(cec_raw)
  )[1]
  
  if (is.na(sample_col)) {
    stop(
      cec_filename,
      ": cannot identify Sample/Plot column. Available columns: ",
      paste(names(cec_raw), collapse = ", ")
    )
  }
  
  # Map standard names to the Dutch column names in the CSV
  cec_column_map <- c(
    H = "H...cmolg.kg.1.",
    Ca = "CALCIUM..cmolg.kg.1.",
    K = "KALIUM..cmolg.kg.1.",
    Mg = "MAGNESIUM..cmolg.kg.1.",
    Na = "NATRIUM..cmolg.kg.1.",
    Al = "Aluminium..cmolg.kg.1.",
    Fe = "Ijzer..cmolg.kg.1.",
    Mn = "Mangaan..cmolg.kg.1.",
    CEC = "CEC..cmolg.kg.1.",
    Base_saturation = "Basesaturation...."
  )
  
  available_map <- cec_column_map[
    cec_column_map %in% names(cec_raw)
  ]
  
  if (length(available_map) < 2) {
    stop(
      cec_filename,
      ": fewer than two expected CEC variables found. Available columns: ",
      paste(names(cec_raw), collapse = ", ")
    )
  }
  
  missing_cols <- setdiff(
    unname(cec_column_map),
    names(cec_raw)
  )
  
  if (length(missing_cols) > 0) {
    warning(
      "CEC columns not found and skipped: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Create a clean dataset with standard variable names
  cec_data <- data.frame(
    sample = as.integer(
      readr::parse_number(as.character(cec_raw[[sample_col]]))
    )
  )
  
  for (variable in names(available_map)) {
    source_col <- unname(available_map[variable])
    
    cec_data[[variable]] <- readr::parse_number(
      as.character(cec_raw[[source_col]])
    )
  }
  
  # Average duplicate rows for each sample, if present
  cec_value_names <- setdiff(names(cec_data), "sample")
  
  cec_data <- cec_data |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(cec_value_names),
        ~ if (all(is.na(.x))) {
          NA_real_
        } else {
          mean(.x, na.rm = TRUE)
        }
      ),
      .groups = "drop"
    )
  
  # Prefix variables for the integrated analysis
  cec_names <- names(cec_data)
  cec_names[cec_names != "sample"] <- paste0(
    "CEC_",
    cec_names[cec_names != "sample"]
  )
  names(cec_data) <- cec_names
  
  message("CEC loaded from: ", cec_filename)
  message(
    "CEC variables included: ",
    paste(setdiff(names(cec_data), "sample"), collapse = ", ")
  )
  
} else {
  
  message(
    "CEC file not found in csv_files/. Expected one of: ",
    paste(cec_candidates, collapse = ", ")
  )
  
}

# ---- FTIR peak-area changes (source file is in outputs/, not csv_files/) -----
ftir_path <- file.path(outputs_dir, "ftir_peak_changes_feb2026.csv")
ftir <- NULL
if (file.exists(ftir_path)) {
  ftir_raw <- readr::read_csv(
    ftir_path, show_col_types = FALSE, na = c("", "NA", "#VALUE!")
  )

  if (!"sample" %in% names(ftir_raw)) {
    stop("FTIR file must contain a `sample` column with values such as soil1 ... soil21.")
  }

  # Convert IDs like soil10 to numeric plot IDs (10). Numeric IDs are also accepted.
  ftir <- ftir_raw %>%
    mutate(
      sample = as.character(sample),
      sample = if_else(str_detect(sample, "^soil\\d+$"),
                       str_remove(sample, "^soil"), sample),
      sample = as.integer(sample)
    )

  if (anyNA(ftir$sample)) {
    stop("Could not convert one or more FTIR sample IDs to plot numbers 1-21.")
  }
  if (anyDuplicated(ftir$sample)) {
    stop("FTIR file must have exactly one row per plot; duplicate sample IDs found.")
  }
  if (any(!ftir$sample %in% 1:21)) {
    stop("FTIR sample IDs must correspond to plot numbers 1-21.")
  }

  # Prefer existing change/delta variables in the source table. Calculate a
  # change from *_before / *_after only when no equivalent change variable
  # already exists, avoiding duplicate representations of the same FTIR result.
  change_cols_existing <- names(ftir)[
    str_detect(names(ftir), regex("delta|change_pct|pct_change|percent_change", ignore_case = TRUE))
  ]
  normalise_band_name <- function(x) {
    key <- x %>%
      str_remove(regex("^FTIR_", ignore_case = TRUE)) %>%
      str_remove_all(regex("delta|change_pct|pct_change|percent_change|change|percent|pct", ignore_case = TRUE)) %>%
      str_remove_all("[^A-Za-z0-9]") %>%
      str_to_lower()
    # Harmonise common shorthand/suffix variants so one physical band cannot
    # appear twice under names such as delta_COO_total and COO_change_pct.
    key <- str_remove(key, "(peak)?area$")
    key <- str_remove(key, "^(peak|band)")
    key[key %in% c("coototal", "cooarea", "coocarboxylate", "carboxylate")] <- "coo"
    key[key %in% c("aliphatictotal", "aliphatic", "ch", "chstretch", "chstretching")] <- "ch"
    key[key %in% c("carbonatebend", "carbonatebending", "co3bend")] <- "carbonatebend"
    key[key %in% c("carbonatestretch", "carbonatestretching", "co3stretch")] <- "carbonatestretch"
    key
  }
  existing_change_keys <- normalise_band_name(change_cols_existing)

  before_cols <- names(ftir)[str_ends(names(ftir), "_before")]
  paired_bands <- str_remove(before_cols, "_before")
  paired_bands <- paired_bands[paste0(paired_bands, "_after") %in% names(ftir)]

  for (band in paired_bands) {
    band_key <- normalise_band_name(band)
    has_existing_delta <- any(
      str_detect(str_to_lower(change_cols_existing), "delta") &
        existing_change_keys == band_key
    )
    has_existing_pct <- any(
      str_detect(str_to_lower(change_cols_existing), "change_pct|pct_change|percent_change") &
        existing_change_keys == band_key
    )

    before_name <- paste0(band, "_before")
    after_name <- paste0(band, "_after")
    before_val <- suppressWarnings(as.numeric(ftir[[before_name]]))
    after_val <- suppressWarnings(as.numeric(ftir[[after_name]]))

    if (!has_existing_delta) {
      ftir[[paste0("FTIR_delta_", make.names(band))]] <- after_val - before_val
    }
    # Percent changes are unstable when the initial band area is near zero.
    # Carbonate bend is excluded from newly calculated percentage changes.
    if (!has_existing_pct && band != "Carbonate_bend") {
      ftir[[paste0("FTIR_pct_change_", make.names(band))]] <- ifelse(
        is.finite(before_val) & abs(before_val) > 1e-8,
        100 * (after_val - before_val) / before_val,
        NA_real_
      )
    }
  }

  # Keep ONE FTIR change metric per underlying band. If both an absolute
  # delta and a percentage-change column exist for the same band, prefer the
  # absolute delta. If no delta exists, retain the percentage-change metric.
  # This avoids feeding duplicate representations of the same band into the
  # correlation matrix and PCA.
  change_cols <- names(ftir)[
    str_detect(names(ftir), regex("delta|change_pct|pct_change|percent_change", ignore_case = TRUE))
  ]
  change_cols <- union(change_cols, intersect(
    c("COO_change_pct", "carbonyl_change_pct", "aromatic_change_pct", "CH_change_pct"),
    names(ftir)
  ))

  if (length(change_cols)) {
    change_keys <- normalise_band_name(change_cols)
    change_priority <- case_when(
      str_detect(str_to_lower(change_cols), "delta") ~ 1L,
      str_detect(str_to_lower(change_cols), "change_pct|pct_change|percent_change") ~ 2L,
      TRUE ~ 3L
    )
    selected_change_cols <- tibble(
      column = change_cols,
      band_key = change_keys,
      priority = change_priority
    ) %>%
      arrange(band_key, priority, column) %>%
      distinct(band_key, .keep_all = TRUE) %>%
      pull(column)
  } else {
    selected_change_cols <- character(0)
  }

  ftir_out <- ftir %>% dplyr::select(sample, all_of(selected_change_cols))
  # Prefix source variables to make their origin clear and avoid collisions.
  names(ftir_out) <- ifelse(
    names(ftir_out) == "sample" | str_starts(names(ftir_out), "FTIR_"),
    names(ftir_out), paste0("FTIR_", names(ftir_out))
  )
  ftir <- ftir_out

  write_csv(ftir, file.path(out_dir, "ftir_peak_changes_calculated.csv"))
  message("FTIR loaded from: ", ftir_path, " (", nrow(ftir),
          " plots; existing delta/change variables retained and duplicate changes avoided).")
} else {
  message("FTIR file not found; skipping: ", ftir_path)
}



# ---- Merge one row per plot; assert no row multiplication --------------------
integrated <- metadata %>%
  full_join(cn_wide, by = "sample") %>%
  full_join(agg_data, by = "sample") %>%
  full_join(mineral_data, by = "sample") %>%
  full_join(tessier_data, by = "sample")

if (!is.null(cec_data)) {
  integrated <- integrated %>%
    full_join(cec_data, by = "sample")
}

if (!is.null(ftir)) {
  integrated <- integrated %>%
    full_join(ftir, by = "sample")
}

if (anyDuplicated(integrated$sample)) {
  stop("Integrated data has duplicate sample rows after joining. Check source keys.")
}

integrated <- integrated %>%
  mutate(
    Basalt_type_exp = if_else(
      sample == 9L,
      TRUE,
      Basalt_type_exp
    )
  )

write_csv(
  integrated,
  file.path(out_dir, "integrated_sample_level_data.csv")
)


# ---- Correlation analysis ----------------------------------------------------
# Spearman correlations use pairwise-complete observations. Report n and raw/BH
# adjusted p-values for each pair. Metadata/dose columns are excluded.
correlate_vars <- function(dat, group_name) {
  numeric_data <- dat %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c("sample", names(dat)[str_detect(names(dat), "^n_agg_reps")])) ) %>%
    dplyr::select(where(~sum(is.finite(.x)) >= 3 && is.finite(sd(.x, na.rm = TRUE)) && sd(.x, na.rm = TRUE) > 0))

  if (ncol(numeric_data) < 2) {
    message("Not enough numeric variables for correlations in ", group_name)
    return(invisible(NULL))
  }

  vars <- names(numeric_data)
  cor_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
                    dimnames = list(vars, vars))
  p_mat <- cor_mat
  n_mat <- cor_mat

  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      x <- numeric_data[[i]]
      y <- numeric_data[[j]]
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      n_mat[i, j] <- n
      if (i == j && n >= 3) {
        cor_mat[i, j] <- 1
        p_mat[i, j] <- 0
      } else if (i < j && n >= 3) {
        test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        cor_mat[i, j] <- cor_mat[j, i] <- unname(test$estimate)
        p_mat[i, j] <- p_mat[j, i] <- test$p.value
      }
    }
  }

  pairs <- expand_grid(var1 = vars, var2 = vars) %>%
    filter(match(var1, vars) < match(var2, vars)) %>%
    mutate(
      rho = map2_dbl(var1, var2, ~cor_mat[.x, .y]),
      p_value = map2_dbl(var1, var2, ~p_mat[.x, .y]),
      n = map2_dbl(var1, var2, ~n_mat[.x, .y]),
      p_adj_BH = p.adjust(p_value, method = "BH"),
      abs_rho = abs(rho)
    ) %>%
    arrange(p_adj_BH, desc(abs_rho))

  write_csv(pairs, file.path(out_dir, paste0("spearman_pairs_", group_name, ".csv")))
  write.csv(cor_mat, file.path(out_dir, paste0("spearman_rho_", group_name, ".csv")))
  write.csv(p_mat, file.path(out_dir, paste0("spearman_p_", group_name, ".csv")))
  write.csv(n_mat, file.path(out_dir, paste0("spearman_n_", group_name, ".csv")))

  png(file.path(out_dir, paste0("spearman_heatmap_", group_name, ".png")),
      width = 2600, height = 2400, res = 220)
  corrplot(cor_mat, method = "color", type = "upper", order = "hclust",
           tl.col = "black", tl.cex = 0.55, tl.srt = 45, na.label = " ")
  dev.off()

  invisible(list(cor = cor_mat, p = p_mat, n = n_mat, pairs = pairs, variables = vars))
}

all_results <- correlate_vars(integrated, "all_plots")
# Overlapping experiment flags are intentional: controls and Huhnerberg_50 plots
# can be included in both subsets as defined in treatment_names.csv.
basalt_results <- correlate_vars(filter(integrated, Basalt_type_exp), "basalt_type_exp")
dose_results <- correlate_vars(filter(integrated, Dose_resp_exp), "dose_response_exp")

# ---- Scatterplots between variables from different analyses ------------------
# Only plot cross-domain pairs (e.g. FTIR vs aggregate stability, or FTIR vs
# oxalate/DCB variables). Do not plot pairs within the same analysis/domain.
# PCA variable selection uses cross-analysis pairs with |Spearman rho| >= 0.5.
# Scatterplots are restricted to cross-analysis pairs with |Spearman rho| >= 0.7.
# (The negative cut-offs indicate that these thresholds refer to rho, not p-values.)

classify_analysis <- function(variable) {
  dplyr::case_when(
    stringr::str_starts(variable, "FTIR_") ~ "FTIR",
    
    stringr::str_starts(variable, "Tessier_") ~ "Tessier",
    
    stringr::str_starts(variable, "CEC_") ~ "CEC",
    
    stringr::str_detect(variable, "^(MWD|GMD|D)_") ~
      "Aggregate_stability",
    
    stringr::str_detect(
      variable,
      "^(ox_|dcb_)|_DCB_minus_ox_|_ox_over_DCB$"
    ) ~ "Oxalate_DCB",
    
    stringr::str_starts(variable, "C_") ~ "C_fractions",
    
    TRUE ~ "Other"
  )
}


if (!is.null(all_results) && nrow(all_results$pairs) > 0) {
  treatment_levels <- sort(unique(stats::na.omit(as.character(integrated$Treatment))))
  treatment_palette <- stats::setNames(
    grDevices::hcl.colors(max(1, length(treatment_levels)), palette = "Dynamic"),
    treatment_levels
  )

  cross_analysis_pairs_all <- all_results$pairs %>%
    dplyr::mutate(
      analysis_1 = classify_analysis(var1),
      analysis_2 = classify_analysis(var2)
    ) %>%
    dplyr::filter(
      is.finite(rho), is.finite(p_value),
      analysis_1 != analysis_2,
      analysis_1 != "Other", analysis_2 != "Other"
    ) %>%
    dplyr::arrange(p_value, dplyr::desc(abs_rho))

  # Save all cross-analysis pairs for reference, plus the subset used for plots.
  readr::write_csv(
    cross_analysis_pairs_all,
    file.path(out_dir, "spearman_pairs_all_cross_analysis.csv")
  )

  # Variables entering PCA: any variable involved in at least one cross-analysis
  # association with an absolute Spearman rho of 0.5 or greater.
  selected_cross_pca_vars <- cross_analysis_pairs_all %>%
    dplyr::filter(abs_rho >= 0.5) %>%
    dplyr::select(var1, var2) %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    intersect(all_results$variables)
  readr::write_csv(
    tibble::tibble(variable = selected_cross_pca_vars),
    file.path(out_dir, "PCA_variables_cross_analysis_abs_rho_ge_0.5.csv")
  )

  # Scatterplots only for stronger cross-analysis relationships: |rho| >= 0.7.
  cross_analysis_pairs <- cross_analysis_pairs_all %>%
    dplyr::filter(abs_rho >= 0.7) %>%
    dplyr::arrange(p_value, dplyr::desc(abs_rho))

  readr::write_csv(
    cross_analysis_pairs,
    file.path(out_dir, "spearman_pairs_cross_analysis_scatterplots_abs_rho_ge_0.7.csv")
  )

  if (nrow(cross_analysis_pairs) > 0) {
    for (k in seq_len(nrow(cross_analysis_pairs))) {
      pair <- cross_analysis_pairs[k, ]
      xvar <- pair$var1[[1]]
      yvar <- pair$var2[[1]]
      plot_data <- integrated %>%
        dplyr::select(Treatment, dplyr::all_of(c(xvar, yvar))) %>%
        dplyr::filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))

      p <- ggplot(plot_data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
        geom_point(aes(color = Treatment), size = 2.5, alpha = 0.85, na.rm = TRUE) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7,
                    na.rm = TRUE) +
        labs(
          title = paste0(yvar, " vs ", xvar),
          subtitle = paste0(
            pair$analysis_1[[1]], " vs ", pair$analysis_2[[1]],
            " | Spearman rho = ", round(pair$rho[[1]], 2),
            "; raw p = ", signif(pair$p_value[[1]], 2),
            "; BH-adjusted p = ", signif(pair$p_adj_BH[[1]], 2),
            "; paired n = ", pair$n[[1]]
          ),
          x = xvar, y = yvar, color = "Treatment"
        ) +
        scale_colour_manual(values = treatment_palette, na.translate = FALSE) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom")

      ggsave(
        file.path(out_dir, paste0(
          "scatter_cross_analysis_", sprintf("%03d", k), "_",
          make.names(xvar), "_vs_", make.names(yvar), ".png"
        )),
        p, width = 7, height = 5, dpi = 300
      )
    }
    message("Created ", nrow(cross_analysis_pairs),
            " scatterplots between different analysis domains.")
  } else {
    message("No cross-analysis variable pairs met the plotting criteria.")
  }
}

# ---- Two subset-specific PCAs using strong cross-analysis variables ---------
# Each PCA uses only variables involved in at least one cross-analysis pair with
# |Spearman rho| >= 0.5, and is fitted independently within its experiment.
# The experiment flags intentionally overlap: controls and Huhnerberg 50 t/ha
# plots may occur in both analyses when marked that way in treatment_names.csv.
run_subset_pca <- function(data_subset, subset_name, plot_title, variable_set) {
  safe_name <- make.names(subset_name)
  available_vars <- intersect(variable_set, names(data_subset))
  if (nrow(data_subset) < 3 || length(available_vars) < 2) {
    message("PCA skipped for ", subset_name,
            ": fewer than 3 plots or fewer than 2 selected cross-analysis variables.")
    return(invisible(NULL))
  }

  pca_input <- data_subset %>% dplyr::select(dplyr::all_of(available_vars))
  pca_missing <- tibble(
    variable = names(pca_input),
    n_missing_before_imputation = vapply(pca_input, function(x) sum(!is.finite(x)), integer(1)),
    n_observed = vapply(pca_input, function(x) sum(is.finite(x)), integer(1)),
    median_imputed_value = vapply(pca_input, function(x) {
      if (all(!is.finite(x))) NA_real_ else median(x[is.finite(x)])
    }, numeric(1))
  )

  usable_vars <- pca_missing %>%
    filter(n_observed >= 3, is.finite(median_imputed_value)) %>%
    pull(variable)
  pca_matrix <- as.data.frame(pca_input %>% dplyr::select(all_of(usable_vars)))
  for (nm in names(pca_matrix)) {
    x <- pca_matrix[[nm]]
    x[!is.finite(x)] <- median(x[is.finite(x)])
    pca_matrix[[nm]] <- x
  }
  nonconstant <- names(pca_matrix)[vapply(pca_matrix, function(x) {
    is.finite(sd(x)) && sd(x) > 0
  }, logical(1))]
  pca_matrix <- pca_matrix[, nonconstant, drop = FALSE]

  if (ncol(pca_matrix) < 2) {
    message("PCA skipped for ", subset_name, ": fewer than two usable non-constant variables.")
    return(invisible(NULL))
  }

  pca_fit <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
  variance_pct <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  pca_scores <- as_tibble(pca_fit$x, .name_repair = "minimal") %>%
    mutate(
      sample = data_subset$sample,
      Treatment = as.character(data_subset$Treatment),
      sample_label = paste0("P", sample)
    )
  pca_loadings <- as_tibble(pca_fit$rotation, rownames = "variable")
  pca_variance <- tibble(
    PC = paste0("PC", seq_along(variance_pct)),
    variance_explained_pct = variance_pct,
    cumulative_variance_pct = cumsum(variance_pct)
  )
  contribution <- as_tibble(pca_fit$rotation^2, rownames = "variable") %>%
    mutate(across(-variable, ~100 * .x / sum(.x)))

  write_csv(pca_scores, file.path(out_dir, paste0("PCA_scores_", safe_name, ".csv")))
  write_csv(pca_loadings, file.path(out_dir, paste0("PCA_loadings_", safe_name, ".csv")))
  write_csv(contribution, file.path(out_dir, paste0("PCA_variable_contributions_pct_", safe_name, ".csv")))
  write_csv(pca_variance, file.path(out_dir, paste0("PCA_variance_explained_", safe_name, ".csv")))
  write_csv(pca_missing, file.path(out_dir, paste0("PCA_variable_coverage_and_imputation_", safe_name, ".csv")))

  treatment_levels <- sort(unique(stats::na.omit(pca_scores$Treatment)))
  treatment_palette <- stats::setNames(
    grDevices::hcl.colors(length(treatment_levels), palette = "Dynamic"),
    treatment_levels
  )
  pc1_pct <- round(variance_pct[1], 1)
  pc2_pct <- round(variance_pct[2], 1)
  score_range <- max(abs(c(pca_scores$PC1, pca_scores$PC2)), na.rm = TRUE)
  loading_range <- max(abs(c(pca_fit$rotation[, 1], pca_fit$rotation[, 2])), na.rm = TRUE)
  arrow_scale <- if (is.finite(loading_range) && loading_range > 0) 0.72 * score_range / loading_range else 1
  arrow_df <- pca_loadings %>%
    transmute(variable, PC1 = PC1 * arrow_scale, PC2 = PC2 * arrow_scale)

  p_biplot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey78") +
    stat_ellipse(aes(group = Treatment, colour = Treatment), type = "t", level = 0.68,
                 linewidth = 0.7, show.legend = FALSE, na.rm = TRUE) +
    geom_point(aes(colour = Treatment), size = 3, alpha = 0.9, na.rm = TRUE) +
    geom_text(aes(label = sample_label), colour = "black", size = 3,
              vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
    geom_segment(data = arrow_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 inherit.aes = FALSE,
                 arrow = grid::arrow(length = grid::unit(0.18, "cm")),
                 colour = "grey20", linewidth = 0.55) +
    geom_text(data = arrow_df, aes(x = PC1, y = PC2, label = variable),
              inherit.aes = FALSE, colour = "grey10", size = 2.8,
              check_overlap = TRUE, vjust = -0.45) +
    scale_colour_manual(values = treatment_palette, na.translate = FALSE) +
    labs(title = plot_title,
         subtitle = paste0("Same starting variables as all-plots correlation analysis; n = ", nrow(pca_scores)),
         x = paste0("PC1 (", pc1_pct, "% variance)"),
         y = paste0("PC2 (", pc2_pct, "% variance)"), colour = "Treatment") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical")
  ggsave(file.path(out_dir, paste0("PCA_biplot_", safe_name, ".png")),
         p_biplot, width = 13, height = 10, dpi = 300)

  p_scores <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = Treatment)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey78") +
    stat_ellipse(aes(group = Treatment), type = "t", level = 0.68,
                 linewidth = 0.7, show.legend = FALSE, na.rm = TRUE) +
    geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
    geom_text(aes(label = sample_label), colour = "black", size = 3,
              vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = treatment_palette, na.translate = FALSE) +
    labs(title = paste0(plot_title, " — sample scores"),
         x = paste0("PC1 (", pc1_pct, "% variance)"),
         y = paste0("PC2 (", pc2_pct, "% variance)"), colour = "Treatment") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, paste0("PCA_scores_", safe_name, ".png")),
         p_scores, width = 9, height = 7, dpi = 300)

  scree_df <- pca_variance %>% slice_head(n = min(10, n()))
  p_scree <- ggplot(scree_df, aes(x = factor(PC, levels = PC), y = variance_explained_pct)) +
    geom_col(fill = "grey45") +
    geom_point(aes(y = cumulative_variance_pct), size = 2) +
    geom_line(aes(y = cumulative_variance_pct, group = 1), linewidth = 0.6) +
    labs(title = paste0(plot_title, " — variance explained"), x = "Principal component",
         y = "Variance explained (%) / cumulative variance (%)") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, paste0("PCA_scree_plot_", safe_name, ".png")),
         p_scree, width = 8, height = 5, dpi = 300)

  message("PCA complete for ", subset_name, ": ", ncol(pca_matrix), " variables, ",
          nrow(pca_scores), " plots.")
  invisible(list(scores = pca_scores, loadings = pca_loadings,
                 contribution = contribution, variance = pca_variance))
}

if (!is.null(all_results) && exists("selected_cross_pca_vars") && length(selected_cross_pca_vars) >= 2) {
  # Basalt_type_exp is the application-rate 50 t/ha experiment; Dose_resp_exp
  # is the Huhnerberg dose-response series. Metadata flags define membership.
  pca_application_50 <- run_subset_pca(
    filter(integrated, Basalt_type_exp),
    "application_50_t_ha",
    "PCA — basalt types at 50 t/ha (strong cross-analysis variables)",
    selected_cross_pca_vars
  )
  pca_dose_response <- run_subset_pca(
    filter(integrated, Dose_resp_exp),
    "dose_response",
    "PCA — Huhnerberg dose response (strong cross-analysis variables)",
    selected_cross_pca_vars
  )
} else {
  message("Subset PCAs skipped: fewer than two variables occur in cross-analysis pairs with |Spearman rho| >= 0.5.")
}

# ---- Data coverage and reproducibility --------------------------------------
coverage <- tibble(
  n_plot_rows = nrow(integrated),
  
  plots_with_MAOM_C = if ("C_MAOM" %in% names(integrated))
    sum(!is.na(integrated$C_MAOM)) else 0,
  
  plots_with_oPOM_C = if ("C_oPOM" %in% names(integrated))
    sum(!is.na(integrated$C_oPOM)) else 0,
  
  plots_with_any_Tessier = {
    cols <- grep("^Tessier_", names(integrated), value = TRUE)
    if (length(cols))
      sum(rowSums(!is.na(integrated[, cols, drop = FALSE])) > 0)
    else 0
  },
  
  plots_with_any_CEC = {
    cols <- grep("^CEC_", names(integrated), value = TRUE)
    if (length(cols))
      sum(rowSums(!is.na(integrated[, cols, drop = FALSE])) > 0)
    else 0
  },
  
  plots_with_any_FTIR = {
    cols <- grep("^FTIR_", names(integrated), value = TRUE)
    if (length(cols))
      sum(rowSums(!is.na(integrated[, cols, drop = FALSE])) > 0)
    else 0
  }
)

write_csv(coverage, file.path(out_dir, "data_coverage_summary.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

message("Integrated analysis finished. Results in: ", normalizePath(out_dir))
