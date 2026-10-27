# ==============================================================================
# COMBINED CEC & TESSIER ANALYSES — GROUPED TREATMENTS (Control vs Lime vs Basalts)
# ==============================================================================

# 1. LIBRARIES & SETUP
library(tidyverse)
library(here)
library(emmeans)
library(multcomp)
library(multcompView)
library(broom)
library(ggplot2)
library(car)

theme_set(theme_bw())
library(rstatix)

output_dir <- here("outputs", "groupedtmts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# ==============================================================================
# PART 1: CEC & EXCHANGEABLE CATIONS ANALYSIS (GROUPED)
# ==============================================================================

cat("Running CEC Grouped Analysis...\n")

CEC_Feb_26 <- read.csv(
  here("csv_files", "CEC_Feb_26.csv"),
  header = TRUE
)[1:252, ]

treatment_names <- read.csv(
  here("csv_files", "treatment_names.csv")
)

CEC_Feb_26 <- CEC_Feb_26 %>%
  left_join(
    treatment_names %>% dplyr::rename(Sample = sample),
    by = "Sample"
  ) %>%
  na.omit() %>%
  mutate(
    Treatment = as.character(Treatment),
    Grouped_Treatment = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Lime_2" ~ "Lime",
      Treatment %in% c("Bolsdorfer_50", "Eifelgold_50", "Huhnerberg_50") ~ "Basalts",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grouped_Treatment)) %>%
  mutate(
    Grouped_Treatment = factor(
      Grouped_Treatment,
      levels = c("Control", "Lime", "Basalts")
    )
  )

cec_properties_raw <- c(
  "H...cmolg.kg.1.",
  "CALCIUM..cmolg.kg.1.",
  "KALIUM..cmolg.kg.1.",
  "MAGNESIUM..cmolg.kg.1.",
  "NATRIUM..cmolg.kg.1.",
  "Aluminium..cmolg.kg.1.",
  "Ijzer..cmolg.kg.1.",
  "Mangaan..cmolg.kg.1.",
  "CEC..cmolg.kg.1.",
  "Basesaturation...."
)

# Check that all expected columns exist
missing_cec <- setdiff(cec_properties_raw, names(CEC_Feb_26))

if (length(missing_cec) > 0) {
  stop(
    "Missing CEC columns: ",
    paste(missing_cec, collapse = ", ")
  )
}

# Create a named mapping from raw to clean names
cec_properties_clean <- make.names(cec_properties_raw)
names(cec_properties_clean) <- cec_properties_raw

# Labels for plotting; names must match the clean column names
property_labels <- c(
  "H...cmolg.kg.1."        = "H^'+'~(cmol[c]~kg^{-1})",
  "CALCIUM..cmolg.kg.1."   = "Calcium~(cmol[c]~kg^{-1})",
  "KALIUM..cmolg.kg.1."    = "Potassium~(cmol[c]~kg^{-1})",
  "MAGNESIUM..cmolg.kg.1." = "Magnesium~(cmol[c]~kg^{-1})",
  "NATRIUM..cmolg.kg.1."   = "Sodium~(cmol[c]~kg^{-1})",
  "Aluminium..cmolg.kg.1." = "Aluminium~(cmol[c]~kg^{-1})",
  "Ijzer..cmolg.kg.1."     = "Iron~(cmol[c]~kg^{-1})",
  "Mangaan..cmolg.kg.1."   = "Manganese~(cmol[c]~kg^{-1})",
  "CEC..cmolg.kg.1."       = "CEC~(cmol[c]~kg^{-1})",
  "Basesaturation...."     = "Base~saturation~plain('%')"
)

property_labels_clean <- setNames(
  unname(property_labels[cec_properties_raw]),
  unname(cec_properties_clean)
)

CEC_Feb_26_clean <- CEC_Feb_26 %>%
  rename_with(
    ~ unname(cec_properties_clean[.x]),
    all_of(cec_properties_raw)
  )

# Use the cleaned column names consistently
cec_properties <- unname(cec_properties_clean)



# ==============================================================================
# CEC: WELCH ANOVA + GAMES-HOWELL POST-HOC TESTS
# ==============================================================================

run_grouped_test <- function(data, prop) {
  
  d <- data %>%
    dplyr::select(
      dplyr::all_of(c(prop, "Grouped_Treatment"))
    ) %>%
    dplyr::filter(
      is.finite(.data[[prop]]),
      !is.na(Grouped_Treatment)
    ) %>%
    dplyr::mutate(
      Grouped_Treatment = droplevels(Grouped_Treatment)
    )
  
  groups <- levels(d$Grouped_Treatment)
  
  if (length(groups) < 2 ||
      any(table(d$Grouped_Treatment) < 2)) {
    warning("Insufficient replication for CEC property: ", prop)
    return(NULL)
  }
  
  # Overall Welch ANOVA
  welch <- tryCatch(
    stats::oneway.test(
      stats::reformulate("Grouped_Treatment", response = prop),
      data = d,
      var.equal = FALSE
    ),
    error = function(e) NULL
  )
  
  overall <- tibble::tibble(
    Property = prop,
    Test = "Welch ANOVA",
    P_Overall = if (is.null(welch)) NA_real_ else welch$p.value
  )
  
  # Games-Howell pairwise comparisons
  gh <- tryCatch(
    rstatix::games_howell_test(
      d,
      formula = stats::reformulate(
        "Grouped_Treatment",
        response = prop
      )
    ),
    error = function(e) {
      warning("Games-Howell failed for ", prop, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(gh) || nrow(gh) == 0) {
    return(list(
      overall = overall,
      pairs = tibble::tibble(),
      letters = tibble::tibble()
    ))
  }
  
  gh <- gh %>%
    dplyr::mutate(Property = prop)
  
  # Compact letters from adjusted Games-Howell p-values
  pvals <- stats::setNames(
    gh$p.adj,
    paste(gh$group1, gh$group2, sep = "-")
  )
  
  # Missing comparisons are conservatively treated as non-significant
  pvals[is.na(pvals)] <- 1
  
  cld <- multcompView::multcompLetters(pvals)$Letters
  
  letters_df <- tibble::tibble(
    Grouped_Treatment = names(cld),
    Letters = unname(cld),
    Property = prop
  )
  
  list(
    overall = overall,
    pairs = gh,
    letters = letters_df
  )
}

cec_results_list <- purrr::map(
  cec_properties,
  ~ run_grouped_test(CEC_Feb_26_clean, .x)
)

cec_test_results <- purrr::map_dfr(
  cec_results_list,
  ~ if (is.null(.x)) tibble::tibble() else
    dplyr::left_join(
      .x$letters,
      .x$overall,
      by = "Property"
    )
)

cec_pairwise_results <- purrr::map_dfr(
  cec_results_list,
  ~ if (is.null(.x)) tibble::tibble() else .x$pairs
)

cec_overall_results <- purrr::map_dfr(
  cec_results_list,
  ~ if (is.null(.x)) tibble::tibble() else .x$overall
)

# Export statistical results
write.csv(
  cec_overall_results,
  file.path(output_dir, "CEC_grouped_overall_tests.csv"),
  row.names = FALSE
)

write.csv(
  cec_pairwise_results,
  file.path(output_dir, "CEC_grouped_GamesHowell_pairwise.csv"),
  row.names = FALSE
)


cec_summary <- CEC_Feb_26_clean %>%
  group_by(Grouped_Treatment) %>%
  summarise(
    across(
      all_of(cec_properties),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) /
          sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

cec_long <- cec_summary %>%
  pivot_longer(
    cols = -Grouped_Treatment,
    names_to = c("Property", ".value"),
    names_pattern = "^(.*)_(mean|se)$"
  )

cec_plot_data <- cec_long %>%
  left_join(
    cec_test_results %>%
      dplyr::select(Grouped_Treatment, Property, Letters),
    by = c("Grouped_Treatment", "Property")
  ) %>%
  group_by(Property) %>%
  mutate(
    y_position = max(mean + se, na.rm = TRUE) +
      0.12 * max(
        diff(range(c(mean - se, mean + se), na.rm = TRUE)),
        abs(max(mean + se, na.rm = TRUE)) * 0.05,
        0.01
      )
  ) %>%
  ungroup()

p_cec_grouped <- ggplot(
  cec_plot_data,
  aes(x = Grouped_Treatment, y = mean)
) +
  geom_point(
    size = 3, shape = 21,
    fill = "white", colour = "black", stroke = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, linewidth = 0.7
  ) +
  geom_text(
    aes(y = y_position, label = Letters),
    colour = "black", fontface = "bold", size = 4, na.rm = TRUE
  ) +
  facet_wrap(
    ~ Property,
    scales = "free_y",
    ncol = 2,
    labeller = as_labeller(
      property_labels_clean,
      default = label_parsed
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  labs(
    x = "",
    y = "Mean ± SE",
    title = "CEC and Exchangeable Cations — Grouped Treatments"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 10
    ),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(
  file.path(output_dir, "CEC_grouped_plot.png"),
  p_cec_grouped,
  width = 10, height = 14, dpi = 300
)

write.csv(
  cec_summary,
  file.path(output_dir, "CEC_grouped_summary.csv"),
  row.names = FALSE
)

write.csv(
  cec_test_results,
  file.path(output_dir, "CEC_grouped_tests.csv"),
  row.names = FALSE
)

# ==============================================================================
# PART 2: TESSIER EXTRACTIONS ANALYSIS (GROUPED)
# ==============================================================================

cat("Running Tessier Grouped Analysis...\n")

metadata <- read.csv(
  here("csv_files", "treatment_names.csv")
) %>%
  mutate(
    sample = as.numeric(sample),
    Treatment = as.character(Treatment),
    Grouped_Treatment = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Lime_2" ~ "Lime",
      Treatment %in% c(
        "Bolsdorfer_50",
        "Eifelgold_50",
        "Huhnerberg_50"
      ) ~ "Basalts",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grouped_Treatment)) %>%
  dplyr::select(sample, Treatment, Grouped_Treatment)

if (anyDuplicated(metadata$sample)) {
  stop("Duplicate sample IDs in treatment_names.csv.")
}

process_tessier_grouped <- function(df, extract_name) {
  
  # Standardise column names for matching
  names(df) <- trimws(names(df))
  sample_candidates <- c(
    "sample", "Sample", "SAMPLE",
    "Sample.ID", "Sample_ID",
    "SampleID", "sample_id"
  )
  
  sample_col <- intersect(sample_candidates, names(df))
  
  if (length(sample_col) == 0) {
    stop(
      "No sample ID column found in Tessier file ",
      extract_name, ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  names(df)[match(sample_col[1], names(df))] <- "sample"
  
  # Continue with the existing required-column checks and processing
  required <- c(
    "sample",
    "Fe..mg.L.", "Al..mg.L.", "Mn..mg.L.",
    "Mg..mg.L.", "Si..mg.L.", "Ca..mg.L.",
    "K..mg.L.", "Na..mg.L.",
    "Corrected.extract.volume..mL.",
    "Mass.soil..g."
  )
  
  missing_cols <- setdiff(required, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in Tessier file ", extract_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(sample = as.numeric(sample)) %>%
    left_join(metadata, by = "sample") %>%
    filter(!is.na(Grouped_Treatment)) %>%
    mutate(
      Extract = extract_name,
      across(
        c(
          Corrected.extract.volume..mL.,
          Mass.soil..g.
        ),
        as.numeric
      ),
      across(
        c(
          Fe..mg.L., Al..mg.L., Mn..mg.L.,
          Mg..mg.L., Si..mg.L., Ca..mg.L.,
          K..mg.L., Na..mg.L.
        ),
        as.numeric
      ),
      across(
        c(
          Fe..mg.L., Al..mg.L., Mn..mg.L.,
          Mg..mg.L., Si..mg.L., Ca..mg.L.,
          K..mg.L., Na..mg.L.
        ),
        ~ .x * Corrected.extract.volume..mL. /
          (1000 * Mass.soil..g.),
        .names = "{.col}_mg_g"
      )
    ) %>%
    dplyr::select(
      sample, Grouped_Treatment, Extract,
      Fe = Fe..mg.L._mg_g,
      Al = Al..mg.L._mg_g,
      Mn = Mn..mg.L._mg_g,
      Mg = Mg..mg.L._mg_g,
      Si = Si..mg.L._mg_g,
      Ca = Ca..mg.L._mg_g,
      K = K..mg.L._mg_g,
      Na = Na..mg.L._mg_g
    ) %>%
    pivot_longer(
      cols = c(Fe, Al, Mn, Mg, Si, Ca, K, Na),
      names_to = "Element",
      values_to = "Concentration"
    ) %>%
    filter(
      is.finite(Concentration),
      !is.na(Concentration)
    )
}

tessier_combined <- bind_rows(
  process_tessier_grouped(read.csv(here("csv_files", "tessier_feb2026_nh4oac.csv")), "NH4OAc"),
  process_tessier_grouped(read.csv(here("csv_files", "tessier_feb2026_acetic.csv")), "Acetic"),
  process_tessier_grouped(read.csv(here("csv_files", "tessier_feb2026_hydroxylamine.csv")), "Hydroxylamine"),
  process_tessier_grouped(read.csv(here("csv_files", "tessier_feb2026_h2o2.csv")), "H2O2")
) %>%
  mutate(
    Extract = factor(Extract, levels = c("NH4OAc", "Acetic", "Hydroxylamine", "H2O2")),
    Element = factor(Element, levels = c("Ca", "Na", "K", "Mg", "Fe", "Al", "Mn", "Si")),
    Grouped_Treatment = factor(Grouped_Treatment, levels = c("Control", "Lime", "Basalts"))
  )

# Calculate difference from Control
control_means <- tessier_combined %>%
  filter(Grouped_Treatment == "Control") %>%
  group_by(Element, Extract) %>%
  summarise(Control_mean = mean(Concentration, na.rm = TRUE), .groups = "drop")

tessier_diff <- tessier_combined %>%
  left_join(control_means, by = c("Element", "Extract")) %>%
  mutate(Difference_mg_g = Concentration - Control_mean)


# ==============================================================================
# TESSIER: WELCH ANOVA + GAMES-HOWELL POST-HOC TESTS
# ==============================================================================

tessier_results_list <- tessier_combined %>%
  dplyr::group_by(Element, Extract) %>%
  dplyr::group_split() %>%
  purrr::map(function(dat) {
    
    element_name <- as.character(dat$Element[1])
    extract_name <- as.character(dat$Extract[1])
    
    dat <- dat %>%
      dplyr::filter(
        is.finite(Concentration),
        !is.na(Grouped_Treatment)
      ) %>%
      dplyr::mutate(
        Grouped_Treatment = droplevels(Grouped_Treatment)
      )
    
    groups <- levels(dat$Grouped_Treatment)
    
    if (length(groups) < 2 ||
        any(table(dat$Grouped_Treatment) < 2)) {
      warning(
        "Insufficient replication: ",
        element_name, " / ", extract_name
      )
      return(NULL)
    }
    
    # Overall Welch ANOVA
    welch <- tryCatch(
      stats::oneway.test(
        Concentration ~ Grouped_Treatment,
        data = dat,
        var.equal = FALSE
      ),
      error = function(e) NULL
    )
    
    overall <- tibble::tibble(
      Element = element_name,
      Extract = extract_name,
      Test = "Welch ANOVA",
      P_Overall = if (is.null(welch)) NA_real_ else welch$p.value
    )
    
    # Pairwise Games-Howell
    gh <- tryCatch(
      rstatix::games_howell_test(
        dat,
        Concentration ~ Grouped_Treatment
      ),
      error = function(e) {
        warning(
          "Games-Howell failed for ",
          element_name, " / ", extract_name,
          ": ", e$message
        )
        NULL
      }
    )
    
    if (is.null(gh) || nrow(gh) == 0) {
      return(list(
        overall = overall,
        pairs = tibble::tibble(),
        letters = tibble::tibble()
      ))
    }
    
    gh <- gh %>%
      dplyr::mutate(
        Element = element_name,
        Extract = extract_name
      )
    
    pvals <- stats::setNames(
      gh$p.adj,
      paste(gh$group1, gh$group2, sep = "-")
    )
    
    pvals[is.na(pvals)] <- 1
    
    cld <- multcompView::multcompLetters(pvals)$Letters
    
    letters_df <- tibble::tibble(
      Grouped_Treatment = names(cld),
      .group = unname(cld),
      Element = element_name,
      Extract = extract_name
    )
    
    list(
      overall = overall,
      pairs = gh,
      letters = letters_df
    )
  })

tessier_overall_results <- purrr::map_dfr(
  tessier_results_list,
  ~ if (is.null(.x)) tibble::tibble() else .x$overall
)

tessier_pairwise_results <- purrr::map_dfr(
  tessier_results_list,
  ~ if (is.null(.x)) tibble::tibble() else .x$pairs
)

tessier_letters <- purrr::map_dfr(
  tessier_results_list,
  ~ if (is.null(.x)) tibble::tibble() else .x$letters
)

# Add label positions based on the plotted differences from control
tessier_label_positions <- tessier_diff %>%
  dplyr::group_by(Element, Extract) %>%
  dplyr::summarise(
    y = max(Difference_mg_g, na.rm = TRUE) +
      0.12 * max(
        diff(range(Difference_mg_g, na.rm = TRUE)),
        0.01
      ),
    .groups = "drop"
  )

tessier_letters <- tessier_letters %>%
  dplyr::left_join(
    tessier_label_positions,
    by = c("Element", "Extract")
  )

# Export Tessier statistical results
write.csv(
  tessier_overall_results,
  file.path(output_dir, "Tessier_grouped_overall_tests.csv"),
  row.names = FALSE
)

write.csv(
  tessier_pairwise_results,
  file.path(output_dir, "Tessier_grouped_GamesHowell_pairwise.csv"),
  row.names = FALSE
)


tessier_plot_data <- tessier_diff %>%
  left_join(tessier_letters, by = c("Element", "Extract", "Grouped_Treatment"))

extract_colors <- c(
  "NH4OAc"        = "#E69F00",
  "Acetic"        = "#56B4E9",
  "Hydroxylamine" = "#009E73",
  "H2O2"          = "#CC79A7"
)

p_tessier_grouped <- ggplot(tessier_plot_data, aes(x = Grouped_Treatment, y = Difference_mg_g, fill = Extract)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_point(aes(color = Extract), position = position_jitter(width = 0.08), size = 1.5, alpha = 0.8) +
  geom_text(
    data = tessier_letters,
    aes(
      x = Grouped_Treatment,
      y = y,
      label = .group
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 3.5,
    fontface = "bold",
    vjust = 0
  )+
  facet_grid(Element ~ Extract, scales = "free_y") +
  scale_fill_manual(values = extract_colors) +
  scale_color_manual(values = extract_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Treatment Group",
    y = expression("Difference from Control (mg g"^-1*" soil)"),
    title = "Tessier Extract Comparisons — Grouped Treatments (Control vs Lime vs Basalts)",
    subtitle = "Absolute difference from Control; Tukey HSD within each Element × Extract"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.margin = margin(15, 15, 15, 15))

ggsave(file.path(output_dir, "Tessier_grouped_plot.png"), p_tessier_grouped, width = 12, height = 16, dpi = 300)

cat("Analysis complete! All grouped results and plots saved successfully.\n")