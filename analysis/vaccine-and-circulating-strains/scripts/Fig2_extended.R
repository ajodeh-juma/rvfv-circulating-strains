#!/usr/bin/env Rscript
# ==============================================================================
# Extended Figure 2 script
# Adds composition-corrected expected Ts/Tv reference lines and statistical
# comparison of observed vs expected, addressing the reviewer's concern about
# the 2:2:8 combinatorial null model.
#
# New additions vs original script:
#   1. TreeTime equilibrium base frequencies per segment (from log files)
#   2. compute_expected_tstv() — composition-corrected neutral Ts/Tv expectation
#   3. compare_observed_expected() — chi-squared test observed vs expected
#   4. plot_tstv_comparison() — new panel: observed Ts/Tv vs expected per vaccine
#   5. Extended plot_substitutions() annotates each bar chart with both the
#      combinatorial (2:2:8) reference line AND the composition-corrected
#      expected Ts/Tv reference line
#   6. Supplementary table of all Ts, Tv, observed and expected ratios
# ==============================================================================

library(colorspace)
library(dplyr)
library(forcats)
library(scales)
library(ggh4x)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggthemes)
library(glue)
library(lemon)
library(lubridate)
library(paletteer)
library(patchwork)
library(RColorBrewer)
library(rprojroot)
library(tidyverse)
library(tidytext)

# exports
source(sprintf("%s/utils.R", "scripts"))

# ==============================================================================
# 1. TreeTime equilibrium base frequencies
#    Source: riftL.log, riftM.log, riftS-NP.log, riftS-NSS.log
#    These are the GTR equilibrium frequencies (pi_i) estimated by TreeTime
#    for the retained 245-sequence dataset after outlier filtering.
# ==============================================================================
treetime_pi <- list(
  Large   = c(A = 0.2931, C = 0.1956, G = 0.2355, T = 0.2658),  # L segment
  Medium  = c(A = 0.2759, C = 0.2121, G = 0.2452, T = 0.2569),  # M segment
  Small_NP  = c(A = 0.2659, C = 0.2268, G = 0.2659, T = 0.2315), # S-NP
  Small_NSS = c(A = 0.2351, C = 0.2115, G = 0.2493, T = 0.2934)  # S-NSS
)

# ==============================================================================
# 2. Composition-corrected expected Ts/Tv
#
# Under a neutral substitution model with equal rates but empirical base
# frequencies, the expected proportion of transitions and transversions is:
#
#   Expected Ts ∝ 2 * (piA*piG + piC*piT)   [purine-purine + pyrimidine-pyrimidine]
#   Expected Tv ∝ 4 * (piA*piC + piA*piT + piG*piC + piG*piT)
#
# Expected Ts/Tv = (2*(piA*piG + piC*piT)) / (4*(piA*piC + piA*piT + piG*piC + piG*piT))
#
# This gives the Ts/Tv ratio expected when all substitution rates are equal
# but base frequencies are those observed in the RVFV dataset — the
# appropriate composition-corrected null for comparing to the observed ratio.
#
# Note: under equal base frequencies (0.25 each), this simplifies to 0.5,
# matching the combinatorial expectation. Under real RVFV frequencies it
# ranges from 0.248 to 0.253 — nearly identical, confirming that base
# composition alone cannot explain the observed Ts/Tv of 5-10.
# ==============================================================================
compute_expected_tstv <- function(pi, W = NULL) {
  # --- Null 1: equal exchangeability + empirical pi ---
  # Tests whether observed Ts/Tv exceeds a composition-corrected neutral
  # expectation under equal substitution rates between all pairs.
  # This is the appropriate INDEPENDENT test of transition excess.
  ts1 <- 2 * (pi["A"] * pi["G"] + pi["C"] * pi["T"])
  tv1 <- 4 * (pi["A"] * pi["C"] + pi["A"] * pi["T"] +
                pi["G"] * pi["C"] + pi["G"] * pi["T"])
  ratio1 <- as.numeric(ts1 / tv1)
  
  result <- list(
    tstv_expected_equal_rates = ratio1,
    ts_fraction_equal_rates   = as.numeric(ts1 / (ts1 + tv1)),
    tv_fraction_equal_rates   = as.numeric(tv1 / (ts1 + tv1))
  )
  
  # --- Null 2: GTR W_ij * pi_i * pi_j (model-predicted, NOT independent) ---
  # This is the expectation under the FITTED GTR model. W_ij was estimated
  # from the same data, so agreement between observed and this expectation
  # confirms model fit, not neutrality. Reported for transparency only.
  if (!is.null(W)) {
    ts2 <- 2 * (pi["A"] * W["AG"] * pi["G"] + pi["C"] * W["CT"] * pi["T"])
    tv2 <- 2 * (pi["A"] * W["AC"] * pi["C"] + pi["A"] * W["AT"] * pi["T"] +
                  pi["G"] * W["GC"] * pi["C"] + pi["G"] * W["GT"] * pi["T"])
    ratio2 <- as.numeric(ts2 / tv2)
    result$tstv_expected_gtr_model  <- ratio2
    result$ts_fraction_gtr_model    <- as.numeric(ts2 / (ts2 + tv2))
    result$tv_fraction_gtr_model    <- as.numeric(tv2 / (ts2 + tv2))
  }
  
  return(result)
}

# Symmetrized exchangeability rates W_ij from TreeTime log files
# (transition pairs: AG, CT; transversion pairs: AC, AT, GC, GT)
treetime_W <- list(
  Large     = c(AG = 2.9436, CT = 4.7664, AC = 0.1539, AT = 0.2938, GC = 0.0560, GT = 0.0826),
  Medium    = c(AG = 2.9631, CT = 4.1648, AC = 0.2525, AT = 0.4591, GC = 0.0817, GT = 0.1846),
  Small_NP  = c(AG = 2.8828, CT = 4.2778, AC = 0.2900, AT = 0.4072, GC = 0.0804, GT = 0.2371),
  Small_NSS = c(AG = 2.4946, CT = 3.3509, AC = 0.2007, AT = 0.4318, GC = 0.1093, GT = 0.1028)
)

# Pre-compute expected values for each segment (both nulls)
expected_tstv <- mapply(compute_expected_tstv,
                        pi = treetime_pi,
                        W  = treetime_W,
                        SIMPLIFY = FALSE)

cat("\n=== Expected Ts/Tv ratios ===\n")
cat(sprintf("%-12s  %-32s  %-32s\n",
            "Segment",
            "Null 1: equal rates + pi (INDEPENDENT)",
            "Null 2: GTR W_ij + pi (model fit check)"))
for (seg in names(expected_tstv)) {
  cat(sprintf("  %-10s  %-32s  %-32s\n",
              seg,
              sprintf("%.4f (Ts: %.4f, Tv: %.4f)",
                      expected_tstv[[seg]]$tstv_expected_equal_rates,
                      expected_tstv[[seg]]$ts_fraction_equal_rates,
                      expected_tstv[[seg]]$tv_fraction_equal_rates),
              sprintf("%.4f (Ts: %.4f, Tv: %.4f)",
                      expected_tstv[[seg]]$tstv_expected_gtr_model,
                      expected_tstv[[seg]]$ts_fraction_gtr_model,
                      expected_tstv[[seg]]$tv_fraction_gtr_model)))
}

# ==============================================================================
# 3. Statistical comparison: observed Ts/Tv vs composition-corrected expected
#
# Chi-squared goodness-of-fit test:
#   Null: Ts and Tv counts are distributed in proportions predicted by the
#         composition-corrected neutral model.
#   Test: chisq.test(c(obs_ts, obs_tv), p = c(exp_ts_frac, exp_tv_frac))
#   This is a more defensible null than the 2:2:8 combinatorial test because
#   it accounts for unequal base frequencies.
# ==============================================================================
compare_observed_expected <- function(data, segment_name, pi_key) {
  exp_vals <- expected_tstv[[pi_key]]
  
  data %>%
    rowwise() %>%
    mutate(
      # Null 1: composition-corrected equal-rate (INDEPENDENT test)
      exp_tstv_null1   = exp_vals$tstv_expected_equal_rates,
      chisq_p_null1    = tryCatch(
        chisq.test(c(ts, tv),
                   p = c(exp_vals$ts_fraction_equal_rates,
                         exp_vals$tv_fraction_equal_rates))$p.value,
        error = function(e) NA_real_
      ),
      fold_above_null1 = (ts / tv) / exp_vals$tstv_expected_equal_rates,
      
      # Null 2: GTR model-predicted (NOT independent — model fit check only)
      exp_tstv_null2   = exp_vals$tstv_expected_gtr_model,
      chisq_p_null2    = tryCatch(
        chisq.test(c(ts, tv),
                   p = c(exp_vals$ts_fraction_gtr_model,
                         exp_vals$tv_fraction_gtr_model))$p.value,
        error = function(e) NA_real_
      ),
      fold_above_null2 = (ts / tv) / exp_vals$tstv_expected_gtr_model,
      
      # Combinatorial 2:2:8 (original, retained for comparison)
      chisq_p_combinatorial = tryCatch(
        chisq.test(c(ts, tv), p = c(2/12, 8/12))$p.value,
        error = function(e) NA_real_
      ),
      obs_tstv = ts / tv,
      segment  = segment_name
    ) %>%
    ungroup() %>%
    select(segment, reference, ts, tv,
           obs_tstv,
           exp_tstv_null1, fold_above_null1, chisq_p_null1,
           exp_tstv_null2, fold_above_null2, chisq_p_null2,
           chisq_p_combinatorial)
}

# ==============================================================================
# 4. New plot: observed Ts/Tv vs composition-corrected expected
# ==============================================================================
plot_tstv_comparison <- function(obs_data, exp_tstv_val, exp_tstv_gtr,
                                 segment_label, subtitle) {
  plot_data <- obs_data %>%
    select(reference, ts, tv) %>%
    mutate(obs_tstv = ts / tv)
  
  p <- ggplot(plot_data, aes(x = reference, y = obs_tstv, fill = reference)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_hline(aes(yintercept = exp_tstv_val,
                   linetype = "Null 1: equal-rate + pi\n(independent test)"),
               colour = "#D32F2F", linewidth = 0.9) +
    geom_hline(aes(yintercept = exp_tstv_gtr,
                   linetype = "Null 2: GTR W_ij + pi\n(model fit check only)"),
               colour = "#1565C0", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", obs_tstv)),
              vjust = -0.5, size = 3.5,
              fontface = "bold",
              family = dviz_font_family) +
    annotate("text",
             x = Inf, y = exp_tstv_val,
             label = sprintf("Expected = %.3f", exp_tstv_val),
             hjust = 1.1, vjust = -0.5,
             colour = "#D32F2F", size = 3.2,
             fontface = "italic",
             family = dviz_font_family) +
    scale_fill_manual(values = c("MP-12"     = "#3CB7CCFF",
                                 "Smithburn" = "#FF7F0FFF",
                                 "Clone-13"  = "#82853BFF"),
                      guide = "none") +
    scale_linetype_manual(name = NULL,
                          values = c(
                            "Null 1: equal-rate + pi\n(independent test)"  = "dashed",
                            "Null 2: GTR W_ij + pi\n(model fit check only)" = "dotted"
                          )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      subtitle = subtitle,
      x        = "Vaccine reference",
      y        = "Observed Ts/Tv ratio"
    ) +
    theme_dviz_open() +
    theme(
      text             = element_text(family = dviz_font_family),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      plot.subtitle    = element_text(hjust = 0.5, margin = margin(b = 10)),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9)
    ) +
    coord_capped_cart(bottom = capped_horizontal("both"),
                      left   = capped_vertical("both"))
  
  return(p)
}

# ==============================================================================
# 5. Extended plot_substitutions():
#    adds composition-corrected expected Ts/Tv line alongside the existing
#    combinatorial null annotations.
# ==============================================================================
plot_substitutions_extended <- function(data, reference, segment, gene,
                                        exp_vals) {
  # exp_vals: the full list from compute_expected_tstv() for this segment,
  # containing tstv_expected_equal_rates, ts_fraction_equal_rates,
  # tv_fraction_equal_rates, tstv_expected_gtr_model, etc.
  # Passing the full list avoids re-deriving fractions from the ratio
  # (which can silently mis-behave when named scalars flow through mapply).
  exp_tstv_null1 <- as.numeric(exp_vals$tstv_expected_equal_rates)
  ts_frac        <- as.numeric(exp_vals$ts_fraction_equal_rates)
  tv_frac        <- as.numeric(exp_vals$tv_fraction_equal_rates)
  exp_tstv_gtr   <- as.numeric(exp_vals$tstv_expected_gtr_model)
  
  mutations <- c('A->C', 'C->A', 'A->G', 'G->A', 'A->T', 'T->A',
                 'C->G', 'G->C', 'C->T', 'T->C', 'G->T', 'T->G')
  
  ct_tc  <- c('C->T', 'T->C')
  ag_ga  <- c('A->G', 'G->A')
  transversions <- setdiff(mutations, c(ct_tc, ag_ga))
  
  # cross-strain Ts/Tv comparison (all three vaccines)
  strain_stats <- data %>%
    filter(variable %in% c("ts", "tv"),
           reference %in% c("MP-12", "Clone-13", "Smithburn")) %>%
    tidyr::pivot_wider(names_from = variable, values_from = value)
  
  ts_tv_matrix <- as.matrix(strain_stats[, c("ts", "tv")])
  chi_strains  <- chisq.test(ts_tv_matrix)
  p_strains    <- format.pval(chi_strains$p.value, digits = 3)
  
  # within-transition asymmetry (AG/GA vs CT/TC)
  ref_data <- data %>%
    filter(reference == !!reference, variable %in% mutations) %>%
    mutate(group = case_when(
      variable %in% ct_tc        ~ "C->T / T->C",
      variable %in% ag_ga        ~ "A->G / G->A",
      TRUE                       ~ "Other (Transversions)"
    ))
  
  ag_ga_count  <- sum(ref_data$value[ref_data$variable %in% c("A->G", "G->A")])
  ct_tc_count  <- sum(ref_data$value[ref_data$variable %in% c("C->T", "T->C")])
  chi_within   <- chisq.test(c(ag_ga_count, ct_tc_count), p = c(0.5, 0.5))
  p_within     <- format.pval(chi_within$p.value, digits = 3)
  sig_label    <- ifelse(chi_within$p.value < 0.05,
                         "Significant Difference (*)",
                         "No Significance (ns)")
  
  # group-level test vs 2:2:8 combinatorial null
  group_summary <- ref_data %>%
    group_by(group) %>%
    summarise(total_val = sum(value),
              n_types   = n(),
              mean_val  = mean(value))
  obs_counts     <- group_summary$total_val
  expected_probs <- group_summary$n_types / sum(group_summary$n_types)
  stat_test      <- chisq.test(obs_counts, p = expected_probs)
  p_val          <- format.pval(stat_test$p.value, digits = 3)
  
  # group-level test vs composition-corrected expected
  ts_count  <- sum(ref_data$value[ref_data$variable %in% c(ct_tc, ag_ga)])
  tv_count  <- sum(ref_data$value[ref_data$variable %in% transversions])
  
  # fractions already extracted above — no arithmetic on potentially-named scalars
  chi_comp  <- tryCatch(
    chisq.test(c(ts_count, tv_count),
               p = c(ts_frac, tv_frac)),
    error = function(e) {
      warning(sprintf(
        "chisq.test failed for %s/%s: %s", reference, segment, e$message))
      list(p.value = NA_real_)
    }
  )
  p_comp    <- format.pval(chi_comp$p.value, digits = 3)
  obs_ratio <- ts_count / tv_count
  
  max_y    <- max(ref_data$value)
  annot_y  <- max_y * 1.1
  annot_x  <- 12
  
  group_cols <- c("C->T / T->C"         = "#FF7F0FFF",
                  "A->G / G->A"          = "#3CB7CCFF",
                  "Other (Transversions)" = "#82853BFF")
  
  p <- ggplot(ref_data, aes(x = reorder(variable, -value),
                            y = value,
                            fill = group)) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    scale_fill_manual(values = group_cols) +
    labs(
      title    = paste(reference, ""),
      subtitle = glue::glue("Gene: {gene} | Segment: {segment}"),
      x        = "Substitution Type",
      y        = "SNV Count",
      fill     = "Mutation Category"
    ) +
    annotate("text",
             x     = annot_x,
             y     = annot_y,
             label = paste0(
               "Strain comparison (Ts/Tv): p = ", p_strains, "\n",
               "Within-strain (AG/GA vs CT/TC): p = ", p_within,
               " ", sig_label, "\n",
               "Group significance (chi-sq, 2:2:8 null): p = ", p_val, "\n",
               sprintf("Observed Ts/Tv = %.2f | Null 1 (comp-corrected) = %.3f | Null 2 (GTR) = %.3f",
                       obs_ratio, exp_tstv_null1, exp_tstv_gtr), "\n",
               "Chi-sq vs comp.-corrected null: p = ", p_comp
             ),
             hjust  = 1,
             vjust  = 1,
             size   = 3.2,
             fontface = "italic",
             color  = "black") +
    theme_classic() +
    theme_dviz_open() +
    theme(
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.35),
      axis.ticks         = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      plot.subtitle      = element_markdown(size = 11, lineheight = 1.2),
      axis.text.x        = element_blank(),
      legend.position    = c(0.8, 0.5)
    ) +
    coord_capped_cart(bottom = capped_horizontal("both"),
                      left   = capped_vertical("both")) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
  
  return(p)
}

# ==============================================================================
# 6. Load data (same as original script)
# ==============================================================================
prepare_tstv_data <- function(segment, segment_length, alignment_dir) {
  
  dir <- file.path(alignment_dir)
  fns <- file.path(dir,
                   list.files(path = dir,
                              pattern = "\\.parsed.bcftools.stats.csv$",
                              recursive = TRUE))
  
  data_list <- list()
  for (i in seq_along(fns)) {
    fn <- fns[i]
    sample_name <- strsplit(basename(fn), ".", fixed = TRUE)[[1]][3]
    df <- read.csv(fn, sep = ",", header = TRUE)
    df$name <- sample_name
    data_list[[i]] <- df
  }
  d <- do.call(rbind, data_list)
  
  # map accessions to vaccine names
  if (segment == "Medium") {
    d$reference[d$name == "DQ380193"] <- "Smithburn"
    d$reference[d$name == "DQ380208"] <- "MP-12"
    d$reference[d$name == "DQ380213"] <- "Clone-13"
  }
  if (segment == "Large") {
    d$reference[d$name == "DQ375430"] <- "Smithburn"
    d$reference[d$name == "DQ375404"] <- "MP-12"
    d$reference[d$name == "DQ375417"] <- "Clone-13"
  }
  if (segment == "Small") {
    d$reference[d$name == "DQ380157"] <- "Smithburn"
    d$reference[d$name == "DQ380154"] <- "MP-12"
    d$reference[d$name == "DQ380182"] <- "Clone-13"
  }
  
  d <- d %>% dplyr::rename(
    "A->C" = A.C, "A->G" = A.G, "A->T" = A.T,
    "C->A" = C.A, "C->G" = C.G, "C->T" = C.T,
    "G->A" = G.A, "G->C" = G.C, "G->T" = G.T,
    "T->A" = T.A, "T->C" = T.C, "T->G" = T.G
  )
  
  d <- d %>% dplyr::select(
    c(samples, reference, ts, tv, tstv, tsALT, tvALT, tstvALT,
      "A->C", "C->A", "A->G", "G->A", "A->T", "T->A",
      "C->G", "G->C", "C->T", "T->C", "G->T", "T->G")
  )
  
  d.fq <- d %>%
    mutate(
      fAC = `A->C` / segment_length, fCA = `C->A` / segment_length,
      fAG = `A->G` / segment_length, fGA = `G->A` / segment_length,
      fAT = `A->T` / segment_length, fTA = `T->A` / segment_length,
      fCG = `C->G` / segment_length, fGC = `G->C` / segment_length,
      fCT = `C->T` / segment_length, fTC = `T->C` / segment_length,
      fGT = `G->T` / segment_length, fTG = `T->G` / segment_length
    )
  
  write.csv(d, file.path(outDir, paste0(segment, "_tstv.metrics.csv")),
            quote = FALSE, row.names = FALSE)
  
  d.melt   <- reshape2::melt(d)
  dfq.melt <- reshape2::melt(d.fq)
  
  list(d.melt, dfq.melt, d)  # return raw data as third element for new analyses
}

# ==============================================================================
# 7. Load segment data and run extended analyses
# ==============================================================================

inDirLarge  <- sprintf("%s/output-dir/riftL/strain-types",  segmentsDir[[1]])
inDirMedium <- sprintf("%s/output-dir/riftM/strain-types",  segmentsDir[[2]])
inDirNSS    <- sprintf("%s/output-dir/riftS-NSS/strain-types", segmentsDir[[3]])
inDirNP     <- sprintf("%s/output-dir/riftS-NP/strain-types",  segmentsDir[[3]])

tstvLarge  <- prepare_tstv_data("Large",  6276, inDirLarge)
tstvMedium <- prepare_tstv_data("Medium", 3591, inDirMedium)
tstvNSS    <- prepare_tstv_data("Small",  795,  inDirNSS)
tstvNP     <- prepare_tstv_data("Small",  735,  inDirNP)

# ==============================================================================
# 8. Compute observed vs expected comparison for each segment
# ==============================================================================
comp_L   <- compare_observed_expected(tstvLarge[[3]],  "L (RdRp)",    "Large")
comp_M   <- compare_observed_expected(tstvMedium[[3]], "M (Glyco.)",  "Medium")
comp_NSS <- compare_observed_expected(tstvNSS[[3]],    "S-NSS (NSs)", "Small_NSS")
comp_NP  <- compare_observed_expected(tstvNP[[3]],     "S-NP (NP)",   "Small_NP")

tstv_comparison_table <- bind_rows(comp_L, comp_M, comp_NSS, comp_NP) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n=== Observed vs Composition-Corrected Expected Ts/Tv ===\n")
print(tstv_comparison_table)

write.csv(tstv_comparison_table,
          file.path(outDir, "tstv_observed_vs_expected_comparison.csv"),
          row.names = FALSE)

# ==============================================================================
# 9. Generate the new Ts/Tv comparison panels
# ==============================================================================
p_comp_L <- plot_tstv_comparison(
  tstvLarge[[3]],
  exp_tstv_val = expected_tstv$Large$tstv_expected_equal_rates,
  exp_tstv_gtr = expected_tstv$Large$tstv_expected_gtr_model,
  segment_label = "Large", subtitle = "L segment | RdRp"
)
p_comp_M <- plot_tstv_comparison(
  tstvMedium[[3]],
  exp_tstv_val = expected_tstv$Medium$tstv_expected_equal_rates,
  exp_tstv_gtr = expected_tstv$Medium$tstv_expected_gtr_model,
  segment_label = "Medium", subtitle = "M segment | NSm/Gn/Gc"
)
p_comp_NSS <- plot_tstv_comparison(
  tstvNSS[[3]],
  exp_tstv_val = expected_tstv$Small_NSS$tstv_expected_equal_rates,
  exp_tstv_gtr = expected_tstv$Small_NSS$tstv_expected_gtr_model,
  segment_label = "Small", subtitle = "S-NSS segment | NSs"
)
p_comp_NP <- plot_tstv_comparison(
  tstvNP[[3]],
  exp_tstv_val = expected_tstv$Small_NP$tstv_expected_equal_rates,
  exp_tstv_gtr = expected_tstv$Small_NP$tstv_expected_gtr_model,
  segment_label = "Small", subtitle = "S-NP segment | NP"
)

# ==============================================================================
# 10. Generate extended bar chart panels (with composition-corrected annotation)
# ==============================================================================

# Large segment
subsLargeMP12_ext <- plot_substitutions_extended(
  tstvLarge[[1]], "MP-12",     "Large", "RdRp",
  expected_tstv$Large)
subsLargeSmithburn_ext <- plot_substitutions_extended(
  tstvLarge[[1]], "Smithburn", "Large", "RdRp",
  expected_tstv$Large)
subsLargeClone13_ext <- plot_substitutions_extended(
  tstvLarge[[1]], "Clone-13",  "Large", "RdRp",
  expected_tstv$Large)

# Medium segment
subsMediumMP12_ext <- plot_substitutions_extended(
  tstvMedium[[1]], "MP-12",     "Medium", "NSm/Gn/Gc",
  expected_tstv$Medium)
subsMediumSmithburn_ext <- plot_substitutions_extended(
  tstvMedium[[1]], "Smithburn", "Medium", "NSm/Gn/Gc",
  expected_tstv$Medium)
subsMediumClone13_ext <- plot_substitutions_extended(
  tstvMedium[[1]], "Clone-13",  "Medium", "NSm/Gn/Gc",
  expected_tstv$Medium)

# Small-NSS
subsNssMP12_ext <- plot_substitutions_extended(
  tstvNSS[[1]], "MP-12",     "Small", "NSs",
  expected_tstv$Small_NSS)
subsNssSmithburn_ext <- plot_substitutions_extended(
  tstvNSS[[1]], "Smithburn", "Small", "NSs",
  expected_tstv$Small_NSS)
subsNssClone13_ext <- plot_substitutions_extended(
  tstvNSS[[1]], "Clone-13",  "Small", "NSs",
  expected_tstv$Small_NSS)

# Small-NP
subsNpMP12_ext <- plot_substitutions_extended(
  tstvNP[[1]], "MP-12",     "Small", "NP",
  expected_tstv$Small_NP)
subsNpSmithburn_ext <- plot_substitutions_extended(
  tstvNP[[1]], "Smithburn", "Small", "NP",
  expected_tstv$Small_NP)
subsNpClone13_ext <- plot_substitutions_extended(
  tstvNP[[1]], "Clone-13",  "Small", "NP",
  expected_tstv$Small_NP)

# ==============================================================================
# 11. Assemble Figure 2 (extended): bar charts + new Ts/Tv comparison panels
# ==============================================================================
clean_theme <- theme(legend.position = "none",
                     axis.title.y    = element_blank(),
                     axis.title.x    = element_blank())

# Original substitution bar chart panel (extended annotation)
subs_grid <- (
  (subsLargeSmithburn_ext + clean_theme) |
    (subsLargeMP12_ext   + theme(axis.title.y = element_blank()) + clean_theme) |
    (subsLargeClone13_ext + theme(axis.title.y = element_blank()))
) /
  (
    (subsMediumSmithburn_ext + clean_theme) |
      (subsMediumMP12_ext + theme(axis.title.y = element_blank()) + clean_theme) |
      (subsMediumClone13_ext + theme(axis.title.y = element_blank()))
  ) /
  (
    (subsNssSmithburn_ext + clean_theme) |
      (subsNssMP12_ext + theme(axis.title.y = element_blank()) + clean_theme) |
      (subsNssClone13_ext + theme(axis.title.y = element_blank()))
  ) /
  (
    (subsNpSmithburn_ext + clean_theme) |
      (subsNpMP12_ext + theme(axis.title.y = element_blank()) + clean_theme) |
      (subsNpClone13_ext + theme(axis.title.y = element_blank()))
  ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0, 1))

# New Ts/Tv comparison panel (observed vs composition-corrected expected)
tstv_comp_grid <- (p_comp_L | p_comp_M) / (p_comp_NSS | p_comp_NP) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Observed vs. Expected Ts/Tv",
    # subtitle = paste0(
    #   "Red dashed = Null 1: composition-corrected equal-rate (independent test, ~0.25)\n",
    #   "Blue dotted = Null 2: GTR W_ij model-predicted (model fit check only, ~7-11)"
    # ),
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0, 1))

# Save outputs
ggsave(file.path(outDir, "Figure2_extended.pdf"),
       subs_grid,
       width = 18, height = 24, units = "in",
       limitsize = FALSE, dpi = 300, bg = "white", device = cairo_pdf)

ggsave(file.path(outDir, "Figure2_extended.png"),
       subs_grid,
       width = 18, height = 24, units = "in",
       limitsize = FALSE, dpi = 300, bg = "white")

ggsave(file.path(outDir, "Figure2_tstv_comparison.pdf"),
       tstv_comp_grid,
       width = 12, height = 10, units = "in",
       limitsize = FALSE, dpi = 300, bg = "white", device = cairo_pdf)

ggsave(file.path(outDir, "Figure2_tstv_comparison.png"),
       tstv_comp_grid,
       width = 12, height = 10, units = "in",
       limitsize = FALSE, dpi = 300, bg = "white")

cat("\n=== Done ===\n")
cat("Outputs:\n")
cat("  Figure2_extended.pdf/png      — bar chart panels with composition-corrected annotation\n")
cat("  Figure2_tstv_comparison.pdf/png — new panel: observed Ts/Tv vs expected\n")
cat("  tstv_observed_vs_expected_comparison.csv — full statistical comparison table\n")

