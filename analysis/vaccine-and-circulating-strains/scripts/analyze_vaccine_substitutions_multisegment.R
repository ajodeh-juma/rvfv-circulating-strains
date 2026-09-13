# ==============================================================================
# Script Name: analyze_vaccine_substitutions_multisegment.R
#
# v4 CHANGE: added `diagnostic_filter` parameter to analyze_vaccine_substitutions().
# When a full scan is run (target_positions = NULL), the output can have
# 1000+ rows (e.g. 1104 nucleotide positions for M/NC_014396), most of which
# are uninformative "Wild-Type / Field Variant" rows. diagnostic_filter lets
# you restrict what gets WRITTEN to output_csv to only the classification(s)
# you care about (e.g. "Strict Target Vaccine Marker"), while the function
# still RETURNS the full, unfiltered table in R (so FDR correction is always
# computed across the complete set of tested positions BEFORE any filtering,
# and downstream/combined analyses are never silently missing rows).
# ==============================================================================
# Script Name: analyze_vaccine_substitutions_multisegment_v3.R
# Author: Viral Genomics & Bioinformatics
# Description: Evaluates vaccine-specific substitutions across L, M, and S segments.
#              Auto-detects (a) whether a run is SELF-REFERENCED (vaccine_reference
#              == target_accession, from the rvfvcirculatingstrains subworkflow) or
#              EXTERNALLY-REFERENCED (an independent RefSeq genome, from the
#              rvfvmutationalprofiling subworkflow), and (b) whether the input
#              mutations file is AMINO-ACID-level (e.g. "I219T") or
#              NUCLEOTIDE-level (e.g. "C656T"). Classification and position
#              filtering are handled correctly for any combination of the two.
#
# ------------------------------------------------------------------------------
# CHANGE LOG
#
# v1: Self-reference detection added -- vaccine_reference == target_accession
#     makes k_vac trivially 0 under the original k_vac>0 classification
#     criterion; classification made mode-aware (field-divergence-based when
#     self-referenced, original k_vac/Fisher-based when externally-referenced).
#
# v2: Coordinate-system auto-detection added -- self-referenced files record
#     AMINO ACID substitutions (e.g. "I219T"); externally-referenced files
#     record NUCLEOTIDE substitutions (e.g. "C656T"). target_positions are
#     converted between systems via codon arithmetic before filtering, and
#     every position is mapped to amino-acid space for domain classification.
#
# v3 (this version): host_distribution, temporal_distribution, and
#     lineage_distribution now report EXCLUSIVITY rather than only a flat
#     breakdown. For each substitution, if every circulating carrier shares
#     the same host / year / lineage, this is now flagged explicitly
#     ("Exclusive to <value> (n=X)"); otherwise the substitution "spans"
#     multiple categories and the full breakdown is still shown for
#     reference. New boolean columns host_exclusive / temporal_exclusive /
#     lineage_exclusive are added for easy filtering. geographic_distribution
#     (country) is unchanged (full breakdown only), since exclusivity was
#     requested only for host/lineage/year.
#
#     PRIMARY METHOD NOTE: the self-referenced (vaccine-accession-as-reference)
#     scheme is the primary analysis. The externally-referenced (RefSeq
#     nucleotide) scheme is retained as a COMPLEMENTARY / validation approach:
#     it lets vaccine-vs-circulating be tested with a standard two-group
#     Fisher's exact test against an independent basis (rather than a single-
#     genome frequency estimate), it cross-validates that markers identified
#     under self-referencing are not artifacts of that coordinate choice, it
#     places different vaccines on one shared coordinate system for direct
#     cross-vaccine comparison, and it offers nucleotide-level resolution
#     within a codon (useful for primer/probe design) that amino-acid-level
#     calls collapse away.
# ==============================================================================
# exports
source(sprintf("%s/utils.R", "scripts"))


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyverse)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 1. Structural & Functional Domain Classifier for L, M, and S Segments
#    (unchanged; operates purely in amino-acid space)
# ------------------------------------------------------------------------------
get_protein_domain <- function(pos, segment_label) {
  seg <- toupper(trimws(segment_label))
  
  if (seg == "L") {
    protein <- "L protein (RdRp)"
    if (pos >= 1 && pos <= 214) {
      domain <- "Endonuclease Domain (ENDO)"
    } else if (pos >= 215 && pos <= 600) {
      domain <- "N-terminal Domain / Linker (PA-C like)"
    } else if (pos >= 601 && pos <= 1200) {
      domain <- "RdRp Core Domain (Palm/Fingers/Thumb, Motifs A-F)"
    } else if (pos >= 1201 && pos <= 2092) {
      domain <- "C-terminal Cap-Binding / Lariat Domain"
    } else {
      domain <- "Extended C-terminal / Unmapped Region"
    }
    
  } else if (seg == "M") {
    if (pos >= 1 && pos <= 153) {
      protein <- "NSm"
      domain  <- if (pos <= 50) "NSm N-terminal Region" else if (pos <= 120) "NSm Core Domain (14kDa)" else "NSm C-terminal / Cleavage Boundary"
    } else if (pos >= 154 && pos <= 687) {
      protein <- "Gn"
      domain  <- if (pos <= 300) "Gn N-terminal Ectodomain (Domain I)" else if (pos <= 500) "Gn Mid-Ectodomain (Domain II)" else if (pos <= 580) "Gn Stem Region" else "Gn Transmembrane / Cytosolic Tail"
    } else if (pos >= 688 && pos <= 1197) {
      protein <- "Gc"
      domain  <- if (pos >= 730 && pos <= 760) "Gc Fusion Loop" else if (pos <= 850) "Gc Head / Domain II" else if (pos <= 1000) "Gc Domain I / III" else "Gc C-terminal / Transmembrane Domain"
    } else {
      protein <- "M Polyprotein"
      domain  <- "Unmapped M Region"
    }
    
  } else if (seg %in% c("S-NP", "NP")) {
    protein <- "Nucleocapsid (N) Protein"
    if (pos >= 1 && pos <= 30) {
      domain <- "N-terminal Arm (Oligomerization)"
    } else if (pos >= 31 && pos <= 180) {
      domain <- "N-lobe Core (RNA-binding Domain)"
    } else if (pos >= 181 && pos <= 245) {
      domain <- "C-lobe (Oligomerization & Stalk)"
    } else {
      domain <- "Unmapped N Region"
    }
    
  } else if (seg %in% c("S-NSS", "NSS")) {
    protein <- "Non-structural (NSs) Protein"
    if (pos >= 1 && pos <= 80) {
      domain <- "NSs N-terminal Self-association Domain"
    } else if (pos >= 81 && pos <= 210) {
      domain <- "NSs Core Domain (PKR / E3-ligase interaction)"
    } else if (pos >= 211 && pos <= 265) {
      domain <- "NSs C-terminal Tail (TFIIH p62 degradation)"
    } else {
      domain <- "Unmapped NSs Region"
    }
    
  } else {
    protein <- "Unknown"
    domain  <- "Unmapped Segment"
  }
  
  return(list(protein = protein, domain = domain))
}

# ------------------------------------------------------------------------------
# 2. Coordinate-system helpers
# ------------------------------------------------------------------------------

.segment_lengths <- list(
  "L"     = list(aa = 2092, nt = 6276),
  "M"     = list(aa = 1197, nt = 3591),
  "S-NP"  = list(aa = 245,  nt = 735),
  "S-NSS" = list(aa = 265,  nt = 795)
)

.normalize_segment_key <- function(segment_label) {
  seg <- toupper(trimws(segment_label))
  if (seg %in% c("S-NP", "NP")) return("S-NP")
  if (seg %in% c("S-NSS", "NSS")) return("S-NSS")
  return(seg)
}

detect_coordinate_system <- function(muts, segment_label) {
  key <- .normalize_segment_key(segment_label)
  expected <- .segment_lengths[[key]]
  max_pos  <- suppressWarnings(max(muts$positions, na.rm = TRUE))
  
  if (is.null(expected) || !is.finite(max_pos)) {
    return(if (max_pos > 2500) "nucleotide" else "amino_acid")
  }
  
  d_aa <- abs(max_pos - expected$aa)
  d_nt <- abs(max_pos - expected$nt)
  if (d_nt < d_aa) "nucleotide" else "amino_acid"
}

convert_positions <- function(positions, from_system, to_system) {
  if (from_system == to_system) return(sort(unique(positions)))
  
  if (from_system == "amino_acid" && to_system == "nucleotide") {
    nt <- unlist(lapply(positions, function(aa) {
      start <- 3 * (aa - 1) + 1
      start:(start + 2)
    }))
    return(sort(unique(nt)))
  }
  
  if (from_system == "nucleotide" && to_system == "amino_acid") {
    return(sort(unique(ceiling(positions / 3))))
  }
  
  stop(sprintf("Unrecognized coordinate system conversion: %s -> %s", from_system, to_system))
}

map_to_amino_acid_position <- function(positions, coordinate_system) {
  if (coordinate_system == "amino_acid") return(positions)
  return(ceiling(positions / 3))
}

# ------------------------------------------------------------------------------
# 3. Main Analysis Function
# ------------------------------------------------------------------------------
analyze_vaccine_substitutions <- function(
    traits_file,
    mutations_file,
    target_accession,
    vaccine_name              = "Vaccine",
    segment_label              = "M",
    output_csv                 = NULL,
    target_positions           = NULL,
    target_position_system     = "amino_acid",
    delim                      = "; ",
    field_marker_threshold     = 99.0,
    exclude_nontarget_vaccines = TRUE,
    diagnostic_filter          = NULL   # (v4) character vector of
    # diagnostic_category CODES to keep in the WRITTEN CSV only (see v5 note
    # below for why this now filters on `diagnostic_category`, a stable code,
    # rather than the human-readable `potential_diagnostic_value` text).
    # Valid self-referenced codes: "vaccine_residue_not_detected",
    # "vaccine_residue_low_frequency", "invariant_site", "polymorphic_site".
    # Valid external-reference codes: "no_shared_occurrence",
    # "low_field_frequency", "shared_lineage_variant", "field_variant".
    # e.g. diagnostic_filter = c("no_shared_occurrence", "low_field_frequency").
    # Default NULL = no filtering, full table written as before. The
    # function's RETURN VALUE is always the full, unfiltered table.
) {
  
  cat(sprintf("\n======================================================================\n"))
  cat(sprintf(">>> Running Analysis: %s (Accession: %s) | Segment: %s\n", vaccine_name, target_accession, segment_label))
  cat(sprintf("======================================================================\n"))
  
  if (!file.exists(traits_file)) stop(paste("Traits file not found:", traits_file))
  if (!file.exists(mutations_file)) stop(paste("Mutations file not found:", mutations_file))
  
  traits <- read.delim(traits_file, header = TRUE, stringsAsFactors = FALSE)
  muts   <- read.csv(mutations_file, header = TRUE, stringsAsFactors = FALSE)
  
  format_counts <- function(vec, delimiter = delim) {
    clean_vec <- na.omit(as.character(vec))
    clean_vec <- clean_vec[trimws(clean_vec) != "" & clean_vec != "NA"]
    if (length(clean_vec) == 0) return("None")
    tbl <- table(clean_vec)
    sorted_names <- sort(names(tbl))
    formatted <- sprintf("%s (%d)", sorted_names, as.numeric(tbl[sorted_names]))
    return(paste(formatted, collapse = delimiter))
  }
  
  # NEW (v3): exclusivity-aware count formatting for host/lineage/year.
  format_exclusivity <- function(vec, delimiter = delim) {
    clean_vec <- na.omit(as.character(vec))
    clean_vec <- clean_vec[trimws(clean_vec) != "" & clean_vec != "NA"]
    if (length(clean_vec) == 0) return("None")
    uniq_vals <- unique(clean_vec)
    if (length(uniq_vals) == 1) {
      return(sprintf("Exclusive to %s (n=%d)", uniq_vals[1], length(clean_vec)))
    }
    tbl <- table(clean_vec)
    sorted_names <- sort(names(tbl))
    formatted <- sprintf("%s (%d)", sorted_names, as.numeric(tbl[sorted_names]))
    return(sprintf("Not exclusive - spans %d categories: %s", length(uniq_vals), paste(formatted, collapse = delimiter)))
  }
  
  is_exclusive <- function(vec) {
    clean_vec <- na.omit(as.character(vec))
    clean_vec <- clean_vec[trimws(clean_vec) != "" & clean_vec != "NA"]
    if (length(clean_vec) == 0) return(NA)
    length(unique(clean_vec)) == 1
  }
  
  coordinate_system <- detect_coordinate_system(muts, segment_label)
  cat(sprintf("[i] Detected mutations file coordinate system: %s\n", coordinate_system))
  
  all_vaccine_accs <- traits$accession[tolower(traits$strain_type) == "vaccine"]
  
  if (!target_accession %in% traits$accession) {
    warning(sprintf("Target accession '%s' for %s not found in traits file! Skipping...", target_accession, vaccine_name))
    return(NULL)
  }
  
  cat(sprintf("[+] Target Vaccine: %s [%s]\n", vaccine_name, target_accession))
  
  if (exclude_nontarget_vaccines) {
    excluded_accessions <- setdiff(all_vaccine_accs, target_accession)
    cat(sprintf("[-] Excluding %d non-target vaccine accession(s) from table (no statistical effect; cosmetic only)\n", length(excluded_accessions)))
    traits_clean <- traits[!traits$accession %in% excluded_accessions, ]
    muts_clean   <- muts[!muts$accession %in% excluded_accessions, ]
  } else {
    cat("[i] Non-target vaccine accessions retained (exclude_nontarget_vaccines = FALSE).\n")
    traits_clean <- traits
    muts_clean   <- muts
  }
  
  self_referenced <- !(target_accession %in% muts$accession)
  
  if (self_referenced) {
    cat(sprintf("[i] Detected SELF-REFERENCED run (vaccine_reference == target_accession).\n"))
    cat(sprintf("    Substitutions represent 'circulating strain differs from vaccine'.\n"))
    cat(sprintf("    Classification will use field-divergence frequency, not a vaccine-side test.\n"))
  } else {
    cat(sprintf("[i] Detected EXTERNALLY-REFERENCED run (vaccine_reference != target_accession).\n"))
    cat(sprintf("    Using vaccine-vs-circulating comparison (k_vac/Fisher's exact test).\n"))
  }
  
  if (!is.null(target_positions) && length(target_positions) > 0 && !all(is.na(target_positions))) {
    native_target_positions <- convert_positions(target_positions, target_position_system, coordinate_system)
    if (target_position_system != coordinate_system) {
      cat(sprintf("[i] target_positions supplied in %s coordinates; converted to %s coordinates for filtering: %s\n",
                  target_position_system, coordinate_system, paste(native_target_positions, collapse = delim)))
    } else {
      cat(sprintf("[+] Filtered target position(s): %s\n", paste(sort(unique(target_positions)), collapse = delim)))
    }
    muts_clean <- muts_clean[muts_clean$positions %in% native_target_positions, ]
  } else {
    cat("[+] Analyzing ALL available positions in dataset.\n")
  }
  
  vac_accessions  <- target_accession
  circ_accessions <- traits_clean$accession[tolower(traits_clean$strain_type) == "non-vaccine"]
  
  n_vac_total  <- length(vac_accessions)
  n_circ_total <- length(circ_accessions)
  n_total      <- n_vac_total + n_circ_total
  
  if (self_referenced && n_vac_total == 1) {
    cat("[!] NOTE: vaccine group n = 1 (single genome). Fisher's exact test reported\n")
    cat("    for transparency only and NOT used for classification in self-referenced mode.\n")
  }
  
  calc_wilson_ci <- function(k, n) {
    if (is.na(n) || n == 0) return(c(0, 0))
    test_res <- prop.test(k, n, conf.level = 0.95, correct = FALSE)
    return(test_res$conf.int[1:2])
  }
  
  unique_positions <- sort(unique(muts_clean$positions))
  if (length(unique_positions) == 0) {
    warning(sprintf("No mutations found matching criteria for %s on Segment %s.", vaccine_name, segment_label))
    return(NULL)
  }
  
  results_list <- list()
  
  for (pos in unique_positions) {
    pos_muts  <- muts_clean[muts_clean$positions == pos, ]
    snp_label <- pos_muts$snps[1]
    
    mapped_aa_pos <- map_to_amino_acid_position(pos, coordinate_system)
    prot_dom      <- get_protein_domain(mapped_aa_pos, segment_label)
    
    mut_accs      <- pos_muts$accession
    vac_mut_accs  <- intersect(mut_accs, vac_accessions)
    circ_mut_accs <- intersect(mut_accs, circ_accessions)
    
    k_vac  <- length(vac_mut_accs)
    k_circ <- length(circ_mut_accs)
    k_tot  <- k_vac + k_circ
    
    af_vac   <- if (n_vac_total > 0) (k_vac / n_vac_total) * 100 else 0
    af_circ  <- if (n_circ_total > 0) (k_circ / n_circ_total) * 100 else 0
    af_total <- if (n_total > 0) (k_tot / n_total) * 100 else 0
    
    af_vac_str   <- sprintf("%.2f%% (%d/%d)", af_vac, k_vac, n_vac_total)
    af_circ_str  <- sprintf("%.2f%% (%d/%d)", af_circ, k_circ, n_circ_total)
    af_total_str <- sprintf("%.2f%% (%d/%d)", af_total, k_tot, n_total)
    
    vac_ci  <- calc_wilson_ci(k_vac, n_vac_total)
    circ_ci <- calc_wilson_ci(k_circ, n_circ_total)
    tot_ci  <- calc_wilson_ci(k_tot, n_total)
    
    ci_formatted <- sprintf("Vac: [%.1f%%-%.1f%%] | Circ: [%.1f%%-%.1f%%] | Overall: [%.1f%%-%.1f%%]",
                            vac_ci[1]*100, vac_ci[2]*100,
                            circ_ci[1]*100, circ_ci[2]*100,
                            tot_ci[1]*100, tot_ci[2]*100)
    
    cont_table <- matrix(c(k_vac, n_vac_total - k_vac,
                           k_circ, n_circ_total - k_circ), nrow = 2, byrow = TRUE)
    fisher_res <- fisher.test(cont_table)
    
    specificity <- if (n_circ_total > 0) ((n_circ_total - k_circ) / n_circ_total) * 100 else 0
    ppv         <- if (k_tot > 0) (k_vac / k_tot) * 100 else 0
    
    circ_meta <- traits_clean[traits_clean$accession %in% circ_mut_accs, ]
    circ_muts <- pos_muts[pos_muts$accession %in% circ_accessions, ]
    
    geo_str <- format_counts(circ_meta$country)
    
    host_str <- format_exclusivity(circ_meta$host)
    temp_str <- format_exclusivity(circ_meta$year)
    lin_str  <- if ("lineage" %in% colnames(circ_muts)) format_exclusivity(circ_muts$lineage) else "N/A"
    
    host_exclusive     <- is_exclusive(circ_meta$host)
    temporal_exclusive <- is_exclusive(circ_meta$year)
    lineage_exclusive  <- if ("lineage" %in% colnames(circ_muts)) is_exclusive(circ_muts$lineage) else NA
    
    # ------------------------------------------------------------------------
    # v5: Bounded-language classification (addresses reviewer comment that
    # "absent from a finite reference dataset cannot automatically be
    # considered fixed or uniquely vaccine-associated").
    #
    # Zero (or near-zero) counts observed in a FINITE sample do not establish
    # that a residue is truly absent, fixed, or population-specific -- they
    # establish an upper bound on how common it could still be, consistent
    # with what was actually observed. We report that bound explicitly
    # (detection_limit_95pct_upper_bound, Wilson score interval -- for a
    # true zero count this closely matches the standard "rule of three"
    # approximation of ~3/n) rather than asserting certainty in either the
    # classification code or its human-readable label.
    #
    # `diagnostic_category` is a STABLE machine-readable code (used for
    # filtering/joins downstream, e.g. in build_manuscript_tables.R).
    # `potential_diagnostic_value` is now a DYNAMIC, per-row human-readable
    # description that always states the actual observed frequency and its
    # 95% CI upper bound, rather than an unqualified label like "vaccine-
    # specific" or "fixed".
    # ------------------------------------------------------------------------
    if (self_referenced) {
      # circ_ci is the Wilson 95% CI for k_circ/n_circ_total (the proportion
      # of circulating sequences DIFFERING from the vaccine). By the
      # symmetry of the Wilson interval under p -> 1-p, the CI for the
      # proportion of circulating sequences SHARING the vaccine's own
      # residue is exactly (1 - circ_ci[2], 1 - circ_ci[1]); we report the
      # upper bound of that complementary interval.
      detection_limit_95pct_upper_bound <- round((1 - circ_ci[1]) * 100, 2)
      vaccine_residue_field_pct <- 100 - af_circ  # % of field sharing vaccine's residue
      
      if (k_circ == n_circ_total && n_circ_total > 0) {
        diagnostic_category <- "vaccine_residue_not_detected"
        diagnostic_value <- sprintf(
          "Vaccine residue not detected in sampled circulating population (0/%d; 95%% CI upper bound = %.2f%%)",
          n_circ_total, detection_limit_95pct_upper_bound)
      } else if (af_circ >= field_marker_threshold) {
        diagnostic_category <- "vaccine_residue_low_frequency"
        diagnostic_value <- sprintf(
          "Vaccine residue at low frequency in sampled circulating population (%.2f%%; 95%% CI upper bound = %.2f%%)",
          vaccine_residue_field_pct, detection_limit_95pct_upper_bound)
      } else if (k_circ == 0) {
        diagnostic_category <- "invariant_site"
        diagnostic_value <- sprintf(
          "Invariant site (vaccine residue observed in 100%% of sampled circulating population, n=%d)",
          n_circ_total)
      } else {
        diagnostic_category <- "polymorphic_site"
        diagnostic_value <- sprintf(
          "Polymorphic site (vaccine residue observed in %.2f%% of sampled circulating population, %d/%d)",
          vaccine_residue_field_pct, n_circ_total - k_circ, n_circ_total)
      }
    } else {
      # circ_ci here is the Wilson 95% CI for k_circ/n_circ_total directly
      # (the proportion of circulating sequences sharing the vaccine's
      # substitution) -- no complement needed, since external-reference mode
      # already scores both groups against the same independent basis.
      detection_limit_95pct_upper_bound <- round(circ_ci[2] * 100, 2)
      
      if (k_vac > 0 && k_circ == 0) {
        diagnostic_category <- "no_shared_occurrence"
        diagnostic_value <- sprintf(
          "No shared occurrence in sampled circulating cohort (0/%d; 95%% CI upper bound = %.2f%%)",
          n_circ_total, detection_limit_95pct_upper_bound)
      } else if (k_vac > 0 && specificity >= 99.0 && fisher_res$p.value < 0.001) {
        diagnostic_category <- "low_field_frequency"
        diagnostic_value <- sprintf(
          "Low frequency in sampled circulating population (%.2f%%; 95%% CI upper bound = %.2f%%; Fisher's exact p = %.4g)",
          af_circ, detection_limit_95pct_upper_bound, fisher_res$p.value)
      } else if (k_vac > 0 && k_circ > 0) {
        diagnostic_category <- "shared_lineage_variant"
        diagnostic_value <- sprintf(
          "Shared with circulating lineage(s) (observed in %.2f%% of sampled circulating population, %d/%d)",
          af_circ, k_circ, n_circ_total)
      } else {
        diagnostic_category <- "field_variant"
        diagnostic_value <- "Not observed in vaccine accession (field variant only)"
      }
    }
    
    results_list[[length(results_list) + 1]] <- data.frame(
      vaccine_strain                      = vaccine_name,
      target_accession                    = target_accession,
      segment                             = segment_label,
      analysis_mode                       = if (self_referenced) "self_referenced" else "external_reference",
      native_coordinate_system            = coordinate_system,
      protein                             = prot_dom$protein,
      position                            = pos,
      mapped_amino_acid_position          = mapped_aa_pos,
      substitution                        = snp_label,
      functional_domain_structural_domain = prot_dom$domain,
      variant_frequency_vaccine           = af_vac_str,
      variant_frequency_circulating       = af_circ_str,
      variant_frequency_overall           = af_total_str,
      potential_diagnostic_value          = diagnostic_value,
      diagnostic_category                 = diagnostic_category,
      detection_limit_95pct_upper_bound   = detection_limit_95pct_upper_bound,
      number_of_sequences                = sprintf("Vac: %d/%d, Circ: %d/%d, Total: %d/%d", k_vac, n_vac_total, k_circ, n_circ_total, k_tot, n_total),
      geographic_distribution             = geo_str,
      host_distribution                   = host_str,
      host_exclusive                      = host_exclusive,
      temporal_distribution               = temp_str,
      temporal_exclusive                  = temporal_exclusive,
      lineage_distribution                = lin_str,
      lineage_exclusive                   = lineage_exclusive,
      confidence_intervals                = ci_formatted,
      fisher_p_value                      = fisher_res$p.value,
      fisher_test_caveat                  = if (self_referenced) "n_vac=1; not used for classification" else "",
      odds_ratio                          = unname(fisher_res$estimate),
      specificity_pct                     = round(specificity, 2),
      ppv_pct                             = round(ppv, 2),
      stringsAsFactors                    = FALSE
    )
  }
  
  final_table <- do.call(rbind, results_list)
  # FDR correction is ALWAYS computed on the full, unfiltered table -- across
  # every position that was actually tested -- before any diagnostic_filter
  # is applied. Filtering first and correcting after would understate the
  # true multiple-testing burden.
  final_table$fdr_adjusted_p_value <- p.adjust(final_table$fisher_p_value, method = "BH")
  
  # (v4/v5) optionally restrict what gets WRITTEN to disk to specific
  # classification(s), matched against the STABLE `diagnostic_category` code
  # -- not `potential_diagnostic_value`, which is now a dynamic, per-row
  # human-readable string (embeds actual counts/percentages) and therefore
  # cannot be matched exactly across rows. The R object returned by this
  # function is always the full table; only the CSV file is filtered.
  table_to_write <- final_table
  if (!is.null(diagnostic_filter)) {
    unknown_labels <- setdiff(diagnostic_filter, unique(final_table$diagnostic_category))
    if (length(unknown_labels) > 0) {
      warning(sprintf("diagnostic_filter contains code(s) not present in this table's classifications: %s",
                      paste(unknown_labels, collapse = ", ")))
    }
    table_to_write <- final_table[final_table$diagnostic_category %in% diagnostic_filter, ]
    cat(sprintf("[i] diagnostic_filter applied: keeping %d/%d position(s) matching: %s\n",
                nrow(table_to_write), nrow(final_table), paste(diagnostic_filter, collapse = ", ")))
  }
  
  setwd(outDir)
  if (!is.null(output_csv)) {
    write.csv(table_to_write, file = output_csv, row.names = FALSE)
    cat(sprintf("[+] Output saved to: %s (%d row(s))\n", output_csv, nrow(table_to_write)))
  }
  
  return(final_table)
}

# ------------------------------------------------------------------------------
# 4. Run configurations
# ------------------------------------------------------------------------------
run_configs <- tibble(
  vaccine_name = c(
    "MP-12", "Smithburn", "13",
    "MP-12", "Smithburn", "Clone 13",
    "MP-12", "Smithburn", "13",
    "MP-12", "Smithburn", "13"
  ),
  target_accession = c(
    "DQ380208", "DQ380193", "DQ380213",
    "DQ375404", "DQ375430", "DQ375417",
    "DQ380154", "DQ380157", "DQ380182",
    "DQ380154", "DQ380157", "DQ380182"
  ),
  segment = c(
    "M", "M", "M",
    "L", "L", "L",
    "S-NP", "S-NP", "S-NP",
    "S-NSS", "S-NSS", "S-NSS"
  ),
  traits_file = c(
    rep(sprintf("%s/output-dir/riftM/strain-types/riftM.traits.txt", segmentsDir[[2]]), 3),
    rep(sprintf("%s/output-dir/riftL/strain-types/riftL.traits.txt", segmentsDir[[1]]), 3),
    rep(sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.traits.txt", segmentsDir[[3]]), 3),
    rep(sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.traits.txt", segmentsDir[[3]]), 3)
  ),
  mutations_file = c(
    sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
    sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380193.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
    sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380213.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
    
    sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375404.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
    sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375430.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
    sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375417.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
    
    sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380154.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
    sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380157.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
    sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380182.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
    
    sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380154.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
    sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380157.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
    sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380182.mutations.per.strain.singleton.csv", segmentsDir[[3]])
  ),
  target_positions = list(
    NULL, NULL, NULL,
    NULL, NULL, NULL,
    NULL, NULL, NULL,
    NULL, NULL, NULL
  )
)

# 5. Iterate cleanly with pmap_dfr (primary, self-referenced analysis)
master_summary <- pmap_dfr(run_configs, function(vaccine_name, target_accession, segment, traits_file, mutations_file, target_positions) {
  analyze_vaccine_substitutions(
    traits_file            = traits_file,
    mutations_file          = mutations_file,
    target_accession        = target_accession,
    vaccine_name             = vaccine_name,
    segment_label            = segment,
    output_csv               = sprintf("RVFV_%s_%s_vaccine_substitutions_summary.csv", vaccine_name, segment),
    target_positions         = target_positions,
    target_position_system   = "amino_acid"
  )
})

# ------------------------------------------------------------------------------
# 6. COMPLEMENTARY run_configs: externally-referenced (nucleotide-level) scan
# ------------------------------------------------------------------------------
# Each M-segment vaccine is re-scored against the SAME independent RefSeq
# mutations file (NC_014396), since that reference is shared across vaccines
# -- unlike the self-referenced files, which each vaccine has its own copy of.
# traits_file is the same riftM.traits.txt as the self-referenced run (per
# your note that this file is common to both subworkflows). target_positions
# left NULL here to do a full scan; set them (in amino-acid coordinates) if
# you only want specific known candidate sites re-validated -- the function
# converts them to nucleotide coordinates automatically either way.
#
# For L and S segments, substitute the equivalent NC_014397 / NC_014395
# mutations files once available; the row structure is identical.

# run_configs_external <- tibble(
#   vaccine_name = c("MP-12", "Smithburn", "13"),                 # M segment vaccines
#   target_accession = c("DQ380208", "DQ380193", "DQ380213"),
#   segment = c("M", "M", "M"),
#   traits_file = rep(
#     sprintf("%s/output-dir/riftM/strain-types/riftM.traits.txt", segmentsDir[[2]]), 3
#   ),
#   mutations_file = rep(
#     sprintf("%s/output-dir/consensus-reference/riftM/strain-types/riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
#     3
#   ),
#   target_positions = list(NULL, NULL, NULL)   # full scan; or e.g. list(c(232,259), c(219,230,249), NULL)
# )
# 
# # Run it, using diagnostic_filter to keep only the informative rows in the
# # output CSV (a full scan can otherwise run to 1000+ rows per vaccine, e.g.
# # 1104 positions were tested for Smithburn/NC_014396 in a full M-segment scan).
# master_summary_external <- pmap_dfr(run_configs_external, function(vaccine_name, target_accession, segment, traits_file, mutations_file, target_positions) {
#   analyze_vaccine_substitutions(
#     traits_file            = traits_file,
#     mutations_file          = mutations_file,
#     target_accession        = target_accession,
#     vaccine_name             = vaccine_name,
#     segment_label            = segment,
#     output_csv               = sprintf("RVFV_%s_%s_external_reference_summary.csv", vaccine_name, segment),
#     target_positions         = target_positions,
#     target_position_system   = "amino_acid"
#     # diagnostic_filter        = c("no_shared_occurrence", "low_field_frequency")
#   )
# })

# Combine both schemes into one supplementary table if useful downstream
# (master_summary_external here is only the ROW OBJECT returned in R, which
# is always the FULL unfiltered table regardless of diagnostic_filter --
# the CSVs on disk are what got filtered):
# combined_summary <- bind_rows(master_summary, master_summary_external)

# ------------------------------------------------------------------------------
# 4. Run configurations
# ------------------------------------------------------------------------------
# run_configs <- tibble(
#   vaccine_name = c(
#     "MP-12", "Smithburn", "13",
#     "MP-12", "Smithburn", "Clone 13",
#     "MP-12", "Smithburn", "13",
#     "MP-12", "Smithburn", "13"
#   ),
#   target_accession = c(
#     "DQ380208", "DQ380193", "DQ380213",
#     "DQ375404", "DQ375430", "DQ375417",
#     "DQ380154", "DQ380157", "DQ380182",
#     "DQ380154", "DQ380157", "DQ380182"
#   ),
#   segment = c(
#     "M", "M", "M",
#     "L", "L", "L",
#     "S-NP", "S-NP", "S-NP",
#     "S-NSS", "S-NSS", "S-NSS"
#   ),
#   traits_file = c(
#     rep(sprintf("%s/output-dir/riftM/strain-types/riftM.traits.txt", segmentsDir[[2]]), 3),
#     rep(sprintf("%s/output-dir/riftL/strain-types/riftL.traits.txt", segmentsDir[[1]]), 3),
#     rep(sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.traits.txt", segmentsDir[[3]]), 3),
#     rep(sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.traits.txt", segmentsDir[[3]]), 3)
#   ),
#   mutations_file = c(
#     sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
#     sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380193.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
#     sprintf("%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380213.mutations.per.strain.singleton.csv", segmentsDir[[2]]),
#     
#     sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375404.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
#     sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375430.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
#     sprintf("%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375417.mutations.per.strain.singleton.csv", segmentsDir[[1]]),
#     
#     sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380154.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
#     sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380157.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
#     sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380182.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
#     
#     sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380154.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
#     sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380157.mutations.per.strain.singleton.csv", segmentsDir[[3]]),
#     sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380182.mutations.per.strain.singleton.csv", segmentsDir[[3]])
#   ),
#   target_positions = list(
#     NULL, NULL, NULL,
#     NULL, NULL, NULL,
#     NULL, NULL, NULL,
#     NULL, NULL, NULL
#   )
# )
# 
# # 5. Iterate cleanly with pmap_dfr
# master_summary <- pmap_dfr(run_configs, function(vaccine_name, target_accession, segment, traits_file, mutations_file, target_positions) {
#   analyze_vaccine_substitutions(
#     traits_file            = traits_file,
#     mutations_file          = mutations_file,
#     target_accession        = target_accession,
#     vaccine_name             = vaccine_name,
#     segment_label            = segment,
#     output_csv               = sprintf("RVFV_%s_%s_vaccine_substitutions_summary.csv", vaccine_name, segment),
#     target_positions         = target_positions,
#     target_position_system   = "amino_acid"
#     # diagnostic_filter = c("High-Confidence Vaccine-Defining Residue",
#     #                       "Vaccine-Specific Residue (absent in ALL sampled field strains)")
#     # diagnostic_filter = c("High-Confidence Vaccine-Defining Residue",
#     #                       "Vaccine-Specific Residue (absent in ALL sampled field strains)",
#     #                       "Invariant Site (vaccine residue universal in field)",
#     #                       "Polymorphic Site (vaccine residue shared by some field strains)"
#     #                       )
#   )
# })



run_configs_external <- tibble(
  vaccine_name = c(
    "MP-12", "Smithburn", "13",
    "MP-12", "Smithburn", "Clone 13",
    "MP-12", "Smithburn", "13",
    "MP-12", "Smithburn", "13"
  ),
  target_accession = c(
    "DQ380208", "DQ380193", "DQ380213",
    "DQ375404", "DQ375430", "DQ375417",
    "DQ380154", "DQ380157", "DQ380182",
    "DQ380154", "DQ380157", "DQ380182"
  ),
  segment = c(
    "M", "M", "M",
    "L", "L", "L",
    "S-NP", "S-NP", "S-NP",
    "S-NSS", "S-NSS", "S-NSS"
  ),
  traits_file = c(
    rep(sprintf("%s/output-dir/riftM/strain-types/riftM.traits.txt", segmentsDir[[2]]), 3),
    rep(sprintf("%s/output-dir/riftL/strain-types/riftL.traits.txt", segmentsDir[[1]]), 3),
    rep(sprintf("%s/output-dir/riftS-NP/strain-types/riftS-NP.traits.txt", segmentsDir[[3]]), 3),
    rep(sprintf("%s/output-dir/riftS-NSS/strain-types/riftS-NSS.traits.txt", segmentsDir[[3]]), 3)
  ),
  mutations_file = c(
    rep(sprintf("%s/output-dir/mutational-profiling/riftM/strain-types/riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv", segmentsDir[[2]]), 3),
    rep(sprintf("%s/output-dir/mutational-profiling/riftL/strain-types/riftL.strain_type.NC_014397.mutations.per.strain.singleton.csv", segmentsDir[[1]]), 3),
    rep(sprintf("%s/output-dir/mutational-profiling/riftS-NP/strain-types-rev/riftS-NP.strain_type.NC_014395.mutations.per.strain.singleton.csv", segmentsDir[[3]]), 3),
    rep(sprintf("%s/output-dir/mutational-profiling/riftS-NSS/strain-types/riftS-NSS.strain_type.NC_014395.mutations.per.strain.singleton.csv", segmentsDir[[3]]), 3)
  ),
  target_positions = list(
    NULL, NULL, NULL,
    NULL, NULL, NULL,
    NULL, NULL, NULL,
    NULL, NULL, NULL
  )
)

master_summary_external <- pmap_dfr(run_configs_external, function(vaccine_name, target_accession, segment, traits_file, mutations_file, target_positions) {
  analyze_vaccine_substitutions(
    traits_file            = traits_file,
    mutations_file          = mutations_file,
    target_accession        = target_accession,
    vaccine_name             = vaccine_name,
    segment_label            = segment,
    output_csv               = sprintf("RVFV_%s_%s_vaccine_nucleotides_summary.csv", vaccine_name, segment),
    target_positions         = target_positions,
    target_position_system   = "amino_acid"
    # diagnostic_filter = c("Strict Target Vaccine Marker",
    #                       "High-Confidence Vaccine-Defining Residue")
    # diagnostic_filter = c("Strict Target Vaccine Marker",
    #                       "High-Confidence Vaccine-Defining Residue",
    #                       "Non-Specific / Shared Lineage Variant",
    #                       "Wild-Type / Field Variant")
  )
})


