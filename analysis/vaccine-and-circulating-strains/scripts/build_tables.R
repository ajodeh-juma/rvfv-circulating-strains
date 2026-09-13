# # ==============================================================================
# # Script Name: build_manuscript_tables.R
# # Author: Viral Genomics & Bioinformatics
# # Description: Builds the MAIN-TEXT table and SUPPLEMENTARY workbook from the
# #              per-vaccine/per-segment CSVs produced by
# #              analyze_vaccine_substitutions_multisegment_v3.R.
# #
# #              Works generically across M, L, and S (S-NP/S-NSS) segments and
# #              any number of vaccines, as long as every relevant CSV (both the
# #              self-referenced/amino-acid runs AND the externally-referenced/
# #              nucleotide runs) lives in one input directory. Routing between
# #              amino-acid and nucleotide results is done using the
# #              `native_coordinate_system` column already written by v3 -- NOT
# #              by filename -- so this script is robust to whatever naming
# #              convention was used for output_csv.
# #
# # IMPORTANT PRE-REQUISITE: every externally-referenced (nucleotide) CSV must
# # have been generated with target_accession SET TO THE VACCINE ACCESSION
# # (e.g. "DQ380193" for Smithburn), NOT the RefSeq reference genome
# # (e.g. "NC_014396"). If target_accession == the RefSeq accession, that run
# # is trivially self-referenced against a NON-VACCINE strain and its rows will
# # be silently wrong. This script does not detect that mistake (it looks
# # correct at the column level), so re-check target_accession before running.
# #
# # OUTPUTS:
# #   1. main_text_vaccine_markers.xlsx (or .csv if openxlsx unavailable)
# #      One row per confirmed marker (amino-acid level preferred; nucleotide-
# #      only markers, e.g. synonymous substitutions, included and flagged).
# #      Compact columns suitable for direct inclusion in the manuscript.
# #
# #   2. supplementary_tables.xlsx (or two .csv files if openxlsx unavailable)
# #      Sheet 1 "Full_Frequency_Table": every position tested, both coordinate
# #        systems, all columns -- addresses the reviewer's request for a
# #        complete frequency table across all circulating and vaccine
# #        sequences.
# #      Sheet 2 "AA_vs_NT_Cross_Validation": per-codon concordance between the
# #        self-referenced (amino-acid) and externally-referenced (nucleotide)
# #        analyses, showing which markers are confirmed by both methods,
# #        which are nucleotide-only (e.g. synonymous), and which are
# #        amino-acid-only (nucleotide route not fully private).
# # ==============================================================================
# source(sprintf("%s/utils.R", "scripts"))
# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(purrr)
#   library(readr)
#   library(stringr)
#   library(tibble)
# })
# 
# # Optional formatting package -- script still works (writes CSV) if absent.
# .has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
# if (.has_openxlsx) {
#   library(openxlsx)
# } else {
#   message("[!] Package 'openxlsx' not found -- falling back to plain CSV output.")
#   message("    Install with install.packages('openxlsx') for formatted .xlsx workbooks.")
# }
# 
# # ------------------------------------------------------------------------------
# # Classification vocabularies (must match analyze_vaccine_substitutions_multisegment_v3.R)
# # ------------------------------------------------------------------------------
# SELF_REFERENCED_MARKER_CLASSES <- c(
#   "Vaccine-Specific Residue (absent in ALL sampled field strains)",
#   "High-Confidence Vaccine-Defining Residue"
# )
# EXTERNAL_REFERENCE_MARKER_CLASSES <- c(
#   "Strict Target Vaccine Marker",
#   "High-Confidence Vaccine Marker"
# )
# 
# SEGMENT_ORDER <- c("M", "L", "S-NP", "S-NSS")
# 
# ALL_MARKER_CLASSES <- c(SELF_REFERENCED_MARKER_CLASSES, EXTERNAL_REFERENCE_MARKER_CLASSES)
# 
# # ------------------------------------------------------------------------------
# # NEW: Cross-vaccine specificity annotation
# # ------------------------------------------------------------------------------
# # A substitution is "vaccine-specific" here if it is flagged as a marker
# # (in EITHER coordinate system's marker classes) for exactly ONE vaccine
# # strain within a given segment. If the SAME position is independently
# # flagged as a marker for two or more vaccine strains (e.g. both MP-12 and
# # Smithburn), it is "shared" -- it can still indicate vaccine origin, but
# # cannot distinguish WHICH vaccine strain is present.
# #
# # Grouping key is (segment, native_coordinate_system, position) using the
# # NATIVE position, not mapped_amino_acid_position. This matters for
# # nucleotide-level data: two different nucleotide positions can map to the
# # same codon (mapped_amino_acid_position) without being the same
# # substitution at all -- using the native position avoids incorrectly
# # treating those as "shared" just because they fall in the same codon.
# add_vaccine_specificity_flags <- function(combined) {
#   marker_rows <- combined %>%
#     filter(potential_diagnostic_value %in% ALL_MARKER_CLASSES) %>%
#     distinct(segment, native_coordinate_system, position, vaccine_strain)
#   
#   sharing_counts <- marker_rows %>%
#     group_by(segment, native_coordinate_system, position) %>%
#     summarise(
#       n_vaccines_with_marker_here = n_distinct(vaccine_strain),
#       vaccines_sharing_this_position = paste(sort(unique(vaccine_strain)), collapse = "; "),
#       .groups = "drop"
#     )
#   
#   combined %>%
#     left_join(sharing_counts, by = c("segment", "native_coordinate_system", "position")) %>%
#     mutate(
#       n_vaccines_with_marker_here = ifelse(is.na(n_vaccines_with_marker_here), 0L, n_vaccines_with_marker_here),
#       vaccine_specific_marker = potential_diagnostic_value %in% ALL_MARKER_CLASSES &
#         n_vaccines_with_marker_here == 1
#     )
# }
# 
# # ------------------------------------------------------------------------------
# # 1. Load every summary CSV in a directory and combine
# # ------------------------------------------------------------------------------
# load_all_results <- function(input_dir, pattern = "_vaccine_(substitutions|nucleotides)_summary\\.csv$") {
#   files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
#   if (length(files) == 0) stop(sprintf("No matching CSV files found in %s (pattern: %s)", input_dir, pattern))
#   
#   cat(sprintf("[i] Found %d CSV file(s) in %s:\n", length(files), input_dir))
#   for (f in files) cat(sprintf("    - %s\n", basename(f)))
#   
#   combined <- map_dfr(files, function(f) {
#     df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
#     
#     # Guard: a CSV with a header but zero data rows (e.g. a vaccine/segment
#     # run that found no qualifying positions) has the right columns but
#     # nrow(df) == 0. Assigning a length-1 value to a new column of a 0-row
#     # data frame trips strict recycling rules in some tidyverse versions
#     # ("replacement has 1 row, data has 0"), so skip such files explicitly
#     # rather than letting that assignment fail.
#     if (nrow(df) == 0) {
#       warning(sprintf("Skipping %s -- file has a header but zero data rows.", basename(f)))
#       return(NULL)
#     }
#     
#     # Guard: read.csv infers column type from the DATA, ignoring the fact
#     # that a field was quoted in the source CSV. One vaccine in this project
#     # is literally named "13" (target_accession DQ380213) -- in a file where
#     # every vaccine_strain value is "13", read.csv silently returns an
#     # INTEGER column, while every other file (e.g. "MP-12", "Smithburn")
#     # correctly reads as character. dplyr::bind_rows() then fails with
#     # "Can't combine ..$vaccine_strain <integer> and ..$vaccine_strain
#     # <character>". Force all known text columns to character explicitly so
#     # this can't happen regardless of which values happen to look numeric.
#     text_cols <- c("vaccine_strain", "target_accession", "segment", "analysis_mode",
#                    "native_coordinate_system", "protein", "substitution",
#                    "functional_domain_structural_domain", "potential_diagnostic_value",
#                    "number_of_sequences", "geographic_distribution", "host_distribution",
#                    "temporal_distribution", "lineage_distribution", "confidence_intervals",
#                    "fisher_test_caveat", "variant_frequency_vaccine",
#                    "variant_frequency_circulating", "variant_frequency_overall")
#     for (col in intersect(text_cols, colnames(df))) {
#       df[[col]] <- as.character(df[[col]])
#     }
#     
#     # Same class of bug can hit numeric columns (e.g. a column that's all NA
#     # in one file reads as logical there but numeric elsewhere). Force known
#     # numeric columns explicitly too.
#     numeric_cols <- c("position", "mapped_amino_acid_position", "fisher_p_value",
#                       "odds_ratio", "specificity_pct", "ppv_pct", "fdr_adjusted_p_value")
#     for (col in intersect(numeric_cols, colnames(df))) {
#       df[[col]] <- as.numeric(df[[col]])
#     }
#     
#     # Logical columns: same reasoning -- ensure TRUE/FALSE/NA read consistently
#     # as logical, not character, across files.
#     logical_cols <- c("host_exclusive", "temporal_exclusive", "lineage_exclusive")
#     for (col in intersect(logical_cols, colnames(df))) {
#       df[[col]] <- as.logical(df[[col]])
#     }
#     
#     required_cols <- c("vaccine_strain", "segment", "native_coordinate_system",
#                        "mapped_amino_acid_position", "potential_diagnostic_value")
#     missing <- setdiff(required_cols, colnames(df))
#     if (length(missing) > 0) {
#       warning(sprintf("Skipping %s -- missing expected column(s): %s",
#                       basename(f), paste(missing, collapse = ", ")))
#       return(NULL)
#     }
#     df$source_file <- basename(f)
#     df
#   })
#   
#   if (is.null(combined) || nrow(combined) == 0) {
#     stop("No usable rows were loaded -- every matching file was empty or missing required columns. See warnings above.")
#   }
#   
#   # Sanity check: warn if any externally-referenced run has target_accession
#   # equal to a value that never appears as "vaccine" strain_type-derived --
#   # we can't check the traits file from here, but we CAN warn if an
#   # "external_reference" row's target_accession looks like a RefSeq ID
#   # (common NCBI RefSeq prefixes), which is the exact mistake this script
#   # is documented to guard against.
#   suspect <- combined %>%
#     filter(analysis_mode == "external_reference") %>%
#     filter(str_detect(target_accession, "^NC_")) %>%
#     distinct(vaccine_strain, segment, target_accession)
#   if (nrow(suspect) > 0) {
#     cat("\n[!] WARNING: the following externally-referenced rows have a RefSeq-\n")
#     cat("    looking target_accession (starts with 'NC_'). This is very likely\n")
#     cat("    the vaccine_reference genome, NOT the vaccine accession itself --\n")
#     cat("    double check these before trusting the results:\n")
#     print(suspect)
#     cat("\n")
#   }
#   
#   cat(sprintf("\n[+] Combined table: %d rows across %d source file(s)\n",
#               nrow(combined), length(files)))
#   combined
# }
# 
# # ------------------------------------------------------------------------------
# # 2. Cross-validation table: amino-acid vs nucleotide concordance per codon
# # ------------------------------------------------------------------------------
# build_cross_validation_table <- function(combined) {
#   
#   is_aa_marker <- function(vals) any(vals %in% SELF_REFERENCED_MARKER_CLASSES)
#   is_nt_marker <- function(vals) any(vals %in% EXTERNAL_REFERENCE_MARKER_CLASSES)
#   
#   cross <- combined %>%
#     group_by(vaccine_strain, segment, mapped_amino_acid_position) %>%
#     summarise(
#       aa_marker       = is_aa_marker(potential_diagnostic_value[native_coordinate_system == "amino_acid"]),
#       nt_marker       = is_nt_marker(potential_diagnostic_value[native_coordinate_system == "nucleotide"]),
#       aa_substitution = paste(unique(substitution[native_coordinate_system == "amino_acid"]), collapse = "; "),
#       nt_substitutions = paste(unique(substitution[native_coordinate_system == "nucleotide" &
#                                                      potential_diagnostic_value %in% EXTERNAL_REFERENCE_MARKER_CLASSES]), collapse = "; "),
#       aa_classification = paste(unique(potential_diagnostic_value[native_coordinate_system == "amino_acid"]), collapse = "; "),
#       n_nucleotide_markers_in_codon = sum(native_coordinate_system == "nucleotide" &
#                                             potential_diagnostic_value %in% EXTERNAL_REFERENCE_MARKER_CLASSES),
#       .groups = "drop"
#     ) %>%
#     filter(aa_marker | nt_marker) %>%
#     mutate(
#       concordance_status = case_when(
#         aa_marker & nt_marker  ~ "Concordant (marker at both amino-acid and nucleotide levels)",
#         aa_marker & !nt_marker ~ "Amino-acid marker only (nucleotide route not fully private in field)",
#         !aa_marker & nt_marker ~ "Nucleotide marker only (e.g. synonymous substitution -- invisible at amino-acid level)",
#         TRUE ~ NA_character_
#       ),
#       segment = factor(segment, levels = SEGMENT_ORDER)
#     ) %>%
#     arrange(segment, vaccine_strain, mapped_amino_acid_position) %>%
#     mutate(segment = as.character(segment))
#   
#   cross
# }
# 
# # ------------------------------------------------------------------------------
# # 3. Main-text table: compact, one row per confirmed marker
# # ------------------------------------------------------------------------------
# build_main_text_table <- function(combined, cross_validation, require_vaccine_specific = TRUE) {
#   
#   # Prefer the amino-acid-level row (cleaner notation, e.g. "I219T") whenever
#   # a marker is confirmed at the amino-acid level; fall back to the
#   # nucleotide-level row (e.g. synonymous markers, only detectable there).
#   aa_rows <- combined %>%
#     filter(native_coordinate_system == "amino_acid",
#            potential_diagnostic_value %in% SELF_REFERENCED_MARKER_CLASSES)
#   
#   nt_only_codons <- cross_validation %>%
#     filter(concordance_status == "Nucleotide marker only (e.g. synonymous substitution -- invisible at amino-acid level)") %>%
#     select(vaccine_strain, segment, mapped_amino_acid_position)
#   
#   nt_rows <- combined %>%
#     filter(native_coordinate_system == "nucleotide",
#            potential_diagnostic_value %in% EXTERNAL_REFERENCE_MARKER_CLASSES) %>%
#     inner_join(nt_only_codons, by = c("vaccine_strain", "segment", "mapped_amino_acid_position"))
#   
#   main <- bind_rows(
#     aa_rows %>% mutate(marker_source = "amino_acid"),
#     nt_rows %>% mutate(marker_source = "nucleotide_only")
#   ) %>%
#     left_join(
#       cross_validation %>% select(vaccine_strain, segment, mapped_amino_acid_position, concordance_status),
#       by = c("vaccine_strain", "segment", "mapped_amino_acid_position")
#     )
#   
#   n_before <- nrow(main)
#   
#   if (require_vaccine_specific) {
#     main <- main %>% filter(vaccine_specific_marker)
#     cat(sprintf("[i] Main-text table narrowed to vaccine-specific (non-shared) substitutions: %d -> %d rows\n",
#                 n_before, nrow(main)))
#     shared_removed <- n_before - nrow(main)
#     if (shared_removed > 0) {
#       cat(sprintf("    (%d marker(s) were shared across >=2 vaccine strains and excluded from the main-text\n", shared_removed))
#       cat(sprintf("     table; they remain visible, flagged, in the supplementary full frequency table.)\n"))
#     }
#   }
#   
#   main %>%
#     mutate(
#       confirmed_by_independent_reference = concordance_status ==
#         "Concordant (marker at both amino-acid and nucleotide levels)",
#       segment = factor(segment, levels = SEGMENT_ORDER)
#     ) %>%
#     arrange(segment, vaccine_strain, mapped_amino_acid_position) %>%
#     mutate(segment = as.character(segment)) %>%
#     select(
#       vaccine_strain,
#       segment,
#       protein,
#       mapped_amino_acid_position,
#       substitution,
#       functional_domain_structural_domain,
#       variant_frequency_circulating,
#       potential_diagnostic_value,
#       vaccine_specific_marker,
#       n_vaccines_with_marker_here,
#       confirmed_by_independent_reference,
#       fisher_p_value,
#       specificity_pct
#     ) %>%
#     rename(
#       amino_acid_position = mapped_amino_acid_position,
#       functional_domain    = functional_domain_structural_domain,
#       field_frequency      = variant_frequency_circulating,
#       diagnostic_value     = potential_diagnostic_value
#     )
# }
# 
# # ------------------------------------------------------------------------------
# # 4. Supplementary full frequency table: everything, tidied and sorted
# # ------------------------------------------------------------------------------
# build_supplementary_full_table <- function(combined) {
#   combined %>%
#     mutate(segment = factor(segment, levels = SEGMENT_ORDER)) %>%
#     arrange(segment, vaccine_strain, native_coordinate_system, position) %>%
#     mutate(segment = as.character(segment)) %>%
#     select(-source_file)
# }
# 
# # ------------------------------------------------------------------------------
# # 5. Formatted output writers
# # ------------------------------------------------------------------------------
# .write_formatted_sheet <- function(wb, sheet_name, df, highlight_col = NULL, highlight_values = NULL) {
#   addWorksheet(wb, sheet_name)
#   writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold",
#                                                           fgFill = "#D9E1F2",
#                                                           border = "Bottom"))
#   freezePane(wb, sheet_name, firstRow = TRUE)
#   setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
#   
#   if (!is.null(highlight_col) && highlight_col %in% colnames(df)) {
#     col_idx <- which(colnames(df) == highlight_col)
#     hl_style <- createStyle(fgFill = "#C6EFCE")
#     match_rows <- which(df[[highlight_col]] %in% highlight_values) + 1  # +1 for header row
#     if (length(match_rows) > 0) {
#       addStyle(wb, sheet_name, style = hl_style, rows = match_rows, cols = col_idx, gridExpand = FALSE)
#     }
#   }
# }
# 
# write_main_text_table <- function(main_table, output_path) {
#   if (.has_openxlsx) {
#     wb <- createWorkbook()
#     .write_formatted_sheet(wb, "Main_Text_Markers", main_table,
#                            highlight_col = "confirmed_by_independent_reference",
#                            highlight_values = TRUE)
#     saveWorkbook(wb, output_path, overwrite = TRUE)
#     cat(sprintf("[+] Main-text table saved to: %s\n", output_path))
#   } else {
#     csv_path <- sub("\\.xlsx$", ".csv", output_path)
#     write.csv(main_table, csv_path, row.names = FALSE)
#     cat(sprintf("[+] Main-text table saved to: %s (CSV fallback)\n", csv_path))
#   }
# }
# 
# write_supplementary_tables <- function(full_table, cross_validation, output_path) {
#   if (.has_openxlsx) {
#     wb <- createWorkbook()
#     .write_formatted_sheet(wb, "Full_Frequency_Table", full_table,
#                            highlight_col = "potential_diagnostic_value",
#                            highlight_values = c(SELF_REFERENCED_MARKER_CLASSES, EXTERNAL_REFERENCE_MARKER_CLASSES))
#     .write_formatted_sheet(wb, "AA_vs_NT_Cross_Validation", cross_validation,
#                            highlight_col = "concordance_status",
#                            highlight_values = "Concordant (marker at both amino-acid and nucleotide levels)")
#     saveWorkbook(wb, output_path, overwrite = TRUE)
#     cat(sprintf("[+] Supplementary workbook saved to: %s\n", output_path))
#   } else {
#     base_path <- sub("\\.xlsx$", "", output_path)
#     write.csv(full_table, paste0(base_path, "_full_frequency_table.csv"), row.names = FALSE)
#     write.csv(cross_validation, paste0(base_path, "_cross_validation.csv"), row.names = FALSE)
#     cat(sprintf("[+] Supplementary tables saved to: %s_full_frequency_table.csv and %s_cross_validation.csv (CSV fallback)\n",
#                 base_path, base_path))
#   }
# }
# 
# # ------------------------------------------------------------------------------
# # 6. Top-level driver
# # ------------------------------------------------------------------------------
# build_manuscript_tables <- function(
#     input_dir,
#     output_dir = ".",
#     main_text_filename = "main_text_vaccine_markers.xlsx",
#     supplementary_filename = "supplementary_tables.xlsx",
#     require_vaccine_specific = TRUE
# ) {
#   combined <- load_all_results(input_dir)
#   combined <- add_vaccine_specificity_flags(combined)
#   
#   n_shared_positions <- combined %>%
#     filter(potential_diagnostic_value %in% ALL_MARKER_CLASSES, n_vaccines_with_marker_here > 1) %>%
#     distinct(segment, native_coordinate_system, position) %>%
#     nrow()
#   if (n_shared_positions > 0) {
#     cat(sprintf("\n[i] %d position(s) are markers for >=2 vaccine strains (shared, not strain-specific):\n", n_shared_positions))
#     print(
#       combined %>%
#         filter(potential_diagnostic_value %in% ALL_MARKER_CLASSES, n_vaccines_with_marker_here > 1) %>%
#         distinct(segment, native_coordinate_system, position, vaccines_sharing_this_position)
#     )
#     cat("\n")
#   }
#   
#   cross_validation <- build_cross_validation_table(combined)
#   main_table        <- build_main_text_table(combined, cross_validation, require_vaccine_specific = require_vaccine_specific)
#   full_table        <- build_supplementary_full_table(combined)
#   
#   cat(sprintf("\n[i] Main-text table: %d confirmed marker(s) across %d segment(s)\n",
#               nrow(main_table), length(unique(main_table$segment))))
#   cat(sprintf("[i] Supplementary full frequency table: %d row(s)\n", nrow(full_table)))
#   cat(sprintf("[i] Cross-validation table: %d codon(s) flagged as a marker at either level\n",
#               nrow(cross_validation)))
#   cat("[i] Cross-validation breakdown:\n")
#   print(table(cross_validation$concordance_status))
#   
#   write_main_text_table(main_table, file.path(output_dir, main_text_filename))
#   write_supplementary_tables(full_table, cross_validation, file.path(output_dir, supplementary_filename))
#   
#   invisible(list(main_table = main_table, full_table = full_table, cross_validation = cross_validation))
# }

# ==============================================================================
# Script Name: build_manuscript_tables.R
# Author: Viral Genomics & Bioinformatics
# Description: Builds the MAIN-TEXT table and SUPPLEMENTARY workbook from the
#              per-vaccine/per-segment CSVs produced by
#              analyze_vaccine_substitutions_multisegment_v3.R.
#
#              Works generically across M, L, and S (S-NP/S-NSS) segments and
#              any number of vaccines, as long as every relevant CSV (both the
#              self-referenced/amino-acid runs AND the externally-referenced/
#              nucleotide runs) lives in one input directory. Routing between
#              amino-acid and nucleotide results is done using the
#              `native_coordinate_system` column already written by v3 -- NOT
#              by filename -- so this script is robust to whatever naming
#              convention was used for output_csv.
#
# IMPORTANT PRE-REQUISITE: every externally-referenced (nucleotide) CSV must
# have been generated with target_accession SET TO THE VACCINE ACCESSION
# (e.g. "DQ380193" for Smithburn), NOT the RefSeq reference genome
# (e.g. "NC_014396"). If target_accession == the RefSeq accession, that run
# is trivially self-referenced against a NON-VACCINE strain and its rows will
# be silently wrong. This script does not detect that mistake (it looks
# correct at the column level), so re-check target_accession before running.
#
# OUTPUTS:
#   1. main_text_vaccine_markers.xlsx (or .csv if openxlsx unavailable)
#      One row per confirmed marker (amino-acid level preferred; nucleotide-
#      only markers, e.g. synonymous substitutions, included and flagged).
#      Compact columns suitable for direct inclusion in the manuscript.
#
#   2. supplementary_tables.xlsx (or two .csv files if openxlsx unavailable)
#      Sheet 1 "Full_Frequency_Table": every position tested, both coordinate
#        systems, all columns -- addresses the reviewer's request for a
#        complete frequency table across all circulating and vaccine
#        sequences.
#      Sheet 2 "AA_vs_NT_Cross_Validation": per-codon concordance between the
#        self-referenced (amino-acid) and externally-referenced (nucleotide)
#        analyses, showing which markers are confirmed by both methods,
#        which are nucleotide-only (e.g. synonymous), and which are
#        amino-acid-only (nucleotide route not fully private).
# ==============================================================================
source(sprintf("%s/utils.R", "scripts"))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

# Optional formatting package -- script still works (writes CSV) if absent.
.has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
if (.has_openxlsx) {
  library(openxlsx)
} else {
  message("[!] Package 'openxlsx' not found -- falling back to plain CSV output.")
  message("    Install with install.packages('openxlsx') for formatted .xlsx workbooks.")
}

# ------------------------------------------------------------------------------
# Classification vocabularies (must match analyze_vaccine_substitutions_multisegment_v3.R)
#
# v5 NOTE: these now match on `diagnostic_category` (a stable machine-readable
# code) rather than `potential_diagnostic_value` (now a dynamic, per-row
# human-readable string that embeds actual counts/percentages, e.g. "Vaccine
# residue not detected in sampled circulating population (0/240; 95% CI
# upper bound = 1.58%)"). That text is unique per row and can no longer be
# matched exactly across rows -- diagnostic_category is what downstream
# logic (this script) filters and joins on.
# ------------------------------------------------------------------------------
SELF_REFERENCED_MARKER_CATEGORIES <- c(
  "vaccine_residue_not_detected",
  "vaccine_residue_low_frequency"
)
EXTERNAL_REFERENCE_MARKER_CATEGORIES <- c(
  "no_shared_occurrence",
  "low_field_frequency"
)

SEGMENT_ORDER <- c("M", "L", "S-NP", "S-NSS")

ALL_MARKER_CATEGORIES <- c(SELF_REFERENCED_MARKER_CATEGORIES, EXTERNAL_REFERENCE_MARKER_CATEGORIES)

# ------------------------------------------------------------------------------
# NEW: Cross-vaccine specificity annotation
# ------------------------------------------------------------------------------
# A substitution is "vaccine-specific" here if it is flagged as a marker
# (in EITHER coordinate system's marker classes) for exactly ONE vaccine
# strain within a given segment. If the SAME position is independently
# flagged as a marker for two or more vaccine strains (e.g. both MP-12 and
# Smithburn), it is "shared" -- it can still indicate vaccine origin, but
# cannot distinguish WHICH vaccine strain is present.
#
# Grouping key is (segment, native_coordinate_system, position) using the
# NATIVE position, not mapped_amino_acid_position. This matters for
# nucleotide-level data: two different nucleotide positions can map to the
# same codon (mapped_amino_acid_position) without being the same
# substitution at all -- using the native position avoids incorrectly
# treating those as "shared" just because they fall in the same codon.
add_vaccine_specificity_flags <- function(combined) {
  marker_rows <- combined %>%
    filter(diagnostic_category %in% ALL_MARKER_CATEGORIES) %>%
    distinct(segment, native_coordinate_system, position, vaccine_strain)
  
  sharing_counts <- marker_rows %>%
    group_by(segment, native_coordinate_system, position) %>%
    summarise(
      n_vaccines_with_marker_here = n_distinct(vaccine_strain),
      vaccines_sharing_this_position = paste(sort(unique(vaccine_strain)), collapse = "; "),
      .groups = "drop"
    )
  
  combined %>%
    left_join(sharing_counts, by = c("segment", "native_coordinate_system", "position")) %>%
    mutate(
      n_vaccines_with_marker_here = ifelse(is.na(n_vaccines_with_marker_here), 0L, n_vaccines_with_marker_here),
      vaccine_specific_marker = diagnostic_category %in% ALL_MARKER_CATEGORIES &
        n_vaccines_with_marker_here == 1
    )
}

# ------------------------------------------------------------------------------
# 1. Load every summary CSV in a directory and combine
# ------------------------------------------------------------------------------
load_all_results <- function(input_dir, pattern = "_(vaccine_(substitutions|nucleotides)|external_reference)_summary\\.csv$") {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(sprintf("No matching CSV files found in %s (pattern: %s)", input_dir, pattern))
  
  cat(sprintf("[i] Found %d CSV file(s) in %s:\n", length(files), input_dir))
  for (f in files) cat(sprintf("    - %s\n", basename(f)))
  
  combined <- map_dfr(files, function(f) {
    df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
    
    # Guard: a CSV with a header but zero data rows (e.g. a vaccine/segment
    # run that found no qualifying positions) has the right columns but
    # nrow(df) == 0. Assigning a length-1 value to a new column of a 0-row
    # data frame trips strict recycling rules in some tidyverse versions
    # ("replacement has 1 row, data has 0"), so skip such files explicitly
    # rather than letting that assignment fail.
    if (nrow(df) == 0) {
      warning(sprintf("Skipping %s -- file has a header but zero data rows.", basename(f)))
      return(NULL)
    }
    
    # Guard: read.csv infers column type from the DATA, ignoring the fact
    # that a field was quoted in the source CSV. One vaccine in this project
    # is literally named "13" (target_accession DQ380213) -- in a file where
    # every vaccine_strain value is "13", read.csv silently returns an
    # INTEGER column, while every other file (e.g. "MP-12", "Smithburn")
    # correctly reads as character. dplyr::bind_rows() then fails with
    # "Can't combine ..$vaccine_strain <integer> and ..$vaccine_strain
    # <character>". Force all known text columns to character explicitly so
    # this can't happen regardless of which values happen to look numeric.
    text_cols <- c("vaccine_strain", "target_accession", "segment", "analysis_mode",
                   "native_coordinate_system", "protein", "substitution",
                   "functional_domain_structural_domain", "potential_diagnostic_value",
                   "diagnostic_category",
                   "number_of_sequences", "geographic_distribution", "host_distribution",
                   "temporal_distribution", "lineage_distribution", "confidence_intervals",
                   "fisher_test_caveat", "variant_frequency_vaccine",
                   "variant_frequency_circulating", "variant_frequency_overall")
    for (col in intersect(text_cols, colnames(df))) {
      df[[col]] <- as.character(df[[col]])
    }
    
    # Same class of bug can hit numeric columns (e.g. a column that's all NA
    # in one file reads as logical there but numeric elsewhere). Force known
    # numeric columns explicitly too.
    numeric_cols <- c("position", "mapped_amino_acid_position", "fisher_p_value",
                      "odds_ratio", "specificity_pct", "ppv_pct", "fdr_adjusted_p_value",
                      "detection_limit_95pct_upper_bound")
    for (col in intersect(numeric_cols, colnames(df))) {
      df[[col]] <- as.numeric(df[[col]])
    }
    
    # Logical columns: same reasoning -- ensure TRUE/FALSE/NA read consistently
    # as logical, not character, across files.
    logical_cols <- c("host_exclusive", "temporal_exclusive", "lineage_exclusive")
    for (col in intersect(logical_cols, colnames(df))) {
      df[[col]] <- as.logical(df[[col]])
    }
    
    required_cols <- c("vaccine_strain", "segment", "native_coordinate_system",
                       "mapped_amino_acid_position", "potential_diagnostic_value",
                       "diagnostic_category", "detection_limit_95pct_upper_bound")
    missing <- setdiff(required_cols, colnames(df))
    if (length(missing) > 0) {
      warning(sprintf("Skipping %s -- missing expected column(s): %s",
                      basename(f), paste(missing, collapse = ", ")))
      return(NULL)
    }
    df$source_file <- basename(f)
    df
  })
  
  if (is.null(combined) || nrow(combined) == 0) {
    stop("No usable rows were loaded -- every matching file was empty or missing required columns. See warnings above.")
  }
  
  # Sanity check: warn if any externally-referenced run has target_accession
  # equal to a value that never appears as "vaccine" strain_type-derived --
  # we can't check the traits file from here, but we CAN warn if an
  # "external_reference" row's target_accession looks like a RefSeq ID
  # (common NCBI RefSeq prefixes), which is the exact mistake this script
  # is documented to guard against.
  suspect <- combined %>%
    filter(analysis_mode == "external_reference") %>%
    filter(str_detect(target_accession, "^NC_")) %>%
    distinct(vaccine_strain, segment, target_accession)
  if (nrow(suspect) > 0) {
    cat("\n[!] WARNING: the following externally-referenced rows have a RefSeq-\n")
    cat("    looking target_accession (starts with 'NC_'). This is very likely\n")
    cat("    the vaccine_reference genome, NOT the vaccine accession itself --\n")
    cat("    double check these before trusting the results:\n")
    print(suspect)
    cat("\n")
  }
  
  cat(sprintf("\n[+] Combined table: %d rows across %d source file(s)\n",
              nrow(combined), length(files)))
  combined
}

# ------------------------------------------------------------------------------
# 2. Cross-validation table: amino-acid vs nucleotide concordance per codon
# ------------------------------------------------------------------------------
build_cross_validation_table <- function(combined) {
  
  is_aa_marker <- function(vals) any(vals %in% SELF_REFERENCED_MARKER_CATEGORIES)
  is_nt_marker <- function(vals) any(vals %in% EXTERNAL_REFERENCE_MARKER_CATEGORIES)
  
  cross <- combined %>%
    group_by(vaccine_strain, segment, mapped_amino_acid_position) %>%
    summarise(
      aa_marker       = is_aa_marker(diagnostic_category[native_coordinate_system == "amino_acid"]),
      nt_marker       = is_nt_marker(diagnostic_category[native_coordinate_system == "nucleotide"]),
      aa_substitution = paste(unique(substitution[native_coordinate_system == "amino_acid"]), collapse = "; "),
      nt_substitutions = paste(unique(substitution[native_coordinate_system == "nucleotide" &
                                                     diagnostic_category %in% EXTERNAL_REFERENCE_MARKER_CATEGORIES]), collapse = "; "),
      aa_classification = paste(unique(potential_diagnostic_value[native_coordinate_system == "amino_acid"]), collapse = "; "),
      n_nucleotide_markers_in_codon = sum(native_coordinate_system == "nucleotide" &
                                            diagnostic_category %in% EXTERNAL_REFERENCE_MARKER_CATEGORIES),
      .groups = "drop"
    ) %>%
    filter(aa_marker | nt_marker) %>%
    mutate(
      concordance_status = case_when(
        aa_marker & nt_marker  ~ "Concordant (marker at both amino-acid and nucleotide levels)",
        aa_marker & !nt_marker ~ "Amino-acid marker only (nucleotide route not fully private in field)",
        !aa_marker & nt_marker ~ "Nucleotide marker only (e.g. synonymous substitution -- invisible at amino-acid level)",
        TRUE ~ NA_character_
      ),
      segment = factor(segment, levels = SEGMENT_ORDER)
    ) %>%
    arrange(segment, vaccine_strain, mapped_amino_acid_position) %>%
    mutate(segment = as.character(segment))
  
  cross
}

# ------------------------------------------------------------------------------
# 3. Main-text table: compact, one row per confirmed marker
# ------------------------------------------------------------------------------
build_main_text_table <- function(combined, cross_validation, require_vaccine_specific = TRUE) {
  
  # Prefer the amino-acid-level row (cleaner notation, e.g. "I219T") whenever
  # a marker is confirmed at the amino-acid level; fall back to the
  # nucleotide-level row (e.g. synonymous markers, only detectable there).
  aa_rows <- combined %>%
    filter(native_coordinate_system == "amino_acid",
           diagnostic_category %in% SELF_REFERENCED_MARKER_CATEGORIES)
  
  nt_only_codons <- cross_validation %>%
    filter(concordance_status == "Nucleotide marker only (e.g. synonymous substitution -- invisible at amino-acid level)") %>%
    select(vaccine_strain, segment, mapped_amino_acid_position)
  
  nt_rows <- combined %>%
    filter(native_coordinate_system == "nucleotide",
           diagnostic_category %in% EXTERNAL_REFERENCE_MARKER_CATEGORIES) %>%
    inner_join(nt_only_codons, by = c("vaccine_strain", "segment", "mapped_amino_acid_position"))
  
  main <- bind_rows(
    aa_rows %>% mutate(marker_source = "amino_acid"),
    nt_rows %>% mutate(marker_source = "nucleotide_only")
  ) %>%
    left_join(
      cross_validation %>% select(vaccine_strain, segment, mapped_amino_acid_position, concordance_status),
      by = c("vaccine_strain", "segment", "mapped_amino_acid_position")
    )
  
  n_before <- nrow(main)
  
  if (require_vaccine_specific) {
    main <- main %>% filter(vaccine_specific_marker)
    cat(sprintf("[i] Main-text table narrowed to vaccine-specific (non-shared) substitutions: %d -> %d rows\n",
                n_before, nrow(main)))
    shared_removed <- n_before - nrow(main)
    if (shared_removed > 0) {
      cat(sprintf("    (%d marker(s) were shared across >=2 vaccine strains and excluded from the main-text\n", shared_removed))
      cat(sprintf("     table; they remain visible, flagged, in the supplementary full frequency table.)\n"))
    }
  }
  
  main %>%
    mutate(
      confirmed_by_independent_reference = concordance_status ==
        "Concordant (marker at both amino-acid and nucleotide levels)",
      segment = factor(segment, levels = SEGMENT_ORDER)
    ) %>%
    arrange(segment, vaccine_strain, mapped_amino_acid_position) %>%
    mutate(segment = as.character(segment)) %>%
    select(
      vaccine_strain,
      segment,
      protein,
      mapped_amino_acid_position,
      substitution,
      functional_domain_structural_domain,
      variant_frequency_circulating,
      potential_diagnostic_value,
      diagnostic_category,
      detection_limit_95pct_upper_bound,
      vaccine_specific_marker,
      n_vaccines_with_marker_here,
      confirmed_by_independent_reference,
      fisher_p_value,
      specificity_pct
    ) %>%
    rename(
      amino_acid_position = mapped_amino_acid_position,
      functional_domain    = functional_domain_structural_domain,
      field_frequency      = variant_frequency_circulating,
      diagnostic_value     = potential_diagnostic_value,
      detection_limit_95pct_upper_bound_pct = detection_limit_95pct_upper_bound
    )
}

# ------------------------------------------------------------------------------
# 4. Supplementary full frequency table: everything, tidied and sorted
# ------------------------------------------------------------------------------
build_supplementary_full_table <- function(combined) {
  combined %>%
    mutate(segment = factor(segment, levels = SEGMENT_ORDER)) %>%
    arrange(segment, vaccine_strain, native_coordinate_system, position) %>%
    mutate(segment = as.character(segment)) %>%
    select(-source_file)
}

# ------------------------------------------------------------------------------
# 5. Formatted output writers
# ------------------------------------------------------------------------------
.write_formatted_sheet <- function(wb, sheet_name, df, highlight_col = NULL, highlight_values = NULL) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold",
                                                          fgFill = "#D9E1F2",
                                                          border = "Bottom"))
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  
  if (!is.null(highlight_col) && highlight_col %in% colnames(df)) {
    col_idx <- which(colnames(df) == highlight_col)
    hl_style <- createStyle(fgFill = "#C6EFCE")
    match_rows <- which(df[[highlight_col]] %in% highlight_values) + 1  # +1 for header row
    if (length(match_rows) > 0) {
      addStyle(wb, sheet_name, style = hl_style, rows = match_rows, cols = col_idx, gridExpand = FALSE)
    }
  }
}

write_main_text_table <- function(main_table, output_path) {
  if (.has_openxlsx) {
    wb <- createWorkbook()
    .write_formatted_sheet(wb, "Main_Text_Markers", main_table,
                           highlight_col = "confirmed_by_independent_reference",
                           highlight_values = TRUE)
    saveWorkbook(wb, output_path, overwrite = TRUE)
    cat(sprintf("[+] Main-text table saved to: %s\n", output_path))
  } else {
    csv_path <- sub("\\.xlsx$", ".csv", output_path)
    write.csv(main_table, csv_path, row.names = FALSE)
    cat(sprintf("[+] Main-text table saved to: %s (CSV fallback)\n", csv_path))
  }
}

write_supplementary_tables <- function(full_table, cross_validation, output_path) {
  if (.has_openxlsx) {
    wb <- createWorkbook()
    .write_formatted_sheet(wb, "Full_Frequency_Table", full_table,
                           highlight_col = "diagnostic_category",
                           highlight_values = c(SELF_REFERENCED_MARKER_CATEGORIES, EXTERNAL_REFERENCE_MARKER_CATEGORIES))
    .write_formatted_sheet(wb, "AA_vs_NT_Cross_Validation", cross_validation,
                           highlight_col = "concordance_status",
                           highlight_values = "Concordant (marker at both amino-acid and nucleotide levels)")
    saveWorkbook(wb, output_path, overwrite = TRUE)
    cat(sprintf("[+] Supplementary workbook saved to: %s\n", output_path))
  } else {
    base_path <- sub("\\.xlsx$", "", output_path)
    write.csv(full_table, paste0(base_path, "_full_frequency_table.csv"), row.names = FALSE)
    write.csv(cross_validation, paste0(base_path, "_cross_validation.csv"), row.names = FALSE)
    cat(sprintf("[+] Supplementary tables saved to: %s_full_frequency_table.csv and %s_cross_validation.csv (CSV fallback)\n",
                base_path, base_path))
  }
}

# ------------------------------------------------------------------------------
# 6. Top-level driver
# ------------------------------------------------------------------------------
build_manuscript_tables <- function(
    input_dir,
    output_dir = ".",
    main_text_filename = "main_text_vaccine_markers.xlsx",
    supplementary_filename = "supplementary_tables.xlsx",
    require_vaccine_specific = TRUE
) {
  combined <- load_all_results(input_dir)
  combined <- add_vaccine_specificity_flags(combined)
  
  n_shared_positions <- combined %>%
    filter(diagnostic_category %in% ALL_MARKER_CATEGORIES, n_vaccines_with_marker_here > 1) %>%
    distinct(segment, native_coordinate_system, position) %>%
    nrow()
  if (n_shared_positions > 0) {
    cat(sprintf("\n[i] %d position(s) are markers for >=2 vaccine strains (shared, not strain-specific):\n", n_shared_positions))
    print(
      combined %>%
        filter(diagnostic_category %in% ALL_MARKER_CATEGORIES, n_vaccines_with_marker_here > 1) %>%
        distinct(segment, native_coordinate_system, position, vaccines_sharing_this_position)
    )
    cat("\n")
  }
  
  cross_validation <- build_cross_validation_table(combined)
  main_table        <- build_main_text_table(combined, cross_validation, require_vaccine_specific = require_vaccine_specific)
  full_table        <- build_supplementary_full_table(combined)
  
  cat(sprintf("\n[i] Main-text table: %d confirmed marker(s) across %d segment(s)\n",
              nrow(main_table), length(unique(main_table$segment))))
  cat(sprintf("[i] Supplementary full frequency table: %d row(s)\n", nrow(full_table)))
  cat(sprintf("[i] Cross-validation table: %d codon(s) flagged as a marker at either level\n",
              nrow(cross_validation)))
  cat("[i] Cross-validation breakdown:\n")
  print(table(cross_validation$concordance_status))
  
  write_main_text_table(main_table, file.path(output_dir, main_text_filename))
  write_supplementary_tables(full_table, cross_validation, file.path(output_dir, supplementary_filename))
  
  invisible(list(main_table = main_table, full_table = full_table, cross_validation = cross_validation))
}

# ------------------------------------------------------------------------------
# 7. EXAMPLE USAGE
#    Point input_dir at the folder containing ALL your output CSVs -- both the
#    self-referenced (amino-acid) and externally-referenced (nucleotide) runs,
#    for every vaccine and every segment (M, L, S-NP, S-NSS). No need to
#    separate them into different folders; routing is automatic via the
#    native_coordinate_system column.
# ------------------------------------------------------------------------------
# result <- build_manuscript_tables(
#   input_dir               = "output-dir/vaccine-marker-summaries",
#   output_dir              = "output-dir/manuscript-tables",
#   main_text_filename      = "main_text_vaccine_markers.xlsx",
#   supplementary_filename  = "supplementary_tables.xlsx"
# )
#
# # Inspect interactively if needed:
# # View(result$main_table)
# # View(result$cross_validation)


# ------------------------------------------------------------------------------
# 7. EXAMPLE USAGE
#    Point input_dir at the folder containing ALL your output CSVs -- both the
#    self-referenced (amino-acid) and externally-referenced (nucleotide) runs,
#    for every vaccine and every segment (M, L, S-NP, S-NSS). No need to
#    separate them into different folders; routing is automatic via the
#    native_coordinate_system column.
# ------------------------------------------------------------------------------
setwd(outDir)

result <- build_manuscript_tables(
  input_dir               = outDir,
  output_dir              = outDir,
  main_text_filename      = "main_text_vaccine_markers.xlsx",
  supplementary_filename  = "supplementary_tables.xlsx"
)

# # Inspect interactively if needed:
# # View(result$main_table)
# # View(result$cross_validation)