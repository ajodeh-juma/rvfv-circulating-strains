#!/usr/bin/env Rscript

library(ape)
library(dplyr)
library(ggnewscale)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(paletteer)
library(phytools)
library(purrr)
library(rprojroot)
library(stringr)
library(patchwork)
library(jsonlite)


# exports
source(sprintf("%s/utils.R", "scripts"))



######################## PLOT SNP DENSITY AND SELECTION PRESSURE ###############


# function to compute SNP density
sliding_window_snp <- function(alignment, vaccine_ref_name, window_size = 200, step = 100) {
  dna <- read.dna(alignment, format = "fasta")
  aln_mat <- as.matrix(dna)
  seq_len <- ncol(aln_mat)
  print(dim(aln_mat))

    # ensure the reference name exists in the alignment
  all_names <- rownames(aln_mat)
  ref_idx <- which(grepl(vaccine_ref_name, all_names))
  
  if (length(ref_idx) == 0) {
    stop(paste("Reference sequence containing", vaccine_ref_name, "not found in alignment."))
  }
  actual_ref_name <- all_names[ref_idx[1]]
  message("using Reference: ", actual_ref_name)

  # initialize storage
  windows <- seq(1, seq_len - window_size, by = step)
  results <- data.frame()

  ref_seq <- aln_mat[vaccine_ref_name, ]

  for (start in windows) {
    end <- start + window_size - 1
    window_aln <- aln_mat[, start:end]

    # calculate distance for this window (proportion of differences)
    # model "raw" calculates the fraction of sites that differ
    dists <- dist.dna(window_aln, model = "raw", pairwise.deletion = TRUE)
    dist_mat <- as.matrix(dists)
    
    # extract distances between all sequences and the vaccine reference
    #     # and convert to percentage
    #     ref_distances <- dist_mat[, actual_ref_name]

    # extract distances to the vaccine strain
    window_dists <- data.frame(
      Position = start + (window_size / 2),
      Mean_SNP_Density = mean(dist_mat[, vaccine_ref_name], na.rm = TRUE) * 100
    )
    results <- rbind(results, window_dists)
  }
  return(results)
}

# sliding_window_snp <- function(alignment_path, vaccine_ref_name,
#                                window_size = 100, step = 50) {
#   # Load the alignment
#   dna <- read.dna(alignment_path, format = "fasta")
#   aln_mat <- as.matrix(dna)
#   seq_len <- ncol(aln_mat)
#   
#   # Ensure the reference name exists in the alignment
#   all_names <- rownames(aln_mat)
#   ref_idx <- which(grepl(vaccine_ref_name, all_names))
#   
#   if (length(ref_idx) == 0) {
#     stop(paste("Reference sequence containing", vaccine_ref_name, "not found in alignment."))
#   }
#   
#   actual_ref_name <- all_names[ref_idx[1]]
#   message("using Reference: ", actual_ref_name)
# 
#   # Initialize storage
#   windows <- seq(1, seq_len - window_size, by = step)
#   results <- list()
#   
#   for (i in seq_along(windows)) {
#     start <- windows[i]
#     end <- start + window_size - 1
#     
#     # Subset the alignment for the current window
#     window_aln <- aln_mat[, start:end]
#     
#     # Calculate distance: 
#     # model = "raw" gives the proportion of sites that differ (p-distance).
#     # pairwise.deletion = TRUE ignores sites with 'n', '-', or '?' in each pair.
#     # pairwise.deletion = TRUE
#     dists <- dist.dna(window_aln, model = "raw")
#     dist_mat <- as.matrix(dists)
#     
#     # Extract distances between all sequences and the vaccine reference
#     # and convert to percentage
#     ref_distances <- dist_mat[, actual_ref_name]
#     
#     # Filter out the self-comparison (distance to itself is 0)
#     ref_distances <- ref_distances[names(ref_distances) != actual_ref_name]
#     
#     results[[i]] <- data.frame(
#       Position = start + (window_size / 2),
#       Mean_SNP_Density = mean(ref_distances, na.rm = TRUE) * 100
#     )
#   }
#   
#   final_results <- bind_rows(results)
#   return(final_results)
# }

################################# Large segment ################################

# snp_data_large_zh548 <- sliding_window_snp(
#   alignment = file.path(segmentsDir[1], "output-dir//riftL/alignment/riftL_filtered.fasta"),
#   vaccine_ref_name = "NC_014397|human|EGY|1977-11-07")

get_summary <- function(data){
  data %>% 
    summarise(
      mean_val = mean(Mean_SNP_Density, na.rm = TRUE),
      sd_val   = sd(Mean_SNP_Density, na.rm = TRUE),
      min_val  = min(Mean_SNP_Density, na.rm = TRUE),
      max_val  = max(Mean_SNP_Density, na.rm = TRUE)
    )
}

get_summary(data = snp_data_large_smithburn)


snp_data_large_smithburn <- sliding_window_snp(
  alignment = file.path(segmentsDir[1], "output-dir/riftL/alignment/riftL_filtered.fasta"),
  vaccine_ref_name = "DQ375430|human|UGA|1944-05-16")
get_summary(data = snp_data_large_smithburn)


snp_data_large_mp12 <- sliding_window_snp(
  alignment = file.path(segmentsDir[1], "output-dir/riftL/alignment/riftL_filtered.fasta"),
  vaccine_ref_name = "DQ375404|human|EGY|1977-11-07")
get_summary(data = snp_data_large_mp12)

snp_data_large_clone13 <- sliding_window_snp(
  alignment = file.path(segmentsDir[1], "output-dir/riftL/alignment/riftL_filtered.fasta"),
  vaccine_ref_name = "DQ375417|human|CAF|1974-01-09")
get_summary(snp_data_large_clone13)


################################# Medium segment ################################

# snp_data_medium_zh548 <- sliding_window_snp(
#   alignment = file.path(segmentsDir[2], "output-dir/riftM/alignment/riftM_filtered.fasta"),
#   vaccine_ref_name = "NC_014396|human|EGY|1977-11-07")

snp_data_medium_smithburn <- sliding_window_snp(
  alignment = file.path(segmentsDir[2], "output-dir/riftM/alignment/riftM_filtered.fasta"),
  vaccine_ref_name = "DQ380193|human|UGA|1944-05-16")
get_summary(snp_data_medium_smithburn)

snp_data_medium_mp12 <- sliding_window_snp(
  alignment = file.path(segmentsDir[2], "output-dir/riftM/alignment/riftM_filtered.fasta"),
  vaccine_ref_name = "DQ380208|human|EGY|1977-11-07")
get_summary(snp_data_medium_mp12)

snp_data_medium_clone13 <- sliding_window_snp(
  alignment = file.path(segmentsDir[2], "output-dir/riftM/alignment/riftM_filtered.fasta"),
  vaccine_ref_name = "DQ380213|human|CAF|1974-01-09")
get_summary(snp_data_medium_clone13)




################################# Small segment - NSs ################################

# snp_data_nss_zh548 <- sliding_window_snp(
#   alignment = file.path(segmentsDir[3], "output-dir/riftS-NSS/alignment/riftS-NSS_filtered.fasta"),
#   vaccine_ref_name = "NC_014395|human|EGY|1977-11-07")

snp_data_nss_smithburn <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NSS/alignment/riftS-NSS_filtered.fasta"),
  vaccine_ref_name = "DQ380157|human|UGA|1944-05-16")
get_summary(snp_data_nss_smithburn)

snp_data_nss_mp12 <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NSS/alignment/riftS-NSS_filtered.fasta"),
  vaccine_ref_name = "DQ380154|human|EGY|1977-11-07")
get_summary(snp_data_nss_mp12)

snp_data_nss_clone13 <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NSS/alignment/riftS-NSS_filtered.fasta"),
  vaccine_ref_name = "DQ380182|human|CAF|1974-01-09")
get_summary(snp_data_nss_clone13)



################################# Small segment - NP ################################

# snp_data_np_zh548 <- sliding_window_snp(
#   alignment = file.path(segmentsDir[3], "output-dir/riftS-NP/alignment/riftS-NP_filtered.fasta"),
#   vaccine_ref_name = "NC_014395|human|EGY|1977-11-07")

snp_data_np_smithburn <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NP/alignment/riftS-NP_filtered.fasta"),
  vaccine_ref_name = "DQ380157|human|UGA|1944-05-16")
get_summary(snp_data_np_smithburn)

snp_data_np_mp12 <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NP/alignment/riftS-NP_filtered.fasta"),
  vaccine_ref_name = "DQ380154|human|EGY|1977-11-07")
get_summary(snp_data_np_mp12)

snp_data_np_clone13 <- sliding_window_snp(
  alignment = file.path(segmentsDir[3], "output-dir/riftS-NP/alignment/riftS-NP_filtered.fasta"),
  vaccine_ref_name = "DQ380182|human|CAF|1974-01-09")
get_summary(snp_data_np_clone13)


# plot snp density alongside gene annotation tracks and selection pressure results

# helper function for hyphy outputs (FEL, FUBAR and SLAC)
# each list should have: Codon_Site, p_value, and Omega (dN/dS)


extract_hyphy_results <- function(file_path, method) {
  data <- fromJSON(file_path)
  
  if (method == "SLAC") {
    sites <- data$MLE$content$`0`$`by-branch`$AVERAGED
    print(sites)
    hyphy_results <- data.frame(
      Codon = 1:nrow(sites),
      Omega = as.numeric(sites[,3]),
      Score = as.numeric(sites[,9]),
      Method = "SLAC") |> 
      mutate(Position = Codon * 3) |>
      select(Codon, Omega, Score, Method, Position)
    
  } else if (method == "FEL") {
    sites <- data$MLE$content$`0`
    hyphy_results <- data.frame(
      Codon = 1:nrow(sites),
      Omega = sites[,6],
      dS = sites[,1],
      dN = sites[,2],
      Score = as.numeric(sites[,5]),
      Method = "FEL") |>
      mutate(Position = Codon * 3,
             Omega = dN / dS) |>
      mutate(Omega = case_when(
        Omega == 'Inf' ~ 10,
        TRUE ~ Omega)) |>
      select(Codon, dS, dN, Omega, Score, Method, Position)
    
  } else if (method == "FUBAR") {
    sites <- data$MLE$content$`0`
    hyphy_results <- data.frame(
      Codon = 1:nrow(sites),
      Omega = NA,
      dS = sites[,1],
      dN = sites[,2],
      dNminusdS = sites[,3],
      Score = sites[,5],
      Method = "FUBAR") |>
      mutate(Position = Codon * 3,
             Omega = dN / dS) |>
      select(Codon, dS, dN, Omega, Score, Method, Position)
  }
  return(hyphy_results)
}



get_significant_sites <- function(hyphy_data){
  
  # get significant sites
  # (Omega) dN/dS > 1 & FEL: p < 0.05 (significant positive), FUBAR: posterior > 0.9
  selection_data <- bind_rows(hyphy_data) |>
    # mutate(Significant = (Method == "FEL" & Score < 0.05) | (Score > 0.9)) |>
    # mutate(Significant = ifelse(Method == "FEL", Score < 0.05, Score > 0.9)) |>
    mutate(Significant = case_when(
      Method == "FEL"   ~ Score < 0.05,
      Method == "FUBAR" ~ Score > 0.9,
      TRUE              ~ FALSE 
    )) |>
    filter(Significant == TRUE) |>
    # filter(Omega >= 1.0) |>
    group_by(Position)
  
  print(selection_data)
  
  consensus_sites <- selection_data |>
    filter(Significant == TRUE & (Omega > 1.0)) |>
    group_by(Position) |>
    filter(n_distinct(Method) >= 1)
  
  print(consensus_sites)
  
  return(list(selection_data, consensus_sites))
  
}



####################### Large segment hyphy results ############################
large_hyphy_files <- 
  file.path(file.path(segmentsDir[1], "output-dir/riftL/hyphy"), 
            list.files(path = file.path(segmentsDir[1], "output-dir/riftL/hyphy"), 
                       pattern = "\\.FEL.json|FUBAR.json"))
names(large_hyphy_files) <- c("FEL", "FUBAR")
large_hyphy_data <- mapply(extract_hyphy_results,
                           file_path = large_hyphy_files,
                           method = names(large_hyphy_files),
                           SIMPLIFY = FALSE)
selection_data_large <- bind_rows(large_hyphy_data)
sel_data_large <- get_significant_sites(hyphy_data = selection_data_large)

####################### Medium segment hyphy results ###########################

medium_hyphy_files <- 
  file.path(file.path(segmentsDir[2], "output-dir/riftM/hyphy"), 
            list.files(path = file.path(segmentsDir[2], "output-dir/riftM/hyphy"), 
                       pattern = "\\.FEL.json|FUBAR.json"))
names(medium_hyphy_files) <- c("FEL", "FUBAR")
medium_hyphy_data <- mapply(extract_hyphy_results,
                            file_path = medium_hyphy_files,
                            method = names(medium_hyphy_files),
                            SIMPLIFY = FALSE)

# consolidate
selection_data_medium <- bind_rows(medium_hyphy_data)
sel_data_medium <- get_significant_sites(hyphy_data = selection_data_medium)


###################### Small segment - NP hyphy results ########################

np_hyphy_files <- 
  file.path(file.path(segmentsDir[3], "output-dir/riftS-NP/hyphy"), 
            list.files(path = file.path(segmentsDir[3], "output-dir/riftS-NP/hyphy"), 
                       pattern = "\\.FEL.json|FUBAR.json"))
names(np_hyphy_files) <- c("FEL", "FUBAR")
np_hyphy_data <- mapply(extract_hyphy_results,
                        file_path = np_hyphy_files,
                        method = names(np_hyphy_files),
                        SIMPLIFY = FALSE)

# consolidate
selection_data_np <- bind_rows(np_hyphy_data)
sel_data_np <- get_significant_sites(selection_data_np)


###################### Small segment - NSs hyphy results ########################

nss_hyphy_files <- 
  file.path(file.path(segmentsDir[3], "output-dir/riftS-NSS/hyphy"), 
            list.files(path = file.path(segmentsDir[3], "output-dir/riftS-NSS/hyphy"), 
                       pattern = "\\.FEL.json|FUBAR.json"))
names(nss_hyphy_files) <- c("FEL", "FUBAR")
nss_hyphy_data <- mapply(extract_hyphy_results,
                         file_path = nss_hyphy_files,
                         method = names(nss_hyphy_files),
                         SIMPLIFY = FALSE)

# consolidate
selection_data_nss <- bind_rows(nss_hyphy_data)
sel_data_nss <- get_significant_sites(selection_data_nss)

# annotations for the different segments
L_segment_gene_annotations <- data.frame(
  gene = as.factor(c("RdRp")),
  start = c(1), 
  end = c(6276),
  colour = c("#CDBA96")
)

M_segment_gene_annotations <- data.frame(
  gene = as.factor(c("NSm", "Gn", "Gc")),
  start = c(1, 480, 2091), 
  end = c(479, 2090, 3591),
  colour = c("#8B8970", "#8B2323", "#458B74")
)

NSS_segment_gene_annotations <- data.frame(
  gene = c("NSs"),
  start = c(1), 
  end = c(795),
  colour = c("#CD9B9B")
)

NP_segment_gene_annotations <- data.frame(
  gene = c("NP"),
  start = c(1), 
  end = c(735),
  colour = c("#CD853F")
)


# function to plot snp density and selection pressure data
plot_snp_density <- function(snp_data, selection_data, significant_sites, 
                             col_area, gene_annotations, reference, segment) {
  
  
  # get mean and sd and define threshold
  mean_snp <- mean(snp_data$Mean_SNP_Density)
  sd_snp <- sd(snp_data$Mean_SNP_Density)
  threshold <- mean_snp + sd_snp
  
  print(paste(segment, ": ", reference))
  hotspots <- snp_data |> 
    filter(Mean_SNP_Density > (mean_snp + sd_snp))
  print(hotspots)
  print("-----------------------------")
  
  
  # identify 
  find_peaks <- function(x) {
    pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 1
    return(pks)
  }
  
  
  peaks <- snp_data[find_peaks(snp_data$Mean_SNP_Density), ]
  peaks_above_threshold <- peaks[peaks$Mean_SNP_Density > threshold, ]
  
  # select relevant columns for the export
  hotspot_export <- peaks_above_threshold |>
    select(Position, Mean_SNP_Density) |>
    rename(Genomic_Coordinate = Position, SNP_Percentage = Mean_SNP_Density)
  
  # get top 3 peaks
  top_3_peaks <- peaks_above_threshold |>
    slice_max(order_by = Mean_SNP_Density, n = 3)
  
  print(paste(segment, ": ", reference))
  print(hotspot_export)
  print("-----------------------------")
  
  # write to output file
  # write.csv(hotspot_export, paste("RVFV_", segment, "_hotspot_coordinates.csv"), row.names = FALSE)
  
  p1 <- ggplot() +
    geom_area(data = snp_data, 
              aes(x = Position, y = Mean_SNP_Density), 
              fill = col_area, alpha = 0.5) +
    geom_line(data = snp_data, 
              aes(x = Position, y = Mean_SNP_Density), color = "black",
              linewidth = 0.3) +
    geom_hline(yintercept = mean_snp, linetype = "dashed", 
               color = "#525252", size = 0.5) +
    # annotate("text", 
    #          x = max(snp_data$Position) * 0.10,
    #          # y = max(snp_data$Mean_SNP_Density) * 0.75,
    #          y = mean(snp_data$Mean_SNP_Density),
    #          label = paste0("Mean: ", round(mean_snp, 3), "\n", 
    #                         "SD: ", round(sd_snp, 3)),
    #          color = "black", 
    #          size = 4, 
    #          fontface = "italic") +
    annotate("label", # Using label for better readability in the corner
             x = max(snp_data$Position), 
             y = max(snp_data$Mean_SNP_Density),
             label = paste0("Mean: ", round(mean_snp, 3), "%\n", 
                            "SD: ", round(sd_snp, 3)),
             hjust = 1.1, vjust = 1.1,
             color = "black", 
             fill = "white",
             alpha = 0.8,
             size = 3.5, 
             fontface = "italic") +
    geom_text_repel(data = top_3_peaks,
                    aes(x = Position, y = Mean_SNP_Density, 
                        label = paste0(round(Mean_SNP_Density, 2), "%")),
                    size = 3.0, 
                    color = "black",
                    fontface = "bold",
                    family = dviz_font_family_bold,
                    box.padding = 0.5,      # space around the text
                    point.padding = 0.3,    # space from the peak point
                    direction = "y",        # nudge them vertically if needed
                    nudge_y = 1,            # push them slightly above the peak
                    segment.color = 'grey50', # add a small line pointing to the peak
                    min.segment.length = 0) +
    # scale_x_continuous(breaks = sort(c(seq(0, max(snp_data$Position), 500), peaks_above_threshold$Position))) +
    geom_rect(data = gene_annotations, 
              aes(xmin = start, xmax = end, ymin = -1, ymax = -0.2, fill = gene), 
              alpha = 0.8) +
    scale_fill_manual(values = gene_annotations$colour,
                      breaks = gene_annotations$gene) +
    labs(# title = paste("Evolutionary Divergence and Selection Pressure"),
         subtitle = glue::glue("Vaccine: {reference} | Segment: {segment}"),
         x = "Position (bp)", 
         y = "SNP Density (%)") +
    geom_vline(data = significant_sites, aes(xintercept = Position),
               color = "#525252", linetype = "dashed", alpha = 0.5) +
    theme_dviz_hgrid() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.3),
      text = element_text(family = dviz_font_family),
      plot.margin = margin(3.5, 1.5, 3.5, 1.5),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()) +
    guides(fill = guide_legend(nrow = 3,
                               override.aes = list(size = 6),
                               title = "Gene", order = 1))
  
  
  p2 <- ggplot(selection_data, 
               aes(x = Position, y = Omega, colour = Significant)) +
    geom_point(aes(color = Method), size = 3, alpha = 0.7) +
    geom_vline(data = significant_sites, aes(xintercept = Position),
               color = "#525252", linetype = "dashed", linewidth = 0.5) +
    geom_text(data = significant_sites, aes(x = Position, y = Omega - 1, 
                                            label = paste("Codon site:", Codon)), 
              angle = 0, hjust = 0.5, size = 3.5, 
              color = "black",
              fontface = "bold",
              family = dviz_font_family_bold) +
    facet_grid(Method ~ ., scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "#525252") +
    scale_color_brewer(palette = "Set1") +
    labs(#title = "Significant Sites Under Selection Pressure",
      #subtitle = "Filtered for FEL (p < 0.05) and FUBAR (Posterior > 0.9)",
      y = expression(italic("dN/dS") ~ (omega)),
      x = "Nucleotide Position (bp)") +
    theme_dviz_hgrid() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.3),
      text = element_text(family = dviz_font_family),
      plot.margin = margin(3.5, 1.5, 3.5, 1.5),
      axis.text.x = element_text(angle = 0, hjust = 1)) +
    guides(fill = guide_legend(nrow = 3,
                               override.aes = list(size = 6),
                               title = "Gene", order = 1))
  
  p <- (p1/p2) + 
    plot_layout(heights = c(1, 2))
  # p <- wrap_elements(p)
  print(p)
  return(p)
}



################# Large #####################

# 
# snp_density_large_zh548 <- plot_snp_density(
#   snp_data = snp_data_large_zh548,
#   selection_data = sel_data_large[[1]],
#   significant_sites = sel_data_large[[2]],
#   col_area = "#8B7D6B",
#   gene_annotations = L_segment_gene_annotations,
#   reference = "ZH-548",
#   segment = "Large")


snp_density_large_smithburn <- plot_snp_density(
  snp_data = snp_data_large_smithburn,
  selection_data = sel_data_large[[1]],
  significant_sites = sel_data_large[[2]],
  col_area = "#8B7D6B",
  gene_annotations = L_segment_gene_annotations,
  reference = "Smithburn",
  segment = "Large")

snp_density_large_mp12 <- plot_snp_density(
  snp_data = snp_data_large_mp12,
  selection_data = sel_data_large[[1]],
  significant_sites = sel_data_large[[2]],
  col_area = "#8B7D6B",
  gene_annotations = L_segment_gene_annotations,
  reference = "MP-12",
  segment = "Large")

snp_density_large_clone13 <- plot_snp_density(
  snp_data = snp_data_large_clone13,
  selection_data = sel_data_large[[1]],
  significant_sites = sel_data_large[[2]],
  col_area = "#8B7D6B",
  gene_annotations = L_segment_gene_annotations,
  reference = "Clone 13",
  segment = "Large")


################# Medium #####################
# snp_density_medium_zh548 <- plot_snp_density(
#   snp_data = snp_data_medium_zh548,
#   selection_data = sel_data_medium[[1]],
#   significant_sites = sel_data_medium[[2]],
#   col_area = "#8B5742",
#   gene_annotations = M_segment_gene_annotations,
#   reference = "ZH-548",
#   segment = "Medium")


snp_density_medium_smithburn <- plot_snp_density(
  snp_data = snp_data_medium_smithburn,
  selection_data = sel_data_medium[[1]],
  significant_sites = sel_data_medium[[2]],
  col_area = "#8B5742",
  gene_annotations = M_segment_gene_annotations,
  reference = "Smithburn",
  segment = "Medium")

ggsave(file.path(outDir, "smithburnM.pdf"),
       snp_density_medium_smithburn, 
       width = 6, height = 8, 
       dpi = 300, device = cairo_pdf)

snp_density_medium_mp12 <- plot_snp_density(
  snp_data = snp_data_medium_mp12,
  selection_data = sel_data_medium[[1]],
  significant_sites = sel_data_medium[[2]],
  col_area = "#8B5742",
  gene_annotations = M_segment_gene_annotations,
  reference = "MP-12",
  segment = "Medium")


snp_density_medium_clone13 <- plot_snp_density(
  snp_data = snp_data_medium_clone13,
  selection_data = sel_data_medium[[1]],
  significant_sites = sel_data_medium[[2]],
  col_area = "#8B5742",
  gene_annotations = M_segment_gene_annotations,
  reference = "Clone 13",
  segment = "Medium")



################# Small - NSs #####################
# snp_density_nss_zh548 <- plot_snp_density(
#   snp_data = snp_data_nss_zh548,
#   selection_data = sel_data_nss[[1]],
#   significant_sites = sel_data_nss[[2]],
#   col_area = "#8B5742",
#   gene_annotations = NSS_segment_gene_annotations,
#   reference = "ZH-548",
#   segment = "Small (NSs)")


snp_density_nss_smithburn <- plot_snp_density(
  snp_data = snp_data_nss_smithburn,
  selection_data = sel_data_nss[[1]],
  significant_sites = sel_data_nss[[2]],
  col_area = "#8B5742",
  gene_annotations = NSS_segment_gene_annotations,
  reference = "Smithburn",
  segment = "Small (NSs)")

snp_density_nss_mp12 <- plot_snp_density(
  snp_data = snp_data_nss_mp12,
  selection_data = sel_data_nss[[1]],
  significant_sites = sel_data_nss[[2]],
  col_area = "#8B5742",
  gene_annotations = NSS_segment_gene_annotations,
  reference = "MP-12",
  segment = "Small (NSs)")

snp_density_nss_clone13 <- plot_snp_density(
  snp_data = snp_data_nss_clone13,
  selection_data = sel_data_nss[[1]],
  significant_sites = sel_data_nss[[2]],
  col_area = "#8B5742",
  gene_annotations = NSS_segment_gene_annotations,
  reference = "Clone 13",
  segment = "Small (NSs)")



################# Small - NP #####################
# snp_density_np_zh548 <- plot_snp_density(
#   snp_data = snp_data_np_zh548,
#   selection_data = sel_data_np[[1]],
#   significant_sites = sel_data_np[[2]],
#   col_area = "#8B5742",
#   gene_annotations = NP_segment_gene_annotations,
#   reference = "ZH-548",
#   segment = "Small (NP)")


snp_density_np_smithburn <- plot_snp_density(
  snp_data = snp_data_np_smithburn,
  selection_data = sel_data_np[[1]],
  significant_sites = sel_data_np[[2]],
  col_area = "#8B5742",
  gene_annotations = NP_segment_gene_annotations,
  reference = "Smithburn",
  segment = "Small (NP)")

snp_density_np_mp12 <- plot_snp_density(
  snp_data = snp_data_np_mp12,
  selection_data = sel_data_np[[1]],
  significant_sites = sel_data_np[[2]],
  col_area = "#8B5742",
  gene_annotations = NP_segment_gene_annotations,
  reference = "MP-12",
  segment = "Small (NP)")

snp_density_np_clone13 <- plot_snp_density(
  snp_data = snp_data_np_clone13,
  selection_data = sel_data_np[[1]],
  significant_sites = sel_data_np[[2]],
  col_area = "#8B5742",
  gene_annotations = NP_segment_gene_annotations,
  reference = "Clone 13",
  segment = "Small (NP)")

###############


clean_panel_theme <- theme(
  legend.position = "none",
  strip.text = element_blank()
)

p_A <- wrap_elements(
  (snp_density_large_smithburn & clean_panel_theme) + 
    theme(axis.title.x = element_blank())
)

p_B <- wrap_elements(
  (snp_density_large_mp12 & theme(axis.title.y = element_blank())) & 
    clean_panel_theme + 
    theme(axis.title.x = element_blank())
)

p_C <- wrap_elements(
  (snp_density_large_clone13 & theme(axis.title.y = element_blank())) + 
    theme(axis.title.x = element_blank())
)

p_D <- wrap_elements(
  (snp_density_medium_smithburn) & 
    clean_panel_theme
)

p_E <- wrap_elements(
  (snp_density_medium_mp12 & theme(axis.title.y = element_blank())) & 
    clean_panel_theme
)

p_F <- wrap_elements(
  (snp_density_medium_clone13 & theme(axis.title.y = element_blank()))
)

figure3_final <- (p_A | p_B | p_C) / (p_D | p_E | p_F) + 
  plot_annotation(tag_levels = list(c("A", "C", "E", "B", "D", "F"))) & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0, 1)
  )


# 1. Define the clean theme
clean_panel_theme <- theme(
  legend.position = "none",
  strip.text = element_blank()
)

# ---------------------------------------------------------
# LARGE SEGMENT (Row 1)
# ---------------------------------------------------------
p_A <- wrap_elements((snp_density_large_smithburn & clean_panel_theme) + theme(axis.title.x = element_blank()))
p_B <- wrap_elements((snp_density_large_mp12 & theme(axis.title.y = element_blank()) & clean_panel_theme) + theme(axis.title.x = element_blank()))
p_C <- wrap_elements((snp_density_large_clone13 & theme(axis.title.y = element_blank())) + theme(axis.title.x = element_blank()))

# ---------------------------------------------------------
# MEDIUM SEGMENT (Row 2)
# ---------------------------------------------------------
p_D <- wrap_elements(snp_density_medium_smithburn & clean_panel_theme + theme(axis.title.x = element_blank()))
p_E <- wrap_elements((snp_density_medium_mp12 & theme(axis.title.y = element_blank()) & clean_panel_theme) + theme(axis.title.x = element_blank()))
p_F <- wrap_elements(snp_density_medium_clone13 & theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()))

# ---------------------------------------------------------
# S-NSS SEGMENT (Row 3)
# ---------------------------------------------------------
p_G <- wrap_elements(snp_density_nss_smithburn & clean_panel_theme + theme(axis.title.x = element_blank()))
p_H <- wrap_elements((snp_density_nss_mp12 & theme(axis.title.y = element_blank()) & clean_panel_theme) + theme(axis.title.x = element_blank()))
p_I <- wrap_elements(snp_density_nss_clone13 & theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()))

# ---------------------------------------------------------
# S-NP SEGMENT (Row 4)
# ---------------------------------------------------------
p_J <- wrap_elements(snp_density_np_smithburn & clean_panel_theme)
p_K <- wrap_elements(snp_density_np_mp12 & theme(axis.title.y = element_blank()) & clean_panel_theme)
p_L <- wrap_elements(snp_density_np_clone13 & theme(axis.title.y = element_blank()))

# ---------------------------------------------------------
# FINAL ASSEMBLY
# ---------------------------------------------------------
# Organized by Rows: L (A-C), M (D-F), S-NSs (G-I), S-NP (J-L)
figure3_final <- (p_A | p_B | p_C) / 
  (p_D | p_E | p_F) / 
  (p_G | p_H | p_I) / 
  (p_J | p_K | p_L) + 
  plot_annotation(tag_levels = 'A') & # Simplified tagging or custom list below
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0, 1)
  )

# If you want to maintain your specific non-sequential tag list:
# plot_annotation(tag_levels = list(c("A", "C", "E", "G", "I", "K", "B", "D", "F", "H", "J", "L")))

# Display the final publication-ready plot
print(figure3_final)



ggsave(file.path(outDir, "Figure3_200_100bp.pdf"), 
       figure3_final, width = 24, height = 21.5, 
       dpi = 300, device = cairo_pdf)

ggsave(file.path(outDir, "Figure3_200_100bp.png"), 
       figure3_final,
       width = 24, height = 21.5, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")



