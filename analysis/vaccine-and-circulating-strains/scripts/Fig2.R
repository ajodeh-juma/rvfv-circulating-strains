#!/usr/bin/env Rscript


library(dplyr)
library(scales)
library(ggnewscale)
library(ggplot2)
library(ggtext)
library(ggthemes)
library(lemon)
library(lubridate)
library(paletteer)
library(patchwork)
library(rprojroot)
library(RColorBrewer)
library(tidyverse)


# exports
source(sprintf("%s/utils.R", "scripts"))



################################################################################
############################## substitution plots ##############################
################################################################################

prepare_tstv_data <- function(segment, segment_length, alignment_dir){
  #'
  #'@param name description
  #'@param name description
  #'@param name description
  
  dir <- file.path(alignment_dir)
  fns <- file.path(dir, 
                   list.files(path = dir, 
                              pattern = "\\.parsed.bcftools.stats.csv$", 
                              recursive = TRUE))
  print(fns)
  
  data_list <- list()
  prefixes <- list()
  for (i in 1:length(fns)) {
    fn = fns[i]
    sample_name <- strsplit(basename(fn), ".", fixed=T)[[1]][3]
    print(strsplit(basename(fn), ".", fixed=T)[[1]][3])
    df <- read.csv(fn, sep=',', header=T)
    df$name <- sample_name
    data_list[[i]] <- df
  }
  d <- do.call(rbind, data_list)
  
  if (segment == 'Medium'){
    d$reference[d$name == 'DQ380193'] <- 'Smithburn'
    d$reference[d$name == 'DQ380208'] <- 'MP-12'
    d$reference[d$name == 'DQ380213'] <- 'Clone-13'
  }
  
  if (segment == 'Large'){
    d$reference[d$name == 'DQ375430'] <- 'Smithburn'
    d$reference[d$name == 'DQ375404'] <- 'MP-12'
    d$reference[d$name == 'DQ375417'] <- 'Clone-13'
    
  }
  
  if (segment == 'Small'){
    d$reference[d$name == 'DQ380157'] <- 'Smithburn'
    d$reference[d$name == 'DQ380154'] <- 'MP-12'
    d$reference[d$name == 'DQ380182'] <- 'Clone-13'
  }
  
  d <- d %>% dplyr::rename(
    'A->C'=A.C, 'A->G'=A.G, 'A->T'=A.T, 'C->A'=C.A, 'C->G'=C.G, 'C->T'=C.T,
    'G->A'=G.A, 'G->C'=G.C, 'G->T'=G.T, 'T->A'=T.A, 'T->C'=T.C, 'T->G'=T.G
  )
  
  
  mutations <- c('A->C', 'C->A', 'A->G', 'G->A', 'A->T', 'T->A',
                 'C->G', 'G->C', 'C->T', 'T->C', 'G->T', 'T->G')
  # cols <- brewer.pal(12, "Paired")
  
  d <- d %>% dplyr::select(c(samples, reference, ts, tv, 
                             tstv, tsALT, tvALT, tstvALT, 
                             'A->C', 'C->A', 'A->G', 'G->A', 'A->T', 'T->A',
                             'C->G', 'G->C', 'C->T', 'T->C', 'G->T', 'T->G'))
  
  d.fq <- d
  print(names(d.fq)[9:20])
  d.fq$sum <- rowSums(d.fq[,names(d.fq)[9:20]])
  
  d.fq <- d.fq %>% 
    mutate(
      fAC = d$`A->C`/segment_length,
      fCA = d$`C->A`/segment_length,
      fAG = d$`A->G`/segment_length,
      fGA = d$`G->A`/segment_length,
      fAT = d$`A->T`/segment_length,
      fTA = d$`T->A`/segment_length,
      fCG = d$`C->G`/segment_length,
      fGC = d$`G->C`/segment_length,
      fCT = d$`C->T`/segment_length,
      fTC = d$`T->C`/segment_length,
      fGT = d$`G->T`/segment_length,
      fTG = d$`T->G`/segment_length
    )

  write.csv(d, file.path(outDir, paste0(segment, "_tstv.metrics.csv")), quote = F,
            row.names = F)
  d.melt <- reshape2::melt(d)
  dfq.melt <- reshape2::melt(d.fq)
  
  res <- list(d.melt, dfq.melt)
  return(res)
}


plot_substitutions <- function(data, reference, segment, gene){
  
  #'@description
  #'A short description...
  #'
  #'@param name description
  #'@param name description
  #'
  
  
  mutations <- c('A->C', 'C->A', 'A->G', 'G->A', 'A->T', 'T->A',
                 'C->G', 'G->C', 'C->T', 'T->C', 'G->T', 'T->G')
  
  # transitions
  ct_tc <- c('C->T', 'T->C')
  ag_ga <- c('A->G', 'G->A')
  
  # transversions
  transversions <- setdiff(mutations, c(ct_tc, ag_ga))
  
  # statistical test for across all the 3 vaccine strains
  strain_stats <- data %>%
    filter(variable %in% c("ts", "tv"),
           reference %in% c("MP-12", "Clone-13", "Smithburn")) %>%
    tidyr::pivot_wider(names_from = variable, values_from = value)
  
  # contingency table for ts/tv
  ts_tv_matrix <- as.matrix(strain_stats[, c("ts", "tv")])
  chi_strains <- chisq.test(ts_tv_matrix)
  p_strains <- format.pval(chi_strains$p.value, digits = 3)
  print(p_strains)

  # process data
  ref_data <- data %>%
    filter(reference == !!reference, variable %in% mutations) %>%
    mutate(group = case_when(
      variable %in% ct_tc ~ "C->T / T->C",
      variable %in% ag_ga ~ "A->G / G->A",
      TRUE ~ "Other (Transversions)"
    ))
  
  # Goodness-of-fit test (Null hypothesis: 50/50 distribution)
  ag_ga_count <- sum(ref_data$value[ref_data$variable %in% c("A->G", "G->A")])
  ct_tc_count <- sum(ref_data$value[ref_data$variable %in% c("C->T", "T->C")])
  chi_within <- chisq.test(c(ag_ga_count, ct_tc_count), p = c(0.5, 0.5))
  p_within <- format.pval(chi_within$p.value, digits = 3)
  sig_label <- ifelse(chi_within$p.value < 0.05, "Significant Difference (*)", "No Significance (ns)")
  print(p_within)
  
  # statistical comparison: Transitions vs. Transversions
  # test if the mean count per mutation type in the groups is significantly different
  group_summary <- ref_data %>%
    group_by(group) %>%
    summarise(total_val = sum(value), n_types = n(), mean_val = mean(value))
  
  # chi-squared test: Observed vs Expected (based on number of available types)
  # expected ratio is 2:2:8 based on the number of mutation types in each category
  obs_counts <- group_summary$total_val
  expected_probs <- group_summary$n_types / sum(group_summary$n_types)
  
  stat_test <- chisq.test(obs_counts, p = expected_probs)
  p_val <- format.pval(stat_test$p.value, digits = 3)
  print(p_val)
  
  max_y <- max(ref_data$value)
  annot_y <- max_y * 1.1
  annot_x <- 12
  
  # highlight the groups using a distinct color palette
  group_cols <- c("C->T / T->C" = "#FF7F0FFF", "A->G / G->A" = "#3CB7CCFF", 
                  "Other (Transversions)" = "#82853BFF")
  
  p <- ggplot(ref_data, aes(x = reorder(variable, -value), y = value, fill = group)) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    scale_fill_manual(values = group_cols) +
    labs(
      title = paste(reference, "SNVs Profile", sep = " "),
      subtitle = glue::glue("Gene: {gene} | Segment: {segment}"),
                            # **Strain Comparison (ts/tv):** p = {p_strains} (sig_label)<br>
                            # **Group Significance (Chi-sq):** p = {p_val}"
      x = "Substitution Type",
      y = "SNV Count",
      fill = "Mutation Category") +
    annotate("text", 
             x = annot_x, 
             y = annot_y, 
             label = paste0("Strain Comparison (ts/tv): p = ", p_strains, "\n",
                            "Within Strain (AG/GA vs CT/TC): p = ", p_within, " ", sig_label, "\n",
                            "Group Significance (Chi-sq): p = ", p_val),
             hjust = 1, # Right-aligned so it doesn't bleed off the edge
             vjust = 1, 
             size = 3.5, 
             fontface = "italic", 
             color = "black") +
    theme_classic() +
    theme_dviz_open() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      plot.subtitle = element_markdown(size = 11, lineheight = 1.2),
      axis.text.x = element_blank(),
      legend.position = c(0.8, 0.5)
    ) +
    coord_capped_cart(bottom=capped_horizontal(c("both")), 
                      left=capped_vertical(c("both"))) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
  
  # summary inset plot showing the average frequency per type
  # p_inset <- ggplot(group_summary, aes(x = group, y = mean_val, fill = group)) +
  #   geom_bar(stat = "identity", width = 0.6) +
  #   scale_fill_manual(values = group_cols) +
  #   labs(title = "Average Count/Type", x = NULL, y = NULL) +
  #   theme_void() +
  #   theme(legend.position = "none", 
  #         plot.title = element_text(size = 8, face = "bold"))
  # 
  # final_plot <- p + inset_element(p_inset, 
  #                                 left = 0.6, bottom = 0.6, 
  #                                 right = 0.95, top = 0.95)
  # 
  return(p)
}


############################### Large segment ##################################

inDirLarge <- sprintf(
  "%s/output-dir/riftL/strain-types",
  segmentsDir[[1]])

tstvLarge <- prepare_tstv_data(segment = "Large", 
                               segment_length = 6276,
                               alignment_dir = inDirLarge)


subsLargeMP12 <- plot_substitutions(data = tstvLarge[[1]], 
                                    reference = "MP-12", 
                                    segment = 'Large', 
                                    gene = "RdRp")
subsLargeMP12

subsLargeMP12 <- subsLargeMP12 + theme(legend.position = "none")
subsLargeMP12


subsLargeClone13 <- plot_substitutions(data = tstvLarge[[1]], 
                                       reference = "Clone-13", 
                                       segment = 'Large', 
                                       gene = "RdRp")
subsLargeClone13

subsLargeSmithburn <- plot_substitutions(data = tstvLarge[[1]], 
                                         reference = "Smithburn", 
                                         segment = 'Large', 
                                         gene = "RdRp")
subsLargeSmithburn <- subsLargeSmithburn + 
  theme(legend.position = "none")
subsLargeSmithburn


plots.large.tstv <- (subsLargeMP12 / subsLargeClone13 / subsLargeSmithburn) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 15))

plots.large.tstv

ggsave(file.path(outDir, "tstv-large.pdf"), plots.large.tstv, 
       width = 12, height = 10, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)


############################## Medium segment ##################################

inDirMedium <- sprintf(
  "%s/output-dir/riftM/strain-types",
  segmentsDir[[2]])

tstvMedium <- prepare_tstv_data(segment = "Medium",
                                segment_length = 3591,
                                alignment_dir = inDirMedium)
subsMediumMP12 <- plot_substitutions(data = tstvMedium[[1]], 
                                     reference = "MP-12", 
                                     segment = 'Medium', 
                                     gene = "NSm/Gn/Gc")
subsMediumMP12 <- subsMediumMP12 + 
  theme(legend.position = "none", 
        axis.text.x = element_blank())

subsMediumMP12

subsMediumClone13 <- plot_substitutions(data = tstvMedium[[1]], 
                                        reference = "Clone-13", 
                                        segment = 'Medium', 
                                        gene = "NSm/Gn/Gc")
subsMediumClone13 <- subsMediumClone13 + theme(
  axis.text.x = element_blank()
)

subsMediumSmithburn <- plot_substitutions(data = tstvMedium[[1]], 
                                          reference = "Smithburn", 
                                          segment = 'Medium', 
                                          gene = "NSm/Gn/Gc")
subsMediumSmithburn <- subsMediumSmithburn + theme(legend.position = "none")
subsMediumSmithburn


plots.medium.tstv <- (subsMediumMP12 / subsMediumClone13 / subsMediumSmithburn) + 
  plot_layout(guides = 'collect', axis_titles = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 15))

ggsave(file.path(outDir, "tstv-medium.pdf"), plots.medium.tstv, 
       width = 15, height = 13, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)


############################ Small segment: NSs ################################

inDirNSS <- sprintf(
  "%s/output-dir/riftS-NSS/strain-types",
  segmentsDir[[3]])

tstvNSS <- prepare_tstv_data(segment = "Small", 
                             segment_length = 795,
                             alignment_dir = inDirNSS)
subsNssMP12 <- plot_substitutions(data = tstvNSS[[1]], 
                                  reference = "MP-12", 
                                  segment = 'Small', 
                                  gene = "NSs")
subsNssMP12 <- subsNssMP12 + theme(legend.position = "none")
subsNssMP12

subsNssClone13 <- plot_substitutions(data = tstvNSS[[1]], 
                                     reference = "Clone-13", 
                                     segment = 'Small', 
                                     gene = "NSs")
subsNssClone13

subsNssSmithburn <- plot_substitutions(data = tstvNSS[[1]], 
                                       reference = "Smithburn", 
                                       segment = 'Small', 
                                       gene = "NSs")
subsNssSmithburn <- subsNssSmithburn + theme(legend.position = "none")
subsNssSmithburn


plots.nss.tstv <- (subsNssMP12 / subsNssClone13 / subsNssSmithburn) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(outDir, "tstv-nss.pdf"), plots.nss.tstv, 
       width = 12, height = 10, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

############################ Small segment: np ################################

inDirNP <- sprintf(
  "%s/output-dir/riftS-NP/strain-types",
  segmentsDir[[3]])

tstvNP <- prepare_tstv_data(segment = "Small", 
                            segment_length = 735,
                            alignment_dir = inDirNP)
subsNpMP12 <- plot_substitutions(data = tstvNP[[1]], 
                                 reference = "MP-12", 
                                 segment = 'Small', 
                                 gene = "NP")
subsNpMP12 <- subsNpMP12 + theme(legend.position = "none")
subsNpMP12


subsNpClone13 <- plot_substitutions(data = tstvNP[[1]], 
                                    reference = "Clone-13", 
                                    segment = 'Small', 
                                    gene = "NP")

subsNpSmithburn <- plot_substitutions(data = tstvNP[[1]], 
                                      reference = "Smithburn", 
                                      segment = 'Small', 
                                      gene = "NP")
subsNpSmithburn <- subsNpSmithburn + theme(legend.position = "none")


plots.np.tstv <- (subsNpMP12 / subsNpClone13 / subsNpSmithburn) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 25))

ggsave(file.path(outDir, "tstv-np.pdf"), 
       plots.np.tstv, 
       width = 12, height = 10, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)





# combine the plots
subsSmithburn <- 
  ((subsLargeSmithburn + 
      labs(title = "Smithburn", tag = "A") + 
      theme(plot.tag = element_text(size = 15)))
   / (subsMediumSmithburn + labs(title = "")) 
   / (subsNpSmithburn + labs(title = "")) 
   / (subsNssSmithburn + labs(title = ""))) +
  plot_layout(axis_titles = "collect")


subsSmithburn


subsMP12 <- 
  ((subsLargeMP12 + 
      labs(title = "MP-12", subtitle = "RdRp", tag = "B") + 
      theme(plot.tag = element_text(size = 15),
            axis.title.y = element_blank(),
            axis.text = element_blank())) 
   / (subsMediumMP12 +
     labs(title = "", subtitle = "NSm/Gn/Gc") +
     theme(axis.title.y = element_blank(),
           axis.text = element_blank())) 
   / (subsNpMP12 +
        labs(title = "", subtitle = "NP") +
        theme(axis.title.y = element_blank(),
              axis.text = element_blank())) 
   / (subsNssMP12 +
        labs(title = "", subtitle = "NSs") +
        theme(axis.title.y = element_blank(),
              axis.text = element_blank()))) +
  plot_layout(axis_titles = "collect")

subsMP12


subsClone13 <- 
  ((subsLargeClone13 + 
      labs(title = "Clone-13", subtitle = "RdRp", tag = "C") + 
      theme(plot.tag = element_text(size = 15),
            axis.title.y = element_blank(),
            axis.text = element_blank())) 
   / (subsMediumClone13 +
        labs(title = "", subtitle = "NSm/Gn/Gc") +
        theme(axis.title.y = element_blank(),
              legend.position = "none",
              axis.text = element_blank())) 
   / (subsNpClone13 +
        labs(title = "", subtitle = "NP") +
        theme(axis.title.y = element_blank(),
              axis.text = element_blank(),
              legend.position = "none")) 
   / (subsNssClone13 +
        labs(title = "", subtitle = "NSs") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none"))) +
  plot_layout(axis_titles = "collect")

subsClone13 



combined.all.tstv.plots <- 
  (subsSmithburn | subsMP12 | subsClone13) +
  plot_layout(guides = 'collect', byrow = F) &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0, 1)
    )

combined.all.tstv.plots

ggsave(file.path(outDir, "Figure2.pdf"), 
       combined.all.tstv.plots,
       width = 18, height = 10, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(outDir, "Figure2.png"), 
       combined.all.tstv.plots,
       width = 18, height = 10, 
       units = "in", limitsize = FALSE,
       dpi = 600, bg="white", device = "png")

