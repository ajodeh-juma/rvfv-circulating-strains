#!/usr/bin/env Rscript

library(colorspace)
library(dplyr)
library(forcats)
library(scales)
library(ggh4x)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggtext)
library(ggthemes)
library(glue)
library(lemon)
library(lubridate)
library(paletteer)
library(patchwork)
library(RColorBrewer)
library(rprojroot)
library(rphylopic)
library(tidyverse)
library(tidytext)

# exports
source(sprintf("%s/utils.R", "scripts"))

#################### glycosylation plots ######################

# Function to print site conservation counts across genomes
print_site_conservation <- function(dataset, segment_name = "Segment") {
  # Total unique genomes evaluated in the dataset
  total_genomes <- dplyr::n_distinct(dataset$Accession)
  
  # Overall conservation per amino acid position
  overall_conservation <- dataset %>%
    dplyr::filter(Result %in% c('+', '++')) %>%
    dplyr::group_by(Position) %>%
    dplyr::summarise(
      Genomes_With_Site = dplyr::n_distinct(Accession),
      Total_Genomes = total_genomes,
      Percentage = round((Genomes_With_Site / total_genomes) * 100, 2),
      Hosts_Present = paste(sort(unique(Host)), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Position)
  
  cat(glue::glue("\n=== Glycosylation Site Conservation for {segment_name} (Total Genomes = {total_genomes}) ===\n\n"))
  print(as.data.frame(overall_conservation), row.names = FALSE)
  
  # Breakdown by Host and Position
  host_totals <- dataset %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(Host_Total = dplyr::n_distinct(Accession), .groups = "drop")
  
  host_conservation <- dataset %>%
    dplyr::filter(Result %in% c('+', '++')) %>%
    dplyr::group_by(Position, Host) %>%
    dplyr::summarise(Genomes_With_Site = dplyr::n_distinct(Accession), .groups = "drop") %>%
    dplyr::left_join(host_totals, by = "Host") %>%
    dplyr::mutate(Percentage = round((Genomes_With_Site / Host_Total) * 100, 2)) %>%
    dplyr::arrange(Position, Host)
  
  cat(glue::glue("\n=== Host-Specific Breakdown for {segment_name} ===\n\n"))
  print(as.data.frame(host_conservation), row.names = FALSE)
  
  return(list(overall = overall_conservation, host = host_conservation))
}


prepare_netglyc_data <- function(segment, indir, gene){
  fns <- file.path(indir, list.files(path = indir, 
                                     pattern = "\\.parsed.netglyc.csv$", 
                                     recursive = TRUE))
  data_list <- list()
  for (i in 1:length(fns)) {
    fn = fns[i]
    sample_name <- str_to_title(strsplit(basename(fn), ".", fixed=T)[[1]][3])
    df <- read.csv(fn, sep=',', stringsAsFactors = F, header = T)
    data_list[[i]] <- df
  }
  data <- bind_rows(data_list)
  data$Host <- str_to_title(data$Host)
  cols <- c("Position", "Potential", "Agreement", "Result")
  data[cols] <- lapply(data[cols], factor)
  
  prefix <- gsub("[^[:alnum:]]", "", gene)
  write.csv(data, file.path(outDir, paste0(segment, "_", prefix, "_netglyc.csv")), 
            quote = F, row.names = F)
  return(data)
}

# Load M segment data only
netglyc.m <- prepare_netglyc_data(segment = "M",
                                  indir = sprintf("%s/output-dir/riftM/netglyc",
                                                  segmentsDir[[2]]),
                                  gene = "NSm/Gn/Gc")

n_host_cols <- length(unique(netglyc.m$Host))

host.cols <- data.frame(Host = sort(unique(netglyc.m$Host)), 
                        Color = paletteer_d("tvthemes::Alexandrite")[1:n_host_cols])

pot.cols <- data.frame(Potential = sort(unique(netglyc.m$Potential)), 
                       Color = divergingx_hcl(24, "Geyser"))

plot_netglyc <- function(dataset, min_genomes, segment, gene){
  
  summary <- dataset %>% 
    dplyr::group_by(Host,Accession) %>% 
    dplyr::summarise(Count=n()) %>% 
    dplyr::group_by(Host) %>% 
    dplyr::summarise(Genomes=n())
  print(summary)
  
  dataset.potential <- dataset[dataset$Result %in% c('+', '++'),] %>%
    dplyr::group_by(Position, Host, Potential) %>%
    dplyr::summarise(Score=mean(Jury), Genomes=n())
  print(dataset.potential)
  
  print(n = 100, dataset |>
          group_by(Potential) |>
          summarise(n=n()))
  
  # Get only events with > 1 genome
  dt <- dataset.potential[dataset.potential$Genomes > min_genomes,]
  
  host.colors <- host.cols[host.cols$Host %in% sort(unique(dt$Host)),]$Color
  
  p <- ggplot(dt, 
              aes(x=Position, y=Genomes, size=Genomes, fill=Host)) +
    geom_point(alpha=1.0, shape=21) +
    scale_size(range = c(1, 10)) +
    scale_fill_manual(name = "Host",
                      values = host.colors,
                      breaks = sort(unique(dt$Host)),
                      na.value = NA,
                      na.translate = F,
                      guide = guide_legend(override.aes = list(linetype = 0,
                                                               shape = 12,
                                                               color = host.colors))) +
    coord_flip() +
    labs(
      subtitle = paste(gene),
      x = "Amino acid position",
      y = "No. of genomes",
      size = "No. of genomes"
    ) +
    theme_dviz_open() +
    theme(
      text = element_text(family = dviz_font_family),
      plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
      panel.grid.minor.x = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3)
    ) +
    guides(fill = guide_legend(override.aes = list(size=12)))
  print(p)
  return(p)
}

plot_netglyc_scores <- function(dataset, min_genomes, subtitle){
  
  summary <- dataset %>%
    dplyr::group_by(Host,Accession) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(Genomes=n())
  
  dataset.potential <- dataset[dataset$Result %in% c('+', '++'),] %>%
    dplyr::group_by(Position, Host, Potential) %>%
    dplyr::summarise(Score=mean(Jury), Genomes=n())
  
  dataset <- dataset %>% mutate(Species = case_when(
    Host == 'Antelope' ~ 'Antidorcas marsupialis',
    Host == "Bat" ~ 'Chiroptera',
    Host == "Buffalo" ~ 'Syncerus caffer',
    Host == 'Cow' ~ 'Bos taurus',
    Host == 'Human' ~ 'Homo sapiens',
    Host == 'Mosquito' ~ 'Aedes aegypti',
    Host == 'Sheep' ~ 'Ovis aries',
    TRUE ~ Host
  ))
  
  dt <- dataset.potential[dataset.potential$Genomes > min_genomes,]
  
  positions <- sort(unique(dt$Position))
  data <- dataset[dataset$Position %in% positions,]
  data$Host <- as.factor(data$Host)
  
  img_dir <- file.path(outDir, "img")
  paths <- list.files(path = img_dir, pattern = "*.png", full.names = TRUE)
  df.paths <- data.frame(Path = paths)
  df.paths$Host <- sub(pattern = "\\..*$", replacement = "", basename(df.paths$Path))
  
  counts_df <- data %>%
    dplyr::group_by(Host,Accession) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(Genomes=n()) %>%
    mutate(genomes = paste0(" (n = ", Genomes, ")"))
  
  df <- left_join(data, df.paths, by="Host")
  df <- left_join(df, counts_df, by="Host")
  
  df$label <- glue("<img src='{df$Path}' width='30' /><br>*{df$Host}{df$genomes}*")
  labels <- unique(df$label)
  names(labels) <- unique(df$Host)
  
  df <- df[df$Genomes > min_genomes & df$Host %in% c('Buffalo', 'Cow', 'Human', 'Sheep'),]
  
  html_labeller <- function(variable, value) {
    return(labels[unique(df$Host)])
  }
  
  host.colors <- host.cols[host.cols$Host %in% sort(unique(df$Host)),]$Color
  
  p <- ggplot(data = df, aes(x = Jury, y = Position, fill = Host)) +
    geom_density_ridges2(
      scale = 1.5,
      bandwidth = 0.1,
      alpha = 0.8
    ) +
    facet_wrap(~ Host, 
               labeller = as_labeller(html_labeller), 
               scales = "free_x") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "#545454") +
    scale_fill_manual(name = "Host",
                      values = host.colors,
                      breaks = names(host.colors)) +
    labs(
      subtitle = subtitle,
      x = "N-glycosylation Score",
      y = "Amino acid position"
    ) +
    theme_ridges() +
    theme_dviz_hgrid() +
    theme(
      text = element_text(family = dviz_font_family),
      plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
      strip.text = element_text(margin = margin(7, 7, 3, 7), colour = "black"),
      strip.text.x = element_markdown(colour = "black"),
      panel.spacing.x = grid::unit(14, "pt"),
      plot.margin = margin(3.5, 1.5, 3.5, 1.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 8),
                                 title = "Host"))
  print(p)
  return(p)
}

# Generate Medium segment plots
res.m <- plot_netglyc(dataset = netglyc.m, 
                      min_genomes = 1,
                      segment = "Medium", 
                      gene = "Glycoproteins (NSm/Gn/Gc)")

res.m.scores <- plot_netglyc_scores(
  dataset = netglyc.m, 
  min_genomes = 1,
  subtitle = "Segment: Medium | Gene: NSm/Gn/Gc")


# Run standalone on Medium segment data
m_conservation <- print_site_conservation(netglyc.m, segment_name = "Medium (NSm/Gn/Gc)")



# Combine res.m.scores and res.m side-by-side
figure5 <- (res.m.scores | res.m) +
  plot_layout(guides = 'collect', 
              widths = c(0.65, 0.35)) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0, 1)
  )

# Export single-segment composite plot
ggsave(file.path(outDir, "Figure5_revised.pdf"), 
       figure5, width = 12, height = 7, 
       dpi = 300, device = cairo_pdf)

ggsave(file.path(outDir, "Figure5_revised.png"), 
       figure5,
       width = 12, height = 7, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")


# #!/usr/bin/env Rscript
# 
# 
# library(colorspace)
# library(dplyr)
# library(forcats)
# library(scales)
# library(ggh4x)
# library(ggnewscale)
# library(ggplot2)
# library(ggrepel)
# library(ggridges)
# library(ggtext)
# library(ggthemes)
# library(glue)
# library(lemon)
# library(lubridate)
# library(paletteer)
# library(patchwork)
# library(RColorBrewer)
# library(rprojroot)
# library(rphylopic)
# library(scales)
# library(tidyverse)
# library(tidytext)
# 
# # exports
# source(sprintf("%s/utils.R", "scripts"))
# 
# 
# #################### glycosylation plots ######################
# 
# prepare_netglyc_data <- function(segment, indir, gene){
#   
#   fns <- file.path(indir, list.files(path = indir, 
#                                      pattern = "\\.parsed.netglyc.csv$", 
#                                      recursive = TRUE))
#   data_list <- list()
#   prefixes <- list()
#   for (i in 1:length(fns)) {
#     fn = fns[i]
#     sample_name <- str_to_title(strsplit(basename(fn), ".", fixed=T)[[1]][3])
#     df <- read.csv(fn, sep=',', stringsAsFactors = F, header = T)
#     # df$Host <- sample_name
#     data_list[[i]] <- df
#   }
#   data <- bind_rows(data_list)
#   data$Host <- str_to_title(data$Host)
#   cols <- c("Position", "Potential", "Agreement", "Result")
#   data[cols] <- lapply(data[cols], factor)
#   
#   prefix <- gsub("[^[:alnum:]]", "", gene)
#   write.csv(data, file.path(outDir, paste0(segment, "_", prefix, "_netglyc.csv")), 
#             quote = F, row.names = F)
#   return(data)
# }
# 
# 
# 
# netglyc.l <- prepare_netglyc_data(segment = "L",
#                                   indir = sprintf("%s/output-dir/riftL/netglyc",
#                                                   segmentsDir[[1]]),
#                                   gene = "RdRp")
# 
# netglyc.m <- prepare_netglyc_data(segment = "M",
#                                   indir = sprintf("%s/output-dir/riftM/netglyc",
#                                                   segmentsDir[[2]]),
#                                   gene = "NSm/Gn/Gc")
# 
# netglyc.np <- prepare_netglyc_data(segment = "NP",
#                                    indir = sprintf("%s/output-dir/riftS-NP/netglyc",
#                                                    segmentsDir[[3]]),
#                                    gene = "NP")
# 
# netglyc.nss <- prepare_netglyc_data(segment = "NSS",
#                                    indir = sprintf("%s/output-dir/riftS-NSS/netglyc",
#                                                    segmentsDir[[3]]),
#                                    gene = "NSS")
# 
# n_host_cols <- max(length(unique(netglyc.l$Host)), 
#                    length(unique(netglyc.m$Host)), 
#                    length(unique(netglyc.np$Host)))
# 
# host.cols <- data.frame(Host=sort(unique((c(unique(netglyc.l$Host),
#                                             unique(netglyc.m$Host),
#                                             unique(netglyc.np$Host))))), 
#                         Color=paletteer_d("tvthemes::Alexandrite")[1:n_host_cols])
# 
# 
# 
# pot.cols <- data.frame(Potential=sort(unique((c(unique(netglyc.l$Potential),
#                                                 unique(netglyc.m$Potential),
#                                                 unique(netglyc.np$Potential))))), 
#                        Color=divergingx_hcl(24, "Geyser"))
# 
# 
# plot_netglyc <- function(dataset, min_genomes, segment, gene){
#   
#   summary <- dataset %>% 
#     dplyr::group_by(Host,Accession) %>% 
#     dplyr::summarise(Count=n()) %>% 
#     dplyr::group_by(Host) %>% 
#     dplyr::summarise(Genomes=n())
#   print(summary)
#   
#   dataset.potential <- dataset[dataset$Result %in% c('+', '++'),] %>%
#     dplyr::group_by(Position, Host, Potential) %>%
#     dplyr::summarise(Score=mean(Jury), Genomes=n())
#   print(dataset.potential)
#   
#   print(n = 100, dataset |>
#           group_by(Potential) |>
#           summarise(n=n()))
#   
#   
#   # Get only events with > 1 genome
#   dt <- dataset.potential[dataset.potential$Genomes > min_genomes,]
#   
#   
#   host.colors <- host.cols[host.cols$Host %in% sort(unique(dt$Host)),]$Color
#   
#   p <- ggplot(dt, 
#               aes(x=Position, y=Genomes, size=Genomes, fill=Host)) +
#     geom_point(alpha=1.0, shape=21) +
#     scale_size(range = c(1, 10)) +
#     scale_fill_manual(name = "Host",
#                       values = host.colors,
#                       breaks = sort(unique(dt$Host)),
#                       na.value = NA,
#                       na.translate = F,
#                       guide = guide_legend(override.aes = list(linetype = 0,
#                                                                shape = 12,
#                                                                color = host.colors))) +
#     coord_flip() +
#     labs(
#       # title = paste("N-glycosylated sites in", segment, "segment", sep = " "),
#       subtitle = paste(gene),
#       x = "Amino acid position",
#       y = "No. of genomes",
#       size = "No. of genomes"
#     ) +
#     theme_dviz_open() +
#     theme(
#       text = element_text(family = dviz_font_family),
#       plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
#       panel.grid.minor.x = element_line(color = "grey80", linewidth = 0.3),
#       panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3)
#       ) +
#     guides(fill = guide_legend(override.aes = list(size=12)))
#   print(p)
#   return(p)
# }
# 
# 
# res.m <- plot_netglyc(dataset = netglyc.m, 
#                       min_genomes = 1,
#                       segment = "Medium", 
#                       gene = "Glycoproteins (NSm/Gn/Gc)")
# 
# res.np <- plot_netglyc(dataset = netglyc.np, 
#                        min_genomes = 1,
#                        segment = "Small", 
#                        gene = "Nucleoprotein (NP)")
# 
# res.l <- plot_netglyc(dataset = netglyc.l,
#                       min_genomes = 1,
#                       segment = "Large", 
#                       gene = "RNA-dependent RNA polymerase (RdRp)")
# 
# 
# plot_netglyc_scores <- function(dataset, min_genomes, subtitle){
#   
#   # summarise the data
#   summary <- dataset %>%
#     dplyr::group_by(Host,Accession) %>%
#     dplyr::summarise(Count=n()) %>%
#     dplyr::group_by(Host) %>%
#     dplyr::summarise(Genomes=n())
#   
#   # select only results with positive outcomes
#   dataset.potential <- dataset[dataset$Result %in% c('+', '++'),] %>%
#     dplyr::group_by(Position, Host, Potential) %>%
#     dplyr::summarise(Score=mean(Jury), Genomes=n())
#   
#   # add species names
#   dataset <- dataset %>% mutate(Species =  case_when(
#     Host == 'Antelope' ~ 'Antidorcas marsupialis',
#     Host == "Bat" ~ 'Chiroptera',
#     Host == "Buffalo" ~ 'Syncerus caffer',
#     Host == 'Cow' ~ 'Bos taurus',
#     Host == 'Human' ~ 'Homo sapiens',
#     Host == 'Mosquito' ~ 'Aedes aegypti',
#     Host == 'Sheep' ~ 'Ovis aries',
#     TRUE ~ Host
#   ))
#   
#   # Get only events with > 1 genome
#   dt <- dataset.potential[dataset.potential$Genomes > min_genomes,]
#   
#   # get unique positions with scores and subset
#   positions <- sort(unique(dt$Position))
#   data <- dataset[dataset$Position %in% positions,]
#   data$Host <- as.factor(data$Host)
#   
#   # create a dataframe with path and host columns
#   img_dir <- file.path(outDir, "img")
#   paths <- list.files(path = img_dir, pattern = "*.png", full.names = TRUE)
#   df.paths <- data.frame(Path = paths)
#   df.paths$Host <- sub(pattern = "\\..*$", replacement = "", basename(df.paths$Path))
#   
#   # count number of genomes per host
#   counts_df <- data %>%
#     dplyr::group_by(Host,Accession) %>%
#     dplyr::summarise(Count=n()) %>%
#     dplyr::group_by(Host) %>%
#     dplyr::summarise(Genomes=n()) %>%
#     mutate(genomes = paste0(" (n = ", Genomes, ")"))
#   
#   
#   # merge dataframes and create a named vector of host and image sources
#   df <- left_join(data, df.paths, by="Host")
#   df <- left_join(df, counts_df, by="Host")
#   
#   df$label <- glue("<img src='{df$Path}' width='30' /><br>*{df$Host}{df$genomes}*")
#   labels <- unique(df$label)
#   names(labels) <- unique(df$Host)
#   
#   df <- df[df$Genomes > min_genomes & df$Host %in% c('Buffalo', 'Cow', 'Human', 'Sheep'),]
#   
#   # custom labeller function to return the HTML strings
#   html_labeller <- function(variable, value) {
#     # 'value' will contain the factor levels (Bat, Cow, Human)
#     # 'image_labels' is the named vector
#     return(labels[unique(df$Host)])
#   }
#   
#   
#   # define host colours
#   host.colors <- host.cols[host.cols$Host %in% sort(unique(df$Host)),]$Color
#   
#   p <- ggplot(data = df, aes(x = Jury, y = Position, fill = Host)) +
#     geom_density_ridges2(#aes(height = after_stat(count)),
#       #stat = "density",
#       scale = 1.5,
#       bandwidth = 0.1,
#       alpha = 0.8
#     ) +
#     facet_wrap(~ Host, 
#                labeller = as_labeller(html_labeller), 
#                scales = "free_x") +
#     geom_vline(xintercept = 0.5, linetype = "dashed", color = "#545454") +
#     scale_fill_manual(name = "Host",
#                       values = host.colors,
#                       breaks = names(host.colors)) +
#     labs(
#       # title = "Distribution of N-glycosylation scores by hosts",
#       subtitle = subtitle,
#       x = "N-glycosylation Score",
#       y = "Amino acid position"
#     ) +
#     theme_ridges() +
#     theme_dviz_hgrid() +
#     theme(
#       text = element_text(family = dviz_font_family),
#       plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
#       strip.text = element_text(margin = margin(7, 7, 3, 7), colour = "black"),
#       strip.text.x = element_markdown(colour = "black"),
#       panel.spacing.x = grid::unit(14, "pt"),
#       plot.margin = margin(3.5, 1.5, 3.5, 1.5)) + 
#     guides(colour = guide_legend(override.aes = list(size = 8),
#                                  title = "Host"))
#   print(p)
#   return(p)
# }
# 
# 
# res.l.scores <- plot_netglyc_scores(
#   dataset = netglyc.l, 
#   min_genomes = 1,
#   subtitle = "Segment: Large | Gene: RdRp")
# 
# res.m.scores <- plot_netglyc_scores(
#   dataset = netglyc.m, 
#   min_genomes = 1,
#   subtitle = "Segment: Medium | Gene: NSm/Gn/Gc")
# 
# res.np.scores <- plot_netglyc_scores(
#   dataset = netglyc.np, 
#   min_genomes = 1,
#   subtitle = "Segment: Small | Gene: NP")
# 
# 
# 
# figure5 <- (((res.l.scores + theme(axis.title = element_blank())) | 
#     (res.l + theme(axis.title = element_blank(),
#                    legend.position = "none"))) / 
#     ((res.m.scores + theme(axis.title.x = element_blank()))  | (res.m + theme(axis.title = element_blank()))) / 
#     ((res.np.scores + theme(axis.title.y = element_blank())) | 
#        (res.np + theme(axis.title.y = element_blank(),
#                        legend.position = "none")))) +
#   plot_layout(guides = 'collect', 
#               widths = c(0.75, 0.4, 0.4), 
#               heights = c(0.75, 0.75, 0.5)) +
#   plot_annotation(tag_levels = "A") & 
#   theme(
#     plot.tag = element_text(size = 18, face = "bold"),
#     plot.tag.position = c(0, 1)
#   )
# 
# 
# ggsave(file.path(outDir, "Figure5.pdf"), 
#        figure5, width = 15, height = 16, 
#        dpi = 300, device = cairo_pdf)
# 
# ggsave(file.path(outDir, "Figure5.png"), 
#        figure5,
#        width = 15, height = 16, 
#        units = "in", limitsize = FALSE,
#        dpi = 300, bg="white", device = "png")
# 
# 
# 
