#!/usr/bin/env Rscript

library(ape)
library(dplyr)
library(ggnewscale)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(jsonlite)
library(lubridate)
library(paletteer)
library(phytools)
library(purrr)
library(rprojroot)
library(stringr)
library(patchwork)
library(tidyverse)
library(treeio)




# exports
source(sprintf("%s/utils.R", "scripts"))


# read metadata
prepare_metadata <- function(traits_fn, lineages_fn, tempest_fn) {
  
  # read lineages file
  lineages <- read.csv(lineages_fn, sep=",", header = TRUE)
  lineages <- lineages |> dplyr::rename("accession"="Query")
  
  # read tempest
  tempest <- read.csv(tempest_fn, sep = "\t")
  tempest$accession <- sapply(strsplit(as.character(tempest$tip), "|", fixed = TRUE), `[`, 1)
  
  # read traits
  traits <- read.csv(traits_fn, sep = "\t")
    
    
  # merge
  list_df <- list(traits, lineages, tempest)
  
  metadata <- list_df |> purrr::reduce(left_join, by = "accession")
  metadata$strain[metadata$strain == '13'] <- 'Clone 13'
  
  # rename
  metadata <- metadata |> rename("date"="date.x", "date2"="date.y")

  
  metadata$strain_type <- str_to_title(metadata$strain_type)
  
  return(metadata)
}


segmentsDir <- sprintf("%s/segments/%s/%s", 
                       find_rstudio_root_file(), 
                       c("L", "M", "S"),
                       "complete/global"
)

metadata_large <- prepare_metadata(
  lineages_fn = file.path(segmentsDir[1], "assignment/output-dir/report/lineages.csv"),
  tempest_fn = file.path(segmentsDir[1], "output-dir/riftL/iqtree/root-to-tip.tsv"),
  traits_fn = file.path(segmentsDir[1], "output-dir/riftL/strain-types/riftL.traits.txt")
)

metadata_medium <- prepare_metadata(
  lineages_fn = file.path(segmentsDir[2], "assignment/output-dir/report/lineages.csv"),
  tempest_fn = file.path(segmentsDir[2], "output-dir/riftM/iqtree/root-to-tip.tsv"),
  traits_fn = file.path(segmentsDir[2], "output-dir/riftM/strain-types/riftM.traits.txt")
)

metadata_np <- prepare_metadata(
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv"),
  tempest_fn = file.path(segmentsDir[3], "output-dir/riftS-NP/iqtree/root-to-tip.tsv"),
  traits_fn = file.path(segmentsDir[3], "output-dir/riftS-NP/strain-types/riftS-NP.traits.txt")
)

metadata_nss <- prepare_metadata(
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv"),
  tempest_fn = file.path(segmentsDir[3], "output-dir/riftS-NSS/iqtree/root-to-tip.tsv"),
  traits_fn = file.path(segmentsDir[3], "output-dir/riftS-NSS/strain-types/riftS-NSS.traits.txt")
)


################################ PLOT TREES ####################################
plot_root_to_tip <- function(metadata, column, subtitle, cols, breaks, 
                             tempest_root) {
  
  # regress genetic distances against time
  fit <- lm(distance ~ date2, data = metadata)
  s_fit <- summary(fit)
  
  r2_val <- round(s_fit$r.squared, 3)
  cor_val <- round(sqrt(r2_val), 3)
  n_val <- nrow(metadata)
  
  # Heuristic Residual Mean Squared (RMS)
  # TempEst typically reports RMS = sum(residuals^2) / (n - 2)
  rms_val <- round(sum(s_fit$residuals^2) / (n_val - 2), 6)
  
  # labelling
  x_min <- min(metadata$date2, na.rm = TRUE)
  x_max <- max(metadata$date2, na.rm = TRUE)
  y_max <- max(metadata$distance, na.rm = TRUE)
  y_min <- min(metadata$distance, na.rm = TRUE)
  
  # get coordinates for the lables
  # x_pos is 2% into the plot; labels are stacked downwards from the top
  x_pos <- x_min + (as.numeric(x_max - x_min) * 0.04)
  v_step <- round((y_max - y_min) * 0.08, 4)
  column_sym <- ensym(column)
  
  # plot
  plot <- ggplot(metadata, aes(date2, distance)) +
    
    geom_point(aes(fill = !!column_sym), size = 6, shape = 21, 
               stroke = 0.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "#000000", 
                fullrange = TRUE, linewidth = 0.8) +
    # Correct X-axis scaling
    # scale_x_date(
    #   name = "Sampling Date (Year)",
    #   date_labels = "%Y",
    #   # breaks = breaks,
    #   breaks = seq(as.Date("1890-01-01"), as.Date("2022-01-01"), by = "10 years"),
    #   #limits = c(as.Date(date_decimal(x_min)), as.Date(date_decimal(x_max))),    # Correctly captures the 1890-2022 span
    #   #expand = expansion(mult = c(0, 0.05))
    # ) +
    # Expand Y-axis by 25% at the top to fit labels
    # scale_y_continuous(name = "Root-to-tip Distance",
    #                    expand = expansion(mult = c(0.05, 0.25))) +
    scale_y_continuous(limits=c(0, y_max + 0.008)) +
    scale_fill_manual(values = cols) +
    labs(title = "",
         subtitle = subtitle,
         x = "Sampling Date (Year)",
         y = "Root-to-tip Distance") +
    
    # annotate
    annotate(geom = "text", x = x_pos, y = y_max + (v_step * 1.5), 
             label = paste0("n = ", n_val), 
             hjust = 0, size = 3.5, family = dviz_font_family) +
    annotate(geom = "text", x = x_pos, y = y_max + (v_step * 0.5), 
             label = paste0("Root (TempEst) = ", tempest_root), 
             hjust = 0, size = 3.5, family = dviz_font_family, color = "#000000") +
    # annotate(geom = "text", x = x_pos, y = y_max + (v_step * 0.5), 
    #          label = paste0("RMS = ", rms_val), 
    #          hjust = 0, size = 5, family = dviz_font_family) +
    annotate(geom = "text", x = x_pos, y = y_max - (v_step * 0.5), 
             label = paste0("Corr.coeff = ", cor_val), 
             hjust = 0, size = 3.5, family = dviz_font_family) +
    annotate(geom = "text", x = x_pos, y = y_max - (v_step * 1.5), 
             label = paste("italic(R)^{2} == ", r2_val), parse = TRUE, 
             hjust = 0, size = 3.5, family = dviz_font_family) +
    theme_dviz_open() +
    theme(
      text = element_text(family = dviz_font_family),
      # plot.margin = margin(3.5, 1.5, 3.5, 1.5),
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.3),
      legend.position = c(0.75, 0.25)
    ) +
    guides(fill = guide_legend(title="Type", ncol = 1, override.aes = list(size = 6)))
  
  return(plot)
}

# Example of how to call the function with a TempEst intercept (e.g., 1890.5)
p1 <- plot_root_to_tip(metadata = metadata_large,
                       column = "strain_type",
                       subtitle = "Temporal signal of RVFV Large segment (2012-2025)",
                       cols = c("#A1ACC8FF",  "#D76E9AFF"),
                       breaks = "10 years",
                       tempest_root = 1891.269)
p1


p2 <- plot_root_to_tip(metadata = metadata_medium,
                       column = "strain_type",
                       subtitle = "Temporal signal of RVFV Medium segment (2012-2025)",
                       cols = c("#A1ACC8FF",  "#D76E9AFF"),
                       breaks = "10 years",
                       tempest_root = 1912.47)
p2


p3 <- plot_root_to_tip(metadata = metadata_np,
                       column = "strain_type",
                       subtitle = "Temporal signal of RVFV S segment (NP) (2012-2025)",
                       cols = c("#A1ACC8FF",  "#D76E9AFF"),
                       breaks = "10 years",
                       tempest_root = 1919.87)
p3

p4 <- plot_root_to_tip(metadata = metadata_nss,
                       column = "strain_type",
                       subtitle = "Temporal signal of RVFV S segment (NSs) (2012-2025)",
                       cols = c("#A1ACC8FF",  "#D76E9AFF"),
                       breaks = "10 years",
                       tempest_root = 1895.32)
p4

# color schemes 
colors <- c("#9370DB", "#87CEFF", "#FF0000", "#551A8B", "#000033", 
            "#FF00B6", "#7AC5CD", "#FFD300", "#009FFF", "#8B1C62", 
            "#00FFBE", "#8B5A2B", "#1F9698", "#FFACFD", "#B1CC71")

colors <- c("tomato", "mediumslateblue", "hotpink4", "cadetblue3", 
            "chartreuse4", "burlywood4",  "greenyellow",
            "orange", "aquamarine1", "limegreen", "hotpink1", "navyblue",
            "darksalmon", "lightseagreen", "yellow")



df <- data.frame(lineage=toupper(letters[1:15]), color=colors)

plot_tree <- function(tree_fn, metadata, segment, cols){
  
  
  # get mrsd
  mrsd <- metadata[rev(order(as.Date(metadata$date2, format = "%Y-%m-%d"))),]$date2[1]
  
  # read tree
  tree <- read.tree(tree_fn)
  print(tree$tip.label)
  
  # match taxa names
  metadata$taxa <- metadata$tip
  metadata <- metadata[match(tree$tip.label, metadata$taxa), ]
  print(all(tree$tip.label == metadata$taxa))
  
  
  # make dataframe for clade nodes
  clades.df <- data.frame(clade=unique(metadata$lineage), node=NA)
  
  # find the most recent common ancestor for each clade
  for (i in 1:length(clades.df$clade)) {
    clades.df$node[i] <- MRCA(
      tree,
      metadata$taxa[metadata$lineage == clades.df$clade[i]]
    )
  }
  
  # plot preliminary tree
  tre <- ggtree(tree,
                as.Date = FALSE,
                linetype=NA) %<+% metadata +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=clade),
                   alpha=1,
                   align="right",
                   extend=0.1,
                   show.legend=FALSE) +
    geom_tree(linewidth=0.5) +
    geom_tiplab(aes(label=accession), size=4)
  tre
  
  tre$data
  
  # order clades dataframe to match the tree
  clades.df <- clades.df[match(tre$data |>
                                 dplyr::filter(isTip == "TRUE") |>
                                 dplyr::arrange(y) |>
                                 dplyr::pull(lineage) |>
                                 unique(),
                               clades.df$clade),]
  # add a column with alternating binary value
  clades.df$highlight <- rep(c(0,1), length.out=length(clades.df$clade))
  clades.df$highlight <- as.factor(clades.df$highlight)
  
  
  
  # define tips to label
  # tips_to_label <- metadata[metadata$strain_type == 'vaccine' & metadata$strain %in% c("Clone 13", "MP-12", "Smithburn"),]$taxa
  # print(tips_to_label)
  
  p <-
    ggtree(midpoint.root(tree), linewidth = 0.3,
           as.Date = FALSE,
           color = "#000000") %<+% metadata +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=as.factor(highlight)),
                   alpha = 1,
                   align = "right",
                   extend = 0.01,
                   show.legend = FALSE
    ) +
    vexpand(.05, direction = 1) +
    geom_tree(linewidth = 0.5) +
    geom_tippoint(aes(color = lineage), alpha = 1.0, size = 6) +
    scale_color_manual(name = "Lineage",
                       values = cols,
                       breaks = sort(unique(metadata$lineage))) +
    geom_tiplab(aes(label = ""),
                # subset = label %in% tips_to_label
                align = TRUE,
                linesize = 0.5,
                geom = "text",
                linetype = "dotted",
                nudge_y = 0.05,
                check_overlap = FALSE,
                offset = 0.001,
                size = 1.5,
                alpha = 1,
                fontface = "bold"
    ) +
    labs(subtitle = paste("Segment: ", segment),
         x = "",
         y = "") +
    
    scale_fill_manual(values = c("#F5F5F5", "#ECECEC")) +
    theme_tree2() +
    theme(
      text = element_text(family = dviz_font_family),
      plot.tag = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 14, colour = "#000000", hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
      axis.title = element_text(size = 14, colour = "#000000"),
      axis.text = element_text(size = 14, colour = "#000000"),
      legend.title = element_text(size = 14, color="#000000"),
      legend.text = element_text(size = 14, color="#000000"),
    ) +
    guides(colour = guide_legend(nrow = length(sort(unique(metadata$lineage))),
                                 override.aes = list(size = 6),
                                 title = "Lineage", order = 1)) +
    new_scale_fill() +
    geom_fruit(
      geom = geom_tile,
      mapping = aes(fill = str_to_title(strain_type)),
      alpha = 1.0,
      width = 0.015,
      color = "white",
      offset = 0.4
    ) +
    scale_fill_manual(
      name = "Type",
      values = c("#A1ACC8FF",  "#D76E9AFF"),
      na.translate = FALSE,
      guide = guide_legend(nrow = 2,
                           override.aes = list(size = 6),
                           order=2
      )
    ) +
    new_scale_fill()
  
  print(p)
  return(p)
  
}


plot_tree_scaled <- function(tree_fn, metadata, segment, cols) {
  

  mrsd <- metadata[rev(order(as.Date(metadata$date, format = "%Y-%m-%d"))),]$date[1]
  print(mrsd)

  tree <- read.nexus(tree_fn) # read.beast is preferred for nexus/beast outputs
  print(tree$data)

  # match taxa names
  metadata$taxa <- metadata$tip
  metadata <- metadata[match(tree$tip.label, metadata$taxa), ]
  print(all(tree$tip.label == metadata$taxa))

  clades.df <- data.frame(clade = unique(metadata$lineage), node = NA)
  for (i in 1:nrow(clades.df)) {
    clades.df$node[i] <- MRCA(tree, metadata$taxa[metadata$lineage == clades.df$clade[i]])
  }
  print(clades.df)

  # plot preliminary tree
  tre <- ggtree(tree,
                mrsd = mrsd,
                as.Date = TRUE,
                linetype=NA) %<+% metadata +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=clade),
                   alpha=1,
                   align="right",
                   extend=0.1,
                   show.legend=FALSE) +
    geom_tree(linewidth=0.5) +
    geom_tiplab(aes(label=accession), size=4)
  tre

  # order clades dataframe to match the tree
  clades.df <- clades.df[match(tre$data |>
                                 dplyr::filter(isTip == "TRUE") |>
                                 dplyr::arrange(y) |>
                                 dplyr::pull(lineage) |>
                                 unique(),
                               clades.df$clade),]
  # add a column with alternating binary value
  clades.df$highlight <- rep(c(0,1), length.out=length(clades.df$clade))
  clades.df$highlight <- as.factor(clades.df$highlight)


  # extract annotations
  strain_type <- metadata |>
    dplyr::select(c("taxa", "strain_type")) %>%
    remove_rownames %>%
    column_to_rownames(var="taxa")
  print(strain_type)

  p <- ggtree(tree, mrsd = mrsd, linewidth = 0.3) %<+% metadata +
    # geom_highlight(data = clades.df,
    #                aes(node = node, fill = clade),
    #                alpha = 0.1, show.legend = FALSE) +
    vexpand(.05, direction = 1) +
    hexpand(.2, direction = 1) +
    geom_tree(linewidth = 0.5) +
    geom_tippoint(aes(color = lineage), alpha = 1.0, size = 6) +
    scale_color_manual(name = "Lineage",
                       values = cols,
                       breaks = sort(unique(metadata$lineage))) +
    geom_tiplab(aes(label = paste(country, "")),
                # subset = label %in% tips_to_label
                align = TRUE,
                linesize = 0.5,
                geom = "text",
                linetype = "dotted",
                nudge_y = 0.05,
                check_overlap = FALSE,
                offset = 3.5,
                size = 1.0,
                alpha = 1,
                fontface = "bold"
    ) +
    # geom_label2(aes(label = label,
    #                 subset = !isTip &              # Ignore tip labels
    #                   !is.na(as.numeric(label)) & # Ensure it's a number
    #                   as.numeric(label) > 80),    # Your threshold
    #             size = 2.5,                        # Smaller size for internal nodes
    #             fill = "white", alpha = 0.7,
    #             label.padding = unit(0.1, "lines"),
    #             family = dviz_font_family,
    #             hjust = 1.2, vjust = -0.5) +
    scale_color_manual(values = cols) +
    scale_x_continuous(labels = function(x) round(x, 0)) +
    # scale_x_date(date_labels = "%Y", breaks="10 year") +
    theme_tree2() +
    labs(subtitle = paste("Segment:", segment), x = "Date") +
    theme(
      text = element_text(family = dviz_font_family),
      plot.tag = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, colour = "#000000", hjust = 0.05),
      plot.subtitle = element_text(size = 14, colour = "#000000", hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
      axis.text = element_text(size = 14, colour = "#000000"),
      legend.title = element_text(size = 14, color="#000000"),
      legend.text = element_text(size = 14, color="#000000"),
    ) +
    guides(colour = guide_legend(nrow = length(sort(unique(metadata$lineage))),
                                 override.aes = list(size = 6),
                                 title = "Lineage", order = 1))
  print(p)

  p <- p + ggnewscale::new_scale_fill()

  p1 <- gheatmap(p + theme(legend.position = "none"),
                 strain_type,
                 offset = 5,
                 width = 0.1,
                 colnames_angle = 0,
                 colnames_position = "top",
                 font.size = 10,
                 family = dviz_font_family,
                 legend_title="Strain Type",
                 colnames_offset_y = 0.5,
                 custom_column_labels = "") +
    scale_fill_manual(
      name = "Strain Type",
      values = c("#A1ACC8FF", "#D76E9AFF")
    )

  print(p1)
  
  return(p1)
}



cols <- df[df$lineage %in% sort(unique(metadata_large$lineage)),]$color
n_colors <- length(unique(sort(metadata_large$strain_type)))
strain_type_custom_colors <- rev(paletteer_d("NatParksPalettes::LakeNakuru")[1:n_colors])

tree_fn <- file.path(segmentsDir[1], "output-dir/riftL/iqtree/riftL.treefile")
treeL <- plot_tree(tree_fn = tree_fn,
                   metadata = metadata_large,
                   segment = "Large",
                   cols = cols)

tree_nexus <- file.path(segmentsDir[1], "output-dir/riftL/treetime/riftL_timetree.nexus")

tree_nexus_L <- plot_tree_scaled(
  tree_fn = tree_nexus,
  metadata = metadata_large,
  segment = "Large",
  cols = cols
)



cols <- df[df$lineage %in% sort(unique(metadata_medium$lineage)),]$color
n_colors <- length(unique(sort(metadata_medium$strain_type)))
strain_type_custom_colors <- rev(paletteer_d("NatParksPalettes::LakeNakuru")[1:n_colors])

tree_fn <- file.path(segmentsDir[2], "output-dir/riftM/iqtree/riftM.treefile")
treeM <- plot_tree(tree_fn = tree_fn,
                   metadata = metadata_medium,
                   segment = "Medium",
                   cols = cols)

tree_nexus <- file.path(segmentsDir[2], "output-dir/riftM/treetime/riftM_timetree.nexus")
tree_nexus_M <- plot_tree_scaled(tree_fn = tree_nexus,
                                 metadata = metadata_medium,
                                 segment = "Medium",
                                 cols = cols)



ggsave(file.path(outDir, "treeM.pdf"),
       treeM, width = 8, height = 10, 
       dpi = 300, device = cairo_pdf)



cols <- df[df$lineage %in% sort(unique(metadata_np$lineage)),]$color
n_colors <- length(unique(sort(metadata_np$strain_type)))
strain_type_custom_colors <- rev(paletteer_d("NatParksPalettes::LakeNakuru")[1:n_colors])

tree_fn <- file.path(segmentsDir[3], "output-dir/riftS-NP/iqtree/riftS-NP.treefile")
treeNP <- plot_tree(tree_fn = tree_fn,
                    metadata = metadata_np,
                    segment = "Small (NP)",
                    cols = cols)

tree_nexus <- file.path(segmentsDir[3], "output-dir/riftS-NP/treetime/riftS-NP_timetree.nexus")
tree_nexus_NP <- plot_tree_scaled(tree_fn = tree_nexus,
                                 metadata = metadata_np,
                                 segment = "Small (NP)",
                                 cols = cols)


cols <- df[df$lineage %in% sort(unique(metadata_nss$lineage)),]$color
n_colors <- length(unique(sort(metadata_nss$strain_type)))
strain_type_custom_colors <- rev(paletteer_d("NatParksPalettes::LakeNakuru")[1:n_colors])

tree_fn <- file.path(segmentsDir[3], "output-dir/riftS-NSS/iqtree/riftS-NSS.treefile")
treeNSS <- plot_tree(tree_fn = tree_fn,
                     metadata = metadata_nss,
                     segment = "Small (NSs)",
                     cols = cols)

tree_nexus <- file.path(segmentsDir[3], "output-dir/riftS-NSS/treetime/riftS-NSS_timetree.nexus")
tree_nexus_NSS <- plot_tree_scaled(tree_fn = tree_nexus,
                                  metadata = metadata_nss,
                                  segment = "Small (NSs)",
                                  cols = cols)

# save plots

figure1 <- ((tree_nexus_L | tree_nexus_M) / (tree_nexus_NP | tree_nexus_NSS) | ((p1 + theme(axis.title.x = element_blank(),
                                                                axis.text.x = element_blank()))
                                                    /(p2 + theme(axis.title.x = element_blank(),
                                                                 axis.text.x = element_blank()))
                                                    /(p3 + theme(axis.title.x = element_blank(),
                                                                 axis.text.x = element_blank()))
                                                    /p4)) + 
  plot_annotation(tag_levels = 'A',
                  title = "") +
  plot_layout(guides = "collect",
              widths = c(2, 1),
              heights = c(1, 0.35)) & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(outDir, "Figure1.pdf"), 
       figure1, width = 16, height = 24, 
       dpi = 300, device = cairo_pdf)

ggsave(file.path(outDir, "Figure1.png"), 
       figure1,
       width = 16, height = 24, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

################################################################################

