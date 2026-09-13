#!/usr/bin/env Rscript

library(colorspace)
library(cowplot)
library(dplyr)
library(extrafont)
library(forcats)
library(scales)
library(ggh4x)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(gt)
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
library(showtext)
library(tidyverse)
library(tidytext)


# exports
source(sprintf("%s/utils.R", "scripts"))


####################### mutation maps #########################

get_host_icons <- function(hosts, img_dir = "icons") {
  if (!dir.exists(img_dir)) dir.create(img_dir)
  
  # Map hosts to UUIDs and Scientific Names
  host_map <- data.frame(
    host = c("Antelope", "Bat", "Buffalo", "Cow", "Human", "Mosquito", "Sheep", "Goat"),
    species = c("Antidorcas marsupialis", "Chiroptera", "Syncerus caffer", "Bos taurus", 
                "Homo sapiens", "Aedes aegypti", "Ovis aries", "Capra hircus"),
    uuid = c("320dcfd5-c738-484b-ac28-dab0ac07139b", "00e0ac1c-9e1e-4a30-aac7-5f768ff915ae", 
             "65c4a9b3-dcde-4f0f-9a1f-8d71e74be9ec", "dc5c561e-e030-444d-ba22-3d427b60e58a", 
             "c089caae-43ef-4e4e-bf26-973dd4cb65c5", "48207547-694e-4696-90d9-b61b2084abcc", 
             "14012740-e356-464b-b84d-fa6c4a0b34a9", "c39a1d31-eb36-4b62-ba4d-32e29e0e5a8d")
  )
  
  host_data <- host_map %>% filter(host %in% str_to_title(hosts))
  
  # Download images if they don't exist
  host_data <- host_data %>%
    mutate(path = file.path(img_dir, paste0(host, ".png"))) %>%
    rowwise() %>%
    mutate(img_status = {
      if (!file.exists(path)) {
        img <- get_phylopic(uuid = uuid)
        save_phylopic(img, path = path)
      }
      path
    }) %>%
    ungroup()
  
  return(host_data)
}



codeMutationTags<- function(count, thre){
  over<- names(count[which(as.numeric(count)>=thre)])
  # print(over)
  over<- strsplit(over, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)
  over<- matrix(unlist(over),ncol=3,byrow=T)
  positions<- over[,2]
  over<- cbind(over[,2], ": ", over[,1], " to ", over[,3])
  over<- apply(over, MARG=1, FUN=function(X) { paste(X, collapse='') })
  over<- data.frame(x_pos=as.numeric(positions), y_pos=thre/2, tag=as.character(over), stringsAsFactors=FALSE)
  if (dim(over)[1] >= 1){
    return(over)
  }
}

plot_mutation_map <- function(fpath, threshold=round(90/100*245), 
                              lineage, aa_size, step, subs, 
                              palette, cols, subtitle, selected_mutations){
  
  mutation_df <- read.csv(fpath, header = TRUE)
  df.mut.s <- mutation_df %>% 
    dplyr::group_by(host,snps) %>% 
    dplyr::summarise(count=n()) %>%
    dplyr::arrange(desc(count)) %>% as.data.frame()
  
  lineage_df <- mutation_df %>% 
    dplyr::select(snps, positions, host) %>% 
    dplyr::rename("mutation" = snps, "location" = positions)
  
  
  tags <- codeMutationTags(table(lineage_df$mutation), thre=threshold)
  
  if (dim(tags)[1] >= 1) {
    tags <- tags[order(tags$x_pos,decreasing=FALSE),]
    table <- table(lineage_df$location)
    table <- data.frame(count=as.numeric(table),location=as.numeric(names(table)))
    to_label <- table[table$location %in% tags$x_pos,]
    to_label <- merge(to_label, tags, by.x = "location", by.y = "x_pos")
    to_label$tag <- paste0(sapply(strsplit(gsub(":", "", to_label$tag), " "), `[`, 2),
                           sapply(strsplit(gsub(":", "", to_label$tag), " "), `[`, 1),
                           sapply(strsplit(gsub(":", "", to_label$tag), " "), `[`, 4))
    
    rows <- dim(to_label)[1]
    
    d <- df.mut.s[df.mut.s$snps %in% to_label$tag,]
    d$position <- gsub("\\D+", "", d$snps)
    d$position <- as.numeric(d$position)
    d$snps <- as.factor(d$snps)
    d$host <- str_to_title(d$host)
    
    d <- d %>% mutate(species =  case_when(
      host == 'Antelope' ~ 'Antidorcas marsupialis',
      host == "Bat" ~ 'Chiroptera',
      host == "Buffalo" ~ 'Syncerus caffer',
      host == 'Cow' ~ 'Bos taurus',
      host == 'Human' ~ 'Homo sapiens',
      host == 'Mosquito' ~ 'Aedes aegypti',
      host == 'Sheep' ~ 'Ovis aries',
      host == 'Goat' ~ 'Capra hircus',
      TRUE ~ host
    )) %>% mutate(uuid = case_when(
      host == 'Antelope' ~ '320dcfd5-c738-484b-ac28-dab0ac07139b',
      host == 'Bat' ~ '00e0ac1c-9e1e-4a30-aac7-5f768ff915ae',
      host == 'Buffalo' ~ '65c4a9b3-dcde-4f0f-9a1f-8d71e74be9ec',
      host == 'Cow' ~ 'dc5c561e-e030-444d-ba22-3d427b60e58a',
      host == 'Human' ~ 'c089caae-43ef-4e4e-bf26-973dd4cb65c5',
      host == 'Mosquito' ~ '48207547-694e-4696-90d9-b61b2084abcc',
      host == 'Sheep' ~ '14012740-e356-464b-b84d-fa6c4a0b34a9',
      host == 'Goat' ~ 'c39a1d31-eb36-4b62-ba4d-32e29e0e5a8d'
    ))
    print(d)
    
    # d$url <- sprintf("https://images.phylopic.org/images/%s/vector.svg", d$uuid)
    # d$label <- paste0("<img src='", d$url, "' width='10' />")
    
    # download the images to local computer
    img_dir <- file.path(outDir, "img")
    if (!dir.exists(img_dir)){dir.create(img_dir)}
    
    for (uuid in unique(d$uuid)){
      dest_file <- paste0(unique(d[d$uuid == uuid, ]$host), ".png")
      dest_path <- file.path(img_dir, dest_file)
      print(dest_path)
      
      if (!file.exists(dest_path)) {
        phylopic_obj <- get_phylopic(uuid = uuid, format="vector", preview = FALSE)
        save_phylopic(img = phylopic_obj, path = dest_path, bg = "white")
      }
    }
    
    # create a dataframe with path and host columns
    paths <- list.files(path = img_dir, pattern = "*.png", full.names = TRUE)
    df.paths <- data.frame(path = paths)
    df.paths$host <- sub(pattern = "\\..*$", replacement = "", basename(df.paths$path))
    
    # merge dataframes
    data <- left_join(d, df.paths, by="host")
    data$label <- paste0("<img src='", data$path, "' width='30' /><br>", data$host)
    
    data$label <- glue("<img src='{data$path}' width='30' /><br>*{data$host}*")
    names(data$label) <- data$species
    data <- data %>% arrange(position, species)
    
    # 
    plots_per_page = min(20, length(to_label$tag))
    
    if ( length(to_label$tag) %% 2 != 0 & length(to_label$tag) == 1  ) {
      nrow <- 1
    }
    
    if ( length(to_label$tag) %% 2 == 0 & length(to_label$tag)/2 <= 1) {
      nrow <- plots_per_page/2
    }
    
    if ( length(to_label$tag) %% 2 == 0 & length(to_label$tag)/2 > 1) {
      nrow <- plots_per_page/2
    }
    
    if ( length(to_label$tag) %% 2 != 0 & length(to_label$tag) != 1) {
      nrow <- plots_per_page/2
    }
    
    if (nrow == 1){
      nrow <- 1
    } else if (nrow > 5){
      nrow <- 4
    } else {
      nrow <- nrow - 0.9
      nrow <- round(nrow, 0)
    }
    
    # plots
    data <- data[data$host %in% c('Cow', 'Human', 'Sheep', 'Buffalo', 'Goat', 'Mosquito'),]
    
    p <- ggplot(data = data, aes(x = species, y = count, colour = snps, size = count)) +
      geom_point(alpha=1) +
      facet_wrap(~ fct_reorder(snps, position), nrow = nrow, scales = "fixed") +
      scale_colour_manual(name="Substitution", 
                          values=cols,
                          breaks = to_label$tag,
                          guide=guide_legend(order=1,
                                             override.aes=list(size=12))) +
      scale_x_discrete(labels = data$label) +
      scale_size_continuous(range = c(1, 10)) +
      labs(
        # title = paste("Substitutions by vaccine strain:", lineage, sep = " "),
        subtitle = paste("Vaccine: ", lineage, "|", subtitle),
        x = "Species",
        y = "No. of genomes",
        size = "No. of genomes"
      ) +
      theme_dviz_hgrid() +
      theme(
        text = element_text(family = dviz_font_family),
        axis.ticks = element_line(colour = "grey80", linewidth = 0.3),
        axis.text.x = element_markdown(),
        strip.text = element_text(margin = margin(7, 7, 3, 7), colour = "black"),
        panel.spacing.x = grid::unit(14, "pt"),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = margin(3.5, 1.5, 3.5, 1.5)) + 
      guides(colour = guide_legend(override.aes = list(size = 12),
                                   title = "Substitution", order = 1),
             fill = guide_legend(override.aes = list(size = 12), order = 2))
    
    # save plot
    prefix <- paste0(subtitle, "-", lineage)
    prefix <- gsub("[^[:alnum:]]", "", prefix)
    outfile <- file.path(outDir, paste0(prefix, "-plot-1.pdf"))
    
    if (nrow == 1 & plots_per_page <= 3) {
      width <- plots_per_page * 8.0
      height <- nrow * 5.0
    } else {
      width <- nrow * 8.0
      height <- nrow * 5.0
    }
    print(p)
    
    # ggsave(outfile, p,
    #        width = width, height = height, units = "in",
    #        limitsize = FALSE,
    #        dpi = 300, bg="white", device = cairo_pdf)
    
    p1 <- ggplot(table, aes(x = location, y=count, size = count)) +
      geom_bar(data=to_label, 
               aes(x=location, y=count, fill=tag), 
               stat="identity", width = 2.5, alpha=0.8) +
      scale_x_continuous(limits = c(0, aa_size),breaks = seq(0, aa_size, by = step)) +
      geom_hline(yintercept=threshold, linetype="dashed", color = "#545454") +
      geom_point(data=to_label, 
                 aes(x = location, y = count,  size = count), 
                 alpha = 1
      ) +
      scale_size_continuous(range = c(1, 8)) +
      geom_label_repel(data = to_label,
                       aes(label=tag, fill=tag),
                       size = 5,
                       colour = "#000000",
                       fontface = "plain",
                       nudge_x = 0.25,
                       nudge_y = 0.25,
                       show.legend = FALSE,
                       max.overlaps = 20,
                       angle = 0) +
      scale_fill_manual(name="Substitution", 
                        values=cols,
                        breaks = to_label$tag,
                        guide=guide_legend(order=1,
                                           override.aes=list(size=12),
                                           keywidth = 5,
                                           keyheight = 5)
      ) +
      labs(
        # title = paste("Substitutions by vaccine strain:", lineage, sep = " "),
        subtitle = paste("Vaccine: ", lineage, "|", subtitle),
        x = "Amino acid position",
        y = "No. of genomes",
        size = "No. of genomes"
      ) +
      theme_dviz_open() +
      theme(
        text = element_text(family = dviz_font_family),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        axis.title = element_markdown(face = "bold", size = 15),
      ) +
      guides(colour = guide_legend(title = "No. of genomes", 
                                   override.aes = list(size = 12),
                                   order = 2),
             fill = guide_legend(title = "Substitution",
                                 override.aes = list(size = 12),
                                 order = 1)
      )
    print(p1)
    # save plot
    # prefix <- paste0(subtitle, "-", lineage)
    # prefix <- gsub("[^[:alnum:]]", "", prefix)
    # outfile <- file.path(outDir, paste0(prefix, "-plot-2.pdf"))
    # 
    # ggsave(outfile, p1,
    #        width = width, height = height, units = "in",
    #        limitsize = FALSE,
    #        dpi = 300, bg="white", device = cairo_pdf)
    
    
    # plot selected mutations
    data_selected_mutations <- data[data$snps %in% c('A172V', 'I1244M', 'H259Y',
                                                     'G1182R', 'S278N', 'D95N',
                                                     'I602V', 'V659A'),]
    
    if ( dim(data_selected_mutations)[1] >= 1 ){
      p2 <- ggplot(data = data_selected_mutations, 
                   aes(x = species, y = count, colour = snps, size = count)) +
        geom_point(alpha=1) +
        facet_wrap(~ fct_reorder(snps, position), nrow = nrow, scales = "fixed") +
        scale_colour_manual(name="Substitution", 
                            values=cols,
                            breaks = to_label$tag,
                            guide=guide_legend(order=1,
                                               override.aes=list(size=12))) +
        scale_x_discrete(labels = data_selected_mutations$label) +
        scale_size_continuous(range = c(1, 10)) +
        labs(
          # title = paste("Substitutions by vaccine strain:", lineage, sep = " "),
          subtitle = paste("Vaccine: ", lineage, "|", subtitle),
          x = "Species",
          y = "No. of genomes",
          size = "No. of genomes"
        ) +
        theme_dviz_hgrid() +
        theme(
          text = element_text(family = dviz_font_family),
          axis.ticks = element_line(colour = "grey80", linewidth = 0.3),
          axis.text.x = element_markdown(),
          strip.text = element_text(margin = margin(7, 7, 3, 7), colour = "black"),
          panel.spacing.x = grid::unit(14, "pt"),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.margin = margin(3.5, 1.5, 3.5, 1.5)) + 
        guides(colour = guide_legend(override.aes = list(size = 12),
                                     title = "Substitution", order = 1),
               fill = guide_legend(override.aes = list(size = 12), order = 2))
      
      # save plot
      # prefix <- paste0(subtitle, "-", lineage)
      # prefix <- gsub("[^[:alnum:]]", "", prefix)
      # outfile <- file.path(outDir, paste0(prefix, "-plot-3.pdf"))
      # 
      # ggsave(outfile, p2,
      #        width = width, height = height, units = "in",
      #        limitsize = FALSE,
      #        dpi = 300, bg="white", device = cairo_pdf)
      # 
    } else {
      p2 <- NULL
    }
    print(p2)
    
    l <- list(lineage_df, p, p1, p2)
  }
  return(l)
}

############################### Large segment ##################################

largeSingletonMP12File <- sprintf(
  "%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375404.mutations.per.strain.singleton.csv",
  segmentsDir[1])

largeSingletonClone13File <- sprintf(
  "%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375417.mutations.per.strain.singleton.csv",
  segmentsDir[1])

largeSingletonSmithburnFile <- sprintf(
  "%s/output-dir/riftL/strain-types/riftL.strain_type.DQ375430.mutations.per.strain.singleton.csv",
  segmentsDir[1])


# threshold <- round(90/100*282)
# ncol = 2, nrow = 3
largeMutMapMP12 <- plot_mutation_map(fpath = largeSingletonMP12File,
                                     threshold = round(90/100*245),
                                     lineage = "MP-12",
                                     aa_size = 2092, 
                                     step = 200,
                                     palette = "ag_GrnYl",
                                     cols = paletteer_d("ggsci::category20_d3"),
                                     subtitle = "Segment: Large | Gene: RdRp")


# nrow = 1, ncol = 1
largeMutMapClone13 <- plot_mutation_map(fpath = largeSingletonClone13File,
                                        threshold = round(80/100*245),
                                        lineage = "Clone 13",
                                        aa_size = 2092, 
                                        step = 200,
                                        palette = "YlGnBu",
                                        cols = paletteer_d("ggsci::category20_d3"),
                                        subtitle = "Segment: Large | Gene: RdRp")

largeMutMapSmithburn <- plot_mutation_map(fpath = largeSingletonSmithburnFile,
                                          threshold = round(90/100*245),
                                          lineage = "Smithburn",
                                          aa_size = 2092, 
                                          step = 200,
                                          palette = "Batlow",
                                          cols = paletteer_d("ggsci::category20_d3"),
                                          subtitle = "Segment: Large | Gene: RdRp")

############################### Medium segment #################################

mediumSingletonMP12File <- sprintf(
  "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv",
  segmentsDir[[2]])

mediumSingletonClone13File <- sprintf(
  "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380213.mutations.per.strain.singleton.csv",
  segmentsDir[[2]])

mediumSingletonSmithburnFile <- sprintf(
  "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380193.mutations.per.strain.singleton.csv",
  segmentsDir[[2]])



mediumMutMapMP12 <- plot_mutation_map(fpath = mediumSingletonMP12File,
                                      threshold = round(90/100*245),
                                      lineage = "MP-12",
                                      aa_size = 1197, 
                                      step = 200,
                                      palette = "ag_GrnYl",
                                      cols = c(paletteer_d("ggsci::category10_d3"), paletteer_d("peRReo::nicky")),
                                      subtitle = "Glycoproteins (NSm/Gn/Gc)")

mediumMutMapClone13 <- plot_mutation_map(fpath = mediumSingletonClone13File,
                                         threshold = round(90/100*245),
                                         lineage = "Clone 13",
                                         aa_size = 1197, 
                                         step = 200,
                                         palette = "YlGnBu",
                                         cols = c(paletteer_d("ggsci::category10_d3"), paletteer_d("peRReo::nicky")),
                                         subtitle = "Glycoproteins (NSm/Gn/Gc)")

mediumMutMapSmithburn <- plot_mutation_map(fpath = mediumSingletonSmithburnFile,
                                           threshold = round(90/100*245),
                                           lineage = "Smithburn",
                                           aa_size = 1197, 
                                           step = 200,
                                           palette = "Batlow",
                                           cols = c(paletteer_d("ggsci::category10_d3"), paletteer_d("peRReo::nicky")),
                                           subtitle = "Glycoproteins (NSm/Gn/Gc)")

############################# Small segment: NSs ###############################

nssSingletonMP12File <- sprintf(
  "%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380154.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

nssSingletonClone13File <- sprintf(
  "%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380182.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

nssSingletonSmithburnFile <- sprintf(
  "%s/output-dir/riftS-NSS/strain-types/riftS-NSS.strain_type.DQ380157.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

nssMutMapMP12 <- plot_mutation_map(fpath = nssSingletonMP12File,
                                   threshold = round(90/100*245),
                                   lineage = "MP-12",
                                   aa_size = 245, 
                                   step = 50,
                                   palette = "ag_GrnYl",
                                   cols = paletteer_d("Polychrome::dark"),
                                   subtitle = "Non-structural protein (NSs)")

# nssMutMapClone13 <- plot_mutation_map(fpath = nssSingletonClone13File,
#                                       threshold = round(90/100*245),
#                                       lineage = "Clone 13",
#                                       aa_size = 245, 
#                                       step = 50,
#                                       palette = "YlGnBu",
#                                       cols = paletteer_d("Polychrome::dark"),
#                                       subtitle = "Non-structural protein (NSs)")

nssMutMapSmithburn <- plot_mutation_map(fpath = nssSingletonSmithburnFile,
                                        threshold = round(90/100*245),
                                        lineage = "Smithburn",
                                        aa_size = 245, 
                                        step = 50,
                                        palette = "Batlow",
                                        cols = paletteer_d("Polychrome::dark"),
                                        subtitle = "Non-structural protein (NSs)")

############################## Small segment: NP ###############################

npSingletonMP12File <- sprintf(
  "%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380154.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

npSingletonClone13File <- sprintf(
  "%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380182.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

npSingletonSmithburnFile <- sprintf(
  "%s/output-dir/riftS-NP/strain-types/riftS-NP.strain_type.DQ380157.mutations.per.strain.singleton.csv",
  segmentsDir[[3]])

npMutMapMP12 <- plot_mutation_map(fpath = npSingletonMP12File,
                                  threshold = round(90/100*245),
                                  lineage = "MP-12",
                                  aa_size = 265, 
                                  step = 50,
                                  palette = "ag_GrnYl",
                                  cols = rev(paletteer_d("khroma::soil")),
                                  subtitle = "Nucleoprotein (NP)")

npMutMapClone13 <- plot_mutation_map(fpath = npSingletonClone13File,
                                     threshold = round(5/100*245),
                                     lineage = "Clone 13",
                                     aa_size = 265, 
                                     step = 50,
                                     palette = "YlGnBu",
                                     cols = rev(paletteer_d("khroma::soil")),
                                     subtitle = "Nucleoprotein (NP)")

npMutMapSmithburn <- plot_mutation_map(fpath = npSingletonSmithburnFile,
                                       threshold = round(90/100*245),
                                       lineage = "Smithburn",
                                       aa_size = 265, 
                                       step = 50,
                                       palette = "Batlow",
                                       cols = rev(paletteer_d("khroma::soil")),
                                       subtitle = "Nucleoprotein (NP)")

p <- (mediumMutMapSmithburn[[2]] / mediumMutMapMP12[[2]]  / mediumMutMapClone13[[2]]) +
  plot_layout(guides = 'keep', axis_titles = "collect_x", 
              widths = c(0.2, 0.4, 0.2), heights = c(1.5, 2, 0.5)) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0, 1)
  )
p

# ggsave(file.path(outDir, "mutation-plots.pdf"), p, 
#        width = 27, height = 21, units = "in", 
#        limitsize = FALSE,
#        dpi = 300, bg="white", device = cairo_pdf)



get_common_substitutions <- function(file, snps_count, reference){
  
  df <- read.csv(file)
  df$host <- str_to_title(df$host)
  df1 <- df %>%
    dplyr::group_by(host, snps, positions) %>%
    dplyr::summarise(count=n())
  
  df2 <- df %>% dplyr::group_by(host, snps, positions) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::arrange(desc(count))
  
  data_list <- list()
  df3 <- NULL
  for (i in 1:length(unique(df2$host))){
    df3$host <- unique(df2$host)[i]
    df3$total <- df2[df2$host == unique(df2$host)[i],]$count[[1]]
    data_list[[i]] <- df3
  }
  dat <- bind_rows(data_list)
  
  # get abundant snps
  data <- df1[df1$count > snps_count, ]
  data <- data[order(data$positions),]
  data <- dplyr::left_join(data, dat, by="host")
  data$percentage <- data$count/data$total*100
  data <- data |> dplyr::mutate_if(is.numeric, round, digits=1)
  data <- data |> dplyr::mutate(perc = paste0(sprintf("%4.1f", count / total * 100), "%"))
  
  d <- data %>%
    dplyr::group_by(snps, positions) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::filter(count >= 3) %>%
    dplyr::arrange(desc(positions))
  d$vaccine <- reference
  return(d)
}

largeSubsMP12 <- get_common_substitutions(
  file = largeSingletonMP12File,
  snps_count = 10,
  reference = "MP-12"
)

largeSubsClone13 <- get_common_substitutions(
  file = largeSingletonClone13File,
  snps_count = 10,
  reference = "Clone-13"
)

largeSubsSmithburn <- get_common_substitutions(
  file = largeSingletonSmithburnFile,
  snps_count = 10,
  reference = "Smithburn")

dplyr::left_join(largeSubsSmithburn,
                 largeSubsMP12,
                 by="snps") |>
  dplyr::left_join(largeSubsClone13, by="snps") |>
  dplyr::filter(!is.na(positions.y))


mediumSubsMP12 <- get_common_substitutions(
  file = mediumSingletonMP12File,
  snps_count = 10,
  reference = "MP-12"
)

mediumSubsClone13 <- get_common_substitutions(
  file = mediumSingletonClone13File,
  snps_count = 10,
  reference = "Clone-13"
)

mediumSubsSmithburn <- get_common_substitutions(
  file = mediumSingletonSmithburnFile,
  snps_count = 10,
  reference = "Smithburn")

dplyr::left_join(mediumSubsSmithburn,
                 mediumSubsMP12,
                 by="snps") |>
  dplyr::left_join(mediumSubsClone13, by="snps") |>
  dplyr::filter(!is.na(positions.y))


nssSubsMP12 <- get_common_substitutions(
  file = nssSingletonMP12File,
  snps_count = 10,
  reference = "MP-12"
)

nssSubsClone13 <- get_common_substitutions(
  file = nssSingletonClone13File,
  snps_count = 10,
  reference = "Clone-13"
)

nssSubsSmithburn <- get_common_substitutions(
  file = nssSingletonSmithburnFile,
  snps_count = 10,
  reference = "Smithburn")

dplyr::left_join(nssSubsSmithburn,
                 nssSubsMP12,
                 by="snps") |>
  dplyr::left_join(nssSubsClone13, by="snps") |>
  dplyr::filter(!is.na(positions.y))

npSubsMP12 <- get_common_substitutions(
  file = npSingletonMP12File,
  snps_count = 1,
  reference = "MP-12"
)

npSubsClone13 <- get_common_substitutions(
  file = npSingletonClone13File,
  snps_count = 1,
  reference = "Clone-13"
)

npSubsSmithburn <- get_common_substitutions(
  file = npSingletonSmithburnFile,
  snps_count = 1,
  reference = "Smithburn")

dplyr::left_join(npSubsSmithburn,
                 npSubsMP12,
                 by="snps") |>
  dplyr::left_join(npSubsClone13, by="snps") |>
  dplyr::filter(!is.na(positions.y))



plot_mutations <- function(file, snps_count, perc_genomes, segment, 
                           reference, palette, subtitle){
  
  df <- read.csv(file)
  df$host <- str_to_title(df$host)
  
  # count mutations per host
  df1 <- df %>% 
    dplyr::group_by(host, snps, positions) %>% 
    dplyr::summarise(count=n()) 
  
  df2 <- df %>% dplyr::group_by(host, snps, positions) %>% 
    dplyr::summarise(count=n()) %>% 
    dplyr::arrange(desc(count))
  
  data_list <- list()
  df3 <- NULL
  for (i in 1:length(unique(df2$host))){
    df3$host <- unique(df2$host)[i]
    df3$total <- df2[df2$host == unique(df2$host)[i],]$count[[1]]
    data_list[[i]] <- df3
  }
  
  dat <- bind_rows(data_list)
  print(dat)
  
  # get abundant snps
  data <- df1[df1$count > snps_count, ]
  data <- data[order(data$positions),]
  data <- dplyr::left_join(data, dat, by="host")
  data$percentage <- data$count/data$total*100
  data <- data |> dplyr::mutate_if(is.numeric, round, digits=1) 
  data <- data |> dplyr::mutate(perc = paste0(sprintf("%4.1f", count / total * 100), "%"))
  
  print(data %>%
          dplyr::group_by(snps, positions) %>%
          dplyr::summarise(count=n()) %>%
          dplyr::filter(count >= 3) %>% 
          dplyr::arrange(desc(positions)), 
        n=100)
  
  # add additional information
  data$vaccine_strain <- reference
  data$segment <- segment
  
  
  # add images
  # add species names
  data <- data %>% mutate(species =  case_when(
    host == 'Antelope' ~ 'Antidorcas marsupialis',
    host == "Bat" ~ 'Chiroptera',
    host == "Buffalo" ~ 'Syncerus caffer',
    host == 'Cow' ~ 'Bos taurus',
    host == 'Human' ~ 'Homo sapiens',
    host == 'Mosquito' ~ 'Aedes aegypti',
    host == 'Sheep' ~ 'Ovis aries',
    TRUE ~ host
  ))
  print(unique(data$host))
  
  data <- data[data$host %in% c('Buffalo', 'Cow', 'Goat', 'Human', 'Mosquito', 'Sheep'),]
  data$host <- as.factor(data$host)
  
  # create a dataframe with path and host columns
  img_dir <- file.path(outDir, "img")
  paths <- list.files(path = img_dir, pattern = "*.png", full.names = TRUE)
  df.paths <- data.frame(path = paths)
  df.paths$host <- sub(pattern = "\\..*$", replacement = "", basename(df.paths$path))
  
  # merge dataframes and create a named vector of host and image sources
  df <- left_join(data, df.paths, by="host")
  
  df$label <- glue("<img src='{df$path}' width='30' /><br>*{df$host}*")
  labels <- unique(df$label)
  names(labels) <- unique(df$host)
  
  df <- df %>% arrange(host, positions)
  
  # custom labeller function to return the HTML strings
  html_labeller <- function(variable, value) {
    # 'value' will contain the factor levels (Bat, Cow, Human)
    # 'image_labels' is the named vector
    return(labels[unique(df$host)])
  }
  
  
  df <- df[df$percentage >= perc_genomes,]
  
  # reorder(snps, sort(as.numeric(positions)))
  p <- ggplot(df, aes(x = count,
                      y = fct_reorder(snps, positions),
                      fill = snps)) +
    geom_col() +
    geom_label(data=df,
               aes(label = perc),
               hjust = 1,
               nudge_x = -.5,
               size = 4,
               # fontface = "bold",
               family = "Arial Narrow",
               ## turn into white box without outline
               fill = "white", label.size = 0
    ) +
    labs(
      # title = paste("Frequency and Distribution of RVFV SNPs in Hosts"),
      subtitle = subtitle,
      x = "No. of genomes",
      y = "Substitutions"
    ) +
    facet_grid(rows =  ~ fct_reorder(host, sort(as.numeric(positions))), 
               labeller = as_labeller(html_labeller),
               scales = "free_x" ) +
    scale_fill_discrete_divergingx(palette = palette, 
                                   rev = TRUE, alpha=1.0) +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    theme_dviz_open() +
    theme(
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
      plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, margin = margin(b = 15)),
      strip.text.x = element_markdown(margin = margin(7, 7, 3, 7), colour = "black"),
      strip.background = element_rect(fill = "white"),
      panel.spacing.x = grid::unit(14, "pt"),
      plot.margin = margin(3.5, 1.5, 3.5, 1.5),
      legend.position = "none"
    ) +
    coord_capped_cart(bottom=capped_horizontal(c("both")),
                      left=capped_vertical(c("both"))) +
    guides(color = "none")
  print(p)
  
  # save
  n_grid <- 3
  n_mutations <- length(unique(df$snps))
  print(n_mutations)
  
  # # save plot
  # prefix <- paste0(gene, "-", reference)
  # prefix <- gsub("[^[:alnum:]]", "", prefix)
  # outfile <- file.path(outDir, paste0(prefix, "-plot-4.pdf"))
  # 
  # if ( n_mutations <= 5 ) {
  #   width <- n_grid * 2.5
  #   height <- n_mutations * 1.75
  # } else {
  #   width <-  n_grid * 3.15
  #   height <- n_mutations * 0.75
  # }
  
  # ggsave(outfile, p,
  #        width = width, height = height, units = "in",
  #        limitsize = FALSE,
  #        dpi = 300, bg="white", device = cairo_pdf)
  # 
  # save the outputs to list variable
  print(p)
  l <- list(df, p)
  return(l)
  
}



plotMutationsMediumMP12 <- plot_mutations(file = mediumSingletonMP12File,
                                          snps_count = 10,
                                          perc_genomes = 90,
                                          segment = "M",
                                          reference = "MP-12",
                                          subtitle = "Vaccine: MP-12 | Segment: Medium | Gene: Glycoproteins (NSm/Gn/Gc)",
                                          palette = "Tropic"
)
plotMutationsMediumClone13 <- plot_mutations(file = mediumSingletonClone13File,
                                             snps_count = 10,
                                             perc_genomes = 90,
                                             segment = "M",
                                             reference = "Clone-13",
                                             subtitle = "Vaccine: Clone-13 | Segment: Medium | Gene: Glycoproteins (NSm/Gn/Gc)",
                                             palette = "RdGy"
)
plotMutationsMediumSmithburn <- plot_mutations(file = mediumSingletonSmithburnFile,
                                               snps_count = 10,
                                               perc_genomes = 90,
                                               segment = "M",
                                               reference = "Smithburn",
                                               subtitle = "Vaccine: Smithburn | Segment: Medium | Gene: Glycoproteins (NSm/Gn/Gc)",
                                               palette = "Zissou 1"
)


plotMutationsLargeMP12 <- plot_mutations(file = largeSingletonMP12File,
                                         snps_count = 10,
                                         perc_genomes = 90,
                                         segment = "L",
                                         reference = "MP-12",
                                         subtitle = "Vaccine: MP-12 | Segment: Large | Gene: RdRp",
                                         palette = "Tropic"
)
plotMutationsLargeClone13 <- plot_mutations(file = largeSingletonClone13File,
                                            snps_count = 10,
                                            perc_genomes = 90,
                                            segment = "L",
                                            reference = "Clone-13",
                                            subtitle = "Vaccine: Clone-13 | Segment: Large | Gene: RdRp",
                                            palette = "RdGy"
)
plotMutationsLargeSmithburn <- plot_mutations(file = largeSingletonSmithburnFile,
                                              snps_count = 10,
                                              perc_genomes = 90,
                                              segment = "L",
                                              reference = "Smithburn",
                                              subtitle = "Vaccine: Smithburn | Segment: Large | Gene: RdRp",
                                              palette = "Zissou 1"
)



plotMutationsNssMP12 <- plot_mutations(file = nssSingletonMP12File,
                                       snps_count = 10,
                                       perc_genomes = 90,
                                       segment = "S",
                                       reference = "MP-12",
                                       subtitle = "Vaccine: MP-12 | Segment: Small | Gene: Non-structural (NSs)",
                                       palette = "Tropic"
)
plotMutationsNssClone13 <- plot_mutations(file = nssSingletonClone13File,
                                          snps_count = 30,
                                          perc_genomes = 100,
                                          segment = "S",
                                          reference = "Clone-13",
                                          subtitle = "Vaccine: Clone-13 | Segment: Small | Gene: Non-structural (NSs)",
                                          palette = "RdGy"
)
plotMutationsNssSmithburn <- plot_mutations(file = nssSingletonSmithburnFile,
                                            snps_count = 10,
                                            perc_genomes = 90,
                                            segment = "S",
                                            reference = "Smithburn",
                                            subtitle = "Vaccine: Smithburn | Segment: Small | Gene: Non-structural (NSs)",
                                            palette = "Zissou 1"
)

plotMutationsNpMP12 <- plot_mutations(file = npSingletonMP12File,
                                      snps_count = 0,
                                      perc_genomes = 10,
                                      segment = "S",
                                      reference = "MP-12",
                                      subtitle = "Vaccine: MP-12 | Segment: Small | Gene: Nucleoprotein (NP)",
                                      palette = "Tropic"
)
plotMutationsNpClone13 <- plot_mutations(file = npSingletonClone13File,
                                         snps_count = 0,
                                         perc_genomes = 90,
                                         segment = "S",
                                         reference = "Clone-13",
                                         subtitle = "Vaccine: Clone-13 | Segment: Small | Gene: Nucleoprotein (NP)",
                                         palette = "RdGy"
)
plotMutationsNpSmithburn <- plot_mutations(file = npSingletonSmithburnFile,
                                           snps_count = 1,
                                           perc_genomes = 90,
                                           segment = "S",
                                           reference = "Smithburn",
                                           subtitle = "Vaccine: Smithburn | Segment: Small | Gene: Nucleoprotein (NP)",
                                           palette = "Zissou 1"
)




snps_count <- 10
perc_genomes <- 90

prepare_mutation_data <- function(file){
  
  df <- read.csv(file)
  
  df$host <- stringr::str_to_title(df$host)
  
  # count mutations per host
  df1 <- df %>%
    group_by(host, snps, positions) %>%
    summarise(count=n(), .groups="drop")
  
  # determine total genomes per host
  totals <- df1 %>%
    group_by(host) %>%
    summarise(total=max(count))
  
  # merge totals
  data <- left_join(df1, totals, by="host")
  
  data <- data %>%
    mutate(
      percentage = count/total*100
    )
  
  # apply filters used in plotting function
  data <- data %>%
    filter(count > snps_count,
           percentage >= perc_genomes)
  
  return(data)
  
}

smithburn <- prepare_mutation_data(mediumSingletonSmithburnFile)
mp12 <- prepare_mutation_data(mediumSingletonMP12File)
clone13 <- prepare_mutation_data(mediumSingletonClone13File)

annotate_domains <- function(df){
  
  df %>%
    # distinct(snps, positions) |>
    mutate(
      Domain = case_when(
        Position <= 154 ~ "NSm",
        Position >=155 & Position <=560 ~ "Gn",
        Position >=561 & Position <=690 ~ "Gn-Gc junction",
        Position >=691 & Position <=1120 ~ "Gc",
        Position >1120 ~ "Gc cytoplasmic tail"
      )
    )
  
}

# generate_table3 <- function(df, reference){
#   
#   df %>%
#     annotate_domains() %>%
#     group_by(domain) %>%
#     summarise(
#       substitutions = n(),
#       .groups="drop"
#     ) %>%
#     mutate(vaccine_reference = reference)
#   
# }


# annotate epitopes
annotate_epitopes <- function(df){
  
  df %>%
    mutate(
      `Epitope_cluster` = case_when(
        Position >=180 & Position <=260 ~ "Cluster I",
        Position >=261 & Position <=350 ~ "Cluster II",
        Position >=351 & Position <=430 ~ "Cluster III",
        TRUE ~ NA_character_
      ),
      
      `Epitope_region` = case_when(
        Position >=180 & Position <=260 ~ "Gn neutralizing epitope (180–260)",
        Position >=261 & Position <=350 ~ "Gn neutralizing epitope (261–350)",
        Position >=351 & Position <=430 ~ "Gn neutralizing epitope (351–430)",
        TRUE ~ NA_character_
      )
    )
  
}


aa_group <- function(aa){
  
  case_when(
    aa %in% c("A","V","L","I","M","F","W","P","G") ~ "Hydrophobic",
    aa %in% c("S","T","C","Y","N","Q") ~ "Polar",
    aa %in% c("K","R","H") ~ "Positive",
    aa %in% c("D","E") ~ "Negative"
  )
  
}

classify_substitution <- function(df){
  
  df %>%
    mutate(
      refAA = substr(snps,1,1),
      altAA = substr(snps,nchar(snps),nchar(snps)),
      Position = positions,
      Substitution = snps
    ) %>%
    mutate(
      ref_group = aa_group(refAA),
      alt_group = aa_group(altAA),
      `Biochemical change` = paste(ref_group,"→",alt_group),
      Classification = ifelse(ref_group == alt_group,
                              "Conservative",
                              "Radical")
    )
}



functional_interpretation <- function(df){
  
  df %>%
    mutate(
      `Functional interpretation` = case_when(
        
        Domain == "Gn" &
          Classification == "Radical" ~
          "Potential antigenic or receptor-binding modification",
        
        Domain == "Gc" &
          Classification == "Radical" ~
          "Possible effect on membrane fusion machinery",
        
        Domain == "NSm" ~
          "Possible effect on virion assembly and morphogenesis",
        
        TRUE ~
          "Likely structurally tolerated substitution"
      )
    )
  
}



mutation_references <- tribble(
  
  ~Substitution, ~Reference,
  
  "I219T", "Ikegami et al. 2015 J Virol",
  "N230D", "Gerrard et al. 2007 Virology",
  "Y249H", "Bird et al. 2011 PLoS Pathogens",
  "E294K", "Ikegami & Makino 2011 J Virol",
  "V328A", "Kortekaas et al. 2012 Vaccine"
)

add_references <- function(df){
  
  df %>%
    left_join(mutation_references, by="Substitution")
  
}


generate_table3 <- function(df, vaccine){
  
  df %>%
    distinct(snps, positions) %>%
    rename(Substitution = snps,
           Position = positions) %>%
    annotate_domains() %>%
    mutate(`Vaccine strain` = vaccine) %>%
    left_join(mutation_references, by="Substitution") %>%
    select(Substitution,
           Position,
           Domain,
           `Vaccine strain`,
           Reference)
}


generate_table4 <- function(df, vaccine){
  
  df %>%
    distinct(snps, positions) %>%
    rename(Substitution = snps,
           Position = positions) %>%
    annotate_domains() %>%
    annotate_epitopes() %>%
    filter(!is.na(Epitope_region)) %>%
    mutate(`Vaccine strain` = vaccine) %>%
    left_join(mutation_references, by="Substitution") %>%
    select(Substitution,
           Position,
           Domain,
           Epitope_region,
           `Vaccine strain`,
           Reference)
  
}

generate_table5 <- function(df, vaccine){
  
  df %>%
    
    # keep unique substitutions
    distinct(snps, positions) %>%
    
    # rename variables
    rename(
      Substitution = snps,
      Position = positions
    ) %>%
    
    # biochemical classification
    mutate(
      refAA = substr(Substitution,1,1),
      altAA = substr(Substitution,nchar(Substitution),nchar(Substitution))
    ) %>%
    
    mutate(
      ref_group = aa_group(refAA),
      alt_group = aa_group(altAA),
      `Biochemical change` = paste(ref_group,"→",alt_group),
      Classification = ifelse(ref_group == alt_group,
                              "Conservative",
                              "Radical")
    ) %>%
    
    # annotate structural domains
    annotate_domains() %>%
    
    # annotate epitope information
    annotate_epitopes() %>%
    
    # functional interpretation
    mutate(
      `Functional interpretation` = case_when(
        
        Domain == "Gn" & Classification == "Radical" ~
          "Possible alteration of antigenic surface or receptor-binding interface",
        
        Domain == "Gc" & Classification == "Radical" ~
          "Potential effect on membrane fusion conformational changes",
        
        Domain == "NSm" ~
          "May influence virion assembly or budding",
        
        TRUE ~
          "Likely structurally tolerated substitution"
      )
    ) %>%
    
    # add vaccine information
    mutate(`Vaccine strain` = vaccine) %>%
    
    # add literature references
    add_references() %>%
    
    rename(`Epitope region`  = Epitope_region,
           `Epitope cluster` = Epitope_cluster) %>%
    
    # reorder columns
    select(
      Substitution,
      Position,
      Domain,
      `Epitope region`,
      `Epitope cluster`,
      `Vaccine strain`,
      `Biochemical change`,
      Classification,
      `Functional interpretation`,
      Reference
    ) %>%
    arrange(`Vaccine strain`, Position)
  
}


tbl3 <- bind_rows(
  generate_table3(smithburn, "Smithburn"),
  generate_table3(mp12, "MP-12"),
  generate_table3(clone13, "Clone 13")
)


tbl4 <- bind_rows(
  generate_table4(smithburn, "Smithburn"),
  generate_table4(mp12, "MP-12"),
  generate_table4(clone13, "Clone 13")
)


tbl5 <- bind_rows(
  generate_table5(smithburn, "Smithburn"),
  generate_table5(mp12, "MP-12"),
  generate_table5(clone13, "Clone 13")
)


tbl5_final <- gt(tbl5) %>%
  cols_label(
    Substitution = "Substitution",
    Position = "Position",
    Domain = "Domain",
    `Vaccine strain` = "Vaccine strain",
    `Biochemical change` = "Biochemical change",
    Classification = "Classification",
    `Functional interpretation` = "Functional interpretation",
    `Epitope region` = "Epitope region",
    `Epitope cluster` = "Epitope cluster",
    Reference = "Reference"
  )

write.table(tbl5_final, file.path(outDir, "Table5_medium.tsv"), 
            sep = "\t", quote = F,
            row.names = F)




# -----------------------------
# 1. Domain annotation for RdRp
# -----------------------------

annotate_polymerase_domains <- function(df){
  
  df %>%
    mutate(
      Domain = case_when(
        Position <= 900 ~ "RdRp N-terminal domain",
        Position >=901 & Position <=1500 ~ "Core polymerase region",
        Position >=1501 & Position <=2000 ~ "Palm domain",
        Position >=2001 & Position <=2600 ~ "Thumb domain",
        Position >=2601 ~ "C-terminal region"
      )
    )
  
}

# --------------------------------
# 2. Polymerase motif annotation
# --------------------------------

annotate_motifs <- function(df){
  
  df %>%
    mutate(
      Motif = case_when(
        Position >=1000 & Position <=1030 ~ "Motif A",
        Position >=1050 & Position <=1080 ~ "Motif B",
        Position >=1100 & Position <=1140 ~ "Motif C",
        Position >=1180 & Position <=1210 ~ "Motif D",
        Position >=1240 & Position <=1280 ~ "Motif E",
        Position >=1320 & Position <=1350 ~ "Motif F",
        TRUE ~ NA_character_
      )
    )
  
}

# --------------------------------
# 3. Amino acid biochemical classes
# --------------------------------

aa_group <- function(aa){
  
  case_when(
    aa %in% c("A","V","L","I","M","F","W","P","G") ~ "Hydrophobic",
    aa %in% c("S","T","C","Y","N","Q") ~ "Polar",
    aa %in% c("K","R","H") ~ "Positive",
    aa %in% c("D","E") ~ "Negative"
  )
  
}

# --------------------------------
# 4. Functional classification
# --------------------------------

classify_substitution <- function(df){
  
  df %>%
    mutate(
      refAA = substr(Substitution,1,1),
      altAA = substr(Substitution,nchar(Substitution),nchar(Substitution))
    ) %>%
    mutate(
      `Biochemical change` =
        paste(aa_group(refAA),"→",aa_group(altAA)),
      Classification =
        ifelse(aa_group(refAA)==aa_group(altAA),
               "Conservative","Radical")
    )
  
}

# --------------------------------
# 5. Functional interpretation
# --------------------------------

interpret_function <- function(df){
  
  df %>%
    mutate(
      `Functional interpretation` = case_when(
        
        !is.na(Motif) & Classification=="Radical" ~
          "Possible effect on polymerase catalytic activity",
        
        !is.na(Motif) & Classification=="Conservative" ~
          "Mutation within conserved polymerase motif",
        
        Domain=="Palm domain" & Classification=="Radical" ~
          "May influence nucleotide incorporation",
        
        Domain=="Thumb domain" ~
          "Possible effect on template positioning",
        
        TRUE ~
          "Likely structurally tolerated substitution"
      )
    )
  
}

# --------------------------------
# 6. Generate final table
# --------------------------------

generate_table_L_segment <- function(df, vaccine){
  
  df %>%
    distinct(snps,positions) %>%
    rename(
      Substitution = snps,
      Position = positions
    ) %>%
    
    annotate_polymerase_domains() %>%
    annotate_motifs() %>%
    classify_substitution() %>%
    interpret_function() %>%
    
    mutate(`Vaccine strain` = vaccine) %>%
    
    select(
      Substitution,
      Position,
      Domain,
      Motif,
      `Vaccine strain`,
      `Biochemical change`,
      Classification,
      `Functional interpretation`
    ) %>%
    
    arrange(Position)
  
}


smithburn_L <- prepare_mutation_data(largeSingletonSmithburnFile)
mp12_L <- prepare_mutation_data(largeSingletonMP12File)
clone13_L <- prepare_mutation_data(largeSingletonClone13File)


smithburn_L %>% 
  distinct(snps,positions) %>%
  rename(
    Substitution = snps,
    Position = positions
  )

table_L <- bind_rows(
  generate_table_L_segment(smithburn_L,"Smithburn"),
  generate_table_L_segment(mp12_L,"MP-12"),
  generate_table_L_segment(clone13_L,"Clone-13")
)

table_L

write.table(table_L, file.path(outDir, "Table5_large.tsv"), 
            sep = "\t", quote = F,
            row.names = F)




# ------------------------------
# NSs domain annotation
# ------------------------------

annotate_nss_domain <- function(position){
  
  case_when(
    position <= 60 ~ "N-terminal regulatory region",
    position >60 & position <=120 ~ "Central oligomerization region",
    position >120 ~ "C-terminal functional region"
  )
}

# ------------------------------
# Process each mutation file
# ------------------------------

process_nss <- function(df, vaccine){
  
  # unique substitutions
  data <- df %>%
    distinct(snps, positions)
  
  # extract amino acids
  data <- data %>%
    mutate(
      refAA = substr(snps,1,1),
      altAA = substr(snps,nchar(snps),nchar(snps))
    )
  
  # biochemical annotation
  data <- data %>%
    mutate(
      ref_class = aa_group(refAA),
      alt_class = aa_group(altAA),
      `Biochemical change` = paste(ref_class,"→",alt_class),
      Classification = ifelse(ref_class==alt_class,
                              "Conservative","Radical")
    )
  
  # domain annotation
  data <- data %>%
    mutate(
      Domain = annotate_nss_domain(positions)
    )
  
  # functional interpretation
  data <- data %>%
    mutate(
      `Functional interpretation` = case_when(
        
        Domain == "N-terminal regulatory region" &
          Classification == "Radical" ~
          "Possible impact on host transcription suppression",
        
        Domain == "Central oligomerization region" ~
          "May influence NSs multimerization",
        
        Domain == "C-terminal functional region" ~
          "Potential effect on host interaction",
        
        TRUE ~
          "Likely structurally tolerated substitution"
      )
    )
  
  data <- data %>%
    mutate(`Vaccine strain` = vaccine) %>%
    rename(
      Substitution = snps,
      Position = positions) %>%
    select(
    Substitution,
    Position,
    Domain,
    `Vaccine strain`,
    `Biochemical change`,
    Classification,
    `Functional interpretation`
    ) %>%
    arrange(Position)
  
  return(data)
}

# ------------------------------
# Run analysis for all strains
# ------------------------------


smithburn_nss <- prepare_mutation_data(nssSingletonSmithburnFile)
mp12_nss <- prepare_mutation_data(nssSingletonMP12File)
clone13_nss <- prepare_mutation_data(nssSingletonClone13File)

tbl_NSS <- bind_rows(
  process_nss(smithburn_nss, "Smithburn"),
  process_nss(mp12_nss, "MP-12"),
  process_nss(clone13_nss, "Clone 13")
)
tbl_NSS

write.table(tbl_NSS, file.path(outDir, "Table5_nss.tsv"), 
            sep = "\t", quote = F,
            row.names = F)




# NP domain annotation
annotate_np_domain <- function(position){
  
  case_when(
    position <= 35 ~ "N-terminal oligomerization arm",
    position >35 & position <=110 ~ "RNA-binding groove",
    position >110 & position <=180 ~ "Core structural domain",
    position >180 & position <=220 ~ "NP oligomerization interface",
    position >220 ~ "C-terminal tail"
  )
}

# functional interpretation
interpret_np_function <- function(domain, classification){
  
  case_when(
    
    domain == "RNA-binding groove" & classification == "Radical" ~
      "Possible alteration of RNA binding affinity",
    
    domain == "N-terminal oligomerization arm" &
      classification == "Radical" ~
      "May affect NP oligomerization",
    
    domain == "NP oligomerization interface" ~
      "Possible effect on ribonucleoprotein assembly",
    
    domain == "Core structural domain" ~
      "Potential impact on NP structural stability",
    
    domain == "C-terminal tail" ~
      "Possible interaction change with viral polymerase",
    
    TRUE ~
      "Likely structurally tolerated substitution"
  )
}

# process NP mutation file
process_np <- function(df, vaccine){
  
  # unique substitutions
  data <- df %>%
    distinct(snps, positions)
  
  # extract amino acids
  data <- data %>%
    mutate(
      refAA = substr(snps,1,1),
      altAA = substr(snps,nchar(snps),nchar(snps))
    )
  
  # biochemical annotation
  data <- data %>%
    mutate(
      ref_class = aa_group(refAA),
      alt_class = aa_group(altAA),
      `Biochemical change` = paste(ref_class,"→",alt_class),
      Classification = ifelse(ref_class==alt_class,
                              "Conservative","Radical")
    )
  
  data <- data %>%
    mutate(
      Domain = annotate_np_domain(positions)
    )
  
  data <- data %>%
    mutate(
      `Functional interpretation` =
        interpret_np_function(Domain, Classification)
    )
  
  data <- data %>%
    mutate(`Vaccine strain` = vaccine) %>%
    rename(
      Substitution = snps,
      Position = positions
    ) %>%
    select(
      Substitution,
      Position,
      Domain,
      `Vaccine strain`,
      `Biochemical change`,
      Classification,
      `Functional interpretation`
    ) %>%
    arrange(`Vaccine strain`, Position)
  
  return(data)
}



smithburn_np <- prepare_mutation_data(npSingletonSmithburnFile)
mp12_np <- prepare_mutation_data(npSingletonMP12File)
clone13_np <- prepare_mutation_data(npSingletonClone13File)

tbl_NP <- bind_rows(
  process_np(smithburn_np, "Smithburn"),
  process_np(mp12_np, "MP-12"),
  #process_np(clone13_np, "Clone 13")
)
tbl_NP

write.table(tbl_NP, file.path(outDir, "Table5_np.tsv"), 
            sep = "\t", quote = F,
            row.names = F)



clean_plot <- function(p) {
  p + theme(
    axis.title.y = element_blank(), # Remove Y labels for inner columns
    axis.title.x = element_blank(), # Remove X labels for inner rows
    strip.text = element_blank(),   # Remove strip text (except for top row)
    legend.position = "none"        # Remove legends
  )
}


# ROW 1: Large
row1 <- (clean_plot(plotMutationsLargeSmithburn[[2]]) | 
           clean_plot(plotMutationsLargeMP12[[2]]) | 
           clean_plot(plotMutationsLargeClone13[[2]]))

# ROW 2: Medium
row2 <- (clean_plot(plotMutationsMediumSmithburn[[2]]) | 
           clean_plot(plotMutationsMediumMP12[[2]]) | 
           clean_plot(plotMutationsMediumClone13[[2]]))

# ROW 3: Nss (Note the spacer_plot for Clone-13)
spacer_plot <- plot_spacer() 
row3 <- (clean_plot(plotMutationsNssSmithburn[[2]]) | 
           clean_plot(plotMutationsNssMP12[[2]]) | 
           spacer_plot)

# ROW 4: Np
row4 <- (clean_plot(plotMutationsNpSmithburn[[2]]) | 
           clean_plot(plotMutationsNpMP12[[2]]) | 
           clean_plot(plotMutationsNpClone13[[2]]))

# We wrap the assembly in an overarching structure
final_figure <- (row1 / row2 / row3 / row4) + 
  plot_layout(widths = c(1.15, 1.3, 1, 0.75),
              heights = c(0.5, 0.5, 0.5, 0.5)) +
  plot_annotation(
    title = "Frequency and Distribution of RVFV SNPs in Hosts",
    tag_levels = "A"
  ) & 
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0, 1)
  )

final_figure

ggsave(file.path(outDir, "Figure4.pdf"), 
       plot = final_figure, 
       width = 28, 
       height = 25, 
       dpi = 300, 
       device = cairo_pdf)

ggsave(file.path(outDir, "Figure4.tiff"), 
       plot = final_figure, 
       width = 28, 
       height = 25, 
       dpi = 300, 
       device = "tiff")

ggsave(file.path(outDir, "Figure4.png"), 
       final_figure,
       width = 28, height = 25, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")



mutations_df <- tbl5

M_segment_gene_annotations <- data.frame(
  gene = as.factor(c("NSm", "Gn", "Gc")),
  start = c(1, 480, 2091) / 3, 
  end = c(479, 2090, 3591) / 3,
  colour = c("#8B8970", "#8B2323", "#458B74")
)

mutations_df <- mutations_df %>%
  mutate(
    Epitope.cluster = ifelse(is.na(`Epitope cluster`) | `Epitope cluster` == "NA", "Non-Epitopic", `Epitope cluster`),
    Biochem_Class = paste(Classification, ":", `Biochemical change`)
  )

p_medium <- ggplot() +
  geom_rect(data = M_segment_gene_annotations, 
            aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.1, fill = colour), 
            alpha = 0.7, color = "black", size = 0.3) +
  geom_text(data = M_segment_gene_annotations, 
            aes(x = (start + end)/2, y = 0, label = gene), 
            color = "white", fontface = "bold") +
  
  geom_segment(data = mutations_df, 
               aes(x = Position, xend = Position, y = 0.1, yend = 0.7, color = Epitope.cluster),
               size = 0.5, linetype = "dotted") +
  geom_point(data = mutations_df, 
             aes(x = Position, y = 0.7, shape = Classification, color = Epitope.cluster), 
             size = 3.5, stroke = 1) +
  
  geom_text_repel(data = mutations_df, 
                  aes(x = Position, y = 0.7, label = Substitution),
                  size = 3, fontface = "bold", box.padding = 0.4) +
  
  facet_wrap(~`Vaccine strain`, ncol = 1, strip.position = "left") +
  scale_fill_identity(guide = "legend", name = "Gene Region", 
                      labels = c("Gc (Fusion)", "Gn (Attachment)", "NSm (Assembly)")) +
  scale_color_brewer(palette = "Set1", name = "Epitope Cluster") +
  scale_shape_manual(values = c("Conservative" = 21, "Radical" = 24), name = "Mutation Type") +
  labs(# title = "Structural and Antigenic Landscape of RVFV Vaccine Substitutions",
       subtitle = "M-Segment",
       x = "Amino Acid Position (M-segment polyprotein)",
       y = "",
       caption = "Positions normalized from nucleotide annotations (1-3591 nt).") +
  theme_dviz_open() +
  theme(
    text = element_text(family = dviz_font_family),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, colour = "#000000", 
                                 hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )

p_medium


ggsave(file.path(outDir, "RVFV_glycoprotein_Mutation_Profile.pdf"), 
       plot = p_medium, width = 15, height = 12, dpi = 300,
       device = cairo_pdf)


L_segment_gene_annotations <- data.frame(
  gene = as.factor(c("RdRp")),
  start = c(1) / 3, 
  end = c(6276) / 3,
  colour = c("#CDBA96")
)


l_mutations_df <- table_L %>%
  mutate(Motif = ifelse(is.na(Motif) | Motif == "NA", "Non-Motif", Motif))

p_large <- ggplot() +
  geom_rect(data = L_segment_gene_annotations, 
            aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.1, fill = colour), 
            alpha = 0.6, color = "black", size = 0.3) +
  geom_text(data = L_segment_gene_annotations, 
            aes(x = (start + end)/2, y = 0, label = gene), 
            color = "black", fontface = "bold", size = 5) +
    geom_segment(data = l_mutations_df, 
               aes(x = Position, xend = Position, y = 0.1, yend = 0.7, color = Motif),
               alpha = 0.5, linetype = "solid") +
  geom_point(data = l_mutations_df, 
             aes(x = Position, y = 0.7, shape = Classification, color = Motif), 
             size = 4, stroke = 1.2) +
    geom_text_repel(data = l_mutations_df, 
                  aes(x = Position, y = 0.7, label = Substitution),
                  box.padding = 0.6, point.padding = 0.4, segment.size = 0.3,
                  size = 3.5, fontface = "bold", max.overlaps = 15) +
  
  facet_wrap(~`Vaccine strain`, ncol = 1, strip.position = "left") +
  
  scale_fill_identity(guide = "legend", name = "Gene Region", 
                      labels = "RdRp (Polymerase)") +
  scale_color_manual(values = c("Non-Motif" = "grey60", "Motif A" = "#E41A1C", 
                                "Motif E" = "#377EB8", "Motif F" = "#4DAF4A"), 
                     name = "Conserved Motifs") +
  scale_shape_manual(values = c("Conservative" = 21, "Radical" = 24), name = "Mutation Type") +
  
  labs(#title = "Structural Distribution of Mutations in the RVFV RdRp (L Segment)",
       subtitle = "L-segment",
       x = "Amino Acid Position",
       y = "",
       caption = "Positions normalized from 6276 nt reference sequence.") +
  theme_dviz_open() +
  theme(
    text = element_text(family = dviz_font_family),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, colour = "#000000", 
                                 hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )
p_large


ggsave(file.path(outDir, "RVFV_RdRp_Mutation_Profile.pdf"), 
       plot = p_large, width = 15, height = 12, dpi = 300,
       device = cairo_pdf)


NSS_segment_gene_annotations <- data.frame(
  gene = as.factor(c("NSs")),
  start = c(1) / 3, 
  end = c(795) / 3,
  colour = c("#CD9B9B")
)

nss_mutations_df <- tbl_NSS

p_nss <- ggplot() +
  geom_rect(data = NSS_segment_gene_annotations, 
            aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.1, fill = colour), 
            alpha = 0.8, color = "black", size = 0.4) +
  geom_text(data = NSS_segment_gene_annotations, 
            aes(x = (start + end)/2, y = 0, label = gene), 
            color = "white", fontface = "bold", size = 5) +
  
  # Mutation Lollipops
  geom_segment(data = nss_mutations_df, 
               aes(x = Position, xend = Position, y = 0.1, yend = 0.8, color = Domain),
               alpha = 0.6, linetype = "dashed", size = 0.6) +
  geom_point(data = nss_mutations_df, 
             aes(x = Position, y = 0.8, 
                 shape = Classification, 
                 color = Domain), 
             size = 4, stroke = 1.2) +
  
  # Labels: Substitution Labels
  geom_text_repel(data = nss_mutations_df, 
                  aes(x = Position, y = 0.8, label = Substitution),
                  box.padding = 0.5, point.padding = 0.3, segment.size = 0.2,
                  size = 3.5, fontface = "bold", max.overlaps = 20) +
  
  # Faceting by Vaccine Strain
  facet_wrap(~`Vaccine strain`, ncol = 1, strip.position = "left") +
  
  # Custom Scales and Themes
  scale_fill_identity(guide = "legend", name = "Protein", labels = "NSs (IFN Antagonist)") +
  scale_color_brewer(palette = "Dark2", name = "Functional Domain") +
  scale_shape_manual(values = c("Conservative" = 21, "Radical" = 24), name = "Mutation Type") +
  
  labs(#title = "Mutation Profile of the RVFV NSs Protein (Small Segment)",
       subtitle = "S-segment (NSs)",
       x = "Amino Acid Position",
       y = "",
       caption = "Positions normalized from 795 nt (NSs ORF). Note: Clone 13 contains a large natural deletion not shown as point mutations.") +
  theme_dviz_open() +
  theme(
    text = element_text(family = dviz_font_family),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, colour = "#000000", 
                                 hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )
p_nss

ggsave(file.path(outDir, "RVFV_NSs_Mutation_Profile.pdf"), 
       plot = p_nss, width = 15, height = 12, dpi = 300,
       device = cairo_pdf)



NP_segment_gene_annotations <- data.frame(
  gene = c("NP"),
  start = c(1), 
  end = c(735) / 3, # Normalizing to Amino Acid scale
  colour = c("#CD853F")
)



np_mutations_df <- tbl_NP

p_np <- ggplot() +
  geom_rect(data = NP_segment_gene_annotations, 
            aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.1, fill = colour), 
            alpha = 0.8, color = "black", size = 0.4) +
  geom_text(data = NP_segment_gene_annotations, 
            aes(x = (start + end)/2, y = 0, label = gene), 
            color = "white", fontface = "bold", size = 5) +
  
  geom_segment(data = np_mutations_df, 
               aes(x = Position, xend = Position, y = 0.1, yend = 0.8, color = Domain),
               alpha = 0.6, linetype = "solid", size = 0.7) +
  geom_point(data = np_mutations_df, 
             aes(x = Position, y = 0.8, 
                 shape = Classification, 
                 color = Domain), 
             size = 4, stroke = 1.2) +
  
  geom_text_repel(data = np_mutations_df, 
                  aes(x = Position, y = 0.8, label = Substitution),
                  box.padding = 0.6, point.padding = 0.4, segment.size = 0.3,
                  size = 3.5, fontface = "bold") +
  
  facet_wrap(~`Vaccine strain`, ncol = 1, strip.position = "left") +
  
  scale_fill_identity(guide = "legend", name = "Protein", labels = "NP (Nucleoprotein)") +
  scale_color_brewer(palette = "Set1", name = "Structural Domain") +
  scale_shape_manual(values = c("Conservative" = 21, "Radical" = 24), name = "Mutation Type") +
  
  labs(#title = "Mutation Distribution in the RVFV Nucleoprotein (NP)",
       subtitle = "S-segment (NP)",
       x = "Amino Acid Position",
       y = "",
       caption = "Positions normalized from 735 nt (NP ORF).") +
  theme_dviz_open() +
  theme(
    text = element_text(family = dviz_font_family),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, colour = "#000000", 
                                 hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )

print(p_np)




# 5. Combine using Patchwork
combined_figure <- (p_large / p_medium / p_np / p_nss) + 
  plot_layout(guides = "keep") +
  plot_annotation(
    tag_levels = 'A',
    title = "Structural and Antigenic Landscape of RVFV Vaccine Substitutions",
    subtitle = "Mapping of Fixed Substitutions in Smithburn, MP-12, and Clone 13 Strains",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  legend.position = "right")
  ) &
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0, 1)
  )

combined_figure

ggsave(file.path(outDir, "Figure4-v2.pdf"), 
       combined_figure, width = 24, height = 28, 
       dpi = 300, device = cairo_pdf)

ggsave(file.path(outDir, "Figure4-v2.png"), 
       combined_figure,
       width = 24, height = 28, 
       units = "in", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

