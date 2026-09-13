#!/usr/bin/env Rscript

library(rprojroot)


# utility functions
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


# exports
is.waive <- function(x) inherits(x, "waiver")

`%!in%` <- Negate(`%in%`)

dviz_font_family <- "Arial Narrow"
dviz_font_family_bold <- "Arial Narrow"
dviz_font_family_condensed <- "Arial Narrow"
dviz_font_family_bold_condensed <- "Arial Narrow"


theme_dviz_open <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                            rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14) {
  half_line <- font_size / 2
  
  cowplot::theme_half_open(font_size = font_size, font_family = font_family, line_size = line_size,
                           rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large)  %+replace%
    theme(
      plot.margin = margin(half_line/2, 1.5, half_line/2, 1.5),
      complete = TRUE
    )
}

theme_dviz_hgrid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
                             colour = "grey90") {
  half_line <- font_size / 2
  
  cowplot::theme_minimal_hgrid(font_size = font_size, font_family = font_family, line_size = line_size,
                               rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large,
                               colour = colour)  %+replace%
    theme(
      plot.margin = margin(half_line/2, 1.5, half_line/2, 1.5),
      complete = TRUE
    )
}

theme_dviz_vgrid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
                             colour = "grey90") {
  half_line <- font_size / 2
  
  cowplot::theme_minimal_vgrid(font_size = font_size, font_family = font_family, line_size = line_size,
                               rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large,
                               colour = colour)  %+replace%
    theme(
      plot.margin = margin(half_line/2, 1.5, half_line/2, 1.5),
      complete = TRUE
    )
}


# define base directory
workingDir <- sprintf(
  "%s", 
  find_rstudio_root_file()
)

segmentsDir <- sprintf("%s/segments/%s/%s", 
                       find_rstudio_root_file(), 
                       c("L", "M", "S"),
                       "complete/global")


# output directories
outDir <- file.path(workingDir, "results/figures")
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)