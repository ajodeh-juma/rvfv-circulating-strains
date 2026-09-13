#!/usr/bin/env Rscript

###########################################
######## load/install packages ############
###########################################


library(ape)
library(dplyr)
library(readr)
library(phytools)
library(dendextend)
library(magrittr)
library(lubridate)
library(rprojroot)
library(gridExtra)
library(countrycode)
library(stringr)

# exports
source(sprintf("%s/utils.R", "scripts"))

userDir <- path.expand("~")


# read metadata
prepare_metadata <- function(metadata_fn, lineages_fn, ncbi_accession_file) {
  
  # read metadata
  meta <- read.csv(metadata_fn, sep=",", header = TRUE)
  meta$taxa <- paste(meta$accession, meta$host,  meta$country, meta$date, sep="|")
  meta$date <- as.Date(as.character(meta$date), format = "%Y-%m-%d")
  meta$year <- sapply(strsplit(as.character(meta$date), "-", fixed = TRUE), `[`, 1)
  meta$sample_date <- decimal_date(meta$date)
  
  # fix strain labelled 13 as Clone 13
  meta$strain[meta$strain == '13'] <- 'Clone 13'
  
  # read lineages file
  lineages <- read.csv(lineages_fn, sep=",", header = TRUE)
  
  # merge
  metadata <- left_join(meta, lineages, by=c("accession"="Query"))
  
  ncbi <- read.csv(ncbi_accession_file, header = T, sep = ",") |> 
    dplyr::rename("ncbi_accession"="accession")
  
  ncbi <- ncbi |>
    dplyr::select(c("ncbi_accession", "strain"))
  
  ncbi$strain <- gsub("Rwanda[0-9]{4}", "", ncbi$strain)
  
  # merge and select columns

  # accession number, 
  # strain,
  # country, 
  # host, 
  # sampling year, 
  # segment availability, 
  # lineage, 
  # inclusion/exclusion status 
  
  meta <- left_join(metadata, ncbi, by="strain") |> 
    distinct(accession, .keep_all = TRUE) |>
    select(c("ncbi_accession", "strain", "country", "year", "Lineage", "host"))
  
  meta$country <- countrycode(meta$country, origin = "iso3c", 
                              destination = "country.name")
  meta <- meta |> 
    rename("Accession"="ncbi_accession", "Country"="country",
           "Strain"="strain", "Year"="year", "Host"="host") |>
    mutate(Host = str_to_title(Host))
   
  return(meta)
}


metadata_large <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[1], "shared-strains/data/RVFV-L.csv"),
  lineages_fn = file.path(segmentsDir[1], "assignment/output-dir/report/lineages.csv"),
  ncbi_accession_file = file.path(userDir, "projects/RVFV/RVFV-accessions/L-segment/sequences.csv")
)

metadata_medium <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[2], "shared-strains/data/RVFV-M.csv"),
  lineages_fn = file.path(segmentsDir[2], "assignment/output-dir/report/lineages.csv"),
  ncbi_accession_file = file.path(userDir, "projects/RVFV/RVFV-accessions/M-segment/sequences.csv")
)

metadata_np <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[3], "shared-strains/data/RVFV-S.csv"),
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv"),
  ncbi_accession_file = file.path(userDir, "projects/RVFV/RVFV-accessions/S-segment/sequences.csv")
)

metadata_nss <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[3], "shared-strains/data/RVFV-S.csv"),
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv"),
  ncbi_accession_file = file.path(userDir, "projects/RVFV/RVFV-accessions/S-segment/sequences.csv")
)







# Load data and rename Accession columns
df_L <- metadata_large |> rename(L = Accession, Lineage_L = Lineage)
df_M <- metadata_medium |> rename(M = Accession, Lineage_M = Lineage)
df_S <- metadata_np |> rename(S = Accession, Lineage_S = Lineage)

# Join datasets on Strain and common metadata columns
metadata_cols <- c("Strain", "Country", "Year", "Lineage", "Host")



combined_df <- df_L %>%
  full_join(df_M, by = "Strain", suffix = c("_L", "_M")) %>%
  full_join(df_S, by = "Strain") %>%
  rename(
    Country_S = Country, 
    Year_S    = Year, 
    Host_S    = Host
  ) %>%
  mutate(
    Country = coalesce(Country_L, Country_M, Country_S),
    Year    = coalesce(Year_L,    Year_M,    Year_S),
    Host    = coalesce(Host_L,    Host_M,    Host_S)
  ) %>%
  select(
    Strain, 
    L, M, S, 
    Lineage_L, Lineage_M, Lineage_S, 
    Country, Year, Host
  )

write_csv(combined_df, "combined_segments_metadata.csv")







# # library(lubridate)
# # library(dplyr)
# # library(stringr)
# # library(countrycode)
# # library(colorspace)
# # library(tidyverse)
# # library(ggrepel)
# # library(ggplot2)
# # library(patchwork)
# 
# # get path of the script
# if (rstudioapi::isAvailable()){
#   if (require('rstudioapi') != TRUE){
#     install.packages('rstudioapi', repos = r)
#   } else {
#     library(rstudioapi)
#   }
#   wd <- dirname(getActiveDocumentContext()$path)
# } else {
#   wd <- getwd()
# }
# 
# 
# 
# userDir <- path.expand("~")
# # output directories
# figuresDir <- file.path(wd, "figures")
# if (!dir.exists(figuresDir)) dir.create(figuresDir, recursive = TRUE)
# 
# setwd(wd)
# 
# # colours
# cols_countries <- c("#EEEE00", "#5D5E36", "#732039", "#17664D", "#30CED8",
#                     "#4E11A8", "#E88817", "#FF85A9", "#666666", "#91793C",
#                     "#DC143C", "#EE00EE", "#F7D2BC", "#FF69B4", "#AFEEEE",
#                     "#9A3324", "#6495ED", "#12B385")
# 
# cols_lineages <- c("#9370DB", "#30CED8", "#FF0000", "#82A3FF", "#000033",
#                    "#567391", "#FF00B6", "#E88817", "#009FFF", "#FE5000",
#                    "#00FFBE", "#8B5A2B", "#36312E", "#BB85AB", "#B1CC71")
# 
# cols_hosts <- c("#655454", "#B3C4DD", "#D6C499", "#947059", "#BAA844", 
#                 "#ffb302", "#606024")
# 
# 
# country.codes <-  c("BDI", "BFA", "CAF", "EGY", "GIN", "KEN", 
#                     "MDG", "MRT", "MYT", "NAM", "RWA", "SAU", 
#                     "SDN", "SEN", "TZA", "UGA", "ZAF", "ZWE")
# countrynames <- countrycode(country.codes, origin = 'iso3c', 
#                             destination = 'country.name')
# hosts <- c("buffalo", "cow", "goat", "human", "mosquito", "sheep", 
#            "springbok")
# hosts <- str_to_title(hosts)
# df_countries <- data.frame(country = country.codes,
#                            countryname = countrynames,
#                            color = cols_countries)
# df_lineages <- data.frame(Lineage = toupper(letters[1:15]),
#                           color = cols_lineages)
# df_hosts <- data.frame(host=hosts,
#                        color=cols_hosts)
# 
# names(cols_countries) <- as.character(df_countries$countryname)
# names(cols_lineages) <- as.character(df_lineages$Lineage)
# names(cols_hosts) <- as.character(df_hosts$host)
# 
# # function to merge datasets
# prepare_metadata <- function(metadata_file, lineages_file, ncbi_accession_file, prefix, outdir) {
#   
#   # read metadata
#   meta <- read.csv(metadata_file)
#   meta$date <- sapply(strsplit(as.character(meta$taxa), "|", fixed = TRUE), `[`, 4)
#   meta$accession <- sapply(strsplit(as.character(meta$taxa), "|", fixed = TRUE), `[`, 1)
#   meta$date <- as.Date(as.character(meta$date), format = "%Y-%m-%d")
#   meta$year <- sapply(strsplit(as.character(meta$date), "-", fixed = TRUE), `[`, 1)
#   meta$sample_date <- decimal_date(meta$date)
#   
#   # read lineages file and merge with metadata
#   lineages <- read.csv(lineages_file, header = T)
#   m <- left_join(meta,lineages, by=c("accession"="Query"))
#   
#   # # read tempest data
#   # tempest_data <- read.csv(tempest_file, header = T, sep = "\t")
#   # tempest_data$date2 <- sapply(strsplit(as.character(tempest_data$tip), "|", fixed = TRUE), `[`, 4)
#   # tempest_data$date2 <- as.Date(as.character(tempest_data$date2), format = "%Y-%m-%d")
#   # tempest_data$accession <- sapply(strsplit(tempest_data$tip, "|", fixed = TRUE), `[`, 1)
#   # 
#   
#   # merge dataframes
#   meta <- m %>% dplyr::select(c(taxa, accession, strain, Lineage, country, location, host,year,date))
#   meta$countryname <- countrycode(meta$country, origin = 'iso3c', destination = 'country.name')
#   
#   # fix host names
#   meta$host[meta$host == "bovine"] <- "cow"
#   meta$strain[meta$strain == 'Ken 9800523'] <- 'Ken 199800523'
#   
#   # read ncbi updated accessions file
#   ncbi <- read.csv(ncbi_accession_file, header = T, sep = ",") |> 
#     dplyr::rename("ncbi_accession"="accession") 
#   ncbi <- ncbi |>
#     dplyr::select(c("ncbi_accession", "strain"))
#   ncbi$strain <- gsub("Rwanda[0-9]{4}", "", ncbi$strain)
# 
#   meta <- left_join(meta, ncbi, by="strain")
#   
#   meta$Lineage[meta$accession == 'KY196500'] <- 'H' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'DQ380167'] <- 'G' # wrongly assigned as E
#   meta$Lineage[meta$accession == 'DQ380165'] <- 'G' # wrongly assigned as E
#   meta$Lineage[meta$accession == 'DQ380168'] <- 'G' # wrongly assigned as E
#   meta$Lineage[meta$accession == 'DQ380161'] <- 'G' # wrongly assigned as E
#   meta$Lineage[meta$accession == 'DQ380160'] <- 'G' # wrongly assigned as E
#   
#   
#   meta$Lineage[meta$accession == 'DQ375434'] <- 'J' # wrongly assigned as I
#   meta$Lineage[meta$accession == 'MG659874'] <- 'H' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'MG659825'] <- 'H' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'MG659886'] <- 'H' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'KY126678'] <- 'H' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'DQ375425'] <- 'D' # wrongly assigned as C
#   meta$Lineage[meta$accession == 'DQ375433'] <- 'O' # wrongly assigned as L
#   
#   # convert host common names to scientific names
#   meta$sciname <- NA
#   
#   meta <- meta %>%
#     mutate(sciname = case_when(
#       host == 'cow' ~ 'Bos taurus',
#       host == 'human' ~ 'Homo sapiens',
#       host == 'mosquito' ~ 'Aedes',
#       host == 'sheep' ~ 'Ovis aries',
#       host == 'bat' ~ 'Chiroptera',
#       host == 'goat' ~ 'Capra',
#       host == 'buffalo' ~ 'Bubalina',
#       host == 'unknown' ~ 'Aedes',
#       host == 'springbok' ~ 'Antidorcas',
#       TRUE ~ host)
#     )
#   
#   meta$uid <- NA
#   meta <- meta %>%
#     mutate(uid = case_when(
#       sciname == 'Bos taurus' ~ 'ea92020a-298a-4977-b6d3-dccbe816bb5e',
#       sciname == 'Homo sapiens' ~ 'c089caae-43ef-4e4e-bf26-973dd4cb65c5',
#       sciname == 'Aedes' ~ '37fd37ec-5ac1-4429-bee9-ce66a3a2eca4',
#       sciname == 'Ovis aries' ~ 'a9297cbd-10ca-457b-b5d5-a0b038720df7',
#       sciname == 'Chiroptera' ~ '384a1c02-9a04-4688-a1e5-eb15966707e3',
#       sciname == 'Capra' ~ '52c28a94-15e2-482d-ac4d-f75ce7b09b95',
#       sciname == 'Bubalina' ~ '5ede126f-6a18-4a08-898d-e48828e7dcaa',
#       sciname == 'Antidorcas' ~ '320dcfd5-c738-484b-ac28-dab0ac07139b',
#       TRUE ~ sciname)
#     )
#   meta$host <- str_to_title(meta$host)
#   meta <- meta[order(meta$ncbi_accession),]
#   
#   meta.local <- meta |> dplyr::filter(str_detect(ncbi_accession, "^PP|OR"))
#   meta.ncbi <- meta |> dplyr::filter(!str_detect(ncbi_accession, "^PP|OR"))
#   
#   meta.local$source <- 'Study'
#   meta.ncbi <- meta.ncbi[meta.ncbi$accession == meta.ncbi$ncbi_accession,]
#   meta.ncbi$source <- 'NCBI'
#   meta <- rbind(meta.ncbi, meta.local)
#   # print(meta[meta$countryname == 'Rwanda',])
# 
#   # summarise the dataset
# 
#   # how many human and livestock samples were included,
#   dataset_by_host <- meta |>
#     dplyr::group_by(host) |>
#     dplyr::summarise(n=n())
#   df1 <- as.data.frame(dataset_by_host)
#   
#   write.table(df1,
#             file.path(outdir, paste0(prefix, "-global-by-host.tsv")),
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE)
#   
#   # how many were from NCBI
#   dataset_ncbi <- dim(meta.ncbi)[1]
#   print(dataset_ncbi)
# 
#   # your own collection of outbreak samples
#   dataset_local <- dim(meta.local)[1]
#   print(dataset_local)
# 
#   # and the breakdown of samples by country.
#   dataset_by_country <- meta |>
#     dplyr::group_by(countryname) |>
#     dplyr::summarise(n=n()) |>
#     dplyr::rename("Country"="countryname")
#   df2 <- as.data.frame(dataset_by_country)
#   write.table(df2,
#               file.path(outdir, paste0(prefix, "-global-by-country.tsv")),
#               sep = "\t",
#               quote = FALSE,
#               row.names = FALSE)
#   
#   
#   # and the breakdown of samples by year.
#   dataset_by_year <- meta |>
#     dplyr::group_by(year) |>
#     dplyr::summarise(n=n()) |>
#     dplyr::rename("Year"="year")
#   df <- as.data.frame(dataset_by_year)
#   
#   p <- ggplot(df, aes(x=Year, y=n)) +
#     geom_bar(data=df, aes(x=Year, y=n, fill=Year), 
#              stat="identity", alpha=1) +
#     xlab("Year") + 
#     ylab("No. of genomes") +
#     ggtitle("", subtitle = prefix) +
#     scale_fill_discrete_divergingx(name="No. of Sequences",
#                                        palette = "Temps",
#                                        rev = TRUE,
#                                        alpha=0.5) +
#     theme_classic(base_size = 30) +
#     theme(
#       text=element_text(family="Helvetica"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.line = element_line(colour = "#000000", linewidth = 0.6),
#       axis.ticks.y = element_line(colour = "#000000", linewidth = 0.6),
#       axis.ticks.x = element_blank(),
#       plot.subtitle = element_text(colour = "#000000", face = "bold"),
#       plot.tag = element_text(colour = "#000000", face = "bold"),
#       axis.text.x = element_text(colour = "#000000", face = "bold", angle = 90),
#       axis.text.y = element_text(colour = "#000000", face = "bold"),
#       axis.title = element_text(colour = "#000000", face = "bold"),
#       legend.position = "none",
#       legend.key = element_blank())
#   
#   print(p)
# 
#   # and the breakdown of samples by country local
#   dataset_by_country_local <- meta.local |>
#     dplyr::group_by(countryname) |>
#     dplyr::summarise(n=n()) |>
#     dplyr::rename("Country"="countryname")
#   df3 <- as.data.frame(dataset_by_country_local)
#   write.table(df3,
#               file.path(outdir, paste0(prefix, "-outbreak-by-country.tsv")),
#               sep = "\t",
#               quote = FALSE,
#               row.names = FALSE)
#   
#   
#   to_write <- meta |> 
#     dplyr::select(c('ncbi_accession', 'Lineage',
#                     'host', 'countryname', 'location',
#                     'year', 'source')) |>
#     dplyr::rename("Accession"="ncbi_accession",
#                   'Lineage'='Lineage',
#                   "Host"="host",
#                   "Country"="countryname",
#                   "Location"="location",
#                   "Year"="year",
#                   "Source"="source")
#   print(dim(to_write[to_write$Lineage=='C',]))
#   write.table(to_write,
#               file.path(outdir, paste0(prefix, "-all-metadata.tsv")),
#               sep = "\t",
#               quote = FALSE,
#               row.names = FALSE)
#   
#   return(p)
# }
# 
# 
# ############################ small segment #####################################
# 
# metadata_fn <- file.path(wd, "S-global.geolocations.csv")
# lineages_fn <- file.path(wd, "S-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/S-segment/sequences.csv")
# 
# meta.small <- prepare_metadata(metadata_file = metadata_fn,
#                                 lineages_file = lineages_fn,
#                                 ncbi_accession_file = ncbi_fn,
#                                prefix = "S-segment",
#                                outdir = dirname(ncbi_fn))
# 
# write.table(unique(sort(meta.small$countryname)), 
#             file.path(dirname(ncbi_fn), "countries.txt"),
#             row.names = F
#             )
# 
# ############################ lineage C small segment #####################################
# 
# metadata_fn <- file.path(userDir, "projects/RVFV/continuous/segments/S/complete/global/lineage-sequences/nss/C.geolocations.csv")
# lineages_fn <- file.path(wd, "S-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/S-segment/sequences.csv")
# 
# meta.small.c <- prepare_metadata(metadata_file = metadata_fn,
#                                lineages_file = lineages_fn,
#                                ncbi_accession_file = ncbi_fn,
#                                prefix = "NSSLineageC",
#                                outdir = dirname(ncbi_fn))
# 
# # ms <- dplyr::left_join(meta.small, meta.small.c, by="taxa") |>
# #   mutate(Lineage.x = coalesce(Lineage.x, Lineage.y))
# 
# # small <- ms |> 
# #   dplyr::select(c('ncbi_accession.x', 'Lineage.x',
# #                   'host.x', 'countryname.x', 'location.x',
# #                   'date.x')) |>
# #   dplyr::rename("Accession"="ncbi_accession.x",
# #                 'Lineage'='Lineage.x',
# #                 "Host"="host.x",
# #                 "Country"="countryname.x",
# #                 "location"="location.x")
# # 
# # write.table(small,
# #             file.path(dirname(ncbi_fn), paste0("S-segment-all-metadata.tsv")),
# #             sep = "\t",
# #             quote = FALSE,
# #             row.names = FALSE)
# 
# ##################################################################################
# 
# 
# 
# 
# ############################ Medium segment #####################################
# 
# metadata_fn <- file.path(wd, "M-global.geolocations.csv")
# lineages_fn <- file.path(wd, "M-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/M-segment/sequences.csv")
# 
# meta.medium <- prepare_metadata(metadata_file = metadata_fn,
#                                 lineages_file = lineages_fn,
#                                 ncbi_accession_file = ncbi_fn,
#                                 prefix = "M-segment",
#                                 outdir = dirname(ncbi_fn))
# 
# 
# ############################ lineage C medium segment #####################################
# 
# metadata_fn <- file.path(userDir, "projects/RVFV/continuous/segments/M/complete/global/lineage-sequences/C.geolocations.csv")
# lineages_fn <- file.path(wd, "M-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/M-segment/sequences.csv")
# 
# meta.medium.c <- prepare_metadata(metadata_file = metadata_fn,
#                                  lineages_file = lineages_fn,
#                                  ncbi_accession_file = ncbi_fn,
#                                  prefix = "MediumLineageC",
#                                  outdir = dirname(ncbi_fn))
# 
# # mm <- dplyr::left_join(meta.medium, meta.medium.c, by="taxa") |>
# #   mutate(Lineage.x = coalesce(Lineage.x, Lineage.y))
# # 
# # medium <- mm |> 
# #   dplyr::select(c('ncbi_accession.x', 'Lineage.x',
# #                   'host.x', 'countryname.x', 'location.x',
# #                   'date.x')) |>
# #   dplyr::rename("Accession"="ncbi_accession.x",
# #                 'Lineage'='Lineage.x',
# #                 "Host"="host.x",
# #                 "Country"="countryname.x",
# #                 "location"="location.x")
# # medium
# # write.table(medium,
# #             file.path(dirname(ncbi_fn), paste0("M-segment-all-metadata.tsv")),
# #             sep = "\t",
# #             quote = FALSE,
# #             row.names = FALSE)
# 
# ##################################################################################
# 
# 
# ############################ Large segment #####################################
# 
# metadata_fn <- file.path(wd, "L-global.geolocations.csv")
# lineages_fn <- file.path(wd, "L-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/L-segment/sequences.csv")
# 
# meta.large <- prepare_metadata(metadata_file = metadata_fn,
#                                 lineages_file = lineages_fn,
#                                 ncbi_accession_file = ncbi_fn,
#                                 prefix = "L-segment",
#                                 outdir = dirname(ncbi_fn))
# 
# ############################ lineage C small segment #####################################
# 
# metadata_fn <- file.path(userDir, "projects/RVFV/continuous/segments/L/complete/global/lineage-sequences/C.geolocations.csv")
# lineages_fn <- file.path(wd, "L-global.lineages.csv")
# ncbi_fn <- file.path(userDir, "projects/RVFV/RVFV-accessions/L-segment/sequences.csv")
# 
# meta.large.c <- prepare_metadata(metadata_file = metadata_fn,
#                                  lineages_file = lineages_fn,
#                                  ncbi_accession_file = ncbi_fn,
#                                  prefix = "LargeLineageC",
#                                  outdir = dirname(ncbi_fn))
# 
# ##################################################################################
# # ml <- dplyr::left_join(meta.large, meta.large.c, by="taxa") |>
# #   mutate(Lineage.x = coalesce(Lineage.x, Lineage.y))
# # 
# # large <- ml |> 
# #   dplyr::select(c('ncbi_accession.x', 'Lineage.x',
# #                   'host.x', 'countryname.x', 'location.x',
# #                   'date.x')) |>
# #   dplyr::rename("Accession"="ncbi_accession.x",
# #                 'Lineage'='Lineage.x',
# #                 "Host"="host.x",
# #                 "Country"="countryname.x",
# #                 "location"="location.x")
# # large
# # write.table(large,
# #             file.path(dirname(ncbi_fn), paste0("L-segment-all-metadata.tsv")),
# #             sep = "\t",
# #             quote = FALSE,
# #             row.names = FALSE)
# 
# 
# ######### combine seraphim results ##############
# 
# read_seraphim <- function(file, segment){
#   data <- read.csv(file, header = TRUE, sep = ";")
#   data$segment <- segment
#   data <- data |>
#     dplyr::rename(
#       "Segment"= "segment",
#       "EnvironmentalFactor" = "environmental.factor",
#       "RegressionCoefficient" = "regression.coefficient",
#       "Qstatistic" = "Q.statistic",
#       "p(Q) > 0" = "p.Q....0"
#     )
#   data$k <- gsub("k", "", sapply(strsplit(data$EnvironmentalFactor, "_"), `[`, 2))
#   data$EnvironmentalFactor <- paste0(sapply(strsplit(data$EnvironmentalFactor, "_"), `[`, 1),
#                                      " (", sapply(strsplit(data$EnvironmentalFactor, "_"), `[`, 3), ")")
#   
#   data
#   data <- data |> dplyr::select(Segment,
#                                 k, EnvironmentalFactor,
#                                 RegressionCoefficient,
#                                 Qstatistic, "p(Q) > 0", BF)
#   return(data)
# }
# L <- read_seraphim(file = "~/projects/RVFV/data_and_scripts/RVFV-L_1.csv",
#                    segment = "L")
# M <- read_seraphim(file = "~/projects/RVFV/data_and_scripts/RVFV-M_1.csv",
#                    segment = "M")
# S <- read_seraphim(file = "~/projects/RVFV/data_and_scripts/RVFV-N_1.csv",
#                    segment = "S")
# 
# combined <- bind_rows(L, M, S)
# 
# write.table(combined,
#             file.path("~/projects/RVFV/data_and_scripts", "combined_seraphim_1.tsv"),
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE)
# 
# 
# seqs_by_year <- ((((meta.large + theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank())) 
#   / (meta.medium + theme(
#     axis.title.x = element_blank())) / 
#     (meta.small + theme(
#       axis.title.y = element_blank()
#     ))) +
#     plot_layout(guides = "keep",
#                 widths = c(2, 2, 2),
#                 heights = c(1, 1, 1.5),
#                 byrow = T)) | 
#   (((meta.large.c + theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())) 
#    / (meta.medium.c + theme(
#      axis.title.x = element_blank(),
#      axis.title.y = element_blank())) / 
#      (meta.small.c + theme(
#        axis.title.y = element_blank())))) +
#   plot_layout(guides = "keep",
#               widths = c(1, 1, 1),
#               heights = c(1, 1, 1.5),
#               byrow = T)) +
# 
# 
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 30, colour = "#000000", face = "bold"),
#         plot.subtitle = element_text(size = 26, colour = "#000000", face = "bold"),
#         strip.text.x = element_text(size=25, colour = "#000000", face="bold"),
#         legend.title = element_text(size = 24, colour = "#000000", face = "bold"))
# 
# 
# # seqs_by_year <- ((meta.large + theme(
# #                      axis.title.x = element_blank())) 
# #   / (meta.medium + theme(
# #                          axis.title.x = element_blank())) / 
# #     meta.small) +
# #   plot_layout(guides = "keep",
# #               widths = c(2, 2, 2),
# #               heights = c(1, 1, 1.5),
# #               byrow = T) + 
# #   plot_annotation(tag_levels = "A") &
# #   theme(plot.tag = element_text(size = 30, colour = "#000000", face = "bold"),
# #         plot.subtitle = element_text(size = 26, colour = "#000000", face = "bold"),
# #         strip.text.x = element_text(size=25, colour = "#000000", face="bold"),
# #         legend.title = element_text(size = 24, colour = "#000000", face = "bold"))
#   
# 
# ggsave(file.path(figuresDir, "FigureS1.pdf"),
#        seqs_by_year,
#        width = 30, height = 21, units = "in", 
#        limitsize = FALSE,
#        dpi = 300, bg="white", device = cairo_pdf)
# ggsave(file.path(figuresDir, "FigureS1.png"),
#        seqs_by_year,
#        width = 30, height = 21, units = "in", 
#        limitsize = FALSE,
#        dpi = 300, bg="white")
