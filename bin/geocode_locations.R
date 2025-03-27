#!/usr/bin/env Rscript

library(argparse)
library(tidygeocoder)
library(dplyr)
library(countrycode)
library(R.utils)

library(ggplot2)
library(maps)
# library(sf)



usage <- function() {
  usage.text <- '\nUsage: geocode_locations.R --infile <path to the geolocations file having the columns "taxa", "country" and "location" --prefix <prefix to the output filename>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--infile", default=NULL, help="path to the metadata csv file containing the columns: taxa, country and location")
parser$add_argument("--prefix", default=NULL, help="prefix to the output filename")
args <- parser$parse_args()


if (is.null(args$infile)) {
  parser$print_help()
  stop("Please provide the metadata csv file having column labels taxa, country and location", call.=FALSE)
}
if (!file.exists(args$infile)) {
  parser$print_help()
  stop(paste("The following metadata csv file don't exist:", 
             paste(args$infile, sep='', collapse=' '), sep=' '), call.=FALSE)
} else {
  in_file <- args$infile
}

prefix <- args$prefix

# fix country location data
geolocations <- read.csv(in_file)
geolocations$location[geolocations$location == 'Grande Island'] <- 'Grande Comore'
geolocations$location[geolocations$location == 'Gatooma'] <- 'Kadoma'
geolocations$location[geolocations$location == 'Beatrice'] <- 'Sinoia'
geolocations$location[geolocations$location == 'Fada-Ngourma'] <- 'Gourma'
geolocations$location[geolocations$location == 'Sinoia'] <- 'Chinhoyi'
geolocations$location[geolocations$location == 'Bela-Bela'] <- 'Warmbaths'
geolocations$location[geolocations$location == 'Tanga'] <- ''
geolocations$location[geolocations$location == 'Mombassa'] <- 'Mombasa'
geolocations$location[geolocations$location == 'Sada'] <- 'Sada, Mayotte'
geolocations$location[geolocations$location == 'Magadi'] <- 'Magadi, Kenya'
geolocations$location[geolocations$location == 'Embu'] <- 'Embu, Kenya'
geolocations$location[geolocations$location == 'Meru Central'] <- 'Meru, Kenya'
geolocations$location[geolocations$location == 'Meru South'] <- 'Meru, Kenya'
geolocations$location[geolocations$location == 'Gezira'] <- 'Wad Madani'
geolocations$location[geolocations$location == 'Natal'] <- 'KwaZulu-Natal'
geolocations$location[geolocations$location %in% c('St Joseph', 'Salisbury')] <- ''
geolocations$location[geolocations$location == 'St Joseph'] <- ''
geolocations$location[geolocations$location == 'Thika'] <- 'Thika, Kenya'
geolocations$location[geolocations$location == 'Bura'] <- 'Bura, Kenya'

geolocations <- geolocations %>% 
  dplyr::mutate(location = ifelse(location %in% "", countrycode(country, origin="iso3c", destination='country.name'), location))


# geocode
geocode_tbl <- tidygeocoder::geocode(geolocations,
                      address = "location",
                      method = "osm") %>% 
  dplyr::select(taxa, lat, long) %>% 
  dplyr::rename("traits" = taxa)

# write to output file
outdir <- dirname(R.utils::getAbsolutePath(in_file))


# plot
tbl <- dplyr::left_join(geocode_tbl, geolocations,
                        by=c("traits"="taxa")) %>%
  dplyr::select(traits, lat, long, country, location)

# xlim and ylim need to be changed. These were specific for Africa and the Middle East
p <- ggplot(tbl, aes(long, lat), color = "grey99") +
  borders("world", xlim = c(-17.520278, 55.47519), ylim = c(-34.65302, 37.34698)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = location)) +
  theme_void()

write.table(geocode_tbl, 
          file = file.path(outdir, paste0(prefix, "_geocoded.txt")),
          row.names = F,
          sep = "\t",
          quote = F)

ggsave(file.path(outdir, paste0(prefix, "_geocoded.pdf")),
       p, width = 25, height = 20, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "pdf")
