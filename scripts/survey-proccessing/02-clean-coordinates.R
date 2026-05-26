### Script for the cleaning of Project Lac coordinates


#### 1. load packages and data ----

my_packages <- c('tidyverse', 'janitor', 'mapview', 'readxl', 'skimr', 'sf', 
                 'rnaturalearth', 'rnaturalearthdata', 'rgeos')
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , 'Package'])]   
if(length(not_installed)) install.packages(not_installed)  
lapply(my_packages, function(x) library(x, character.only = T))

# load in path to data dump 
dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'

#### 2. Read in project lac cleaned data from Tim Alexander ----

# read in project lac data
lac <- tibble(readRDS(file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_03052022.RDS")))

# convert factor columns to character
i <- sapply(lac, is.factor)
lac[i] <- lapply(lac[i], as.character)

#### 3. Look at the properties of the georeferenced locations ----

# Most are in WGS84 but some are not
lac %>% 
  dplyr::select(Coord_format, Coordinates_start_N, Coordinates_start_E) %>% 
  unique() %>% 
  .$Coord_format %>% 
  table

# check the xy
lac_xy <- lac %>% 
  select(Coord_format, Coordinates_start_N, Coordinates_start_E, Lake) %>% 
  unique()


# plot using map view and interact
lac_xy %>% 
  na.omit() %>% 
  st_as_sf(coords= c("Coordinates_start_E", "Coordinates_start_N"), crs=4326) %>% 
  mapview(zcol = 'Coord_format')

# On inspection of the data only the following appear to have potential issues with 
# projections of the data
# The following record types appear to be slightly erroneous in their location
# 'Converted_from_lambert' -> these are in the location Leman, Annency and Bourget

# plot a map of those suspicious records
lac_xy %>% 
  filter(Coord_format == 'Converted_from_lambert') %>% 
  distinct(Lake)

lac_xy %>% 
  na.omit() %>% 
  filter(Coord_format == 'Converted_from_lambert') %>% 
  st_as_sf(coords= c("Coordinates_start_E", "Coordinates_start_N"), crs=4326) %>% 
  mapview(zcol = 'Coord_format')

lac_xy %>% 
  na.omit() %>% 
  filter(Coord_format == 'Lambert') %>% 
  st_as_sf(coords= c("Coordinates_start_E", "Coordinates_start_N"), crs=4326) %>% 
  mapview(zcol = 'Coord_format')


#### 4. Provide corrections to latitudes and longitudes ----

## ACTUALLY THIS IS NO LONGER NEEDED AS WE FOCUS IN ON SWISS LAKES







