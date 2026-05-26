## Script to clean the geographic coordinates of project lac records
# Aim of script is to 
# process and clean the non-standard data that are not included in Tim Alexanders script

#### 0. Set WD to data and read in libraries ----

dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/project-lac/Data files/'

#### 1.  Read in data files and libraries ----

# Load data from *.CSV file
fish <- read.csv(paste0(dd, "FISH_PL.csv"),
                 sep = ",",
                 header = TRUE,
                 na.strings = "NA",
                 strip.white=TRUE,
                 stringsAsFactors = TRUE)

# Load PL database extract from *.CSV file
setting <- read.csv(paste0(dd, "SETTING_ALL.csv"),
                    sep = ",",
                    header = TRUE,
                    na.strings = "NA",
                    strip.white=TRUE,
                    stringsAsFactors = TRUE)

setting <- subset(setting, Project != "Projet Lac")

