### script to read in project lac data and estimate thermal niche limits

library(data.table)

# HERE YOU CAN SET YOUR OWN DIRECTORY SO THE SCRIPT RUNS ON ANY WYSS PROJECT COMPUTER
user <- 'cw21p621'
dd <- paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')

#### 1. read in project lac data ----

lac <- readRDS(paste0(dd, 'project-lac/Data files/Lake_climate_fish_data_18052022.RDS'))

# quickly investigate the species names again
sp_summary <- lac[Final_identification != ''][,c('Taxa_latin', 'Taxa_latin_aggreg', 'Final_identification')][order(Taxa_latin)][,.N, by = list(Taxa_latin, Taxa_latin_aggreg, Final_identification)]
View(sp_summary)
sp_summary <- lac[Final_identification != ''][,c('Taxa_latin', 'Taxa_latin_aggreg', 'Final_identification')][order(Taxa_latin)][,.N, by = list(Taxa_latin, Taxa_latin_aggreg, Final_identification)][Taxa_latin != Final_identification]
View(sp_summary)

writexl::write_xlsx(sp_summary, path = paste0(dd, 'project-lac/Data files/taxa_id_investigation.xlsx'))


## decisions
# in consulatation with tim alexander we identified the corrected species names

