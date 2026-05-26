#### clean up the taxonomic names

user <- 'cw21p621'
dd <- paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')
lac <- tibble(readRDS(file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_depth_climate_18052022.RDS")))

#### 1. read in the data and the identification key ----

## Match species names of project lac to the manual corrections 
taxa_corrected <- readxl::read_xlsx(path = paste0(dd, 'project-lac/Data files/taxa_id_investigation-EDITED_TA.xlsx'))
lac_taxa <- left_join(lac, 
                      taxa_corrected[c('Taxa_latin', 'Taxa_latin_aggreg', 'Final_identification', 'FINAL')], 
                      by = c('Taxa_latin', 'Taxa_latin_aggreg', 'Final_identification'))

# correcting 979 records with iffy taxonomic identifications
table(is.na(lac_taxa$FINAL))

# get those to fix
lac_taxa_unclean                            <- lac_taxa[lac_taxa$Final_identification != '',]
lac_taxa_unclean_to_modify                  <- lac_taxa_unclean[lac_taxa_unclean$Taxa_latin != lac_taxa_unclean$Final_identification,]
lac_taxa_unclean_to_modify$Taxa_latin_FINAL <- lac_taxa_unclean_to_modify$FINAL

# create columns for the subsets that dont need modifying
lac_taxa_unclean_correct <- lac_taxa_unclean[lac_taxa_unclean$Taxa_latin == lac_taxa_unclean$Final_identification,]
lac_taxa_correct         <- lac_taxa[lac_taxa$Final_identification == '',]
lac_taxa_allcorrect      <- rbind(lac_taxa_correct, lac_taxa_unclean_correct)
lac_taxa_allcorrect$Taxa_latin_FINAL <- lac_taxa_allcorrect$Taxa_latin

# bind together corrected and already correct species names
lac_final <- rbind(lac_taxa_allcorrect, lac_taxa_unclean_to_modify)
nrow(lac_final) == nrow(lac_taxa)

# rename some remaining errors
lac_final$Taxa_latin_FINAL <- ifelse(lac_final$Taxa_latin_FINAL == 'Coregonus_sp_Bodenbalchen', 'Coregonus_litoralis', lac_final$Taxa_latin_FINAL)  
lac_final$Taxa_latin_FINAL <- ifelse(lac_final$Taxa_latin_FINAL == 'Coregonus_sp_Sarnerfelchen', 'Coregonus_sarnensis', lac_final$Taxa_latin_FINAL)        
lac_final$Taxa_latin_FINAL <- ifelse(lac_final$Taxa_latin_FINAL == 'Coregonus_sp_Schwebbalchen', 'Coregonus_intermundia', lac_final$Taxa_latin_FINAL) ## unsure which this should be
lac_final$Taxa_latin_FINAL <- ifelse(lac_final$Taxa_latin_FINAL == 'Coregonus_sp_Zugerbalchen', 'Coregonus_helveticus', lac_final$Taxa_latin_FINAL)        
# here I cross-checked with the lakes in the data for the species and the lakes reported in project lac
#lac_final %>% filter(Taxa_latin_FINAL == 'Coregonus_sp_Zugerbalchen') %>% pull(Lake)

# clean hyphenated names
lac_final$Taxa_latin_FINAL <- gsub('-', '_', lac_final$Taxa_latin_FINAL)

# Missing values are only for Coregonus_suidteri who occure only in Lake Hallwil
lac_final <- subset(lac_final, Taxa_latin_FINAL!="Coregonus_suidteri")

## save file with cleaned names
saveRDS(lac_final,   file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_depth_climate_cleantaxa_20052022.RDS"))
write.csv(lac_final, file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_depth_climate_cleantaxa_20052022.csv"), row.names = F)


