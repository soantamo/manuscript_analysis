#######refining errors in species' identity that were noticed upon examination of the
# data frame, following the table 16 in the Projet Lac "Synthesis report"

pacman::p_load(tidyverse, readxl)

species_occ_lakes <- read_xlsx('data/species_occurrences_lakes.xlsx')

# C. brienzii should not be present in lake Thun, add to Coregonus sp.

species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Coregonus_brienzii") & Lake %in% c("Thun"),
                          "Coregonus_sp", as.character(Species)))

# Salvelinus umbla in Lake Thun should be Salvelinus sp. 

species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Salvelinus_umbla") & Lake %in% c("Thun"),
                          "Salvelinus_sp", as.character(Species)))


# Squalius squalus in Biel and Neuchatel to Squalius cephalus

species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Squalius_squalus") & Lake %in% c("Biel", 
                                                                           "Neuchatel"),
                          "Squalius_cephalus", as.character(Species)))


#######Pooling species to be able to include more samples
# Profundal Cottus gobio are pooled

species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c(
    "Cottus_gobio_Profundal_Thun",
    "Cottus_gobio_Profundal_Lucerne", "Cottus_gobio_Profundal_Walen"
  ),
  "Cottus_sp_Profundal", as.character(Species)
  ))


# Salvelinus is pooled into "profundal" and "limnetic" habitats

# profundal

species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c( "Salvelinus_sp_Profundal_dwarf_Thun", 
                                          "Salvelinus_sp_Profundal_dwarf_VWS", 
                                          "Salvelinus_sp_Profundal_extreme_Thun",
                                          "Salvelinus_sp_Profundal_Walen_I", 
                                          "Salvelinus_sp_Profundal_Walen_II", 
                                          "Salvelinus_profundus"),
                          "Salvelinus_sp_Profundal", as.character(Species)
  ))


# limnetic 
species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Salvelinus_sp_Limnetic_Thun", "Salvelinus_sp_Limnetic_VWS"),
                          "Salvelinus_sp_Limnetic", as.character(Species)))


# Coregonus species are pooled to new groups based on ecomorphs,
# following De-Dayne et al. 2022

# Albeli
species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Coregonus_albellus", "Coregonus_candidus", "Coregonus_confusus",
                                         "Coregonus_heglingus", "Coregonus_zugensis"),
                          "Coregonus_sp_albeli", as.character(Species)))
# Balchen
species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Coregonus_alpinus", "Coregonus_arenicolus", "Coregonus_duplex",
                                         "Coregonus_helveticus", "Coregonus_palaea"),
                          "Coregonus_sp_balchen", as.character(Species)))

# Felchen
species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Coregonus_brienzii", "Coregonus_fatioi", "Coregonus_intermundia",
                                         "Coregonus_litoralis", "Coregonus_macrophthalmus", "Coregonus_zuerichensis"),
                          "Coregonus_sp_felchen", as.character(Species)))

# Large pelagic
species_occ_lakes <- species_occ_lakes |>
  mutate(Species = ifelse(Species %in% c("Coregonus_acrinasus", "Coregonus_wartmanni"),
                          "Coregonus_sp_large_pelagic", as.character(Species)))

# Benthic profundal
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species == "Coregonus_profundus",
                          "Coregonus_sp_benthic_profundal", as.character(Species)))
# Pelagic profundal
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species == "Coregonus_nobilis",
                          "Coregonus_sp_benthic_profundal", as.character(Species))) # here gets combined with benthic form due to low sample size

# "Coregonus_sarnensis": could not be assigned to Albeli or Felchen


# Add Cottus gobio assignments based on ecotypes

# Po profundal
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species == "Cottus_sp_Po_profundal",
                          "Cottus_sp_Profundal", as.character(Species)))

# Littoral - the unknown are only in sarnen and biel, and max depth of 60m. 
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Cottus_gobio_Aare_littoral", "Cottus_gobio_unknownlineage"),
                          "Cottus_gobio_littoral", as.character(Species)))

# according to projet lac report, keep Rhine cottus lineage seperate, these are closer related to river lineage.

# checking table with Barbara
species_occ_lakes %>% filter(grepl('Barbatula', Species)) %>% select(Lake, Species) %>% unique

# Update Barbatula based on Calegari 2025
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Barbatula_sp_Lineage_I", "Barbatula_sp_Lineage_II") & 
                            !Lake %in% c('Constance', 'Constance', 'Geneva'), 
                          "Barbatula_ommata", 
                          as.character(Species)))

species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Barbatula_sp_Lineage_I", "Barbatula_sp_Lineage_II", "Barbatula_quignardi") & 
                            Lake %in% c('Constance', 'Constance', 'Geneva'), 
                          "Barbatula_affinisfluvicola", 
                          as.character(Species)))


# Transfer phoxinus septimania to sp in poschiavo and chalain
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Phoxinus_septimaniae"),
                          "Phoxinus_sp", as.character(Species)))

species_occ_lakes %>% filter(grepl('Phoxinus', Species)) %>% select(Lake, Species) %>% unique

# Transfer phoxinus septimania to sp in poschiavo and chalain
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Salaria_fluviatilis_French", "Salaria_fluviatilis_Italian"),
                          "Salaria_fluviatilis", as.character(Species)))

# Combine Esox into both species to better recover niche
species_occ_lakes <- species_occ_lakes |> 
  mutate(Species = ifelse(Species %in% c("Esox_lucius", "Esox_cisalpinus"),
                          "Esox_spp", as.character(Species)))


# check species list
species <- species_occ_lakes |> 
  distinct(Species) %>% 
  pull(Species) %>% 
  sort

# save as edited version
write.csv(species_occ_lakes, file = 'data/species_occurrences_lakes_May2025.csv', row.names = F)


### ---- Save as a new file the occurrence of the species in the lakes

species_occ_lakes %>% 
  select(Species, Lake) %>% 
  unique() %>% 
  arrange(Species) %>% 
  write_csv(., file = 'data/new_species_per_lake_categories.csv')



