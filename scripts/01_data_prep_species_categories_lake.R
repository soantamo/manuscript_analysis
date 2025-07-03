
# information about categorization in each lake for the subgroups
species_endemism <- read_excel("data/species_category_per_lake.xlsx") |> 
  rename(endemism = detail_category)

# C. brienzii should not be present in lake Thun, add to Coregonus sp.

species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Coregonus_brienzii") & fLake %in% c("Thun"),
                          "Coregonus_sp", as.character(species)))

# Salvelinus umbla in Lake Thun should be Salvelinus sp. 

species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Salvelinus_umbla") & fLake %in% c("Thun"),
                          "Salvelinus_sp", as.character(species)))


# Squalius squalus in Biel and Neuchatel to Squalius cephalus

species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Squalius_squalus") & fLake %in% c("Biel", 
                                                                           "Neuchatel"),
                          "Squalius_cephalus", as.character(species)))


#######Pooling species to be able to include more samples
# Profundal Cottus gobio are pooled

species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c(
    "Cottus_gobio_Profundal_Thun",
    "Cottus_gobio_Profundal_Lucerne", "Cottus_gobio_Profundal_Walen"
  ),
  "Cottus_sp_Profundal", as.character(species)
  ))


# Salvelinus is pooled into "profundal" and "limnetic" habitats

# profundal

species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c( "Salvelinus_sp_Profundal_dwarf_Thun", 
                                          "Salvelinus_sp_Profundal_dwarf_VWS", 
                                          "Salvelinus_sp_Profundal_extreme_Thun",
                                          "Salvelinus_sp_Profundal_Walen_I", 
                                          "Salvelinus_sp_Profundal_Walen_II", 
                                          "Salvelinus_profundus"),
                          "Salvelinus_sp_Profundal", as.character(species)
  ))


# limnetic 
species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Salvelinus_sp_Limnetic_Thun", "Salvelinus_sp_Limnetic_VWS"),
                          "Salvelinus_sp_Limnetic", as.character(species)))


# Coregonus species are pooled to new groups based on ecomorphs,
# following De-Dayne et al. 2022

# Albeli
species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Coregonus_albellus", "Coregonus_candidus", "Coregonus_confusus",
                                         "Coregonus_heglingus", "Coregonus_zugensis"),
                          "Coregonus_sp_albeli", as.character(species)))
# Balchen
species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Coregonus_alpinus", "Coregonus_arenicolus", "Coregonus_duplex",
                                         "Coregonus_helveticus", "Coregonus_palaea"),
                          "Coregonus_sp_balchen", as.character(species)))

# Felchen
species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Coregonus_brienzii", "Coregonus_fatioi", "Coregonus_intermundia",
                                         "Coregonus_litoralis", "Coregonus_macrophthalmus", "Coregonus_zuerichensis"),
                          "Coregonus_sp_felchen", as.character(species)))

# Large pelagic
species_endemism <- species_endemism |>
  mutate(species = ifelse(species %in% c("Coregonus_acrinasus", "Coregonus_wartmanni"),
                          "Coregonus_sp_large_pelagic", as.character(species)))

# Benthic profundal
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species == "Coregonus_profundus",
                          "Coregonus_sp_benthic_profundal", as.character(species)))
# Pelagic profundal
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species == "Coregonus_nobilis",
                          "Coregonus_sp_benthic_profundal", as.character(species))) # here gets combined with benthic form due to low sample size

# "Coregonus_sarnensis": could not be assigned to Albeli or Felchen


# Add Cottus gobio assignments based on ecotypes

# Po profundal
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species == "Cottus_sp_Po_profundal",
                          "Cottus_sp_Profundal", as.character(species)))

# Littoral - the unknown are only in sarnen and biel, and max depth of 60m. 
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Cottus_gobio_Aare_littoral", "Cottus_gobio_unknownlineage"),
                          "Cottus_gobio_littoral", as.character(species)))

# according to projet lac report, keep Rhine cottus lineage seperate, these are closer related to river lineage.

# checking table with Barbara
species_endemism %>% filter(grepl('Barbatula', species)) %>% select(fLake, species) %>% unique

# Update Barbatula based on Calegari 2025
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Barbatula_sp_Lineage_I", "Barbatula_sp_Lineage_II") & 
                            !fLake %in% c('Constance', 'Constance', 'Geneva'), 
                          "Barbatula_ommata", 
                          as.character(species)))

species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Barbatula_sp_Lineage_I", "Barbatula_sp_Lineage_II", "Barbatula_quignardi") & 
                            fLake %in% c('Constance', 'Constance', 'Geneva'), 
                          "Barbatula_affinisfluvicola", 
                          as.character(species)))


# Transfer phoxinus septimania to sp in poschiavo and chalain
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Phoxinus_septimaniae"),
                          "Phoxinus_sp", as.character(species)))

species_endemism %>% filter(grepl('Phoxinus', species)) %>% select(fLake, species) %>% unique

# Transfer phoxinus septimania to sp in poschiavo and chalain
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Salaria_fluviatilis_French", "Salaria_fluviatilis_Italian"),
                          "Salaria_fluviatilis", as.character(species)))

# Combine Esox into both species to better recover niche
species_endemism <- species_endemism |> 
  mutate(species = ifelse(species %in% c("Esox_lucius", "Esox_cisalpinus"),
                          "Esox_spp", as.character(species)))


# check species list
species <- species_endemism |> 
  distinct(species) %>% 
  pull(species) %>% 
  sort

# save as edited version
write.csv(species_endemism, file = 'data/species_category_per_lake_May2025.csv', row.names = F)

