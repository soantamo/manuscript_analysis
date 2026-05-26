############################################################
#############    PROJET LAC - LOAD DATA    #################
############################################################
                  ## Tim Alexander ##

#### 0. Aims of the scripts ----

### Aims ###
# Load fish data, species and lake metadata
# Load libraries
# Rename variables
# Exclude data: French lakes, winter samples, unsuccessful actions, crayfish (later)
# Manage single->batteries: combine pelagic nets into battery, sum net areas for the "new" batteries, average coords for each action 
# combine electrofishing methods
# Calc relative position of the fish in watercolumn
# Aggregate species
# Aggregate habitats
# Correct false/incomplete net mesh labels 
# Remove 70mm mesh and 5mm mesh from CEN pelagic and adjust net areas accordingly

#### 1. Load libraries, set working directory, and read in  fish and settings data ---- 

# load libraries
library(vegan)
library(reshape) # required for cast
library(plyr)
library(stringr)    # required for str_split_fixed
library(ggplot2)

# Set plotting parameters
par(mar=c(5,5,2,2))

# Set working directory
#setwd('D:/Projet Lac/Data files/')
setwd('C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/project-lac/Data files/')

# Load data from *.CSV file
fish <- read.csv("FISH_PL.csv",
                 sep = ",",
                 header = TRUE,
                 na.strings = c("NA"),
                 strip.white=TRUE,
                 stringsAsFactors = TRUE)

# Load PL database extract from *.CSV file
setting <- read.csv("SETTING_ALL.csv",
										sep = ",",
										header = TRUE,
										na.strings = c("NA"),
										strip.white=TRUE,
										stringsAsFactors = TRUE)

setting <- subset(setting, Project == "Projet Lac")

#### 2. Correct inaccurate Bourget data ----

# Fish and net details for Bourget actions Bourget_0045 and Bourget_0043 seem to have been mixed up so have been switched
tmp <- which(setting$Fishec_Action == "Bourget_0045")
setting[which(setting$Fishec_Action == "Bourget_0043"),]$Fishec_Action <- "Bourget_0045"
setting[tmp,]$Fishec_Action <- "Bourget_0043"
setting
setting$Observation <- ifelse(setting$Fishec_Action == "Bourget_0043" | setting$Fishec_Action == "Bourget_0045", 
                              "Fish and net details for Bourget actions Bourget_0045 and Bourget_0043 seem to have been mixed up so have been switched", 
                              as.character(setting$Observation))


#### 3. Merge fish and settings data -----

lac <- merge(fish, setting, by = "Fishec_Action", all.x = T)

names(lac) <- sub("Observation.x", "dbo_fish_Observation", names(lac))
names(lac) <- sub("Observation.y", "dbo_setting_Observation", names(lac))
names(lac) <- sub("Operator.x", "dbo_fish_Operator", names(lac))
names(lac) <- sub("Operator.y", "dbo_setting_Operator", names(lac))

names(lac) <- sub("Catch_altitudem.x", "dbo_fish_Catch_altitudem", names(lac))
names(lac) <- sub("Fishec_Action", "dbo_setting_Fishec_Action", names(lac))
lac$dbo_setting_Catch_altitudem <- lac$dbo_fish_Catch_altitudem
lac$dbo_fish_Fishec_Action <- lac$dbo_setting_Fishec_Action

# save(lac, file = "D:/Projet Lac/Data files/DataRAW_250618.Rdata")

lac$Weightg <- as.numeric(as.character(lac$Weightg))
lac$Eff_Lot <- as.numeric(as.character(lac$Eff_Lot))


#### 4. Load species, lake and fish data ----

# Load species information
species_info <- read.csv("Species_info.csv",
										sep = ",",
										header = TRUE, 
										na.strings = c("NA"), 
										strip.white=TRUE,
										stringsAsFactors = TRUE)

# Load lake information
lake_info <- read.csv("lake_info.csv",
										sep = ",",
										header = TRUE, 
										na.strings = c("NA"), 
										strip.white=TRUE,
										stringsAsFactors = TRUE)

fish_data <- read.csv("Fish_data.csv",
										sep = ",",
										header = TRUE, 
										na.strings = c("NA"), 
										strip.white=TRUE,
										stringsAsFactors = TRUE)

#### 5. Clean empty columns ----

# Deleted empty columns
# summary(lac)
lac$Waterbody_type <- NULL
lac$Altitude_m <- NULL
lac$Electric_pass <- NULL
lac$Substrate_type_max <- NULL
lac$Substrate_type_min <- NULL
lac$Created_dt.x <- NULL
lac$Last_mod_user.x <- NULL
lac$Last_mod_dt.x  <- NULL
lac$behaviour <- NULL
lac$Foto_field_numbers <- NULL
lac$Folder_field_photos <- NULL
lac$Fish_source <- NULL
lac$Fish_origin <- NULL
lac$Generation <- NULL
lac$Electro_pass <- NULL
lac$Secci.m. <- NULL
lac$Temp <- NULL
lac$Created_dt.y <- NULL
lac$Last_mod_dt.y <- NULL
lac$Last_mod_user.y <- NULL
lac$Last_mod_user.y <- NULL
lac$Location4 <- NULL
lac$Mesh_type <- NULL
lac$Mesh_measure <- NULL
lac$Depth_median <- NULL
lac$Depth_mean <- NULL
lac$Distance._shore1.m. <- NULL
lac$Distance._shore2.m. <- NULL
lac$Coordinates_end_N <- NULL
lac$Coordinates_end_E <- NULL
lac$Fishing_time <- NULL
lac$Substrate_type_general <- NULL
lac$Substrate_type_measure <- NULL
lac$Substrate_type_min <- NULL
lac$Substrate_type_max <- NULL
lac$Electric_pass <- NULL
lac$Catch_altitudem.y <- NULL
lac$Net_depth_IN <- NULL
lac$Net_depth_OUT <- NULL
lac$Conductivity_mS <- NULL

# Deleted columns that are not useful 
lac$Taxa_Code <- NULL
lac$Sex <- NULL
lac$Maturity <- NULL
lac$Local_name <- NULL
lac$Altitude_m <- NULL


# rename database variables
colnames(lac)[which(colnames(lac) == "Location1")] <- "Lake"

colnames(lac)[which(colnames(lac) == "Type_Fishing")] <- "Type_fishing"
colnames(lac)[which(colnames(lac) == "Date_Fishing")] <- "Date_fishing"
colnames(lac)[which(colnames(lac) == "Date_Setting")] <- "Date_setting"
colnames(lac)[which(colnames(lac) == "Time_Picking")] <- "Time_picking"
colnames(lac)[which(colnames(lac) == "Time_Setting")] <- "Time_setting"
colnames(lac)[which(colnames(lac) == "Type_Lot")] <- "Type_lot"
colnames(lac)[which(colnames(lac) == "Net_Mesh")] <- "Net_mesh"

colnames(lac)[which(colnames(lac) == "Fishing_Quality")] <- "Fishing_quality"
colnames(lac)[which(colnames(lac) == "Width_Elec.m.")] <- "Width_elec"
colnames(lac)[which(colnames(lac) == "Length_Elec.m.")] <- "Length_elec"
colnames(lac)[which(colnames(lac) == "Depth_Range.m.")] <- "Depth_range"
colnames(lac)[which(colnames(lac) == "Depth_Min.m.")] <- "Depth_min"
colnames(lac)[which(colnames(lac) == "Depth_Max.m.")] <- "Depth_max"
colnames(lac)[which(colnames(lac) == "Taxa_Latin")] <- "Taxa_latin"
colnames(lac)[which(colnames(lac) == "Num_Net")] <- "Num_net"
colnames(lac)[which(colnames(lac) == "Eff_Lot")] <- "Num_indiv"
colnames(lac)[which(colnames(lac) == "Surf_distancem")] <- "Surf_dist"
colnames(lac)[which(colnames(lac) == "Area_net.m2.")] <- "Net_area"
colnames(lac)[which(colnames(lac) == "Length_Min")] <- "Length_min"
colnames(lac)[which(colnames(lac) == "Length_Max")] <- "Length_max"

colnames(lac)[which(colnames(lac) == "dbo_setting_Fishec_Action")] <- "Fishec_action"
colnames(lac)[which(colnames(lac) == "dbo_setting_Operator")] <- "Setting_operator"
colnames(lac)[which(colnames(lac) == "dbo_setting_Observation")] <- "Setting_observation"
colnames(lac)[which(colnames(lac) == "dbo_fish_Observation")] <- "Fish_observation"
colnames(lac)[which(colnames(lac) == "dbo_fish_Operator")] <- "Fish_operator"
colnames(lac)[which(colnames(lac) == "dbo_setting_Catch_altitudem")] <- "Setting_altitudem"
colnames(lac)[which(colnames(lac) == "dbo_fish_Catch_altitudem")] <- "Floor_dist"

colnames(lac)[which(colnames(lac) == "Fishec_Num")] <- "Fishec_num"

lac$Lake <- factor(lac$Lake)


#### 6. Clean and rename dataset names ----

# Identify dataset
lac$Dataset <- with(lac, ifelse( Lake == "Annecy" | Lake ==   "Bourget", "NoVERTorElectro",
                       ifelse(Lake == "Biel" | Lake ==   "Sarnen", "PLv2",
                              ifelse(Lake ==  "Como" | Lake == "Mezzola" | Lake ==  "Varese" | Lake == "Idro" | Lake == "Iseo", "Graia",
                                     "PLv1"))))
lac$Dataset <- with(lac, ifelse( Lake == "Maggiore" &  Setting_operator == "GRAIA", "Graia", as.character(lac$Dataset)))
lac$Dataset <- with(lac, ifelse( Lake == "Maggiore" &  Setting_operator == "CNR", "Graia", as.character(lac$Dataset)))
lac$Dataset <- with(lac, ifelse( Lake == "Maggiore" &  Setting_operator == "CNR-GRAIA", "Graia", as.character(lac$Dataset)))
lac$Dataset <- with(lac, ifelse( Lake == "Garda" &  grepl("Garda_additional", lac$Fishec_action), "Graia", as.character(lac$Dataset)))
lac$Dataset <- with(lac, ifelse( Lake == "Garda" &  Setting_operator == "GRAIA", "Graia", as.character(lac$Dataset)))
lac$Dataset <- factor(lac$Dataset)


#### 7. Exclude lugano winter and leman 2010 excursions and non-standard types of fishing ----


# Exclude additional data
lac <- subset(lac, Lake != "Lugano_winter")
lac <- subset(lac, Lake != "Leman_2010")

lac$Lake <- factor(lac$Lake)

# Excluded other types of fishing
lac <- subset(lac, Type_fishing == "CEN_benthic_net" |
                Type_fishing == "CEN_benthic_CH" |
                Type_fishing == "CEN_pelagic_net" |
                Type_fishing == "CEN_pelagic_net_IT" |
                Type_fishing == "electric_fishing_boat" |
                Type_fishing == "electric_fishing_foot" |
                Type_fishing == "electrofishing_IT" |
                Type_fishing == "Vertical_benthic_net" |
                Type_fishing == "Vertical_pelagic_battery" |
                Type_fishing == "Vertical_pelagic_net" |
                Type_fishing == "Vertical_benthic_battery_CH" |
                Type_fishing == "Vertical_pelagic_battery_CH")
lac$Type_fishing <- factor(lac$Type_fishing)

# Correct drifted vnet fishing quality
lac$Fishing_quality[which(lac$Fishec_action == "Leman_0153")] <- "drift"

# Use only successful fishing actions
lac <- subset(lac, Fishing_quality == "OK")

#### 8. Manage other aspects of data ----


# Correct depth range
lac$Depth_range <- ifelse(lac$Type_fishing == "CEN_benthic_net" & lac$Depth_range == "001-003", "000-003", as.character(lac$Depth_range))
# To date this was listed as 001-003 which makes sense, but for calculations 000-003 is easier to handle

# Exclude just plain wrong actions
lac <- lac[-which(lac$Fishec_num == 76029),]       # Stickleback in bodensee caught at 209m [caught during net deployment or retrieval]

# Excl perch at 306m in Geneva 
# lac$Num_indiv <- ifelse(lac$Fishec_num == 35588, 0, as.numeric(lac$Num_indiv))    -> 35588 is a perch in shallow-set net in Brienz
                  

# Preliminary work to identify vert batteries
lac$Num_net <- with(lac, 
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 5 & Depth_max <= 8, "8m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 8 & Depth_max <= 15, "15m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 15 & Depth_max <= 20, "20m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 20 & Depth_max <= 35, "35m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 35 & Depth_max <= 50, "50m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 50 & Depth_max <= 100, "100m",			 			 			 			 
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 100 & Depth_max <= 200, "200m",
ifelse(Type_fishing == "Vertical_pelagic_net" & Num_net == "" & Depth_max > 200 & Depth_max <= 300, "300m",
			 as.character(lac$Num_net))))))))))
			


#### 9. Create labels to identify single pelagic nets set together as battery ----
lac$Net_action <- do.call(paste, c(lac[c("Lake","Date_setting", "Num_net")], sep = "_"))
lac$Action <- ifelse(lac$Type_fishing == "Vertical_pelagic_net", 
										 as.character(lac$Net_action), 
										 as.character(lac$Fishec_action))
lac$Action <- factor(lac$Action)


#### 10. Correct label for net set later than others in its battery ----
lac$Action[which(lac$Action == "Poschiavo_12/08/2012_20m")] <- "Poschiavo_15/08/2012_50m"
lac$Action[which(lac$Action == "Geneva_26/09/2012_100m")] <- "Geneva_19/09/2012_100m"

#### 11. Sum net area for single -> battery vnets ----

setting_info <- subset(lac, select = c(Fishec_action, Action, Net_area))                 # Extract columns containing setting information 
unique_fishecaction <- setting_info[!duplicated(setting_info[,c("Fishec_action")]),]     # Remove duplicate entries for each Fishec_action
Newnetarea <- tapply(unique_fishecaction$Net_area, unique_fishecaction$Action, sum)      # Sum Fishec_action net areas within 'Actions'
Newnetarea <- data.frame(Newnetarea)
setting_info <- setting_info[!duplicated(setting_info[,c("Action")]),]                   # Remove duplicate fishec_actions for each Action
setting_info <- data.frame(setting_info, row.names = setting_info$Action)                # Make Actions into row labels
setting_info <- merge(setting_info, Newnetarea, by = 0)                                  # Merge setting info and summed net area dataframes by rownames
setting_info <- setting_info[-1]                                                         # Remove duplicate column
lac$Newnetarea <- setting_info$Newnetarea[match(lac$Action, setting_info$Action)] #ckd   # Return new net areas for each Action into the original dataframe
lac$Oldnetarea <- lac$Net_area                                                           # keep copy of old net areas from database



#### 12. Action-level coordinates (centroid) for single -> battery vnets ----

setting_info <- subset(lac, select = c(Action, Fishec_action,  Coordinates_start_N, Coordinates_start_E))
unique_fishecaction <- setting_info[!duplicated(setting_info[,c("Fishec_action")]),]
ActionCoords_N <- with(unique_fishecaction, tapply(Coordinates_start_N, Action, mean))
ActionCoords_E <- with(unique_fishecaction, tapply(Coordinates_start_E, Action, mean))
Newcoords <- data.frame(cbind(ActionCoords_N, ActionCoords_E)) 
setting_info <- setting_info[!duplicated(setting_info[,c("Action")]),] 
setting_info <- data.frame(setting_info, row.names = setting_info$Action)
setting_info <- merge(setting_info, Newcoords, by = 0)
setting_info <- setting_info[-1]
lac$ActionCoords_lat <- setting_info$ActionCoords_N[match(lac$Action, setting_info$Action)]    #ckd
lac$ActionCoords_long <- setting_info$ActionCoords_E[match(lac$Action, setting_info$Action)] 


#### 13. Mean depth for single -> battery vnets ----
setting_info <- subset(lac, select = c(Fishec_action, Action, Depth_max))
unique_fishecaction <- setting_info[!duplicated(setting_info[,c("Fishec_action")]),]
lac$Fishec_depthmax <- lac$Depth_max    # create copy of original depths for use with analyses of vnet single data 
Depthmax_action <- tapply(unique_fishecaction$Depth_max, unique_fishecaction$Action, mean)
Depthmax_action <- data.frame(Depthmax_action)
setting_info <- setting_info[!duplicated(setting_info[,c("Action")]),] 
setting_info <- data.frame(setting_info, row.names = setting_info$Action)
setting_info <- merge(setting_info, Depthmax_action, by = 0)
setting_info <- setting_info[-1]
lac$Depthmax_action <- setting_info$Depthmax_action[match(lac$Action, setting_info$Action)]   #ckd



#### 14. Change Type_fishing for pelagic vert singles to batteries and combine electrofishing methods ----

lac$Type_fishing_old <- lac$Type_fishing
lac$Type_fishing <- ifelse(lac$Type_fishing == "Vertical_pelagic_net", "Vertical_pelagic_battery",
										ifelse(lac$Type_fishing == "electric_fishing_foot", "electrofishing",
										ifelse(lac$Type_fishing == "electric_fishing_boat", "electrofishing",
										as.character(lac$Type_fishing))))
lac$Type_fishing <- factor(lac$Type_fishing)


#### 15. Count number of fishec actions contributing to action ----

setting_info <- subset(lac, select = c(Fishec_action, Action))                            # Extract columns containing setting information 
unique_fishecaction <- setting_info[!duplicated(setting_info[,c("Fishec_action")]),]     # Remove duplicate entries for each Fishec_action
unique_fishecaction$value <- 1
numnets <- with(unique_fishecaction, tapply(value, Action, sum))      # Sum Fishec_action net areas within 'Actions'
numnets <- data.frame(numnets)
setting_info <- setting_info[!duplicated(setting_info[,c("Action")]),]                   # Remove duplicate fishec_actions for each Action
setting_info <- data.frame(setting_info, row.names = setting_info$Action)                # Make Actions into row labels
setting_info <- merge(setting_info, numnets, by = 0)                                  # Merge setting info and summed net area dataframes by rownames
setting_info <- setting_info[-1]                                                         # Remove duplicate column
lac$Numnets <- setting_info$numnets[match(lac$Action, setting_info$Action)] #ckd   # Return new net areas for each Action into the original dataframe
lac$Numnets <- ifelse(lac$Type_fishing_old == "Vertical_pelagic_net", as.numeric(lac$Numnets),
  						 ifelse(lac$Type_fishing_old == "Vertical_pelagic_battery", 7, 
  						 		NA))


#### 16. Change type_fishing to benthic for X and Y VERT nets ----
lac$Type_fishing[grep("^X", lac$Num_net)] <- "Vertical_benthic_net"
lac$Type_fishing[grep("^Y", lac$Num_net)] <- "Vertical_benthic_net"

# Create short names
lac$Typefish <- ifelse(lac$Type_fishing == "CEN_benthic_net", "CENben",
                       ifelse(lac$Type_fishing == "CEN_pelagic_net", "CENpel",
                       ifelse(lac$Type_fishing == "CEN_pelagic_net_IT", "CENpelIT",
											 ifelse(lac$Type_fishing == "Vertical_benthic_net", "VERTshal",
											 ifelse(lac$Type_fishing == "Vertical_pelagic_battery", "VERTdeep",
											 ifelse(lac$Type_fishing == "CEN_benthic_CH", "CENbenCH",
											 ifelse(lac$Type_fishing == "Vertical_benthic_battery_CH", "VERTshalCH",
											 ifelse(lac$Type_fishing == "Vertical_pelagic_battery_CH", "VERTdeepCH",
                       ifelse(lac$Type_fishing == "electrofishing", "electro",
                              ifelse(lac$Type_fishing == "electrofishing_IT", "electroIT",
                                     as.character(lac$Type_fishing)))))))))))
lac$Typefish <- factor(lac$Typefish)

# Where floor_dist < 0, set as 0 (extra net below bar?)... also surfdist
lac$Floor_dist <- ifelse(lac$Floor_dist < 0, 0, lac$Floor_dist)
lac$Surf_dist <- ifelse(lac$Surf_dist < 0, 0, lac$Surf_dist)

# Create depth mid
lac$Depth_mid <- with(lac, ifelse(Typefish == "electro" | 
                                    Typefish == "electroIT" | 
                                    Typefish == "CENben" | 
                                    Typefish == "VERTshalCH" | 
                                    Typefish == "CENpel" | 
                                    Typefish == "CENpelIT" | 
                                    Typefish == "VERTshal" | 
                                    Typefish == "CENbenCH",
												apply(lac[,c("Depth_min", "Depth_max")], 1, mean), 
												   ifelse(Typefish == "VERTdeep" | 
                                    Typefish == "VERTdeepCH", Surf_dist, NA)))

# Calculate floor_dist for missing entries (e.g. Maggiore)
lac$Floor_dist <- ifelse(lac$Typefish == "VERTdeep" & is.na(lac$Floor_dist) == T, (lac$Depth_max - lac$Surf_dist), lac$Floor_dist)

# Calculate surf_dist for missing entries (e.g. Constance)
lac$Surf_dist <- ifelse(lac$Typefish == "VERTdeep" & is.na(lac$Surf_dist) == T, (lac$Depth_max - lac$Floor_dist), lac$Surf_dist)

# Where surf_dist < 0, set as 0 
lac$Surf_dist <- ifelse(lac$Surf_dist < 0, 0, lac$Surf_dist)

### Relative vertical position in watercolumn [only meaningful for vertical pelagic nets]	
lac$Posit_watercol <- round(lac$Floor_dist / (lac$Surf_dist + lac$Floor_dist),3)

### Maximum surveyed depth of lakes
lake_maxdepth <- tapply(lac$Depth_max[-which(is.na(lac$Depth_max)==TRUE)], lac$Lake[-which(is.na(lac$Depth_max)==TRUE)], max)
lac$Maxdepth_lake <- lake_maxdepth[match(lac$Lake, row.names(lake_maxdepth))]    #ckd # add max lake depth to lac 

### "Horizontal" position of action in lake relative to max lake depth 
lac$Relpos_lakedepth <-  round(with(lac, Fishec_depthmax/Maxdepth_lake), 3)


### Protocol labels
lac$Protocol <- with(lac, ifelse(Type_fishing == "CEN_benthic_net" | 
                                   Type_fishing == "CEN_pelagic_net" | 
                                   Type_fishing == "CEN_pelagic_net_IT" | 
                                   Type_fishing ==  "CEN_benthic_CH", 
                                 "CEN",
								          ifelse(Type_fishing == "Vertical_benthic_net" | 
								                   Type_fishing == "Vertical_pelagic_battery" | 
								                   Type_fishing == "Vertical_benthic_battery_CH" | 
								                   Type_fishing == "Vertical_pelagic_battery_CH", 
								                 "VERT", 
								          ifelse(Type_fishing == "electrofishing" |
								                   Type_fishing == "electrofishing_IT", 
								                 "electro", 
								                 "NO_METHOD"))))
lac$Protocol <- factor(lac$Protocol)


### RESET CEN depth ranges
lac$Depth_range_OLD <- lac$Depth_range

lac$Depth_range  <-     with(lac, ifelse(Protocol == "VERT", NA,
                                  ifelse(Typefish == "CENpel" & Depth_mid >= 0 & Depth_mid <= 6, "000-006", 
                                  ifelse(Typefish == "CENpelIT", Typefish, 
                                  ifelse(Typefish == "CENben" & Depth_mid >= 0 & Depth_mid <= 3, "000-003", 
                                  ifelse(Typefish == "CENbenCH" & Depth_mid >= 0 & Depth_mid <= 3, "000-003", 
                                  ifelse(Typefish == "CENben" & Depth_mid > 3 & Depth_mid <= 6, "003-006",        
                                  ifelse(Typefish == "CENbenCH" & Depth_mid > 3 & Depth_mid <= 6, "003-006",        
                                  ifelse(Protocol == "CEN" & Depth_mid > 6 & Depth_mid <= 12, "006-012", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 12 & Depth_mid <= 20, "012-020", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 20 & Depth_mid <= 35, "020-035", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 35 & Depth_mid <= 50, "035-050", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 50 & Depth_mid <= 75, "050-075", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 75 & Depth_mid <= 100, "075-100", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 100 & Depth_mid <= 125, "100-125", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 125 & Depth_mid <= 150, "125-150",
                                  ifelse(Protocol == "CEN" & Depth_mid > 150 & Depth_mid <= 175, "150-175", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 175 & Depth_mid <= 200, "175-200", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 200 & Depth_mid <= 225, "200-225", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 225 & Depth_mid <= 250, "225-250", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 250 & Depth_mid <= 275, "250-275", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 275 & Depth_mid <= 300, "275-300", 
                                  ifelse(Protocol == "CEN" & Depth_mid > 300 & Depth_mid <= 325, "300-325",
                                  ifelse(Protocol == "CEN" & Depth_mid > 325 & Depth_mid <= 350, "325-350",
                                  ifelse(Protocol == "CEN" & Depth_mid > 350 & Depth_mid <= 375, "350-375",
                                  ifelse(Protocol == "CEN" & Depth_mid > 375 & Depth_mid <= 400, "375-400",
                                  as.character("ERROR")))))))))))))))))))))))))))
lac$Depth_range <- factor(lac$Depth_range)

#### 17. Identify fish habitat ----

lac$FishHabitat <- factor(with(lac, ifelse(Typefish == "VERTdeep" & Floor_dist <= 3, "Benthic",
                                    ifelse(Typefish == "VERTdeep" & Floor_dist > 3, "Pelagic",
                                    ifelse(Typefish == "VERTdeepCH" & Floor_dist <= 3, "Benthic",
                                    ifelse(Typefish == "VERTdeepCH" & Floor_dist > 3, "Pelagic",
                                    ifelse(Typefish == "CENpel" | Typefish == "CENpelIT", "Pelagic",
                                    ifelse(Typefish == "CENben" | Typefish == "CENbenCH", "Benthic",
                                    ifelse(Typefish == "electro" | Typefish == "electroIT" | Typefish == "VERTshal" | Typefish == "VERTshalCH", "Littoral", 
                                    NA)))))))))

#### 18. Deal with erroneous net mesh and net areas ----

# Correct false/incomplete net mesh labels 
lac$Net_mesh  <-       ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0" & lac$Net_area == 14, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m",
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0" & lac$Net_area == 21, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1.5m",			 #fixed?
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0" & lac$Net_area == 28, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m",
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0" & lac$Net_area == 42, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x3m",
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0" & lac$Net_area == 70, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x5m",			 #fixed?
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m" & lac$Net_area == 28, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m",			 #fixed?	    
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m" & lac$Net_area == 14, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m",				 #fixed?		    			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_" & lac$Net_area == 42, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x3m",			 #fixed?
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70" & lac$Net_area == 20, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70" & lac$Net_area == 60, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x3m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70" & lac$Net_area == 100, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x5m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.01m_+_2mX70" & lac$Net_area == 20, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.01m_+_6mX70" & lac$Net_area == 20, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.02m_+_2mX70" & lac$Net_area == 40, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70",			 
											 ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.02m_+_6mX70" & lac$Net_area == 40, "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70",			 
											               as.character(lac$Net_mesh))))))))))))))))

# Correct net area
lac$Net_area <- with(lac, ifelse(Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70" & Net_area == 28, 40,
                                        lac$Net_area))
													

# Remove fish caught in 70mm mesh (mesh was not used consistently - only for catching specimens for museum)
#lac <- subset(lac, Meshmm != 70) 	
lac <- lac[-which(lac$Meshmm == 70),] 	

# Correct net areas for 70mm removal
lac$Newnetarea <- ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x1m_+_6mX70", 14,
									ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x2m_+_6mX70", 28,
									ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x3m_+_6mX70", 42,
									ifelse(lac$Net_mesh == "10.0/15.0/20.0/30.0/40.0/50.0/60.0x5m_+_6mX70", 70,			 
									lac$Newnetarea))))							 
																 
# Remove fish caught in 5mm mesh in CEN pelagic (5mm mesh is not meant to be included in the CEN protocol)
lac <- lac[-which(lac$Type_fishing == "CEN_pelagic_net" & lac$Meshmm == "5"),]

# Correct net areas for 5mm removal
lac$Newnetarea <- ifelse(lac$Num_net == "BP2A", 165,
									ifelse(lac$Num_net == "BP2B", 165,
									ifelse(lac$Num_net == "BP3A", 165,
									ifelse(lac$Num_net == "BP3B", 165,
									lac$Newnetarea))))

### Save old Net_area column and assign battery-aggregated net areas (Newnetarea) into original Net_area column
lac$Old_Net_area <- lac$Net_area
lac$Net_area <- lac$Newnetarea

#### 19. Identify mid depth of fishing action ----

lac$Depthmid <- apply(lac[,c("Depth_min", "Depth_max")], 1, mean)

#### 20. Remove fish that were caught while lifting ----

lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught while lifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught during setting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_during_lifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_during_net_lifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_during_net_lifting_caught_alive_net_pochet")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_during_relifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_during_setting?")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "CAUGHT_WHEN_LIFTING?")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "caught_while_lifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "doubts_on_depth?")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "Perch at almost 200m -> caught while lifting")) == T),]   
lac <- lac[which(is.na(str_extract(lac$Fish_observation, "while_lifting")) == T),]   

lac <- lac[which(is.na(str_extract(lac$Fish_observation, "TJA: Burbot in pelagic zone is very unlikely")) == T),]   


#### 21. Manage species data ----

# Remove crayfish
lac <- subset(lac, lac$Taxa_latin != "Orconectes_limosus") 
lac <- subset(lac, lac$Taxa_latin != "Pacifastacus_leniusculus")

lac$Taxa_latin <- factor(lac$Taxa_latin)

# Insert label in Taxa_latin to identify actions without fish
lac$Taxa_latin <- ifelse(lac$Fish_observation == "Action_without_fish", "NO_FISH", as.character(lac$Taxa_latin))
lac$Taxa_latin <- ifelse(lac$Fish_observation == "action_without_fish", "NO_FISH", as.character(lac$Taxa_latin))
lac$Taxa_latin <- factor(lac$Taxa_latin)
lac$Weightg[which(lac$Taxa_latin == "NO_FISH")] <- 0 
lac$Num_indiv[which(lac$Taxa_latin == "NO_FISH")] <- 0 


# Aggregate species
lac$Taxa_latin_aggreg <- lac$Taxa_latin

lac$Taxa_latin_aggreg <- ifelse(grepl("Coregonus", lac$Taxa_latin_aggreg), "Coregonus_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Rutilus", lac$Taxa_latin_aggreg), "Rutilus_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Salvelinus", lac$Taxa_latin_aggreg), "Salvelinus_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg), "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Alosa", lac$Taxa_latin_aggreg), "Alosa_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Salmo", lac$Taxa_latin_aggreg), "Salmo_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Phoxinus", lac$Taxa_latin_aggreg), "Phoxinus_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Gasterosteus", lac$Taxa_latin_aggreg), "Gasterosteus_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Alburnus", lac$Taxa_latin_aggreg), "Alburnus_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Esox", lac$Taxa_latin_aggreg), "Esox_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Gasterosteus", lac$Taxa_latin_aggreg), "Gasterosteus_sp", as.character(lac$Taxa_latin_aggreg))
# lac$Taxa_latin_aggreg <- ifelse(grepl("Squalius", lac$Taxa_latin_aggreg), "Squalius_sp", as.character(lac$Taxa_latin_aggreg))

lac$Taxa_latin_aggreg <- with(lac, ifelse(Taxa_latin == "Salvelinus_namaycush", "Salvelinus_namaycush", as.character(lac$Taxa_latin_aggreg)))

lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Bourget", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Geneva", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Chalain", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Saint-Point", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Remoray", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Rousses", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Brenet", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Neuchatel", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Morat", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Biel", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))
lac$Taxa_latin_aggreg <- ifelse(grepl("Scardinius", lac$Taxa_latin_aggreg) & lac$Lake == "Zurich_Untersee", "Scardinius_sp", as.character(lac$Taxa_latin_aggreg))

lac$Taxa_latin_aggreg <- factor(lac$Taxa_latin_aggreg)


#### 22. Rename habitat variables ----

### Create new variable with pelagic habitats
lac$Pel_habcat <- ifelse(lac$Habitat_1 == "CMIN", "CMIN",
								  ifelse(lac$Habitat_1 == "CMED", "CMED",
									ifelse(lac$Habitat_1 == "CMAX", "CMAX",
									ifelse(lac$Habitat_1 == "TINF", "TINF",
									ifelse(lac$Habitat_1 == "TSUP", "TSUP",
									"")))))
lac$Pel_habcat <- factor(lac$Pel_habcat)


### Manage littoral habitat categories

# Create littoral habitat column 
Action_habitats <- subset(lac, select = c(Action, Habitat_1))
Action_habitats <- Action_habitats[!duplicated(Action_habitats[,c("Action")]),]
Action_habitats$Habitat_1 <- ifelse(Action_habitats$Habitat_1 == "TINF", "",
														 ifelse(Action_habitats$Habitat_1 == "TMED", "",
														 ifelse(Action_habitats$Habitat_1 == "TSUP", "",
														 ifelse(Action_habitats$Habitat_1 == "CMAX", "",
														 ifelse(Action_habitats$Habitat_1 == "CMED", "",
														 ifelse(Action_habitats$Habitat_1 == "CMIN", "",
														 			as.character(Action_habitats$Habitat_1)))))))
lac$Litt_habcat <- Action_habitats$Habitat_1[match(lac$Action, Action_habitats$Action)]
lac$Litt_habcat <- factor(lac$Litt_habcat)

# Aggregate similar categories 
lac$Litt_habcat_aggr <-   ifelse(lac$Litt_habcat == "HLE", "HEL",
                          ifelse(lac$Litt_habcat == "HLD", "HEL",
													ifelse(lac$Litt_habcat == "HYI", "HYD",
												  ifelse(lac$Litt_habcat == "FNO", "SED",
												  ifelse(lac$Litt_habcat == "FNM", "SED",
										 			ifelse(lac$Litt_habcat == "BLS", "BLO",
												  ifelse(lac$Litt_habcat == "DAL", "BED",
												  ifelse(lac$Litt_habcat == "GAL", "COB",
												  ifelse(lac$Litt_habcat == "GLS", "COB",
												 	as.character(lac$Litt_habcat))))))))))

lac$Litt_habcat <- factor(lac$Litt_habcat)
lac$Habitat_1 <- factor(lac$Habitat_1)


#### 23. Re-assign pelagic habitat categories ----
lac$DSUP <- lake_info$DSUP[match(lac$Lake, lake_info$Lake)]    # identify lower limit of pelagic categories
lac$DMIN <- lake_info$DMIN[match(lac$Lake, lake_info$Lake)] 
lac$DMED <- lake_info$DMED[match(lac$Lake, lake_info$Lake)] 
lac$DMAX <- lake_info$DMAX[match(lac$Lake, lake_info$Lake)] 

lac$Lakemaxdepth <- lake_info$Depth_max[match(lac$Lake, lake_info$Lake)] 
lac$Depthmax_action <- ifelse(lac$Depthmax_action > lac$Lakemaxdepth, lac$Lakemaxdepth, lac$Depthmax_action)

lac$Pel_habcat <-   with(lac, ifelse(Typefish == "VERTdeep" | Typefish == "VERTdeepCH",
															        ifelse(Depthmax_action <= DSUP, "DINF",
																			ifelse(Depthmax_action > DSUP & Depthmax_action <= DMIN, "DSUP",
																			ifelse(Depthmax_action > DMIN & Depthmax_action <= DMED, "DMIN",
										  								ifelse(Depthmax_action > DMED & Depthmax_action <= DMAX, "DMED",
										  								ifelse(Depthmax_action > DMAX & Depthmax_action <= Lakemaxdepth, "DMAX",
										  						    "ERROR"))))),""))
lac$Pel_habcat <- factor(lac$Pel_habcat)



#### 24. Identify basins ----


basins <- read.csv("Basins_LucerneZurichConstance.csv",
										sep = ",",
										header = TRUE, 
										na.strings = "NA", 
										strip.white=TRUE,
										stringsAsFactors = TRUE)

lac$Basin <- basins$Basin[match(lac$Fishec_action, basins$Fishec_action)]
lac$Basin <- as.character(lac$Basin)
lac$Basin[which(is.na(lac$Basin) == TRUE)] <- "Not_determined"
lac$Basin[which(lac$Lake == "Zurich" & lac$Basin == NA)] <- "ZurichUntersee"
lac$Basin <- with(lac, ifelse(Lake == "Zurich" & Basin == "Not_determined", "ZurichUntersee",
                       ifelse(Lake == "Constance" & Basin == "Not_determined", "ConstanceObersee",
											 ifelse(Lake == "Lucerne" & Basin == "Not_determined", "Lucerne_transitional",
															as.character(Basin)))))
lac$Basin <- factor(with(lac, ifelse(Basin == "Not_determined", as.character(Lake), as.character(Basin))))
lac$Basin2 <- factor(lac$Basin)
lac$Basin <- with(lac, ifelse(Lake == "Lucerne", "Lucerne",  as.character(Basin)))
lac$Basin <- factor(lac$Basin)



# Add other Constance Basins
cons_basins <- read.csv("Basins_Constance_detail.csv",
										sep = ",",
										header = TRUE, 
										na.strings = "NA", 
										strip.white=TRUE,
										stringsAsFactors = TRUE)

Bas2 <- cons_basins$Basin2[match(lac$Action, cons_basins$Action)]
lac$Basin2 <- ifelse(lac$Lake == "Constance", as.character(Bas2), as.character(lac$Basin2))
lac$Basin2 <- factor(lac$Basin2)


#### Conor Adjustment 
lac$Basin3 <- as.character(lac$Basin2)
lac$Basin3[which(lac$Lake == 'Constance')] <- as.character(lac$Basin[which(lac$Lake == 'Constance')])

#### 25. Phosphorous data ----
# Take best available phosphous data - 5yr mean where available, otherwise Phos_mod
lake_info$Phos5yrMean[is.na(lake_info$Phos5yrMean)] <- "NA"           # use lake_info to keep alldata NAs as NAs
lake_info$Phos_best <- ifelse(lake_info$Phos5yrMean == "NA", lake_info$Phos_mod, lake_info$Phos5yrMean)




#### 26. Cleanup redundant dataframes ----

rm(Newnetarea,
	 cons_basins,
	 lac2,
	 numnets,
	 Bas2,
	unique_fishecaction,
	Action_habitats,
	setting_info,
	ActionCoords_N,
	ActionCoords_E,
	Newcoords,
	Depthmax_action,
	lake_maxdepth,
	basins)



#### 27. Calc length mean/mid for lots ----

lac$Length_mid <- ifelse(lac$Type_lot == "L", apply(lac[,c("Length_min", "Length_max")], 1, mean), NA)
lac$Length_mid <- ifelse(is.na(lac$Lengthmm) == F, lac$Lengthmm, lac$Length_mid)   #  replace where data available


#### 28. Abundance and weight depending on soak time and net area ----

# Backup num_indiv
lac$Num_indiv_raw <- lac$Num_indiv
lac$Weightg_raw <- lac$Weightg

# Calculate fishing time                                                                                 #
library(lubridate)

# Parse and extract month and year from Date_Setting
Date_setting_pars <- parse_date_time(as.character(lac$Date_setting), "dmy", tz = "Europe/Zurich")
Year_setting <- year(Date_setting_pars) 
Month_setting <- month(Date_setting_pars)

# Parse and extract month and year from Date_Fishing
Date_fishing_pars <- parse_date_time(as.character(lac$Date_fishing), "dmy", tz = "Europe/Zurich")
Year_fishing <- year(Date_fishing_pars) 
Month_fishing <- month(Date_fishing_pars)

# Parse and extract hour from Time_Setting
TimeSetting_pars <- parse_date_time(as.character(lac$Time_setting), "HM", tz = "Europe/Zurich")
Hour_setting <- hour(TimeSetting_pars) 

# Parse and extract hour from Time_Picking
TimePicking_pars <- parse_date_time(as.character(lac$Time_picking), "HM", tz = "Europe/Zurich")
Hour_picking <- hour(TimePicking_pars) 


### Calculate fishing time between setting and retrieving
# Specify that this variable is a date/time (join date strings, specify format and timezone)
DateTime_SetNet_pars <- parse_date_time(as.character(paste(lac$Date_setting, lac$Time_setting)), 
                                        "%d%m%Y %H%M", 
                                        tz = "Europe/Zurich")             

DateTime_RetrieveNet_pars <- parse_date_time(as.character(paste(lac$Date_fishing, lac$Time_picking)), 
                                             "%d%m%Y %H%M", 
                                             tz = "Europe/Zurich")

# Calculate difference between setting and retrieving and return in minutes (can also return "hours", "days")
Time_fishing_mins <- difftime(DateTime_RetrieveNet_pars, DateTime_SetNet_pars, 
                              tz = "Europe/Zurich", 
                              units = "mins")  

Time_fishing_hrs <- difftime(DateTime_RetrieveNet_pars, DateTime_SetNet_pars, 
                             tz = "Europe/Zurich", 
                             units = "hours")  

##### Attach new variables and name factors
lac <- cbind(lac, 
             Time_fishing_hrs)

lac$Time_fishing_hrs <- as.numeric(lac$Time_fishing_hrs)

rm(         Date_fishing_pars, 
            Year_fishing, 
            Month_fishing, 
            Date_setting_pars, 
            Year_setting, 
            Month_setting,
            Hour_setting,
            Hour_picking,
            Time_fishing_mins,
            Time_fishing_hrs,
            TimeSetting_pars, 
            TimePicking_pars,
            DateTime_SetNet_pars,
            DateTime_RetrieveNet_pars)


lac$Time_fishing_hrs[which(is.na(lac$Time_fishing_hrs) == T)] <- 14                                      #
lac$Time_fishing_hrs <- with(lac, ifelse(Lake == "Morat" & Typefish == "VERTshal",                       #
																				 14.28,     # Mean soaktime for benthic vnets across other lakes #
																				 lac$Time_fishing_hrs))                                          #
lac$Weightg_soak <-  ifelse(lac$Type_fishing != "electrofishing", lac$Weightg / lac$Time_fishing_hrs * 14, lac$Weightg)                                                 #
lac$Num_indiv_soak <-  ifelse(lac$Type_fishing != "electrofishing", lac$Num_indiv / lac$Time_fishing_hrs * 14, lac$Num_indiv)                                              #

lac$Num_indiv <- lac$Num_indiv_raw
lac$Weightg <- lac$Weightg_raw

#### 29. Save final object ----

actions <- lac[!duplicated(lac[,c("Action")]),]     

### Save data for fast re-loading
save(lac, lake_info, species_info, fish_data, actions, 
     file = "PLDB_final_03052022.Rdata")

# some quicck summaries of the data to view
# datasets are:
# view(dfSummary(lac))
# view(dfSummary(fish_data))
# view(dfSummary(lake_info))
# view(dfSummary(species_info))
# view(dfSummary(actions))

#### 30. descritions of the project lac dataset columns ----

## ACTION ID
# Fishec_action

## LAKE ID
lac %>% dplyr::select(Lake, Location2, Basin, Basin2, Basin3) %>% unique() # here I made basin 3 which should cover our needs
# Lake
# Basin3 # contains the intra-lake basins

## HABITAT CATEGORIES
# FishHabitat      # classified as Littoral, Benthic or Pelagic
# Pel_habcat       # not sure onthe definition of these categories
# Litt_habcat      # these habitat categories are for the shallow benthic nets and electric fishing
# Litt_habcat_aggr # these habitat categories are for the shallow benthic nets and electric fishing aggregated when similar
# Habitat_1        # unaggregated habitat codes

## LOCATION OF FISHING
# Coordinates_start_N # these likely needed to be cleaned but perhaps not depending on the lakes we use
# Coordinates_start_E # these likely needed to be cleaned but perhaps not depending on the lakes we use
# Coord_format        # settings for coordinates (needed for correct conversions)

## SPECIES IDENTITY 
# Fishec_num           # unique identifier for all observations within the fishec group.
# Taxa_latin           # latin name
# Taxa_latin_aggreg    # aggregation when hard to ID: Coregonus, Salvelinus, Alosa, Salmo, Phoxinus, Gasterosteus, Scardinius (in some lakes).
# Final_identification # 80% empty strings but some where errors are corrected? How does this compare to Taxa_latin
# Comments_Identity    # Comments on how the final ID is decided
## Might be worth making a column called 'Taxa_latin_final' which replaces the corrections

## SPECIES PROPERTIES
# Lengthmm       # measured length
# Length_mid     # this length is for 'lots' (i.e., aggregations of individuals)
# Weightg_raw    # unadjusted weight in grams
# Num_indiv_raw  # unadjusted number of individuals
# Weightg_soak   # weight in grams adjusted by soak time and net area
# Num_indiv_soak # number of individuals adjusted by soak time and area

## FISHING METHODS
# Protocol      # CEN VERT or ELETRO
# Type_fishing  # Type_fishing kept 12 types of project lac data
# Typefish      # shortened Type_fishing kept 12 types of project lac data
# Meshmm        # recording of mm depth
# Net_area      # area of the net in m2
# Numnets       # Sum Fishec_action net areas within 'Actions'

## FISHING DEPTH RECORDS
# Surf_dist     # distance from the lake surface
# Floor_dist    # distance from the lake bottom
# Depth_range   # range of depths fished in a given protocol
# Depth_min     # minimum depth of fishing action net
# Depth_max     # maximum depth of fishing action net
# Depthmid      # mean of minimum and maximum depth
# Depth_pelagic      # depth of specimin caught during pelagic net setting
# Setting_altitudem  # depth of net setting

## TIME AND DATE OF FISHING
# Date_setting     # date of laying nets
# Date_fishing     # date of net retrieval
# Time_setting     # time of net setting
# Time_picking     # time of net retrieval
# Time_fishing_hrs # duration of fishing

#### 31. Create a concise version of the project lac dataset with key columns and their descriptions above ----
## select important columns above
lac_concise <- lac %>% dplyr::select(
  ## ACTION ID
  Fishec_action,
  ## LAKE ID
  Lake,Basin3, 
  ## HABITAT CATEGORIES
  FishHabitat, Pel_habcat, Litt_habcat, Litt_habcat_aggr, Habitat_1, 
  ## LOCATION OF FISHING
  Coordinates_start_N, Coordinates_start_E, Coord_format, 
  ## SPECIES IDENTITY
  Fishec_num, Taxa_latin, Taxa_latin_aggreg, Final_identification, Comments_Identity, 
  ## SPECIES PROPERTIES
  Lengthmm, Length_mid, Weightg_raw, Num_indiv_raw, Weightg_soak, Num_indiv_soak, 
  ## FISHING METHODS
  Protocol, Type_fishing, Typefish, Meshmm, Net_area, Numnets, 
  ## FISHING DEPTH RECORDS
  Surf_dist, Floor_dist, Depth_range, Depth_min, Depth_max, Depthmid, Depth_pelagic, Setting_altitudem, 
  ## TIME AND DATE OF FISHING
  Date_setting,Date_fishing, Time_setting, Time_picking, Time_fishing_hrs)

# convert factor columns to character
i <- sapply(lac_concise, is.factor)
lac_concise[i] <- lapply(lac_concise[i], as.character)


#### 32. Filter to the focal lakes for out analysis ----

# list lakes in this data
sort(unique(lac_concise$Lake)) 

# get the lakes of interest with matching temperature data
lac_names <- c("Annecy", "Aulnes", "Biel", "Bonlieu", "Bourget", "Brenet", "Bret", "Brienz", "Chalain",
               "Como", "Constance", "Garda", "Geneva", "Hallwil", "Idro", "Iseo", "Joux", "Lucerne",
               "Lugano", "Maggiore", "Mezzola", "Morat", "Neuchatel", "Poschiavo", "Remoray", "Rousses", "Saint-Point",
               "Sarnen", "Sils", "Thun", "Varese", "Walen", "Zug", "Zurich") 

# get the potential lake names
lake_t <- sort(list.files('C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/simstrat1D-currentday/lake_temperature/rawData_perLake', pattern = '.dat', recursive = T))

lakes_with_t <- c("Biel", "Brienz", "Constance", "Geneva", "Hallwil", "Joux", "Lucerne",
                  "Lugano", "Maggiore", "Morat", "Neuchatel", "Poschiavo", 
                  "Sarnen", "Sils", "Thun", "Walen", "Zug", "Zurich")

# filter the lakes that also have temperature data
lac_filter_t <- lac_concise %>% dplyr::filter(Lake %in% lakes_with_t)


#### 33. remove survey occasions that are far out of the season ----

library(tidyverse)

# remove survey occasions that are far out of the season
lugano_remove <- lac_filter_t %>% 
  mutate(yday = lubridate::yday(as.Date(.$Date_setting))) %>% 
  filter(Lake == 'Lugano', yday < 275) %>% 
  pull(Fishec_action) %>% unique()

neuchatel_remove <- lac_filter_t %>% 
  mutate(yday = lubridate::yday(as.Date(.$Date_setting))) %>% 
  filter(Lake == 'Neuchatel', yday > 314) %>% 
  pull(Fishec_action) %>% unique()


# remove the actions that are outliers
lac_filter_t <- lac_filter_t %>% filter(!Fishec_action %in% c(lugano_remove, neuchatel_remove))

#### 35. add day month year and yearday to data ----

lac_filter_t <- lac_filter_t %>% 
  mutate(day = lubridate::day(as.Date(.$Date_setting, format = '%d/%m/%Y')),
         month = lubridate::month(as.Date(.$Date_setting, format = '%d/%m/%Y')),
         year = lubridate::year(as.Date(.$Date_setting, format = '%d/%m/%Y')),
         yday = lubridate::yday(as.Date(.$Date_setting, format = '%d/%m/%Y')))

#### 34. save final objects ----

# save the smaller version of the dataset with filtered lakes
# save output to a consistent place in the data dump
saveRDS(lac_filter_t, 
        file = paste0('C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/project-lac/Output files/main dataset/', 
                      "PLDB_final_short_03052022.RDS"))
write.csv(lac_filter_t, 
        file = paste0('C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/project-lac/Output files/main dataset/', 
                      "PLDB_final_short_03052022.csv"), 
        row.names = F)

