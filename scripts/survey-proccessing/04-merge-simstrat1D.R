# dataset locations 
dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'
lac <- readRDS(paste0(dd, 'project-lac/Output files/main dataset/', "PLDB_final_short_depth_11052022.RDS"))
Alldata_simstrat <- readRDS(paste0(dd, "simstrat1D-currentday/processed-data/Alldata_simstrat.RDS"))


#### Merge the simstrat extractions to project lac data ---- 

#Merge temperature data with project lac by lake, date and depth
# lac <- tibble(readRDS(file = paste0(dd, "project-lac/Data files/PLDB_final_short_depth_11052022.RDS")))

#Temperature data of Lake Lugano are also split but not yet in project lac data
#No temp data for lake Sils as model starts only in 2014
lac$LakeBasin<-paste(lac$Lake,lac$Basin3)
table(lac$LakeBasin)

lac <- lac %>%
  mutate(LakeBasin = if_else(LakeBasin == 'Biel Biel', 'LakeBiel', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Brienz Brienz', 'LakeBrienz', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Constance ConstanceObersee', 'UpperLakeConstance', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Constance ConstanceUntersee', 'LowerLakeConstance', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Geneva Geneva', 'LakeGeneva', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Hallwil Hallwil', 'LakeHallwil', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Joux Joux', 'LacdeJoux', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Lucerne LucerneGersau', 'LakeLucerne-Gersauer-andTreibbecken', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Lucerne LucerneUrnersee', 'LakeLucerne-Urnersee', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Lucerne LucerneVitznau', 'LakeLucerne-KreuztrichterandVitznauerbecken', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Lugano Lugano', 'UpperLakeLugano', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Maggiore Maggiore', 'LakeMaggiore', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Morat Morat', 'LakeMurten', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Neuchatel Neuchatel', 'LakeNeuchatel', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Poschiavo Poschiavo', 'LagodiPoschiavo', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Sarnen Sarnen', 'LakeSarnen', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Sils Sils', 'LakeSils', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Thun Thun', 'LakeThun', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Walen Walen', 'Walensee', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Zug Zug', 'LakeZug', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Zurich ZurichObersee', 'UpperLakeZurich', LakeBasin),
         LakeBasin = if_else(LakeBasin == 'Zurich ZurichUntersee', 'LowerLakeZurich', LakeBasin))



#Depth_sample round this column
lac$Depth <- as.character(round(lac$Depth_sample, digits=0))
lac$Date_setting <- as.character(lac$Date_setting)
lac$Date_fishing <- as.character(lac$Date_fishing)


#### 8. Final data file ready for analyses ---- 
library(data.table)
#convert to datatable for joining
lac_dt <- data.table(lac)

# merge
class(lac_dt$Date_setting) == class(Alldata_simstrat$Date_setting)
class(lac_dt$LakeBasin) == class(Alldata_simstrat$LakeBasin)
class(lac_dt$Depth) == class(Alldata_simstrat$Depth)
lac_climate <- merge(lac_dt, Alldata_simstrat, by = c('Date_setting', 'LakeBasin', 'Depth'), all.x=T)

# check the fields are a similar type
lac_dt[LakeBasin == 'LakeHallwil'] %>% pull(Date_setting) %>% unique
Alldata_simstrat[LakeBasin == 'LakeHallwil'] %>% pull(Date_setting) %>% unique
lac_dt[LakeBasin == 'LakeHallwil'] %>% pull(LakeBasin) %>% unique
Alldata_simstrat[LakeBasin == 'LakeHallwil'] %>% pull(LakeBasin) %>% unique
lac_dt[LakeBasin == 'LakeHallwil'] %>% pull(Depth) %>% unique
Alldata_simstrat[LakeBasin == 'LakeHallwil'] %>% pull(Depth) %>% unique
lac_climate[LakeBasin == 'LakeHallwil'] %>% pull(Depth) %>% unique
lac_climate[LakeBasin == 'LakeHallwil'] %>% pull(LakeBasin) %>% unique
lac_climate[LakeBasin == 'LakeHallwil'] %>% pull(Date_setting) %>% unique
lac_climate[LakeBasin == 'LakeHallwil'] %>% pull(max_temp_per_day ) %>% unique


## -- check ranges of temperature for each lake and that all lakes have valid data
lake_temp   <- lac_climate %>% 
  filter(!is.infinite(min_temp_per_day), !is.infinite(max_temp_per_day)) %>% 
  group_by(LakeBasin) %>% 
  do(min_temp = quantile(.$min_temp_per_day, na.rm = T, 0.05), 
     max_temp = quantile(.$max_temp_per_day, na.rm = T, 0.95)) %>% 
  unnest()



# sils and hallwil are still missing - do the date ranges not match up? 
lac_dt %>% 
  filter(LakeBasin %in% c('LakeHallwil', 'LakeSils')) %>% 
  group_by(LakeBasin) %>% 
  do(min_date = min(unique(as.Date(.$Date_setting))), 
     max_date = max(unique(as.Date(.$Date_setting)))) %>% 
  unnest()

Alldata_simstrat %>% 
  filter(LakeBasin %in% c('LakeHallwil', 'LakeSils')) %>% 
  group_by(LakeBasin) %>% 
  do(min_date = min(unique(as.Date(.$Date_setting))), 
     max_date = max(unique(as.Date(.$Date_setting)))) %>% 
  unnest()


## ok so Sil and Hallwil dont have the right dates.

# remove from further analysis
lac_climate <- lac_climate %>% filter(!LakeBasin %in% c('LakeHallwil', 'LakeSils')) 
  
# save file with depth information
saveRDS(lac_climate, 
        file = paste0(dd, 'project-lac/Output files/main dataset/', "PLDB_final_short_depth_climate_18052022.RDS"))

write.csv(lac_climate, 
          file = paste0(dd, 'project-lac/Output files/main dataset/', "PLDB_final_short_depth_climate_18052022.csv"), 
          row.names = F)

