### Script to extract the correct depth information for all records in project lac ----
#### 1. load in libraries ----
library(tidyverse)
library(summarytools)

#### 2. load in data ----

# dataset locations 
dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'

# load in the project lac data
lac <- tibble(readRDS(file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_03052022.RDS")))

#### 3. Preliminary investigations of depth categories ----

# get the relevant types of fishing
sort(table(lac$Type_fishing))

#### 4. Assign depth column estimates: Electric fishing ----

depth_electro <- c('electrofishing')

lac_electro <- lac %>% 
  filter(Type_fishing %in% depth_electro) %>% 
  mutate(Depth_sample = 1) 

plot(lac_electro$Depth_max, lac_electro$Depth_sample)
plot(lac_electro$Depth_min, lac_electro$Depth_sample)


#### 5. Assign depth column estimates: Benthic nets ----

# filter the four types of depth extractions
depth_min_max_type <- c('CEN_benthic_net',
                        'CEN_pelagic_net', 
                        'CEN_benthic_CH', 
                        'CEN_pelagic_net_IT')
# 'Special_net', 
# 'Pic_pelagic_net', 
# 'Special_net_winter', 
# 'Pic_benthic_net', 
# 'Ole_diving', 
# 'M15', 
# 'Special_net_profundus_20mm', 
# 'pelagic_net', 
# 'Eel_trap')

# General comment that could be better if depth_min should be -1.5 meters (the height of the net) because it is the min within the setting area

## extract depths for CEN sampling scheme
lac_CEN <- lac %>% 
  filter(Type_fishing %in% depth_min_max_type) %>% 
  mutate(Depth_sample = (.$Depth_max + (.$Depth_min-1.5)) / 2) %>% 
  mutate(Depth_sample = ifelse(.$Depth_sample < 0, 0, .$Depth_sample))

#### Vertical_benthic_net & Vertical_benthic_battery_CH
# All vertical nets extend through the watercolumn from the lake 
# floor to the lake surface. The columns you describe provide the depth of the 
# lake where the net was set. There is no further information on the height of the fish in these 
# "shallow"-set nets, so fish caugh in one of these nets were therefore caught between 0 to depth_max_m.
lac %>% 
  filter(Type_fishing %in% c('Vertical_benthic_net', 'Vertical_benthic_battery_CH')) %>% 
  do(Depth_range_2 = .$Depth_max - .$Depth_min, 
     Depth_max   = .$Depth_max, 
     Depth_min   = .$Depth_min, 
     Depth_range = .$Depth_range) %>% 
  unnest() -> lac_vertBent_hists
hist(lac_vertBent_hists$Depth_range)
hist(lac_vertBent_hists$Depth_range_2)
hist(lac_vertBent_hists$Depth_max)
hist(lac_vertBent_hists$Depth_min)
plot(lac_vertBent_hists$Depth_min, lac_vertBent_hists$Depth_max)


## treat these seperately as removing those with too large a depth range
lac_vertBent <- lac %>% 
  filter(Type_fishing %in% c('Vertical_benthic_net', 'Vertical_benthic_battery_CH')) %>% 
  mutate(Depth_range_2 = .$Depth_max - (.$Depth_min-1.5)) %>% 
  filter(Depth_range_2 < 5) %>% 
  select(-Depth_range_2) %>% 
  mutate(Depth_sample = (.$Depth_max + (.$Depth_min-1.5)) / 2) %>% 
  mutate(Depth_sample = ifelse(.$Depth_sample < 0, 0, .$Depth_sample))
hist(lac_vertBent$Depth_sample)



#### 6. Assign depth column estimates: Pelagic nets ----

pelagic_distance_type <- c('Vertical_pelagic_battery', 
                           'Vertical_pelagic_battery_CH', 
                           'Vertical_pelagic_net')

lac_pelagic <- lac %>% 
  filter(Type_fishing %in% pelagic_distance_type) %>% 
  mutate(Depth_sample = .$Surf_dist)

hist(lac_pelagic$Depth_sample)


#### 7. Combine back together all datasets ----

# combine back the data with the depth samples
lac_recombined <- bind_rows(lac_electro, lac_CEN, lac_vertBent, lac_pelagic)

#### 8. Clean date column ----

library(lubridate)
class(lac_recombined$Date_fishing)
class(lac_recombined$Date_setting)
unique(lac_recombined$Date_fishing)
unique(lac_recombined$Date_setting)

# ok seems to work correctly
test_date <- data.frame(unique(lac_recombined$Date_fishing), 
           unique(as.Date(lac_recombined$Date_fishing, format = c("%d/%m/%Y"))))


# reset the dates
lac_recombined %>% 
  mutate(Date_fishing = as.Date(.$Date_fishing, format = c("%d/%m/%Y")), 
         Date_setting = as.Date(.$Date_setting, format = c("%d/%m/%Y"))) -> lac_recombined_2

#### 9. Set the dates ----

# save file with depth information
saveRDS(lac_recombined_2, file = paste0(dd, 'project-lac/Output files/main dataset/', "PLDB_final_short_depth_11052022.RDS"))
write.csv(lac_recombined_2, file = paste0(dd, 'project-lac/Output files/main dataset/', "PLDB_final_short_depth_11052022.csv"), 
          row.names = F)


















