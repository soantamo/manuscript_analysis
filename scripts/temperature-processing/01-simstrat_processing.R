
#Data annotation from Simstrat model output for Lake analyses
#Coded by D. Josi 12.5.2022

#### 1. set up working directory and libraries ----

# generalise working directory
# HERE YOU CAN SET YOUR OWN DIRECTORY SO THE SCRIPT RUNS ON ANY WYSS PROJECT COMPUTER
user <- 'cw21p621'
user <- 'djosi'
dd <- paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')
setwd(paste0(dd, 'simstrat1D-currentday/lake_temperature/rawData_perLake/Test'))

library(tidyverse)
library(dplyr)
library(zoo)
library(readr)
library(readxl)
library(data.table)

#### 2. read in files ----

#Function to read all files
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

#Annotate all files with filename
tbl_with_sources <-
  list.files(pattern = "*.dat", 
             full.names = T) %>% 
  map_df(~read_plus(.))

#### 3. Shorten to only dates that exist in the lac data ----

##Specify the start of the lake temperature modelling
#Somehow there are two timezones in the dataset CET and CEST. I did not resolve this

colnames(tbl_with_sources)[1] <- "Date"
as_datetime <- function(x) as.POSIXct("1981-01-01") + as.difftime(x, units = "days")
tbl_with_sources$Date <- as_datetime(tbl_with_sources$Date)

max(tbl_with_sources$Date)
min(tbl_with_sources$Date)
median(tbl_with_sources$Date)


##select timeframe within Project Lac based on file 'PLDB_final_short_depth_11052022.RDS'
# load in the project lac data
lac <- tibble(readRDS(file = paste0(dd, "project-lac/Output files/main dataset/PLDB_final_short_depth_11052022.RDS")))

min_lac_date <- min(lac$Date_setting)
#"2011-08-18"
max_lac_date <- max(lac$Date_setting)
#"2017-09-21"

# filter the data to only be between the date range
tbl_with_sources <- filter(tbl_with_sources, Date >= (min_lac_date-365), Date <= (max_lac_date+35))

#min and max are in different time zones
max(tbl_with_sources$Date)
min(tbl_with_sources$Date)

##remove time values and remain with only dates
tbl_with_sources$Date <- format(as.POSIXct(tbl_with_sources$Date,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d')


#### 4. Change table structure to long version -----

#gc()
#memory.limit(9999999999)
final <- gather(tbl_with_sources,
                key = "Depth",
                value = "Temperature",
                -Date,-filename)
#gc()

final$filename <- substring(final$filename,3)
final$filename <- substring(final$filename,1, nchar(final$filename)-6)

final$Depth <- substring(final$Depth,2) #only works because 0 also has a "-"

final$Depth <- as.numeric(final$Depth)
final$Depth <- as.character(final$Depth)


#### 5. Calculate the average value per Day/Depth/Lake and summarise for 7 and 30 day periods before ----

# data.table is muy rapido!
# library(data.table)

# # summarise mean
# final <- final %>% 
#   group_by(Depth, filename, Date) %>% 
#   summarize(mean_temp_per_day = mean(Temperature, na.rm = T))
final_dt <- data.table(final)
final_dt <- final_dt[,.(mean_temp_per_day = mean(Temperature, na.rm = T), 
                        min_temp_per_day  = min(Temperature, na.rm = T), 
                        max_temp_per_day  = max(Temperature, na.rm = T)),
                     .(Depth, filename, Date)]



##calculate the average/min/max per date over the last 7 and 30 days
#careful dates have to be sorted for the roll functions to work correctly
#sort according to 1) Location 2) Depth 3) Date
final_dt  <- final_dt[order(filename, Depth, Date)]
final_agg <- tibble(final_dt)

# final <- final %>% 
#   arrange(Date) %>%
#   arrange(Depth) %>% 
#   arrange(filename)

# check the structure of NAs
final_nas_only <- final_agg[is.na(final_agg$mean_temp_per_day),]
final_nas_only %>% 
  select(Depth, filename) %>% 
  distinct() %>% 
  group_by(filename) %>% 
  do(min_depth = min(as.numeric(as.character(.$Depth)), na.rm = T)) %>% 
  unnest() %>% 
  View()

# OK so we make the dataframe shorter by removing NAs that are beyond the maximum depth
final_agg <- final_agg[!is.na(final_agg$mean_temp_per_day),]


# calculate rolling properties across the dates
## testing!
# sample_2 <- sample(1:nrow(final_agg), 1)
# final_agg %>% 
#   slice(sample_2:(sample_2+50)) %>% 
#   group_by(Depth, filename) %>% 
#   mutate(mean_last_7days  = rollmeanr(mean_temp_per_day,  k=7,  fill=NA),
#          max_last_7days   = rollmaxr(mean_temp_per_day,   k=7,  fill=NA),
#          min_last_7days   = rollapplyr(mean_temp_per_day, width=7,  FUN = min, fill=NA),
#          mean_last_30days = rollmeanr(mean_temp_per_day,  k=30, fill=NA),
#          max_last_30days  = rollmaxr(mean_temp_per_day,   k=30, fill=NA),
#          min_last_30days  = rollapplyr(mean_temp_per_day, width=30, FUN = min, fill=NA)) %>% 
#   View

#gc()
#memory.limit(9999999999)
final_roll <- final_agg %>% 
  group_by(Depth, filename) %>%
  mutate(mean_last_7days  = rollmeanr(mean_temp_per_day,  k=7,  fill=NA),
         max_last_7days   = rollmaxr(max_temp_per_day,   k=7,  fill=NA),
         min_last_7days   = rollapplyr(min_temp_per_day, width=7,  FUN = min, fill=NA),
         mean_last_30days = rollmeanr(mean_temp_per_day,  k=30, fill=NA),
         max_last_30days  = rollmaxr(max_temp_per_day,   k=30, fill=NA),
         min_last_30days  = rollapplyr(min_temp_per_day, width=30, FUN = min, fill=NA)) %>% 
  ungroup()
#gc()


final_roll %>% 
  slice(sample(1:nrow(.), 10000)) %>% 
  select(mean_temp_per_day:min_last_30days) %>% 
  PerformanceAnalytics::chart.Correlation(.)

#LacdeJoux (misses  n=373), I don't know why
#LakeSarnen and (misses n=373), I don't know why
#LakeHallwil here the modelling only starts at the 10.08.2012




#### 6.  Merge annual average temperature per depth with general information of the lakes ----

annual <- read_csv("annualAvgDepthTempCHlakes.csv")

annual <- gather(annual,
                 key = "filename",
                 value = "mean_temp_per_year",
                 -depth)

colnames(annual)[1] <- "Depth"

#make all values positive
annual$Depth<-abs(annual$Depth)
annual$Depth<-as.character(annual$Depth)

# Add general information of the lakes
info_lakes <- read_excel("info_lakes.xlsx")
colnames(info_lakes)[1] <- "filename"
colnames(info_lakes)[3] <- "Depth_max_lake"

info_lakes<-info_lakes %>%
  mutate(filename = if_else(filename == 'LakeNeuchâtel', 'LakeNeuchatel', filename))


Datainfo<-merge(info_lakes, annual, by=c("filename"), all.y=T)
Datainfo$Depth<-as.character(Datainfo$Depth)




#merge info files with temperature data
#This is the final data file for the temperatures
#takes a long time to process

Alldata_simstrat <- merge(data.table(Datainfo %>% 
                                       mutate(Depth = as.character(.$Depth))), 
                          data.table(final_roll %>% 
                                       mutate(Depth = as.character(.$Depth))),
                          by = c("Depth", "filename"), 
                          all.y=T)

# complete cases for speed
Alldata_simstrat <- Alldata_simstrat[complete.cases(Alldata_simstrat),]

Alldata_simstrat  <- Alldata_simstrat[order(filename, Depth, Date)]

# Alldata_simstrat <- Alldata_simstrat %>% 
#   arrange(Date)%>%
#   arrange(Depth)%>% 
#   arrange(filename)

table(Alldata_simstrat$filename)

# make date compatible with project lac data
names(Alldata_simstrat)[names(Alldata_simstrat) == 'filename'] <- 'LakeBasin'
names(Alldata_simstrat)[names(Alldata_simstrat) == 'Date'] <- 'Date_setting'

# Date to character
Alldata_simstrat$Date_setting <- as.character(Alldata_simstrat$Date_setting)

## save output
saveRDS(Alldata_simstrat, paste0(dd, "simstrat1D-currentday/processed-data/Alldata_simstrat.RDS"))


#### End of script ----