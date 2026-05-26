### script to aggregate for all relevant lakes the data for the different GCMs


#### 1. load in libraries and data paths ----
library(tidyverse)
library(data.table)
library(parallel)

# dataset locations 
# HERE YOU CAN SET YOUR OWN DIRECTORY SO THE SCRIPT RUNS ON ANY WYSS PROJECT COMPUTER
user <- 'cw21p621'
dd <- paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')

#### 2. Read in file for aggregation ----

# vector to correct lake names
match_lakes <- c('bielersee' = 'Biel', 
                 'bodensee' = 'Constance', # here there are two options which makes it a bit more confusing 
                 'brienzersee' = 'Brienz', 
                 'lacleman' = 'Geneva',
                 'lacdejoux' = 'Joux', 
                 'lacneuchatel' = 'Neuchatel', 
                 'murtensee' = 'Murten', 
                 #'thunersee' = ???, 
                 'vierwaldstaettersee' = 'Lucerne', # four basins available
                 'zuerichsee' = 'Zurich'#,           # here upper and lower available
                 #'zugersee' = ???
)


# list of the files from the ch2018 simstrat water temperature
all_ch2018 <- list.files(paste0(dd, 'simstrat1D-ch2018/watertemperature/WaterTemperature'), full.names = T)
all_ch2018_short <- list.files(paste0(dd, 'simstrat1D-ch2018/watertemperature/WaterTemperature'), full.names = F)

# filter to those that are relevant for our analysis
all_ch2018_list <- lapply(lapply(match_lakes, function(x) grepl(x, all_ch2018)), function(x) all_ch2018[x])
all_ch2018_short_list <- lapply(lapply(match_lakes, function(x) grepl(x, all_ch2018_short)), function(x) all_ch2018_short[x])

# make dataframe out of example set of data
df_file_names <- lapply(all_ch2018_short_list, function(x){ 
  z <-  data.frame(str_split(x, '_', simplify = T))
  #z[[8]] <- NULL
  colnames(z) <- c('simstrat', 'lake', 'model1', 'model2', 'model3', 'resolution', 'rcp', 'ignore')
  return(z %>% dplyr::select(-ignore) %>% mutate(filename = x, full_filename = paste0(dd, 'simstrat1D-ch2018/watertemperature/WaterTemperature/',x)))})


# loop over each lake and merge all data files and aggregate by RCP
lapply(df_file_names[1:length(df_file_names)], function(x){
  
  # split by the rcp
  lapply(split(x, paste0(x$rcp, x$lake)), function(x1){
    
    # apply read and matrix mean over rcp
    listed_fread <- lapply(x1$full_filename, function(files) as.matrix(fread(files, skip = 1, header = T)))
    matrix_mean <- Reduce("+", listed_fread) / length(listed_fread)
    
    # create directory structure
    output_dir <- paste0(dd, 'simstrat1D-ch2018/aggregated-watertemperature/')
    dir.create(output_dir, recursive = T)
    file_name <- unique(paste0(x1$lake, '_', x1$rcp, '_mean-temperature.csv'))
    fwrite(matrix_mean, file = paste0(output_dir, file_name))
    
  })
  
})



## potential to explore the variation between GCMs

x = df_file_names[1:length(df_file_names)][[1]]
x1 = split(x, paste0(x$rcp, x$lake))[[1]]
# loop over each lake and merge all data files and aggregate by RCP
lapply(df_file_names[1:length(df_file_names)], function(x){
  
  # split by the rcp
  lapply(split(x, paste0(x$rcp, x$lake)), function(x1){
    
    # apply read and matrix mean over rcp
    listed_fread <- lapply(x1$full_filename, function(files) as.matrix(fread(files, skip = 1, header = T)))
    all_files <- rbindlist(lapply(listed_fread, function(x){
                                                colnames(x)[1] <- 'date' 
                                                gather(data.table(x), key = 'depth', value = 'temperature', -date)}))
    
    all_files <- na.omit(all_files)
    all_files <- all_files[ , list(mean_temp = mean(temperature), max_temp = max(temperature), min_temp = min(temperature)), 
                      by = .(date, depth)]
    all_files$range_temp <- all_files$max_temp-all_files$min_temp
    
    # apply transformation for space saving
    all_files[,c('mean_temp', 'max_temp', 'min_temp', 'range_temp')] <- round(all_files[,c('mean_temp', 'max_temp', 'min_temp', 'range_temp')]*100)
    
    # create directory structure
    output_dir <- paste0(dd, 'simstrat1D-ch2018/aggregated-watertemperature-all-properties/')
    dir.create(output_dir, recursive = T)
    file_name <- unique(paste0(x1$lake, '_', x1$rcp, '_mean-temperature.csv'))
    fwrite(all_files, file = paste0(output_dir, file_name))
    
  })
  
})


## testing process
# test_lake <- df_file_names[[1]]
#
# test_lake_rcp26 <- split(test_lake, paste0(test_lake$rcp, test_lake$lake))[[2]]
#
# read in all files per RCPs 
# listed_fread <- lapply(test_lake_rcp26$full_filename, function(x) as.matrix(fread(x, skip = 1, header = T)))
# matrix_mean <- Reduce("+", listed_fread) / length(listed_fread)
# 
# fwrite(matrix_mean, )


