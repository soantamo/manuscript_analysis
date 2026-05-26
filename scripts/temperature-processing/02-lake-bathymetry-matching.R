#### Script aims to link together the simstrat climate data with the 3D depth properties of the lakes

#### 1. Read in libraries ----

library(tidyverse)
library(sf)
library(terra)
library(mapview)
library(data.table)
library(parallel)

# dataset locations 
# HERE YOU CAN SET YOUR OWN DIRECTORY SO THE SCRIPT RUNS ON ANY WYSS PROJECT COMPUTER
user <- 'cw21p621'
dd <- paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')

#### 2. Handle depth data provided by 3d bathymetry ----

bathy_files <- list.files(paste0(dd, 'swiss-BATHY3D/raw-xyz'), full.names = T, pattern = '.xyz')
bathy_files_short <- list.files(paste0(dd, 'swiss-BATHY3D/raw-xyz'), full.names = F, pattern = '.xyz')

# assign altitude to turn bathymetry at altitude into depth from 0m surface
# values taken from simstrat.eawag.ch/lakes
altitude_corrections <- c('bielersee' = 429, 
                          'bodensee' = 395, 
                          'brienzersee' = 564, 
                          'lacdejoux' = 1004, 
                          'lacleman'  = 372,
                          'lacneuchatel' = 429, 
                          'murtensee' = 429, 
                          'thunersee' = 558, 
                          'vierwaldstaettersee' = 434, 
                          'zuerichsee' = 406, 
                          'zugersee' = 417)

# set temporary file location to my hard-drive
terraOptions(temp = 'D:/temp')
for(i in 5){
  
  # read in first file
  file_i <- bathy_files[i]
  print(file_i)
  
  # list the xyz files
  file_i_xyz <- list.files(file_i, full.names = T)
  
  # loop through 10 at a time to save memory
  file_i_xyz_splits <- split(file_i_xyz, ceiling(seq_along(file_i_xyz)/10))
  
  # LOOK THROUGH ALL THE SPLITS OF THE FILES
  for(z in 1:length(file_i_xyz_splits)){
    
    print((z/length(file_i_xyz_splits)*100))
  
  # loop through and read table
  file_i_xyz_read <- lapply(file_i_xyz_splits[[z]], function(x) read.table(x, header = T) )
  
  # bind table 
  file_i_xyz_bind <- bind_rows(file_i_xyz_read)
  
  # aggregate to 10mx10m
  file_i_xyz_bind$x <- plyr::round_any(file_i_xyz_bind$x, 10)
  file_i_xyz_bind$y <- plyr::round_any(file_i_xyz_bind$y, 10)
  
  # remove duplicated
  file_i_xyz_bind <- data.frame(unique(data.table(file_i_xyz_bind), by = c("x", "y")))

  # round to nearest meter
  file_i_xyz_bind[,3] <- round(file_i_xyz_bind[,3])
  
  # convert to rast
  file_i_xyz_rast <- rast(file_i_xyz_bind)
  
  # get the altitude correction
  alt <- altitude_corrections[which(sapply(names(altitude_corrections), function(x) grepl(x, file_i)))]
  file_i_xyz_rast <- file_i_xyz_rast - alt
  file_i_xyz_rast[file_i_xyz_rast>0] = NA
  
  # convert back to xyz
  xyz <- data.frame(crds(file_i_xyz_rast), na.omit(values(file_i_xyz_rast)))
  
  # clean up 
  rm(file_i_xyz_rast, file_i_xyz_bind, file_i_xyz_read, file_i_xyz)
  gc()
  
  # set the maximum depth per 10m x 10m grid cell
  names(xyz)[3] <- 'max_depth'
  
  # remove NAs
  xyz <- xyz[!is.infinite(xyz[,3]),]
  xyz <- xyz[!is.na(xyz[,3]),]
  
  # create new column for each potential depth
  for(col in (1+ncol(xyz)):(abs(min(xyz$max_depth))+1+ncol(xyz))){
    xyz[,col] <- NA
    names(xyz)[col] <- paste0('depth_', col-4)
    }
  
  # pivot longer
  xyz <- gather(xyz, key = 'depth', value = 'temperature', names(xyz)[-c(1:3)])
  
  # convert depth to numeric
  xyz$depth <- -as.numeric(gsub('depth_', '', xyz$depth))
  
  # remove depths that are not > or = to the max depth
  xyz <- xyz[xyz$depth >= xyz$max_depth,]
  
  # remove superfluous columns
  xyz$temperature <- NULL
  xyz$max_depth <- NULL
  
  # now save file of xyz that can be link to depth - temperature files
  ##### TAKE CARE HERE BECAUSE IF NOT RUNNING FOR THE FIRST TIME, IF THE FILE EXISTS, IT WILL APPEND THE FILES
  if(!file.exists(paste0(dd, 'swiss-BATHY3D/clean-xyz/', gsub('.xyz', '.csv', bathy_files_short[i])))){
  dir.create(paste0(dd, 'swiss-BATHY3D/clean-xyz/'))
  fwrite(xyz, file = paste0(dd, 'swiss-BATHY3D/clean-xyz/', gsub('.xyz', '.csv', bathy_files_short[i])))
  }else{
      fwrite(xyz, 
           file = paste0(dd, 'swiss-BATHY3D/clean-xyz/', gsub('.xyz', '.csv', bathy_files_short[i])), 
           append = T)
  }
  rm(xyz)
  gc()
  source('scripts/clear_tmp_files.R')
  }
  
}


### 3. estimate volume at each depth ----

# read in cleaned data above
clean_bathy       <- list.files(paste0(dd, 'swiss-BATHY3D/clean-xyz/'), full.names = T)
clean_bathy_short <- list.files(paste0(dd, 'swiss-BATHY3D/clean-xyz/'), full.names = F)

for(i in 1:length(clean_bathy)){
bathy_read <- fread(clean_bathy[i])
bathy_volume_m3 <- bathy_read[,list(volume_m3 = .N*100),by = list(depth)]
dir.create(paste0(dd, 'swiss-BATHY3D/clean-volume/'))
fwrite(bathy_volume_m3, paste0(dd, 'swiss-BATHY3D/clean-volume/', clean_bathy_short[i]))
}




