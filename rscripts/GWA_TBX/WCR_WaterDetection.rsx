#  Copyright (c) 2017, GeoVille Information Systems GmbH
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.
#
#
# Date created: 06.05.2017
# Date last modified: 09.05.2017
#
#
# __author__ = "Christina Ludwig"
# __version__ = "1.0"


# INPUT PARAMTERS --------------------------------------------------------

# -> FOR QGIS Processing modules
##Water Detection=name
##Water Cycle Regime=group
##Directory_containing_indices=folder
##Directory_containing_TWI=folder
##Output_Directory=folder
##Directory_containing_SAR=folder
##Start_Date= optional string
##End_Date= optional string
##Minimum_water_probability=number 45
##Tile_size_in_meter= optional number 1800
#Global_threshold_based_on = selection Mean;Median
##Minimum_mapping_unit = number 3
##Minimum_AOI_coverage = number 40
##Plot_water_probability= Boolean False
##Plot_certainty_indicator= Boolean False

sink(file.path(Output_Directory, "log.txt"))

# -> Test parameters for execution in R
debug <- F
starttime <- proc.time()

if (debug) {
  .libPaths("C:\\Users\\ludwig\\.qgis2\\processing\\rlibs")
  Directory_containing_indices= "T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\WCR\\indices"
  Directory_containing_TWI= "T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\TWI"
  #Directory_containing_SAR= "T:/Processing/2687_GW_A/02_Interim_Products/SAR/site_98/EQUI7_AF010M/E014N068T1/Seasonal"#
  Directory_containing_SAR= "" # "T:/Processing/2687_GW_A/02_Interim_Products/SAR/site_98/EQUI7_AF010M"
  Output_Directory="T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\test"
  Plot_water_probability <- T
  Minimum_water_probability = 45
  Minimum_AOI_coverage <- 40
  Start_Date <- " "
  End_Date <- " "
  Tile_size_in_meter = 1800
  Minimum_mapping_unit = 3
  Plot_certainty_indicator = T
}
Increase_sensitivity=FALSE

#sink(file.path(Output_Directory, "log.txt"))

library(raster)
library(rgdal)
library(rpanel)
library(stringr)
library(GWAutils)

Start_Date <- str_trim(Start_Date)
End_Date <- str_trim(End_Date)

Minimum_AOI_coverage <- Minimum_AOI_coverage / 100.
Minimum_water_probability <- Minimum_water_probability * 10.

plot_opt = FALSE
debug_mode = FALSE

mmu_wet <- Minimum_mapping_unit
mmu_water <- Minimum_mapping_unit
Nodata_value = -32768
Nodata_value_mask = 255
Lower_than = FALSE
NDVIthresh = 8500

# Index combination for water and wetness detection
indices_water <- c( "mNDWI", "NDWI")

# Indices that need to be reversed so that water / wetness has high values
reverse <- c("NDVI")

#default
pythonPath = ""


# Set raster settings ----------------------------------------------------


# CHECK INPUT PARAMETERS ------------------------------------------------

# Index directory
if (!file.exists(Directory_containing_indices)) {
  cat(Directory_containing_indices)
  rp.messagebox("Directory containing indices does not exist!")
  stop()
} else {
  cat("\nDirectory containing indices: ", Directory_containing_indices)
}

# TWI directory
if (!file.exists(Directory_containing_TWI)) {
  cat(Directory_containing_TWI)
  rp.messagebox("Directory containing TWI does not exist!")
  stop()
} else {
  cat("\nDirectory containing TWI: ", Directory_containing_TWI)
}

# SAR directory
if (!file.exists(Directory_containing_SAR)) {
  cat("\nDirectory containing SAR does not exist! Processing will be done without SAR data.")
  SAR <- FALSE
} else {
  cat("\nDirectory containing SAR: ", Directory_containing_SAR)
  SAR <- TRUE
}


# Get dates of input scenes
# ---------------------------

search_pattern <- paste0('*', '_NDVI.vrt')
NDVIfiles <- list.files(Directory_containing_indices, pattern = glob2rx(search_pattern), full.names = TRUE,
                        include.dirs = TRUE, recursive = TRUE)
if (length(NDVIfiles)==0) {
  cat("ERROR:No NDVI files found.")
  stop()
}

Title <- paste(substr(basename(NDVIfiles[1]), 1,2), "_wat",Minimum_water_probability/10., "_win", Tile_size_in_meter, "_mmu", Minimum_mapping_unit, sep="")

# Output directory
if (!file.exists(Output_Directory)) {
  rp.messagebox("Output directory does not exist!")
  stop("Output directory does not exist!")
} else {
  cat("\nOutput directory: ", Output_Directory,"\n")
  Output_Directory_site <- file.path(Output_Directory, Title)
  if (!file.exists(Output_Directory_site)) {
    dir.create(Output_Directory_site, recursive=TRUE)
  }
}

# Check resolution
resolutions <- c()
for (f in NDVIfiles) {
  ndvi <- raster(f)
  resolutions <- c(resolutions, res(ndvi))
}
if (length(unique(resolutions)) > 1) {
  rp.messagebox("Invalid input data: Input scenes have different resolutions. Lansat and Sentinel indices must be stored in separate folders.")
  stop("\nInput files have different resolutions.")
}

# Get tile size from resolution -----------------------------------------------

Tile_size_water = as.numeric(Tile_size_in_meter) /  resolutions[1]
Tile_size_water2 = (as.numeric(Tile_size_in_meter)*0.66) /  resolutions[1]

if (Tile_size_in_meter < 1000) {
  rp.messagebox("Invalid input parameter: Tile size must be greater than 1000m.")
  stop()
} else if ((Tile_size_water > ndvi@ncols) || (Tile_size_water > ndvi@nrows)) {
  rp.messagebox("Invalid input parameter: Tile size must not be greater than image size. Reduce the value for 'Tile Size in meters'.")
  stop()
}

# DATES --------------------------
dates <- list()
for (i in 1:length(NDVIfiles)) {
  dates[[i]] <- as.Date(substr(basename(NDVIfiles[i]),5,12), "%Y%m%d")
}

if (Start_Date == "") {
  Start_Date = dates[[1]]
} else {
  Start_Date <- as.Date(Start_Date, "%Y%m%d")
  if (is.na(Start_Date)) {
    Start_Date = dates[[1]]
    rp.messagebox("Invalid input parameter: Format of 'Start date' is not valid. Earliest scene", Start_Date, "is choosen as start date instead.")
  }
}

if (End_Date == "") {
  End_Date = tail(dates, n=1)[[1]]
} else {
  End_Date <- as.Date(End_Date, "%Y%m%d")
  if (is.na(End_Date)) {
    End_Date = tail(dates, n=1)[[1]]
    rp.messagebox("Invalid input parameter: Format of 'End date' is not valid. Earliest scene", End_Date, "is choosen as end date instead.")
  }
}

if (End_Date < Start_Date) {
  rp.messagebox("Invalid input parameter: 'Start date' must be earlier than 'End date'.")
}


# Filter dates by start and end date
filteredDates <- list()
for (d in 1:length(dates)) {
  if ((Start_Date <= dates[[d]]) && (dates[[d]] <= End_Date)) {
    filteredDates <- c(filteredDates, dates[d])
  }
}

if (length(filteredDates) == 0) {
  rp.messagebox("Invalid input data: No scenes available within given time period. Adjust start and end date parameters.")
  stop()
} else {
  dates <- filteredDates
}


# Topographic wetness index
# ----------------------------------------------------------------------------

# Raster for resampling
ras <- raster(NDVIfiles[1])

# Search TWI mask files
search_pattern <- paste0(substr(Title,1,2),'*TWIbinary_resampled.tif$')
pathTWImask <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
                          include.dirs = TRUE, recursive = TRUE)

# Check if twi mask matches resolution and extent of sample index
if (length(pathTWImask) >= 1) {
  pathTWImask <- pathTWImask[1]
  twimask <-raster(pathTWImask)
  if ( (twimask@extent != ras@extent) || (res(twimask) != res(ras)) ) {
    remove(twimask)
    pathTWImask <- c()
  }
}

if (length(pathTWImask)==0) {
  search_pattern <- paste0('*TWIbinary.tif')
  pathTWImask <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
                            include.dirs = TRUE, recursive = TRUE)
  if (length(pathTWImask)==0) {
    stop("No TWI mask file found.")
  } else if (length(pathTWImask)>1) {
    stop("Too many TWI mask files found.")
  }else {
    print("Resampling TWI mask ...")
    twimask <- raster(pathTWImask)
    twimask <- projectRaster(twimask, ras, method="ngb")
    pathTWImask <- file.path(Directory_containing_TWI, paste0(substr(Title,1,2), "_TWIbinary_resampled.tif"))
    writeRaster(twimask, pathTWImask, format = "GTiff", datatype='INT1U', overwrite = TRUE)
  }
}



# 2.  Compute WATER and WETNESS masks
# -----------------------------------

main <- function(d) {
  
  rasterOptions(maxmemory=1000000000)
  removeTmpFiles(0.1)
  
  twimask <- raster(pathTWImask)
  
  # Get year and month
  year <- as.integer(format(d, "%Y"))
  month <- as.integer(format(d, "%m"))
  datestring <- format(d, "%Y%m%d")
  prettydatestring <- format(d, "%Y-%m-%d")
  
  if (SAR) {
    SARmonth <- as.integer(month / 3) * 3
    if (SARmonth == 0) {
      SARmonth <- 12
      SARyear <- year - 1
    } else {
      SARyear <- year
    }
  }
  
  outfile_name <- paste0(Title, "_d", datestring)
  
  # OPTICAL BASED WATER MASK
  # ------------------------
  
  # Water probability
  pred_water <- NA
  pred_water <- aggregateIndices(indices_water, Directory_containing_indices, datestring)
  
  if (is.null(pred_water)) {
    cat("WARNING: No indices for scene",prettydatestring, "\n")
    return()
  }
  
  validFraction <- cellStats(!is.na(pred_water), "sum") / (pred_water@ncols * pred_water@nrows)
  if (validFraction < Minimum_AOI_coverage) {
    cat("WARNING: Not enough valid pixels for scene", prettydatestring, ". Scene is skipped.\n")
    return()
  }
  
  # Create water mask
  watermask <- NA
  watermask <- getWaterMask(pred_water, Tile_size_water, Lower_than, Minimum_water_probability, datestring, Output_Directory_site)
  
  if (!(typeof(watermask[[1]])=="S4")) {
    watermask <- getWaterMask(pred_water, Tile_size_water2, Lower_than, Minimum_water_probability, datestring, Output_Directory_site)
    if (!(typeof(watermask[[1]])=="S4")) {
      cat("WARNING: No water mask derived for", prettydatestring,"because no valid water threshold was found.\n")
      return()
    }
  }
  
  # Quality indicator
  if (Plot_certainty_indicator) {
    qualityWater <- abs(pred_water - watermask[[2]]) / 10.
    outfile_qualityWater <- file.path(Output_Directory_site, paste0(outfile_name, "_water_certainty.tif"))
    rf <- writeRaster(qualityWater, filename = outfile_qualityWater, format = "GTiff", datatype='INT2S', overwrite = TRUE)
  }
  
  
  # RADAR BASED WATER MASK
  # ----------------------
  
  if (SAR) {
    search_pattern <- paste0('*', sprintf("%04d",SARyear), sprintf("%02d",SARmonth), '*SFRQWATER*.tif$')
    SARfiles <- list.files(Directory_containing_SAR, pattern = glob2rx(search_pattern), full.names = TRUE,
                           include.dirs = TRUE, recursive = TRUE)
    
    if (length(SARfiles)!=0) {
      
      if (length(SARfiles) > 1) {
        SARimgs <- c()
        for (f in SARfiles) {
          img <- raster(f)
          img <- projectRaster(img, watermask[[3]], method="bilinear")
          SARimgs <- c(SARimgs, img)
        }
      
        rasters1.mosaicargs <- SARimgs
        rasters1.mosaicargs$fun <- mean
        waterFreq_SAR <- do.call(raster::mosaic, rasters1.mosaicargs)
      } else {
        waterFreq_SAR <- raster(SARfiles[1])
        waterFreq_SAR <- projectRaster(waterFreq_SAR, watermask[[3]], method="bilinear")
      }

      waterFreq_SAR[waterFreq_SAR==0] <- 0.1
      # Get joined water mask
      watermaskSAR <- getWaterMask_SAR(waterFreq_SAR, 100, FALSE, 0, "SAR", Output_Directory_site)
      if (!(typeof(watermaskSAR[[1]])=="S4")) {
        watermaskSAR <- getWaterMask_SAR(waterFreq_SAR, 50, FALSE, 0, "SAR", Output_Directory_site)
      }

      if ((typeof(watermaskSAR[[1]])=="S4")) {
        # Set pixels outside of potential water body mask to 0 (dry)
        watermaskSAR[[1]][watermask[[3]]==0] <- 0
        # Join all masks
        watermask[[1]][watermaskSAR[[1]]==1] <- 1
      }

    } else {
      cat("No SAR water frequency found for ", year, " - ", month)
    }
  }
  
  # Apply TWI mask
  if (!is.null(twimask)) {
    twimask_sub <- crop(twimask, pred_water)
    watermask[[1]][is.na(watermask[[1]])] <- 255
    try(watermask[[1]] <- overlay(twimask_sub, watermask[[1]], fun=function(x, y) ifelse((x==1) & (y!=255), 0, y)), silent=F)
    watermask[[1]][watermask[[1]]==255] <- NA
  }
  
  # SAVE TO FILE
  # ------------
  
  if (Plot_water_probability) {
    outfile_predWater <- file.path(Output_Directory_site, paste0(outfile_name, "_water_probability.tif"))
    rf <- writeRaster(pred_water / 10., filename = outfile_predWater, format = "GTiff", datatype='INT2U', overwrite = TRUE)
  }
  
  outfile_waterMask <- file.path(Output_Directory_site, paste0(outfile_name, '_watermaskRaw.tif'))
  rf <- writeRaster(watermask[[1]], filename = outfile_waterMask, format = "GTiff", datatype='INT1U', overwrite = TRUE)
  
  # Sieve mask using gdal
  outfile_sieveRaw = file.path(Output_Directory_site, paste0(outfile_name, '_watermask_sieveRaw.tif'))
  cmd = paste(paste0(pythonPath, "gdal_sieve"), "-st", mmu_water,"-8", outfile_waterMask, outfile_sieveRaw, sep=" ")
  sucess <- try(system(cmd))
  
  # Convert
  outfile_sieve = file.path(Output_Directory_site, paste0(outfile_name, '_watermask.tif'))
  cmd = paste('gdal_translate -ot Byte -of GTiff -co COMPRESS=LZW',  outfile_sieveRaw, outfile_sieve, sep=" ")
  system(cmd)
  
  # Delete sieved file
  file.remove(outfile_sieveRaw)
  file.remove(outfile_waterMask)
  
}

for (d in dates) {
  tryCatch(main(d), error = function(e) {cat(message(e))})
}

# Over all quality
if (Plot_certainty_indicator) {
  search_pattern <- paste0('*water_certainty.tif$')
  qualityWaterFiles <- list.files(Output_Directory_site, pattern = glob2rx(search_pattern), full.names = TRUE,
                                  include.dirs = TRUE, recursive = TRUE)
  qualityRas <- c()
  for (f in qualityWaterFiles) {
    ras <- raster(f)
    qualityRas <- c(qualityRas, ras)
  }
  
  qualityRas <- brick(qualityRas)
  valid <- sum(!is.na(qualityRas))
  certainty <- sum(qualityRas, na.rm=T) / valid
  
  outfile_quality_water <- file.path(Output_Directory_site, paste0(Title, "_water_overallCertainty.tif"))
  rf <- writeRaster(qualityRas, filename = outfile_quality_water, format = "GTiff", datatype='INT2S', overwrite = TRUE)
}

# Benchmarking
endtime = Sys.time()
delta = proc.time() - starttime
deltaStr = paste(floor(delta[3] / 60.)," min, ", round((delta[3] %% 60.),2)," sec"," ")
logfile = file.path(Output_Directory, paste0("Benchmark_WaterDetection_",substr(Sys.time(),1,10),"-",substr(Sys.time(),12,13),substr(Sys.time(),15,16),".txt"))
sink(logfile)
cat("Execution time: ",deltaStr, "\n\n")
cat("Directory_containing_indices:", Directory_containing_indices,"\n")
cat("Directory_containing_TWI", Directory_containing_TWI,"\n")
cat("Directory_containing_SAR", Directory_containing_SAR,"\n")
cat("Output_Directory" ,Output_Directory,"\n")
cat("Plot_water_probability" ,Plot_water_probability,"\n")
cat("Minimum_water_probability" ,Minimum_water_probability,"\n" )
cat("Minimum_AOI_coverage", Minimum_AOI_coverage,"\n")
cat("Start_Date", Start_Date,"\n")
cat("End_Date" ,End_Date,"\n")
cat("Tile_size_in_meter", Tile_size_in_meter,"\n")
cat("Minimum_mapping_unit" ,Minimum_mapping_unit,"\n")
cat("Plot_certainty_indicator", Plot_certainty_indicator,"\n")
sink()

rp.messagebox("Water detection was successful!")
removeTmpFiles()
sink()
