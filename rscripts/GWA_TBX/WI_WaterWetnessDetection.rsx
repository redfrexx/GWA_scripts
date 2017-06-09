#  Copyright (c) 2017, GeoVille Information Systems GmbH
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.
#
#
# Date created: 06.05.2017
# Date last modified: 09.06.2017
#
#
# __author__ = "Christina Ludwig"
# __version__ = "1.0"


# INPUT PARAMTERS

# -> FOR QGIS Processing modules
##Water and Wetness Detection=name
##Wetland Inventory=group
##Directory_containing_indices=folder
##Directory_containing_TWI=folder
##Output_Directory=folder
##Directory_containing_SAR=folder
##Start_Date=optional string
##End_Date=optional string
##Minimum_water_probability=number 45
##Minimum_wetness_probability=number 50
##Tile_size_in_meter=optional number 1800
#Threshold_calculation=selection Median;Mean
##Minimum_mapping_unit=number 3
##Minimum_AOI_coverage=number 40
#Increase_sensitivity=Boolean False
##Plot_water_and_wetness_probability= Boolean False
##Plot_certainty_indicator= Boolean False

debug <- F
starttime <- proc.time()

# -> Test parameters for execution in R
if (debug) {
Directory_containing_indices= "T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\WI\\indices"
Directory_containing_TWI= "T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\TWI"
Directory_containing_SAR= "" # "T:/Processing/2687_GW_A/02_Interim_Products/SAR/site_98/EQUI7_AF010M"
Output_Directory="T:\\Processing\\2687_GW_A\\03_Products\\GWA-TOOLBOX\\02_InterimProducts\\test"
Plot_water_and_wetness_probability <- T
Minimum_water_probability = 50
Minimum_wetness_probability = 50
Minimum_AOI_coverage <- 40
Start_Date <- " "
End_Date <- " "
Increase_sensitivity=FALSE
Tile_size_in_meter = 1800
Minimum_mapping_unit = 3
Global_threshold_based_on = 0
Plot_certainty_indicator = T
.libPaths("C:\\Users\\ludwig\\.qgis2\\processing\\rlibs")
}

#if("raster" %in% rownames(installed.packages()) == FALSE) {install.packages("raster")}
#if("rgdal" %in% rownames(installed.packages()) == FALSE) {install.packages("rgdal")}
#if("pastecs" %in% rownames(installed.packages()) == FALSE) {install.packages("pastecs")}
#if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}
#if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel")}
#if("foreach" %in% rownames(installed.packages()) == FALSE) {install.packages("foreach")}
#install.packages("O:/2687_GW_A/04_CODE/R/compiledPackages/GWAutils.zip", repos=NULL, type="win.binary")

library(raster)
library(rgdal)
library(GWAutils)
library(rpanel)
library(stringr)

Start_Date <- str_trim(Start_Date)
End_Date <- str_trim(End_Date)

Minimum_AOI_coverage <- Minimum_AOI_coverage / 100.
Minimum_water_probability <- Minimum_water_probability * 10.
Minimum_wetness_probability <- Minimum_wetness_probability * 10.
Increase_sensitivity = F

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
indices_wet <- c("NDVI", "NDMI", "ND-NIR-SWIR2")

# Indices that need to be reversed so that water / wetness has high values
reverse <- c("NDVI")

#default
pythonPath = ""


# Set raster settings ----------------------------------------------------

if (!debug) {
sink(file.path(Output_Directory, "log.txt"))
}

removeTmpFiles(h=0.25)

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

Title <- paste(substr(basename(NDVIfiles[1]), 1,2), "_wat",Minimum_water_probability/10., "_wet", Minimum_wetness_probability/10., "_win", Tile_size_in_meter, "_mmu", Minimum_mapping_unit, sep="")

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
Tile_size_wet = as.numeric(Tile_size_in_meter) /  resolutions[1]
Tile_size_wet2 = (as.numeric(Tile_size_in_meter)*0.66) /  resolutions[1]

if (Tile_size_in_meter < 1000) {
rp.messagebox("Invalid input parameter: Tile size must be greater than 1000m.")
stop()
} else if ((Tile_size_wet > ndvi@ncols) || (Tile_size_wet > ndvi@ncols) || (Tile_size_water > ndvi@nrows) || (Tile_size_water > ndvi@nrows)) {
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
stop()
}
}

if (End_Date == "") {
End_Date = tail(dates, n=1)[[1]]
} else {
End_Date <- as.Date(End_Date, "%Y%m%d")
if (is.na(End_Date)) {
End_Date = tail(dates, n=1)[[1]]
rp.messagebox("Invalid input parameter: Format of 'End date' is not valid. Earliest scene", End_Date, "is choosen as end date instead.")
stop()
}
}

if (End_Date < Start_Date) {
rp.messagebox("Invalid input parameter: 'Start date' must be earlier than 'End date'.")
stop()
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


# TWI
search_pattern <- paste0(substr(Title,1,2),'*TWI_resampled.tif$')
pathTWI <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

# Check if existing resampled TWI matches the spectral indices
if (length(pathTWI) >= 1) {
pathTWI <- pathTWI[1]
twi <-raster(pathTWI)
if ( (twi@extent != ras@extent) || (res(twi) != res(ras)) ) {
pathTWI <- c()
remove(twi)
}
}

# Resample TWI
if (length(pathTWI)==0) {
search_pattern <- paste0('*TWI.tif')
pathTWI <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)
if (length(pathTWI)==0) {
stop("No TWI file found.")
} else if (length(pathTWI)>1) {
stop("Too many TWI files found. ")
}else {
print("Resampling TWI ...")
twi <- raster(pathTWI)
twi <- projectRaster(twi, ras, method="bilinear")
twi <- sigmoid01(twi)
pathTWI <- file.path(Directory_containing_TWI, paste0(substr(Title,1,2), "_TWI_resampled.tif"))
writeRaster(twi, pathTWI ,format = "GTiff", datatype='FLT4S', overwrite = TRUE)
}
}

#remove(ras)

# 2.  Compute WATER and WETNESS masks
# -----------------------------------

main <- function(d) {

rasterOptions(maxmemory=100000000)
removeTmpFiles(0.25)

twi <- raster(pathTWI)
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

# RADAR BASED WATER MASK
# ----------------------

if (SAR) {
search_pattern <- paste0('*', sprintf("%04d",year), sprintf("%02d",SARmonth), '*SFRQWATER*.tif$')
SARfiles <- list.files(Directory_containing_SAR, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

if (length(SARfiles)!=0) {

SARimgs <- c()
for (f in SARfiles) {
img <- raster(f)
img <- projectRaster(img, watermask[[3]], method="bilinear")
SARimgs <- c(SARimgs, img)
}

rasters1.mosaicargs <- SARimgs
rasters1.mosaicargs$fun <- mean
waterFreq_SAR <- do.call(mosaic, rasters1.mosaicargs)

waterFreq_SAR <- projectRaster(waterFreq_SAR, watermask[[3]], method="bilinear")
waterFreq_SAR[waterFreq_SAR==0] <- 0.1
# Get joined water mask
watermaskSAR <- getWaterMask_SAR(waterFreq_SAR, 100, FALSE, 20, "SAR", Output_Directory_site)
if (!(typeof(watermaskSAR[[1]])=="S4")) {
watermaskSAR <- getWaterMask_SAR(waterFreq_SAR, 25, FALSE, 20, "SAR", Output_Directory_site)
if ((typeof(watermaskSAR[[1]])=="S4")) {
# Set pixels outside of potential water body mask to 0 (dry)
watermaskSAR[[1]][watermask[[3]]==0] <- 0
# Join all masks
watermask[[1]][watermaskSAR[[1]]==1] <- 1
}
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

# Quality indicator
if (Plot_certainty_indicator) {
qualityWater <- abs(pred_water - watermask[[2]]) / 10.
outfile_qualityWater <- file.path(Output_Directory_site, paste0(outfile_name, "_water_certainty.tif"))
rf <- writeRaster(qualityWater, filename = outfile_qualityWater, format = "GTiff", datatype='INT2S', overwrite = TRUE)
}

# SAVE TO FILE
# ------------

if (Plot_water_and_wetness_probability){
outfile_predWater <- file.path(Output_Directory_site, paste0(outfile_name, "_water_probability.tif"))
rf <- writeRaster(pred_water / 10., filename = outfile_predWater, format = "GTiff", datatype='INT2U', overwrite = TRUE)
}

outfile_waterMask <- file.path(Output_Directory_site, paste0(outfile_name, '_watermaskRaw.tif'))
rf <- writeRaster(watermask[[1]], filename = outfile_waterMask, format = "GTiff", datatype='INT1U', overwrite = TRUE)

# Sieve mask using gdal
outfile_sieveRaw = file.path(Output_Directory_site, paste0(outfile_name, '_watermask_sieveRaw.tif'))
cmd = paste("gdal_sieve", "-st", mmu_water,"-8", outfile_waterMask, outfile_sieveRaw, sep=" ")
sucess <- try(system(cmd))

# Convert
outfile_sieve = file.path(Output_Directory_site, paste0(outfile_name, '_watermask.tif'))
cmd = paste('gdal_translate -ot Byte -of GTiff -co COMPRESS=LZW',  outfile_sieveRaw, outfile_sieve, sep=" ")
system(cmd)

# Delete sieved file
file.remove(outfile_sieveRaw)
file.remove(outfile_waterMask)

# WETNESS MASK
# -------------------------------------------------------------------------

# Water mask
watermask_sieve <- raster(outfile_sieve)

# NDVI mask
search_pattern <- paste0('*',datestring,'*NDVI.vrt')
pathNDVI <- list.files(Directory_containing_indices, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)
NDVI <- raster(pathNDVI)

# OPTICAL BASED WETNESS MASK
# --------------------------

# Wetness probability
pred_wet <- NA
pred_wet <- aggregateIndices(indices_wet, Directory_containing_indices,datestring)

if (is.null(pred_wet)) {
cat("WARNING: No indices for wetness detection found for scene",prettydatestring, "\n")
return()
}



if (Increase_sensitivity == T) {
# Clip out detected water
pred_wet <- crop(pred_wet, watermask_sieve)
pred_wet <- mask(pred_wet, watermask_sieve, maskvalue=1)
if (!is.null(twi)) {   # Crop twi to indices
twi_sub <- crop(twi, pred_wet)
pred_wet <- sigmoid(pred_wet)
pred_wet <- pred_wet * (twi_sub * 0.75 + 0.25)
}
} else {
if (!is.null(twi)) {   # Crop twi to indices
twi_sub <- crop(twi, pred_wet)
pred_wet <- sigmoid(pred_wet)
pred_wet <- pred_wet * (twi_sub * 0.75 + 0.25)
}
# Clip out detected water
pred_wet <- crop(pred_wet, watermask_sieve)
pred_wet <- mask(pred_wet, watermask_sieve, maskvalue=1)
}

# Save to file
if (Plot_water_and_wetness_probability) {
outfile_predWet <- file.path(Output_Directory_site, paste0(outfile_name, "_wetness_probability.tif"))
rf <- writeRaster(pred_wet / 10., filename = outfile_predWet, format = "GTiff", datatype='INT2S', overwrite = TRUE)
}


# Create wetness mask
wetmask <- NA
wetmask <- getWetMask(pred_wet, Tile_size_wet, Lower_than, Minimum_wetness_probability,date , Output_Directory_site)

# Run again using smaller tile size if no threshold was found
if (!(typeof(wetmask[[1]])=="S4")) {
wetmask <- getWetMask(pred_wet, Tile_size_wet2, Lower_than, Minimum_wetness_probability, date, Output_Directory_site)
if (!(typeof(wetmask[[1]])=="S4")) {
cat("WARNING: No wetness mask derived for", prettydatestring,"because no valid wetness threshold was found.\n")
return()
}
}


# RADAR BASED WETNESS MASK
# ------------------------
if (SAR) {
search_pattern <- paste0('*', sprintf("%04d",year), sprintf("%02d",SARmonth), '*SFRQWET*.tif$')
SARfiles <- list.files(Directory_containing_SAR, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

if (length(SARfiles) != 0) {

SARimgs <- c()
for (f in SARfiles) {
img <- raster(f)
img <- projectRaster(img, watermask_sieve, method="bilinear")
SARimgs <- c(SARimgs, img)
}

rasters1.mosaicargs <- SARimgs
rasters1.mosaicargs$fun <- mean
waterFreq_SAR <- do.call(mosaic, rasters1.mosaicargs)

wetFreq_SAR <- mask(wetFreq_SAR, watermask_sieve, maskvalue=1)

# Merge SAR and OPTICAL Frequencies
wetFreq_SAR <- wetFreq_SAR * (pred_wet/1000.)

# Get SAR water mask
wetmaskSAR <- getWetMask_SAR(wetFreq_SAR, 100, FALSE, 30, "SAR", Output_Directory_site)

# Run again using smaller tile size if no threshold was found
if (!(typeof(wetmaskSAR[[1]])=="S4")) {
wetmaskSAR <- getWetMask_SAR(wetFreq_SAR, 25, FALSE, 30, "SAR", Output_Directory_site)
if ((typeof(wetmaskSAR[[1]])=="S4")) {
# Join masks
wetmask[[1]][wetmaskSAR[[1]]==1] <- 1
}
}

} else {
cat("No SAR wetness frequency found for ", year, " - ", month)
}
}

# Apply NDVI mask
wetmask[[1]][NDVI>NDVIthresh] <- 0

# Quality indicator
if (Plot_certainty_indicator) {
qualityWet <- abs(pred_wet - wetmask[[2]]) / 10.
outfile_qualityWet<- file.path(Output_Directory_site, paste0(outfile_name, "_wetness_certainty.tif"))
rf <- writeRaster(qualityWet, filename = outfile_qualityWet, format = "GTiff", datatype='INT2S', overwrite = TRUE)
}

# Save to file
if (Plot_water_and_wetness_probability) {
outfile_predWet <- file.path(Output_Directory_site, paste0(outfile_name, "_wetness_probability.tif"))
rf <- writeRaster(pred_wet / 10., filename = outfile_predWet, format = "GTiff", datatype='INT2S', overwrite = TRUE)
}

# Export mask
outfile_wetMask <- file.path(Output_Directory_site, paste0( outfile_name, "_", paste0(indices_wet, collapse=""), '_wetmaskRaw.tif'))
rf <- writeRaster(wetmask[[1]], filename = outfile_wetMask, format = "GTiff", datatype='INT1U', overwrite = TRUE)

# Sieve mask using gdal
outfile_sieveRaw = file.path(Output_Directory_site, paste0(outfile_name, "_", paste0(indices_wet, collapse=""), '_wetmask_sieveRaw.tif'))
cmd = paste("gdal_sieve", "-st", mmu_wet,"-8", outfile_wetMask, outfile_sieveRaw, sep=" ")
sucess <- try(system(cmd))

# Convert to change data type
outfile_sieve = file.path(Output_Directory_site, paste0(outfile_name, '_wetmask.tif'))
cmd = paste('gdal_translate -ot Byte -of GTiff -co COMPRESS=LZW',  outfile_sieveRaw, outfile_sieve, sep=" ")
system(cmd)

# Delete sieved file
file.remove(outfile_sieveRaw)
file.remove(outfile_wetMask)
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
rf <- writeRaster(certainty, filename = outfile_quality_water, format = "GTiff", datatype='INT2S', overwrite = TRUE)

# Wetness certainty
search_pattern <- paste0('*wetness_certainty.tif$')
qualityWetFiles <- list.files(Output_Directory_site, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)
qualityRas <- c()
validObs <- c()
for (f in qualityWetFiles) {
ras <- raster(f)
qualityRas <- c(qualityRas, ras)
}

qualityRas <- brick(qualityRas)
valid <- sum(!is.na(qualityRas))
certainty <- sum(qualityRas, na.rm=T) / valid

outfile_quality_wet <- file.path(Output_Directory_site, paste0(Title, "_wetness_overallCertainty.tif"))
rf <- writeRaster(certainty, filename = outfile_quality_wet, format = "GTiff", datatype='INT2S', overwrite = TRUE)

}

rp.messagebox("Water and wetness detection was successful!")
removeTmpFiles()

