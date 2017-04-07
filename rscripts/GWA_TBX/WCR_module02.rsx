# Copyright (c) 2016, GeoVille Information Systems GmbH
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.
#
# """
# Random Forest model to identify relevant spectral indices for water and wetness detection
#
#
# Keyword arguments:
#   - pathIN:           input directory of spectral indices
#   - pathSHP:          path to shapefile with with training data
#   - tileID:           ID of test site
#   - season:           season for which indices should be used
#   - classCol:         name of column in shapefile that indicates class
#
# Output:
#   - None
#
#
# Dependencies:
#
# Date created: 6.04.2016
# Date last modified: 9.01.2017
#
# """
#
# __author__ = "Christina Ludwig"
# __version__ = "1.0"

# INPUT PARAMTERS
# ===========================================================

##Directory_containing_indices=folder
#Shapefile_with_AOI=file
##Directory_containing_TWI=folder
##Output_Directory=folder

pythonPath = ""
#Directory_containing_indices="I:/temp/QGIStool/test2/indices"
#Shapefile_with_AOI="I:/2687_GW_A/02_Interim_Products/Toolbox/sampleData/site77.shp"
#Directory_containing_TWI="X:/02_Interim_Products/run_09/site98/TWI"
#Output_Directory=NULL


# install.packages("raster")
# install.packages("maptools")
# install.packages("stats")
# install.packages("randomForest")
# install.packages("caret")
# install.packages("corrplot")
# install.packages("stringr")
# install.packages("cluster")
# install.packages("glmnet")
# install.packages("rgdal")
# install.packages("stabs")
# install.packages("glmnet")
# install.packages("tcltk")
# install.packages("pastecs")
# install.packages("e1071")

library(raster)
library(maptools)
library(stats)
library(stringr)
library(rgdal)
library(pastecs)
library(plyr)

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script.dir)

# FUNCTIONS

# Rescale to same scale
rescale <- function (index) {
maxVal <- cellStats(index, "max", na.rm=T) # quantile(index, probs=c(0.99)) # cellStats(index, "max", na.rm=T) # quantile(index, probs=c(0.99)) # cellStats(index, "min", na.rm=T) #
minVal <- cellStats(index, "min", na.rm=T) #quantile(index, probs=c(0.01)) #  #  # cellStats(index, "max", na.rm=T) # # index@data@min
index[index>maxVal] <- maxVal
index[index<minVal] <- minVal
index.scaled <- (index-minVal)/(maxVal-minVal) * 1000.
#index.sigmoid <- 1. - 1./(1. + exp(1.*index.scaled))
return(index.scaled)
}

rescale_perc <- function (index) {
maxVal <- quantile(index, probs=c(0.99)) # cellStats(index, "max", na.rm=T) #  # cellStats(index, "max", na.rm=T) # quantile(index, probs=c(0.99)) # cellStats(index, "min", na.rm=T) #
minVal <- quantile(index, probs=c(0.01)) # cellStats(index, "min", na.rm=T) #  # # index@data@min
index[index>maxVal] <- maxVal
index[index<minVal] <- minVal
index.scaled <- (index-minVal)/(maxVal-minVal) * 1000.
#index.sigmoid <- 1. - 1./(1. + exp(1.*index.scaled))
return(index.scaled)
}

rescale01 <- function (index) {
maxVal <- cellStats(index, "max", na.rm=T) #quantile(index, probs=c(0.99)) # cellStats(index, "min", na.rm=T) # #
minVal <- cellStats(index, "min", na.rm=T) #quantile(index, probs=c(0.01)) # cellStats(index, "max", na.rm=T) #  #
index[index>maxVal] <- maxVal
index[index<minVal] <- minVal
index.scaled <- (index-minVal)/(maxVal-minVal)
#index.sigmoid <- 1. - 1./(1. + exp(1.*index.scaled))
return(index.scaled)
}

# Tile-based thresholding
tileThreshold <- function(img, splitTileSize, thresLimit, plotname, na.rm = TRUE) {

if (plot_opt) {
pdf(file=file.path(Output_Directory, paste0(plotname, "_tileSize", splitTileSize, ".pdf")), onefile = T)
}

# compute statistics (STDV) for each image tile
nTilesY <- floor(nrow(img) / splitTileSize)
aspectratio <- max(1, round(nrow(img)/ncol(img)), na.rm = TRUE)
nTilesX <- floor(ncol(img) / splitTileSize)
#thr <- vector("numeric", nTilesY * nTilesX)
tilestats <- as.data.frame(matrix(NA, nr=nTilesY*nTilesX, 4))
names(tilestats) <- c('i','j','SDV','fracNA')
for (i in 1:nTilesY) {
cat("row",i,"...\n")
rowmin <- max(1,(i-1)*splitTileSize+1)
m <- getValues(img, row = rowmin, nrows = splitTileSize, format = "matrix")
m[m == 0] <- NA
for (j in 1:nTilesX) {
colmin <- max(1,(j-1)*splitTileSize/aspectratio+1)
colmax <- min(ncol(img),(j-1)*splitTileSize/aspectratio+splitTileSize/aspectratio)
v <- as.vector(m[,colmin:colmax])
k <- (i-1)*nTilesX+j
if (any(!is.na(v))) {
tilestats[k,1] <- i
tilestats[k,2] <- j
tilestats[k,3] <- sd(v, na.rm = TRUE)
tilestats[k,4] <- sum(is.na(v)) / sum(is.finite(v))
}
}
}

# find optimal threshold searching for two-peak-distributions
# (1) search for SDV > 95% percentile
id_good <- which(is.finite(rowSums(tilestats)))
if ( length(id_good) > 0 ){
tilestats <- tilestats[id_good, ]
} else {
cat('WARNING: No valid image values found -> RETURN TO MAIN PROGRAM\n')
dev.off()
return(1)
}
id <- which((tilestats$SDV > quantile(tilestats$SDV, 0.95)) & (tilestats$fracNA < 0.1))
ii <- tilestats$i[id]
jj <- tilestats$j[id]
water_thres <- numeric(length = length(id)) * NA
frac_peaks <- numeric(length = length(id)) * NA
for (ix in 1:length(id)){
i <- ii[ix]
j <- jj[ix]
rowmin <- max(1,(i-1)*splitTileSize+1)
m <- getValues(img, row = rowmin, nrows = splitTileSize, format = "matrix")
m[m == 0] <- NA
colmin <- max(1,(j-1)*splitTileSize/aspectratio+1)
colmax <- min(ncol(img),(j-1)*splitTileSize/aspectratio+splitTileSize/aspectratio)
v_db <- as.vector(m[,colmin:colmax])
dens <- density(v_db, na.rm = TRUE)
if ( plot_opt ){
par(fig=c(0, 1, 0, 0.6))
sub_img <- matrix(data = v_db, nrow = splitTileSize, ncol = splitTileSize, byrow = FALSE)
image(sub_img, col = grey(c(1:255)/255), asp = 1)
#par(fig=c(0, 1, 0, 0.7), new=TRUE)
par(fig=c(0.1, 0.9, 0.5, 1), new=TRUE)
q2_98 <- quantile(v_db, c(0.02,0.98), na.rm = TRUE)
plot(dens, xlim = c(q2_98[1], q2_98[2]), main="")
}

# Determine local minimum of the two peak distribution to identify threshold
tp <- turnpoints(dens$y)
id_pits <- which(tp$pits == TRUE)

# Use only those TPs which are LE than predefined limit (thresLimit)
good_TP <- which(dens$x[id_pits] >= thresLimit)
if ( !length(good_TP) ){
cat('No valid minimum found for tile ', i, '/', j, '\n')
next
}
id_pits <- id_pits[good_TP]
n_pits <- length(id_pits)

if (n_pits == 1){
water_thres[ix] <- dens$x[id_pits]
frac_peaks[ix] <- sum(dens$y[1:id_pits]) / sum(dens$y[id_pits:length(dens$y)])
if ( plot_opt ){
abline(v = dens$x[id_pits], col = 'red')
text(q2_98[2]*0.9, max(dens$y)*0.9, labels = paste0('Thres: ', format(water_thres[ix], digits = 4)))
#par(fig=c(0, 1, 0, 0.6), new=TRUE)
#sub_img <- matrix(data = v_db, nrow = splitTileSize, ncol = splitTileSize, byrow = FALSE)
#image(sub_img, col = grey(c(1:255)/255), asp = 1)
}
} else if (n_pits > 1) {
# Determine the highest pit
tmp_dens_y <- dens$y
max_min_id <- which.max(dens$y[id_pits])

# Check whether broad pit or isolated 'mini-pit' (check 20 values left and right to the pit)
wd <- 30
loLim <- ifelse(id_pits[max_min_id] - wd > 0, id_pits[max_min_id] - wd, 1)
upLim <- ifelse(id_pits[max_min_id] + wd <= length(dens$y), id_pits[max_min_id] + wd, length(dens$y))
pos_min <- which.min(dens$y[loLim:upLim]) - wd
if (pos_min != 1){
tmp_dens_y[id_pits[max_min_id]] <- 0.
max_min_id <- which.max(tmp_dens_y[id_pits])
loLim <- ifelse(id_pits[max_min_id] - wd > 0, id_pits[max_min_id] - wd, 1)
upLim <- ifelse(id_pits[max_min_id] + wd <= length(dens$y), id_pits[max_min_id] + wd, length(dens$y))
pos_min <- which.min(tmp_dens_y[loLim:upLim]) + wd
if (pos_min != 1){
cat('WARNING: Unable to determine appropriate local minimum.\n')
next
}
}
water_thres[ix] <- dens$x[id_pits[max_min_id]]
frac_peaks[ix] <- sum(dens$y[1:id_pits[max_min_id]]) / sum(dens$y[id_pits[max_min_id]:length(dens$y)])
if ( plot_opt ){
abline(v = dens$x[id_pits[max_min_id]], col = 'red')
text(q2_98[2]*0.9, max(dens$y)*0.9, labels = paste0('Thres: ', format(water_thres[ix], digits = 4)))
#par(fig=c(0, 1, 0, 0.6), new=TRUE)
#sub_img <- matrix(data = v_db, nrow = splitTileSize, ncol = splitTileSize, byrow = FALSE)
#image(sub_img, col = grey(c(1:255)/255), asp = 1)
}
}
if ( debug_mode ){ browser() }
}

# Filter tiny peaks
id_good_peaks <- which(frac_peaks > 0.15)
if ( !length(id_good_peaks) ){
cat('WARNING: No appropriate threshold found using tile-based method.\n')
dev.off()
return(1)
}
water_thres <- water_thres[id_good_peaks]

if (plot_opt) {
if (length(water_thres) != 0) {
# Histogram of all thresholds
df <- data.frame(thresholds=water_thres)
mean_thres <- mean(df$thresholds)
median_thres <- median(df$thresholds)
sd_thres <- sd(df$thresholds)
min_thres <- mean_thres - sd_thres
max_thres <- median_thres + sd_thres

dev.off()

h <- ggplot(data=df, aes(thresholds))
h + geom_density() +
geom_vline(xintercept=mean_thres, color="red")+ geom_vline(xintercept=median_thres, color="blue") +
geom_vline(xintercept=min_thres, color="green")+ geom_vline(xintercept=max_thres, color="green") +
theme(plot.margin = unit(c(1.5,1,1,1), "cm")) +
labs(title="Histogram of local thresholds", x="Thresholds") +
annotate("text", label = paste0("Mean: ", round(mean_thres,2), "\n Median: ", round(median_thres,2),"\n StD: ", round(sd_thres,2)), x = mean_thres-50, y = max(density(water_thres)$y)-0.001, size = 4, colour = "black")

ggsave(file.path(Output_Directory, paste0(plotname, "_histogram_thresholds_", splitTileSize, ".png")), plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 25, height = 15, units = c("cm"), dpi = 300)
}
}

return(water_thres)
}

getWaterMask <- function(Layer, splitTileSize, threshold, thresLimit_water,d, Output_Directory) {

water_thres <- NA

plotname <- d
try(water_thres <- tileThreshold(Layer, splitTileSize, thresLimitWater, plotname))

if ((length(water_thres) == 1) | sum(is.na(water_thres)) > 0 ) {
return(NA)
}

# use (mean - 1 * STDV) as threshold
theta <- median(water_thres, na.rm = TRUE)#  sd(water_thres, na.rm = TRUE)
#theta <- median(water_thres, na.rm = TRUE)# - sd(water_thres, na.rm = TRUE)

# Create watermask
if (!is.finite(theta)) {
cat("No threshold could be computed.\n")
cat("Returning empty raster.\n")
watermask <- raster(Layer)
return(1)
} else {
if (Lower_than) {
watermask <- Layer < theta
cat("Automatically computed threshold:", theta, "\n")
#write(paste(d, theta, sep=";"), file=file.path(Output_Directory, "thresholds.txt"), append=T, sep=";")
} else {
watermask <- Layer > theta
cat("Automatically computed threshold:", theta, "\n")
#write(paste(d, theta, sep=";"), file=file.path(Output_Directory, "thresholds.txt"), append=T, sep=";")
}
}
return(list(watermask, theta))
}

# Sigmoid
sigmoid <- function (index) {
index.scaled <- index / 100. -5
index.sigmoid <- (1. - 1./(1. + exp(1. * index.scaled))) * 1000.
return(index.sigmoid)
}

# Sigmoid
sigmoid01 <- function (index) {
maxVal <- cellStats(index, "max", na.rm=T) # quantile(index, probs=c(0.99)) # cellStats(index, "max") #
minVal <- cellStats(index, "min", na.rm=T) # quantile(index, probs=c(0.01)) # cellStats(index, "min") # #  # # index@data@min
index[index>maxVal] <- maxVal
index[index<minVal] <- minVal
index.scaled <- (index-minVal)/(maxVal-minVal) * 10.0 - 5.0
index.sigmoid <- (1. - 1./(1. + exp(1. * index.scaled)))
return(index.sigmoid)
}


aggregateIndices_water <- function(indices_water, pathIndices, d) {

Nodata_value <- -32768
pred <- NA
bands <- c()

# Read in indices
for (idx in indices_water) {

# Load  file
rasterFile = NULL
search_pattern <- paste0(d, '_', idx, '.tif', sep="")
rasterFile <- list.files(pathIndices, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

if (length(rasterFile) == 0) {
return(NULL)
}

cat("Loading: ", basename(rasterFile), "\n")
index_band <- raster(rasterFile)
NAvalue(index_band) <- Nodata_value

# Weight with TWI
index_band <- index_band * 0.1

# Reverse value range if necessary
if (idx %in% reverse) {
index_band <- 1000. - index_band
}

# Sum up with previous indices
bands <- c(bands, index_band)

}

pred <- Reduce("+", bands) / length(indices_water)
return(pred)

}



#rootPath <- dirname(Directory_containing_indices)
#inputFolder <- basename(Directory_containing_indices)

if (!file.exists(Output_Directory)) {
Output_Directory <- dirname(Directory_containing_indices)
}

Output_Directory <- file.path(Output_Directory, "watermasks")

if (!file.exists(Output_Directory)) {
dir.create(Output_Directory)
}

reverse <- c("ND068A")

# Water index combination
indices_water <- c("ND0311", "ND038A")

tileID <- ""

# TWI
# Example file for resampling
search_pattern <- paste0('*', '_ND0311.tif')
rasterFiles <- list.files(Directory_containing_indices, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)
if (length(rasterFiles) == 0) {
stop("No NDVI found.")
}

# Get dates of scenes
dates <- list.files(Directory_containing_indices, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

# One exemplatory raster to reproject twi
ras <- raster(rasterFiles[1])

# Search for TWI
search_pattern <- paste0('*TWIbinary_resampled.tif')
pathTWI <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)

if (length(pathTWI)==0) {
search_pattern <- paste0('*TWIbinary.tif')
pathTWI <- list.files(Directory_containing_TWI, pattern = glob2rx(search_pattern), full.names = TRUE,
include.dirs = TRUE, recursive = TRUE)
if (length(pathTWI)==0) {
stop("TWI not found.")
}
twi <- raster(pathTWI)
twi <- projectRaster(twi, ras)
twi <- resample(twi, ras, method="ngb")
twi <- crop(twi, ras)
#twi <- rescale01(twi)
writeRaster(twi, file.path(Directory_containing_TWI, paste0(tileID, "_TWIbinary_resampled.tif")))
} else {
print(pathTWI)
twi <- raster(pathTWI)
}

remove(ras)


# 2.  Compute WATER masks

# SETTINGS for Otsu Thresholding
plot_opt = FALSE
Tile_size_water = 100
Tile_size_water2 = 50
mmu_water <- 3
Nodata_value = -32768
Nodata_value_mask = 255
Lower_than = FALSE
thresLimitWater = 500
debug_mode = FALSE

for (f in dates) {

date <- substr(basename(f),1,8)

# delete pred_water if already exists
pred_water <- NA
pred_water <- aggregateIndices_water(indices_water, Directory_containing_indices, date)

if (is.null(pred_water)) {
print(paste("No indices for  ",date, sep=" ") )
next()
}

validFraction <- cellStats(!is.na(pred_water), "sum") / (pred_water@ncols * pred_water@nrows)

if (validFraction < 0.2) {
print("Not enough observations in composite.")
next()
}

# Save to file
#NAvalue(pred_water) <- Nodata_value

outfile_name <- date #paste(tileID, "_",date, sep="_")
print(outfile_name)
outfile_predWater <- file.path(Output_Directory, paste0(outfile_name, "_predWater.tif"))
rf <- writeRaster(pred_water, filename = outfile_predWater, format = "GTiff", datatype='INT2U', overwrite = TRUE)

twi_sub <- crop(twi, pred_water)

# Get watermask
watermask <- NA
watermask <- getWaterMask(pred_water, Tile_size_water, Lower_than, thresLimitWater, date, Output_Directory)

if (!(typeof(watermask[[1]])=="S4")) {
cat("Try smaller tile size ...")
watermask <- getWaterMask(pred_water, Tile_size_water2, Lower_than, thresLimitWater, date, Output_Directory)
if (!(typeof(watermask[[1]])=="S4")) {
next()
}
}

# Cut out TWI mask
watermask2 <- watermask[[1]]
watermask2[is.na(watermask2)] <- 255
try(watermask2 <- overlay(twi_sub, watermask2, fun=function(x, y) ifelse((x==1) & (y!=255), 0, y)), silent=F)
watermask2[watermask2==255] <- NA
watermask[[1]] <- watermask2

# Plot
plot(watermask[[1]])

# Save to file
#NAvalue(watermask) <- Nodata_value_mask
outfile_waterMask <- file.path(Output_Directory, paste0(outfile_name, '_watermask.tif'))
rf <- writeRaster(watermask[[1]], filename = outfile_waterMask, format = "GTiff", datatype='INT1U', overwrite = TRUE)

# Sieve mask using gdal
outfile_sieveRaw = file.path(Output_Directory, paste0(outfile_name, '_watermask_sieveRaw.tif'))
cmd = paste(paste0(pythonPath, "gdal_sieve"), "-st", mmu_water,"-8", outfile_waterMask, outfile_sieveRaw, sep=" ")
system(cmd)

# Convert
outfile_sieve = file.path(Output_Directory, paste0(outfile_name, '_watermask_sieve.tif'))
cmd = paste('gdal_translate -ot Byte -of GTiff -co COMPRESS=LZW',  outfile_sieveRaw, outfile_sieve, sep=" ")
system(cmd)

# Delete sieved file
file.remove(outfile_sieveRaw)
file.remove(outfile_waterMask)

}
