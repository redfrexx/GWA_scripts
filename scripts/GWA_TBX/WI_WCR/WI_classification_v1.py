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

##Wetland Inventory Classification = name
##Wetland Inventory=group
##ParameterFile|pathIN|Directory containing water and wetness masks|True|False
##OutputDirectory|pathOUT|Output directory
##ParameterBoolean|exportSeasonalFrequencies|Export seasonal water and wetness frequencies
##*ParameterString|startDate|Start date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterString|endDate|End date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True

debug = False

# test paths
if debug:
    pathIN = r"I:\2687_GW_A\02_Interim_Products\Toolbox\02_InterimProducts\LS_waterWetnessMasks_wat450_wet500_win1800_mmu3"
    pathOUT = r"I:\2687_GW_A\02_Interim_Products\Toolbox\02_InterimProducts\LS_waterWetnessMasks_wat450_wet500_win1800_mmu3"
    exportSeasonalFrequencies = True
    startDate = ""
    endDate = ""

import os, sys
import numpy as np
import fnmatch
import gdal
from shutil import copyfile
import datetime as dt

if not debug:
    from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException

# Load additional library
#here = "."
here = os.path.dirname(scriptDescriptionFile)
pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)
	
import RSutils.RSutils as rsu

# Functions ------------------------------------------------------------------------------------------

def calculateFrequency(inFiles, extent):

    masks_fused = []
    for mask_OPT in inFiles:
        jointExt = rsu.getJointExtent(mask_OPT, AOIextent=extent)
        watermask_opt,geotrans = rsu.raster2array(mask_OPT, jointExt)[:2]
        if watermask_opt.shape[1] < extent.ncol or watermask_opt.shape[0] < extent.nrow:
            watermask_opt = rsu.padArray(watermask_opt, geotrans, extent)
        masks_fused.append(watermask_opt)

    masks = np.array(masks_fused).astype("uint8")

    # Number of valid observations
    validobs = np.nansum(masks!=255, axis=0).astype("float32")

    # Number of water/wetness detections
    Sum = np.nansum(masks==1, axis=0).astype("float32")

    # Water/wet frequency
    #Freq = Sum / validobs

    # Replace NAN
    Sum = np.where(validobs==0, np.nan, Sum)
    #Freq = np.where(validobs==0, np.nan, Freq)

    del masks

    return ( Sum, validobs)


# Check input parameters ------------------------------------------------------------------------------

if not os.path.exists(pathIN):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Directory containing water and wetness masks' does not exist.")
    print("Invalid input parameters: 'Directory containing water and wetness masks' does not exist.")

if not os.path.exists(pathOUT):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Output directory' does not exist.")
    print("Invalid input parameters: 'Output directory' does not exist.")

qmlDir = os.path.join(here, 'data', 'qml')

pathOUT_class = os.path.join(pathOUT,"classification_WI")
if not os.path.exists(pathOUT_class):
    os.mkdir(pathOUT_class)

# Seasonal frequencies
if exportSeasonalFrequencies:
    pathOUT_freqs = os.path.join(pathOUT_class,"seasonal_frequencies")
    if not os.path.exists(pathOUT_freqs):
        os.mkdir(pathOUT_freqs)

# WATER occurance and frequency =========================================================

# Search water masks
waterMaskFiles = [os.path.join(pathIN, f) for f in fnmatch.filter(os.listdir(pathIN), "*_watermask.tif")]
wetMaskFiles = [os.path.join(pathIN, f) for f in fnmatch.filter(os.listdir(pathIN), "*wetmask.tif")]

# Check whether masks exist
if len(waterMaskFiles) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No water masks found")
    else:
        print("No water masks found.")
else:
    print ("Found " + str(len(waterMaskFiles)) + " water mask files.\n")

if len(wetMaskFiles) == 0:
    raise GeoAlgorithmExecutionException("No wetness masks found")
else:
    print ("Found " + str(len(wetMaskFiles)) + " wet mask files.\n")

jointExtent = rsu.getJointExtent(waterMaskFiles)
geotrans, proj = rsu.raster2array(waterMaskFiles[0], jointExtent)[1:3]

waterFreqs = []
validPixels = []
wetFreqs = []

seasons = {"winter":[12,1,2],"spring":[3,4,5],"summer":[6,7,8],"fall":[9,10,11]}
for sname in seasons:
    print(sname)

    waterMasks_season = []
    for f in waterMaskFiles:
        f_name = os.path.basename(f)
        dateidx = f_name.find("_d")
        date = dt.datetime.strptime(f_name[dateidx+2:dateidx+10], "%Y%m%d")
        if (date.month in seasons[sname]):
            waterMasks_season.append(f)

    if len(waterMasks_season) == 0:
        validPixels.append(np.zeros((jointExtent.nrow, jointExtent.ncol)))
        continue

    waterSum, validObs = calculateFrequency(waterMasks_season, jointExtent)
    waterFreq = (waterSum / validObs) * 100.
    waterFreqs.append(waterFreq)
    validPixels.append(validObs)

    # wetness frequency
    wetMasks_season = []
    for f in wetMaskFiles:
        f_name = os.path.basename(f)
        dateidx = f_name.find("_d")
        date = dt.datetime.strptime(f_name[dateidx + 2:dateidx + 10], "%Y%m%d")
        if date.month in seasons[sname]:
            wetMasks_season.append(f)

    if len(wetMasks_season) == 0:
        continue

    wetSum = calculateFrequency(wetMasks_season, jointExtent)[0]
    wetFreq = (wetSum / validObs) * 100.
    wetFreqs.append(wetFreq)

    # Export to file
    if exportSeasonalFrequencies:
        dest = os.path.join(pathOUT_freqs, "WI_waterFrequency_"+sname+".tif")
        rsu.array2raster(waterFreq, geotrans, proj, dest, gdal.GDT_Byte, 255)

        # qml file
        outfile_name = os.path.join(pathOUT_freqs, os.path.basename(dest) + '.qml')
        copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

        dest = os.path.join(pathOUT_freqs, "WI_wetFrequency_" + sname + ".tif")
        rsu.array2raster(wetFreq, geotrans, proj, dest, gdal.GDT_Byte, 255)

        #qml file
        outfile_name = os.path.join(pathOUT_freqs, os.path.basename(dest) + '.qml')
        copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)


waterFreqs = np.array(waterFreqs)
waterFreq_all = np.nansum(waterFreqs, axis=0) / np.nansum(validPixels != 0, axis=0)

wetFreqs = np.array(wetFreqs)
wetFreq_all = np.nansum(wetFreqs, axis=0) / np.nansum(validPixels != 0, axis=0)

del wetFreqs, waterFreqs


validPixels = np.array(validPixels)
validObs = np.nansum(validPixels, axis=0)

# DRY occurance and frequency ==========================================================
dryFreq_all = 100. - waterFreq_all - wetFreq_all

# Sum to unity check --------------------
#sum2unity = wetFreq_all + waterFreq_all + dryFreq
#check = np.where((sum2unity>1) & ~np.isnan(sum2unity), 1, 0)
#if (np.sum(check) != 0):
#    raise GeoAlgorithmExecutionException("Water and wetness masks do not sum to unity.")

# CLASSIFICATION =========================================================================

# Water Wetness Presence Index
WWPI = waterFreq_all + (0.75 * wetFreq_all)

# Permanent Water
permWater = np.where(waterFreq_all >= 80, 1, 0)

# Temporary Water
tempWater = np.where((dryFreq_all <= 75) & (wetFreq_all<75) & (waterFreq_all<80)  & (waterFreq_all >= wetFreq_all), 2, 0) # & ((waterFreq > 0.25) & (waterFreq <= 0.8))

# Permanent Wet
permWet = np.where(wetFreq_all >= 75, 3, 0)

# temporary Wet
tempWet = np.where((dryFreq_all <= 75) & (wetFreq_all<75) & (waterFreq_all<80) & (waterFreq_all < wetFreq_all), 4,0) #& ((wetFreq > 0.25) & (wetFreq <= 0.75))

# Dry
noWater = np.where(dryFreq_all > 75, 10, 0)

classification = permWater + tempWater + permWet + tempWet + noWater

# Reclassify dry pixels from 10 to 0
classification = np.where(classification == 10, 0, classification)


# Export maps to file ===============================================================================

# WATER FREQUENCY

file_name = "WI_total_NUMOB"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(validPixels, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "WI_water_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(waterFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

file_name = "WI_dry_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(dryFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

# Wet frequency
#file_name = "wetness_occurance"
#dest = os.path.join(pathOUT_class, file_name + '.tif')
#rsu.array2raster(wetSum, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "WI_wetness_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(wetFreq, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

# CLASSIFICATION

file_name = "WI_classification"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(classification, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "classification_waterWetness.qml"), outfile_name)

# WWPI

file_name = "WI_WWPI"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(WWPI, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "WWPI.qml"), outfile_name)

del waterFreq_all, wetFreq_all

