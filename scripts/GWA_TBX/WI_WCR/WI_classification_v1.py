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
__author__ = "Christina Ludwig"
__version__ = "1.0"

##Wetland Inventory Classification = name
##Wetland Inventory=group
##ParameterFile|path_input|Directory containing water and wetness masks|True|False
##OutputDirectory|path_output|Output directory
##ParameterBoolean|exportSeasonalFrequencies|Export seasonal water and wetness frequencies
##*ParameterString|startDate|Start date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterString|endDate|End date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterBoolean|spring|Include spring scenes|True
##*ParameterBoolean|summer|Include summer scenes|True
##*ParameterBoolean|fall|Include fall scenes|True
##*ParameterBoolean|winter|Include winter scenes|True

debug = False

# test paths
if debug:
    path_input = r""
    path_output = r""
    exportSeasonalFrequencies = True
    startDate = ""
    endDate = ""
    spring = True
    summer = True
    fall = True
    winter = True
    here = r""

import os, sys
import numpy as np
import fnmatch
import gdal
from shutil import copyfile
import datetime as dt
import time
from processing.tools import dataobjects

if not debug:
    from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
    here = os.path.dirname(scriptDescriptionFile)

# Load additional library
pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)
	
import RSutils.RSutils as rsu
qmlDir = os.path.join(here, 'data', 'qml')

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

if not os.path.exists(path_input):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Directory containing water and wetness masks' does not exist.")
    print("Invalid input parameters: 'Directory containing water and wetness masks' does not exist.")

if not os.path.exists(path_output):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Output directory' does not exist.")
    print("Invalid input parameters: 'Output directory' does not exist.")


pathOUT_class = os.path.join(path_output, "classification_WI")
if not os.path.exists(pathOUT_class):
    os.mkdir(pathOUT_class)

# Seasonal frequencies
if exportSeasonalFrequencies:
    pathOUT_freqs = os.path.join(pathOUT_class,"seasonal_frequencies")
    if not os.path.exists(pathOUT_freqs):
        os.mkdir(pathOUT_freqs)

# Check start and end dates ---------------------------------------------------------------------------------------------------------------------------------
if startDate == "":
    startDate = None
else:
    try:
        startDate = dt.datetime.strptime(startDate, "%Y%m%d")
    except:
        if not debug:
            raise GeoAlgorithmExecutionException("Invalid input parameter: Format of 'Start date' is not valid.")

if endDate == "":
    endDate = None
else:
    try:
        endDate = dt.datetime.strptime(endDate, "%Y%m%d")
    except:
        if not debug:
            raise GeoAlgorithmExecutionException("Invalid input parameter: Format of 'End date' is not valid.")

if endDate is not None and startDate is not None and endDate < startDate:
    raise GeoAlgorithmExecutionException("Invalid input parameters: 'Start date'  must be earlier than 'End date'.")



# WATER occurance and frequency =========================================================

# Search water masks
waterMaskFiles = [os.path.join(path_input, f) for f in fnmatch.filter(os.listdir(path_input), "*_watermask.tif")]
wetMaskFiles = [os.path.join(path_input, f) for f in fnmatch.filter(os.listdir(path_input), "*wetmask.tif")]

# Check whether masks exist
if len(waterMaskFiles) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No water masks found")
    else:
        print("No water masks found.")
else:
    print ("Found " + str(len(waterMaskFiles)) + " water mask files.\n")

if len(wetMaskFiles) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No wetness masks found")
    else:
        print("No wetness masks found.")
else:
    print ("Found " + str(len(wetMaskFiles)) + " wet mask files.\n")

jointExtent = rsu.getJointExtent(waterMaskFiles)
geotrans, proj = rsu.raster2array(waterMaskFiles[0], jointExtent)[1:3]

waterFreqs = []
validPixels = []
validPixels_wet = []
wetFreqs = []


seasons = [["winter",[12,1,2],winter],["spring",[3,4,5], spring],["summer",[6,7,8], summer],["fall",[9,10,11], fall]]
for season in seasons:

    if not season[2]:
        continue

    print(season[0])

    waterMasks_season = []
    for f in waterMaskFiles:
        f_name = os.path.basename(f)
        dateidx = f_name.find("_d")
        date = dt.datetime.strptime(f_name[dateidx+2:dateidx+10], "%Y%m%d")
        if (date.month in season[1]) and (startDate is None or (startDate <= date)) and (endDate is None or (endDate >= date)):
            waterMasks_season.append(f)

    if len(waterMasks_season) == 0:
        validPixels.append(np.zeros((jointExtent.nrow, jointExtent.ncol)))
        print("No observations for %s " % season[0])
        continue

    waterSum, validObs_water = calculateFrequency(waterMasks_season, jointExtent)
    waterFreq = (waterSum / validObs_water) * 100.
    waterFreqs.append(waterFreq)
    validPixels.append(validObs_water)

    # wetness frequency
    wetMasks_season = []
    for f in wetMaskFiles:
        f_name = os.path.basename(f)
        dateidx = f_name.find("_d")
        date = dt.datetime.strptime(f_name[dateidx + 2:dateidx + 10], "%Y%m%d")
        if date.month in season[1] and (startDate is None or (startDate <= date)) and (endDate is None or (endDate >= date)):
            wetMasks_season.append(f)

    if len(wetMasks_season) == 0:
        continue

    wetSum, validObs_wet = calculateFrequency(wetMasks_season, jointExtent)
    wetFreq = (wetSum / validObs_water) * 100.
    wetFreqs.append(wetFreq)
    validPixels_wet.append(validObs_wet)

    # Export to file
    if exportSeasonalFrequencies:
        dest = os.path.join(pathOUT_freqs, "WI_waterFrequency_"+season[0]+".tif")
        rsu.array2raster(waterFreq, geotrans, proj, dest, gdal.GDT_Byte, 255)

        # qml file
        outfile_name = os.path.join(pathOUT_freqs, os.path.basename(dest)[:-4] + '.qml')
        copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

        dest = os.path.join(pathOUT_freqs, "WI_wetFrequency_" + season[0] + ".tif")
        rsu.array2raster(wetFreq, geotrans, proj, dest, gdal.GDT_Byte, 255)

        #qml file
        outfile_name = os.path.join(pathOUT_freqs, os.path.basename(dest)[:-4] + '.qml')
        copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)


validPixels = np.array(validPixels)
validObs = np.nansum(validPixels, axis=0)

validPixels_wet = np.array(validPixels_wet)
validObs_wet = np.nansum(validPixels_wet, axis=0)

waterFreqs = np.array(waterFreqs)
waterFreq_all = np.nansum(waterFreqs, axis=0) / np.nansum(validPixels != 0, axis=0)

wetFreqs = np.array(wetFreqs)
wetFreq_all = np.nansum(wetFreqs, axis=0) / np.nansum(validPixels != 0, axis=0)

wetFreq_all = np.where((validObs == 0), np.nan, wetFreq_all)
waterFreq_all = np.where(validObs == 0, np.nan, waterFreq_all)

del wetFreqs, waterFreqs, validPixels

if np.nansum(validObs) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No water masks found")
    else:
        print("No water masks found.")
        sys.exit(1)

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
rsu.array2raster(validObs, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "WI_water_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(waterFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)
dataobjects.load(dest, os.path.basename(dest))

file_name = "WI_dry_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(dryFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

# Wet frequency
file_name = "WI_wetness_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(wetFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)
dataobjects.load(dest, os.path.basename(dest))

# CLASSIFICATION
file_name = "WI_classification"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(classification, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "classification_waterWetness.qml"), outfile_name)
dataobjects.load(dest, os.path.basename(dest))

# WWPI
file_name = "WI_WWPI"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(WWPI, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "WWPI.qml"), outfile_name)
dataobjects.load(dest, os.path.basename(dest))

del waterFreq_all, wetFreq_all

progress.setConsoleInfo('Classification successful!')
