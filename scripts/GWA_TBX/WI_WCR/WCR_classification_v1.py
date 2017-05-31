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

##Water Cycle Regime Classfication = name
##Water Cycle Regime=group

##ParameterFile|path_input|Directory containing water masks|True|False
##OutputDirectory|path_output|Output directory
##ParameterBoolean|exportSeasonalFrequencies|Export seasonal water frequencies
##*ParameterString|startDate|Start date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterString|endDate|End date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterBoolean|spring|Include spring scenes|True
##*ParameterBoolean|summer|Include summer scenes|True
##*ParameterBoolean|fall|Include fall scenes|True
##*ParameterBoolean|winter|Include winter scenes|True

debug = False

# test paths
if debug:
    path_input = r"T:\Processing\2687_GW_A\03_Products\GWA-TOOLBOX\02_InterimProducts\WI\SE_wat45_wet50_win1800_mmu3"
    path_output = r"T:\Processing\2687_GW_A\03_Products\GWA-TOOLBOX\02_InterimProducts\WCR"
    exportSeasonalFrequencies = True
    startDate = ""
    endDate = ""
    spring = True
    summer = True
    fall = True
    winter = True
    here = "."

import os, sys
import numpy as np
import fnmatch
from matplotlib import pyplot as plt
import gdal
from shutil import copyfile
from os.path import expanduser
import datetime as dt
import time

if not debug:
    from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
    here = os.path.dirname(scriptDescriptionFile)

# Load additional library
pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)

qmlDir = os.path.join(here, 'data', 'qml')

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

if not os.path.exists(path_input):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Directory containing water masks' does not exist.")
    print("Invalid input parameters: 'Directory containing water masks' does not exist.")

if not os.path.exists(path_output):
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input parameters: 'Output directory' does not exist.")
    print("Invalid input parameters: 'Output directory' does not exist.")


pathOUT_class = os.path.join(path_output, "classification_WCR")
if not os.path.exists(pathOUT_class):
    os.mkdir(pathOUT_class)

# Seasonal frequencies
if exportSeasonalFrequencies:
    pathOUT_freqs = os.path.join(pathOUT_class,"seasonal_frequencies")
    if not os.path.exists(pathOUT_freqs):
        os.mkdir(pathOUT_freqs)

if not debug:
    progress.setText(pathOUT_class)

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

if endDate < startDate:
    raise GeoAlgorithmExecutionException("Invalid input parameters: 'Start date'  must be earlier than 'End date'.")


# WATER occurance and frequency =========================================================

# Search water masks
waterMaskFiles = [os.path.join(path_input, f) for f in fnmatch.filter(os.listdir(path_input), "*_watermask.tif")]

# Check whether masks exist
if len(waterMaskFiles) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No water masks found")
    else:
        print("No water masks found.")
else:
    print ("Found " + str(len(waterMaskFiles)) + " water mask files.\n")
    if not debug:
        progress.setText("Found " + str(len(waterMaskFiles)) + " water mask files.\n")

jointExtent = rsu.getJointExtent(waterMaskFiles)
geotrans, proj = rsu.raster2array(waterMaskFiles[0], jointExtent)[1:3]

waterFreqs = []
validPixels = []

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
        continue

    waterSum, validObs = calculateFrequency(waterMasks_season, jointExtent)
    waterFreq = (waterSum / validObs) * 100.
    waterFreqs.append(waterFreq)
    validPixels.append(validObs)

    # Export to file
    if exportSeasonalFrequencies:
        dest = os.path.join(pathOUT_freqs, "WCR_waterFrequency_"+season[0]+".tif")
        rsu.array2raster(waterFreq, geotrans, proj, dest, gdal.GDT_Byte, 255)

        # qml file
        outfile_name = os.path.join(pathOUT_freqs, os.path.basename(dest) + '.qml')
        copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)


waterFreqs = np.array(waterFreqs)
waterFreq_all = np.nansum(waterFreqs, axis=0) / np.nansum(validPixels != 0, axis=0)

del waterFreqs

validPixels = np.array(validPixels)
validObs = np.nansum(validPixels, axis=0)

waterFreq_all = np.where(validObs == 0, np.nan, waterFreq_all)

if np.nansum(validObs) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No water masks found")
    else:
        print("No water masks found.")
        sys.exit(1)

# CLASSIFICATION

minExtent = np.where(waterFreq_all > 80, 1, 0)
maxExtent = np.where(waterFreq_all > 10, 1, 0)

# Permanent Water
classification = np.array(waterFreq_all)
classification[:] = 0
classification = np.where(waterFreq_all >= 80, 1, classification)

# Temporary Water
classification = np.where((waterFreq_all < 80) & (waterFreq_all >= 25), 2, classification)

# Export maps to file ===============================================================================

# WATER FREQUENCY

file_name = "WCR_total_NUMOB"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(validObs, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "WCR_water_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(waterFreq_all, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)


# CLASSIFICATION

file_name = "WCR_minimum_extent"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(minExtent, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "minimum_extent.qml"), outfile_name)

file_name = "WCR_maximum_extent"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(maxExtent, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "maximum_extent.qml"), outfile_name)

file_name = "WCR_classification"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(classification, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "classification_water.qml"), outfile_name)


