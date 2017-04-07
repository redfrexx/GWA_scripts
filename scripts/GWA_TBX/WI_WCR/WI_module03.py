##ParameterFile|pathIN|Directory containing water and wetness masks|True|False
##OutputDirectory|pathOUT|Output directory

import os
import sys
import fnmatch
from shutil import copyfile
import numpy as np
import gdal

here = os.path.dirname(scriptDescriptionFile)
pyDir = os.path.join(here, 'data', 'python')
qmlDir = os.path.join(here, "data", "qml")
if pyDir not in sys.path:
    sys.path.append(pyDir)

import RSutils.RSutils as rsu


# test paths
#pathIN="I:/temp/QGIStool/test_WI/watermasks"
#pathOUT = "I:/temp/QGIStool/test_WI"

pathOUT_class = os.path.join(pathOUT,"classification_WI")
if not os.path.exists(pathOUT_class):
    os.mkdir(pathOUT_class)


# WATER occurance and frequency =========================================================

waterMaskFiles = [os.path.join(pathIN, f) for f in fnmatch.filter(os.listdir(pathIN), "*_watermask_sieve.tif")]
print ("Found "+ str(len(waterMaskFiles)) +" water mask files.\n")

jointExtent = rsu.getJointExtent(waterMaskFiles)
geotrans, proj = rsu.raster2array(waterMaskFiles[0], jointExtent)[1:3]

waterMasks = [rsu.raster2array(f, jointExtent)[0] for f in waterMaskFiles]
waterMasks = np.array(waterMasks)
waterMasks = np.where(waterMasks==255, np.nan, waterMasks)

valid = ~np.isnan(waterMasks)
validPixels = np.nansum(~np.isnan(waterMasks), axis=0).astype("float32")

waterSum = np.nansum(waterMasks==1, axis=0).astype("float32")
waterFreq = waterSum / validPixels

del waterMasks

# WETNESS occurance and frequency =========================================================

wetMaskFiles = [os.path.join(pathIN, f) for f in fnmatch.filter(os.listdir(pathIN), "*wetmask_sieve.tif")]

if len(wetMaskFiles) == 0:
    print("no masks found.")
else:
    print ("Found "+ str(len(wetMaskFiles)) +" wet mask files.\n")

wetMasks = [rsu.raster2array(f,jointExtent)[0] for f in wetMaskFiles]
wetMasks = np.array(wetMasks)
wetMasks = np.where(wetMasks==255, np.nan, wetMasks)

wetSum = np.nansum(wetMasks==1, axis=0).astype("float32")
wetFreq = wetSum / validPixels

del wetMasks

# DRY occurance and frequency ==========================================================

drySum = (validPixels - wetSum - waterSum)
dryFreq = drySum / validPixels

# CLASSIFICATION =========================================================================

# Water Wetness Presence Index
WWPI = ((waterSum + (0.75*wetSum)) / validPixels ) * 100.0

# Permanent Water
permWater = np.where(waterFreq >= 0.8, 1, 0)

# Temporary Water
tempWater = np.where((dryFreq <= 0.75) & (wetFreq<0.75) & (waterFreq<0.8)  & (waterFreq >= wetFreq), 2, 0) # & ((waterFreq > 0.25) & (waterFreq <= 0.8))

# Permanent Wet
permWet = np.where(wetFreq >= 0.75, 3, 0)

# temporary Wet
tempWet = np.where((dryFreq <= 0.75) & (wetFreq<0.75) & (waterFreq<0.8) & (waterFreq < wetFreq), 4,0) #& ((wetFreq > 0.25) & (wetFreq <= 0.75))

# Dry
noWater = np.where(dryFreq > 0.75, 10, 0)

classification = permWater + tempWater + permWet + tempWet + noWater

# Reclassify dry pixels from 10 to 0
classification = np.where(classification == 10, 0, classification)


# Export maps to file ===============================================================================

# WATER FREQUENCY

file_name = "total_NUMOB"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(validPixels, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "water_occurance"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(waterSum, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "water_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(waterFreq, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)


# Wet frequency
file_name = "wetness_occurance"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(wetSum, geotrans, proj, dest, gdal.GDT_Byte, 255)


file_name = "wetness_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(wetFreq, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

# CLASSIFICATION

file_name = "classification"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(classification, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "classification_waterWetness.qml"), outfile_name)

# WWPI

file_name = "WWPI"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(WWPI, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "WWPI.qml"), outfile_name)

del waterFreq, wetFreq, wetSum, waterSum


