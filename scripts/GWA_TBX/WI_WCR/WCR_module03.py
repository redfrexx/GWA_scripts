##ParameterFile|pathIN|Directory containing water masks|True|False
##OutputDirectory|pathOUT|Output directory

import os
import sys
import numpy as np
import fnmatch
import gdal
from shutil import copyfile

here = os.path.dirname(scriptDescriptionFile)
pyDir = os.path.join(here, 'data', 'python')
qmlDir = os.path.join(here, "data", "qml")
if pyDir not in sys.path:
    sys.path.append(pyDir)

import RSutils.RSutils as rsu

# test paths
#pathIN="I:/temp/QGIStool/site77/watermasks"
#pathOUT = "I:/temp/QGIStool/site77"

pathOUT_class = os.path.join(pathOUT,"classification_WCR")
if not os.path.exists(pathOUT_class):
    os.mkdir(pathOUT_class)


# WETNESS occurance and frequency =========================================================

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

# CLASSIFICATION

minExtent = np.where(waterFreq > 0.8, 1, 0)
maxExtent = np.where(waterFreq > 0.1, 1, 0)

# Permanent Water
classification = np.array(waterFreq)
classification[:] = 0
classification = np.where(waterFreq >= 0.8, 1, classification)

# Temporary Water
classification = np.where((waterFreq < 0.8) & (waterFreq >= 0.25), 2, classification)


# Export maps to file

# Water frequency
# ---------------

# Save to file
file_name = "water_occurance"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(waterSum, geotrans, proj, dest, gdal
                 .GDT_Byte, 255)

file_name = "total_NUMOB"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(validPixels, geotrans, proj, dest, gdal.GDT_Byte, 255)

file_name = "water_frequency"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(waterFreq, geotrans, proj, dest, gdal.GDT_Float32, -9999)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "water_wet_frequency.qml"), outfile_name)

file_name = "minimum_extent"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(minExtent, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "minimum_extent.qml"), outfile_name)

file_name = "maximum_extent"
dest = os.path.join(pathOUT_class, file_name + '.tif')
rsu.array2raster(maxExtent, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "maximum_extent.qml"), outfile_name)

file_name = "classification"
dest = os.path.join(pathOUT_class,  file_name + '.tif')
rsu.array2raster(classification, geotrans, proj, dest, gdal.GDT_Byte, 255)

# qml file
outfile_name = os.path.join(pathOUT_class, file_name+'.qml')
copyfile(os.path.join(qmlDir, "classification_water.qml"), outfile_name)

del waterFreq, waterSum


