#  Copyright (c) 2016, GeoVille Information Systems GmbH
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.
#
# """
# Merge Sentinel Scenes from same date and calculate spectral indices
#
#
# Arguments:
#   - path_imagery:           path to directory with Sentinel imagery
#   - path_AOI:                 path to AOI shapefile
#   - path_output:             path to output directory
#   - siteID:                      siteID if there are more than one feature in the shapefile (optional)
#   - extentBuffer:             buffer around AOI in pixels
#
#
# Output:
#   -  merged scenes and spectral indices as GeoTiff files in output directory (path_output)
#
#
# Dependencies:
#
# Date created: 05.01.2017
# Date last modified: 09.01.2017
#
# """
#
# __author__ = "Christina Ludwig"
# __version__ = "1.0"

##ParameterFile|path_imagery|Directory containing imagery|True|False
##ParameterFile|path_AOI|Area of Interest (AOI) as shapefile |False|False|shp
##OutputDirectory|path_output|Output directory
#ParameterNumber|extentBuffer|Buffer around AOI|0|None|0
##ParameterSelection|sensor|Sensor|Landsat;Sentinel

from osgeo import gdal, ogr
import numpy as np
from matplotlib import pyplot as plt
import logging
import fnmatch
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
import os
import sys

here = os.path.dirname(scriptDescriptionFile)
pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)

import RSutils.RSutils as rsu


# TEST PARAMETER

#path_output = "I:/temp/QGIStool/site77"
#path_imagery = "X:/01_RawData/Imagery/site77/Sentinel"
#path_AOI = "I:/2687_GW_A/02_Interim_Products/Toolbox/sampleData/site77.shp"
#sensor = "Sentinel"
extentBuffer=0

# FUNCTIONS --------------------------------------------------------------------------------------------------
def calculateNDIs(bands, outputDir, outputName, geotrans, proj, rescale=False):
    """Calculate normalized difference indices from input bands"""

    bandnames = ["02","03","04","05","06","07","08","8A","11","12"]
    selectedIndices = ["ND0311", "ND038A"]

    nodata = -9999
    for num1, b1 in enumerate(bands,0):
        for num2, b2 in enumerate(bands,0):
            indexname = "ND" + bandnames[num1] + bandnames[num2]
            if (num1 < num2) and (indexname in selectedIndices):
                dest = os.path.join(outputDir, "ND" + bandnames[num1] + bandnames[num2], outputName + "_ND" + bandnames[num1] + bandnames[num2] + ".tif")
                if not os.path.exists(dest):
                    print "Normalized Difference %s / %s" % (bandnames[num1],bandnames[num2])
                    if not os.path.exists(os.path.dirname(dest)):
                        os.mkdir(os.path.dirname(dest))
                    index = (b1 - b2) / (b1 + b2) * 5000. + 5000.
                    if rescale:
                        index = rescale(index)
                    #index = ndimage.median_filter(index, size=(3,3), cval=np.nan)
                    rsu.array2raster(index, geotrans, proj, dest, gdal.GDT_Int16, nodata)

def allIndicesForSeasonComposite(path_composite, pathOUT, bandNo=[1, 2, 3, 4, 5, 6], rescale=False):

    nodata = -32768

    # Output directory
    #pathOUT = os.path.join(os.path.dirname(path_composite), "indices")
    if not os.path.exists(pathOUT):
        os.mkdir(pathOUT)

    # Number of bands
    comp = gdal.Open(path_composite)
    nbands = comp.RasterCount
    geotrans = comp.GetGeoTransform()
    proj = comp.GetProjection()
    del comp

    # Get file title
    try:
        idx_start = os.path.basename(path_composite).index("_") + 1
        title = os.path.basename(path_composite)[idx_start:-4]
    except:
        title = os.path.basename(path_composite)[:-4]

    # Read in bands
    bands = []
    for b in range(1, nbands+1):
        bands.append(rsu.raster2array(path_composite, bandNo=b)[0])

    # Band minimum value
    minimum = abs(min(0, np.nanmin(bands))) + 10
    bands = [band + minimum for band in bands]

    # Calculate normalized difference indices
    calculateNDIs(bands, pathOUT, title, geotrans, proj, rescale=rescale)

    # Calculate simple ratios
    #calculateSRIs(bands, pathOUT, title, geotrans, proj, rescale=rescale)

    # Bands for other indices
    #bands2 = [bands[i-1] for i in bandNo]
    #calculateIndices(bands2, pathOUT, title, geotrans, proj, rescale=rescale)











attrName = None
siteID = None

if sensor == 1:
    sensor = "Sentinel"
elif sensor == 0:
    sensor = "Landsat"

progress.setText("Checking input parameters ...")
progress.setPercentage(0)

# Change directory
#os.chdir(rootDir)


# AOI
path_AOI = os.path.abspath(path_AOI)
if not os.path.exists(path_AOI):
    print "AOI not found: %s" % path_AOI
    logging.critical("AOI not found: %s" % path_AOI)
    raise GeoAlgorithmExecutionException("AOI not found: %s" % path_AOI)

# Landsat imagery
if not os.path.exists(path_imagery):
    print "Imagery folder not found: %s" % path_imagery
    logging.critical("Input folder not found: %s" % path_imagery)
    raise GeoAlgorithmExecutionException("Input folder not found: %s" % path_imagery)

# Create output directories

path_vrt = os.path.join(path_output, "VRTfiles")
if not os.path.exists(path_vrt):
    os.mkdir(path_vrt)

path_composites = os.path.join(path_output, "composites")
if not os.path.exists(path_composites):
    os.mkdir(path_composites)

path_indices= os.path.join(path_output, "indices")
if not os.path.exists(path_indices):
    os.mkdir(path_indices)

# Get paths to scene directories if scene ID list is provided -----------------------------------------

sceneDirs = rsu.searchSceneDirs_sentinel(path_imagery)
print sceneDirs

# Get extent of AOI -----------------------------------------------------------------------------------
if os.path.exists(path_AOI):
    print("right")
    # use band 5 with 20 m resolution to get shapfile extent in raster coordinates
    if sensor == "Landsat":
        exRaster = os.path.join(sceneDirs[0], fnmatch.filter(os.listdir(sceneDirs[0]), "*band1*.tif")[0])
    elif sensor == "Sentinel":
        exRaster = os.path.join(sceneDirs[0], fnmatch.filter(os.listdir(sceneDirs[0]), "*B05.jp2")[0])
    extentAOI = rsu.convertShapeExtentToRasterExtent(path_AOI, exRaster, attributeName=attrName, attributeValue=siteID, buffer=extentBuffer)

else:
    exBands = []
    for d in sceneDirs:
        # use band 5 with 20 m resolution to get shapfile extent in raster coordinates
        if sensor == "Landsat":
            exBands.append(os.path.join(d, fnmatch.filter(os.listdir(d), "*band1*.tif")[0]))
        elif sensor == "Sentinel":
            exBands.append(os.path.join(d, fnmatch.filter(os.listdir(d), "*B05.jp2")[0]))
    extentAOI = rsu.getJointExtent(exBands, buffer=extentBuffer)

# Get input scenes -------------------------------------------------------------------------------------

progress.setPercentage(10)
progress.setText("\nFound %s scene directories:" % len(sceneDirs))

scenes = []

for sD in sceneDirs:
    progress.setText("Checking %s" % sD)
    try:
        if sensor == "Landsat":
            newScene = rsu.LandsatScene(sD, path_output, extentAOI)
        elif sensor == "Sentinel":
            newScene = rsu.SentinelScene(sD, path_output, extentAOI)
        #newScene.calcCloudCoverage()
        #if newScene.cloudCoverage > 70.0:
         #   print "not enough valid pixels."
         #   print newScene.cloudCoverage
          #  continue
        newScene.createVRT(path_vrt)
        scenes.append(newScene)
        del newScene
    except:
        logging.warning("Scene %s couldn't be read." % os.path.basename(sD)[0:16])
        print "Scene %s couldn't be read." % os.path.basename(sD)[0:16]

if len(scenes) == 0:
    raise GeoAlgorithmExecutionException("No scenes for processing available due to insufficient coverage of AOI!")


# Merge scenes from the same date ----------------------

print "\nMerging scenes ..."
progress.setText("\nMerging scenes ...")

# Find unique dates
dates = [sce.date for sce in scenes]
uniqueDates = set(dates)
numDates = float(len(uniqueDates))

progress.setPercentage(20)

refScenes = []
for i,d in enumerate(iter(uniqueDates)):

    # Update progress bar
    state = 20.+((i+1.)/numDates)*50.
    progress.setPercentage(state)

    # Merge scenes
    dest = os.path.join(path_composites, d.strftime("%Y%m%d") + ".tif")
    if not os.path.exists(dest):
        print "Merging scenes for %s" % d.strftime("%Y / %m / %d")
        progress.setText("Merging scenes for %s" % d.strftime("%Y / %m / %d"))
        scenesAtdate = [sce for sce in scenes if (sce.date == d)]
        if len(scenesAtdate)==0:
            continue
        sceneComposite = rsu.ImageComposite(scenesAtdate, extentAOI, stepSize=1000)
        if sensor == "Sentinel":
            sceneComposite.calc_index(scenesAtdate, method="blue", bands=[1,2,3,8,9,10])
            sceneComposite.export_composite(path_composites, d.strftime("%Y%m%d"), transformed=False)
        else:
            sceneComposite.calc_index(scenesAtdate, method="blue")

# Calculate spectral indices ----------------------

progress.setText("\nCalculating spectral indices ...")

composites = [os.path.join(path_composites, f) for f in fnmatch.filter(os.listdir(path_composites), "????????.tif")]
numComp = float(len(composites))

for i, path_comp in enumerate(composites):

    # Update progress bar
    state = 70.+((i+1.)/numComp)*30.
    progress.setPercentage(state)

    compositeScene = rsu.ImageComposite
    progress.setText("Calculating indices for %s" % os.path.basename(path_comp))
    if sensor == "Sentinel":
        allIndicesForSeasonComposite(path_comp, path_indices, bandNo=[1,2,3,8,9,10])
    elif sensor == "Landsat":
        allIndicesForSeasonComposite(path_comp, path_indices, bandNo=[1,2,3,4,5,6])



