#  Copyright (c) 2017, GeoVille Information Systems GmbH
#  All rights reserved.
# 
#  Redistribution and use in source and binary forms, with or without
#  modification, is prohibited for all commercial applications without
#  licensing by GeoVille GmbH.
# 
# 
# Date created: 06.05.2017
# Date last modified: 09.05.2017
# 
# 
# __author__ = "Christina Ludwig"
# __version__ = "1.0"

##Water and Wetness Indices=name
##Wetland Inventory=group
##ParameterFile|path_imagery|Directory containing imagery|True|False
##ParameterSelection|sensor|Sensor|Sentinel-2;Landsat 5/7/8|Sentinel-2
##ParameterSelection|AOItype|Type of AOI|Shapefile;User defined extent;Joint extent of all scenes
##ParameterFile|path_AOI|Shapefile containing AOI (Area of Interest)
##ParameterExtent|extentCoords|User defined Extent||True
##OutputDirectory|path_output|Output directory
##*ParameterNumber|maxCloudCover|Maximum Cloud Coverage|0|100|100
##*ParameterString|tileID|Tile ID (Sentinel-2) or Path/Row (Landsat)||False|True|False
##*ParameterString|startDate|Start date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterString|endDate|End date (YYYYMMDD) - if left empty all available scenes will be used||False|True|True
##*ParameterSelection|WCRonly|Spectral indices for ... |Wetland Inventory;Water Cycle Regime|Wetland Inventory

# TODO Add max cloud coverage parameter

debug = False

if debug:
    # TEST PARAMETERS
    path_output = r"I:\2687_GW_A\Toolbox\02_InterimProducts\site_02"
    path_imagery = r"T:\Processing\2687_GW_A\01_RawData\Imagery\workshop\otherwetland\zipped"
    path_AOI = r"T:\Processing\2687_GW_A\01_RawData\Ancillary_Data\workshop\Shapefiles_Lake_Wet\Wetland.shp"
    sensor = 0
    maxCloudCover = 100.
    tileID = ""
    WCRonly = False
    extentCoords = "-1851429.54238,-1831170.14094,1554149.79985,1571383.44653"
    projWkt = 'PROJCS["WGS 84 / Pseudo-Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Mercator_1SP"],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"],AUTHORITY["EPSG","3857"]]'
    startDate = ""
    endDate = ""
    AOItype=0
    here = r"C:\Users\ludwig\.qgis2\processing\scripts\GWA_TBX\WI_WCR"

# IMPORTS ----------------------------------------------------------------------------------------------
import os, sys
from osgeo import gdal, osr
import numpy as np
from matplotlib import pyplot as plt
import logging
import fnmatch
import datetime as dt
import zipfile

# Load additional library
if not debug:
    here = os.path.dirname(scriptDescriptionFile)

pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)

import RSutils.RSutils as rsu

if not debug:
    from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
    import qgis


# FUNCTIONS #########################################################################################

def calculate_indices_for_scene(scene, outputDir, bandNo, extentAOI=None, WCRonly=False):

    # Output directory
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    # Make file title
    outputName = scene.sensorCode + "_" + scene.tileID + "_" + dt.datetime.strftime(scene.date, "%Y%m%d")

    # Read in bands
    bands = []
    for b in bandNo:
        bands.append(scene.getBand(b, masked=True))

    # Band minimum value
    minimum = abs(min(0, np.nanmin(bands)))# + 10
    bands = [band + minimum for band in bands]
    bands = np.array(bands)
    invalidPixels = np.nansum(np.isnan(bands), axis=0)
    bands = np.where(invalidPixels>0, np.nan, bands)

    # Assign arrays to spectral bands
    # BLUE = bands[0].astype("float32")
    GREEN = bands[1].astype("float32")
    RED = bands[2].astype("float32")
    NIR = bands[3].astype("float32")
    SWIR1 = bands[4].astype("float32")
    SWIR2 = bands[5].astype("float32")

    NAN = np.where(np.isnan(bands) | (bands == 0), 1, 0)
    mask = np.where(np.nansum(NAN, axis=0) > 0, 1, 0)
    if (mask.shape[1] < extentAOI.ncol) or (mask.shape[0] < extentAOI.nrow):
        mask = rsu.padArray(mask, scene.geotrans, scene.proj, extentAOI)

    nodata = -9999

    del bands

    # NDVI - Normalized Difference Vegetation Index (NIR-RED/NIR+RED)
    dest = os.path.join(outputDir, "NDVI", outputName + "_NDVI.tif")
    # if not os.path.exists(dest):
    if not os.path.exists(os.path.dirname(dest)):
        os.mkdir(os.path.dirname(dest))
    print("NDVI")
    NDVI = ((NIR - RED) / (NIR + RED)) * 5000. + 5000.
    if (NDVI.shape[1] < extentAOI.ncol) or (NDVI.shape[0] < extentAOI.nrow):
        NDVI = rsu.padArray(NDVI, scene.geotrans, scene.proj, extentAOI)
    # Save to NDVI to file
    NDVI = np.where(mask == 1, nodata, NDVI)
    rsu.array2raster(NDVI, scene.geotrans, scene.proj, dest, gdal.GDT_Int16, nodata)
    del NDVI

    # NDWI (McFeeters, 1996) and Normalized Difference Water Index
    dest = os.path.join(outputDir, "NDWI", outputName + "_NDWI.tif")
    # if not os.path.exists(dest):
    print("NDWI")
    if not os.path.exists(os.path.dirname(dest)):
        os.mkdir(os.path.dirname(dest))
    NDWI = ((GREEN - NIR) / (GREEN + NIR)) * 5000. + 5000.
    if (NDWI.shape[1] < extentAOI.ncol) or (NDWI.shape[0] < extentAOI.nrow):
        NDWI = rsu.padArray(NDWI, scene.geotrans, scene.proj, extentAOI)
    NDWI = np.where(mask == 1, nodata, NDWI)
    rsu.array2raster(NDWI,  scene.geotrans, scene.proj, dest, gdal.GDT_Int16, nodata)
    del NDWI

    # mNDWI (Xu ,2006)
    dest = os.path.join(outputDir, "mNDWI", outputName + "_mNDWI.tif")
    # if not os.path.exists(dest):
    print("mNDWI")
    if not os.path.exists(os.path.dirname(dest)):
        os.mkdir(os.path.dirname(dest))
    MNDWI = ((GREEN - SWIR1) / (GREEN + SWIR1)) * 5000. + 5000.
    if (MNDWI.shape[1] < extentAOI.ncol) or (MNDWI.shape[0] < extentAOI.nrow):
        MNDWI = rsu.padArray(MNDWI, scene.geotrans, scene.proj, extentAOI)
    MNDWI = np.where(mask == 1, nodata, MNDWI)
    rsu.array2raster(MNDWI, scene.geotrans, scene.proj, dest, gdal.GDT_Int16, nodata)
    del MNDWI

    if not WCRonly:
        # NDMI (Gao, 1996) and Normalized Difference Moisture Index
        dest = os.path.join(outputDir, "NDMI", outputName + "_NDMI.tif")
        # if not os.path.exists(dest):
        print("NDMI")
        if not os.path.exists(os.path.dirname(dest)):
            os.mkdir(os.path.dirname(dest))
        NDMI = ((NIR - SWIR1) / (NIR + SWIR1)) * 5000. + 5000.
        if (NDMI.shape[1] < extentAOI.ncol) or (NDMI.shape[0] < extentAOI.nrow):
            NDMI = rsu.padArray(NDMI, scene.geotrans, scene.proj, extentAOI)
        NDMI = np.where(mask == 1, nodata, NDMI)
        rsu.array2raster(NDMI,  scene.geotrans, scene.proj, dest, gdal.GDT_Int16, nodata)
        del NDMI

        # Norm Diff GREEN - SWIR2
        # dest = os.path.join(outputDir, "ND-GREEN-SWIR2", outputName + "_ND-GREEN-SWIR2.tif")
        # if not os.path.exists(dest):
        #     print("ND-GREEN-SWIR2")
        #     if not os.path.exists(os.path.dirname(dest)):
        #         os.mkdir(os.path.dirname(dest))
        #     ND_G_SW1 = ((GREEN - SWIR2) / (GREEN + SWIR2)) * 5000. + 5000.
        #     if (ND_G_SW1.shape[1] < extentAOI.ncol) or (ND_G_SW1.shape[0] < extentAOI.nrow):
        #         ND_G_SW1 = rsu.padArray(ND_G_SW1, geotrans, proj, extentAOI)
        #     rsu.array2raster(ND_G_SW1, extentAOI.getGeotrans(),extentAOI.getProj(), dest, gdal.GDT_Int16, nodata)
        #     del ND_G_SW1

        dest = os.path.join(outputDir, "ND-NIR-SWIR2", outputName + "_ND0406.tif")
        # if not os.path.exists(dest):
        print("ND0406")
        if not os.path.exists(os.path.dirname(dest)):
            os.mkdir(os.path.dirname(dest))
        ND0406 = ((NIR - SWIR2) / (NIR + SWIR2)) * 5000. + 5000.
        if (ND0406.shape[1] < extentAOI.ncol) or (ND0406.shape[0] < extentAOI.nrow):
            ND0406 = rsu.padArray(ND0406, scene.geotrans, scene.proj, extentAOI)
        ND0406 = np.where(mask == 1, nodata, ND0406)
        rsu.array2raster(ND0406,  scene.geotrans, scene.proj, dest, gdal.GDT_Int16, nodata)
        del ND0406

        # Normalized Multi-band Drought Index (NMDI)
        # dest = os.path.join(outputDir, "NMDI", outputName + "_NMDI.tif")
        # if not os.path.exists(dest):
        #     print("NMDI")
        #     if not os.path.exists(os.path.dirname(dest)):
        #         os.mkdir(os.path.dirname(dest))
        #     NMDI = ((NIR - (SWIR1 - SWIR2)) / (NIR + (SWIR1 - SWIR2))) * 10000.
        #     if (NMDI.shape[1] < extentAOI.ncol) or (NMDI.shape[0] < extentAOI.nrow):
        #         NMDI = rsu.padArray(NMDI,  scene.geotrans,  scene.proj, extentAOI)
        #     rsu.array2raster(NMDI, extentAOI.getGeotrans(), extentAOI.getProj(), dest, gdal.GDT_Int16, nodata)
        #     del NMDI



# INPUT PARAMETERS ##################################################################################

if not debug:
    progress.setPercentage(0)

# AOI
if AOItype == 0:
    if not path_AOI.endswith(".shp"):
        raise GeoAlgorithmExecutionException("Invalid input parameter: 'Shapefile containing AOI' is not a shapefile (must have ending .shp). Please specify the path to the shapefile containing the AOI or choose a different type of AOI." )
        sys.exit(1)
    if not os.path.exists(path_AOI) :
        raise GeoAlgorithmExecutionException("Invalid input parameter: 'Shapefile containing AOI' was not found. Please specify the correct path to the shapefile or choose a different type of AOI." )
        sys.exit(1)
elif AOItype == 1:
    if extentCoords is None:
        raise GeoAlgorithmExecutionException("Invalid input parameter: No user defined extent given. Please select an extent or choose a different type of AOI." )
        sys.exit(1)

# Imagery
if not os.path.exists(path_imagery):
    print "Invalid input parameter: 'Directory containing imagery' not found: %s" % path_imagery
    raise GeoAlgorithmExecutionException("Invalid input parameter: 'Directory containing imagery' not found: %s" % path_imagery)

# Output directory
if not os.path.exists(path_output):
    print "Invalid input parameter: 'Output directory' not found: %s" % path_imagery
    raise GeoAlgorithmExecutionException("Invalid input parameter: 'Output directory' not found: %s" % path_imagery)

# Sensor
if sensor == 0:
    sensor = "Sentinel"
elif sensor == 1:
    sensor = "Landsat"

# Subset by path/row (landsat) or tileID (Sentinel)
if tileID == "":
    tileID = None

# Create output directories
path_vrt = os.path.join(path_output, "VRTfiles")
if not os.path.exists(path_vrt):
    os.mkdir(path_vrt)

path_indices= os.path.join(path_output, "indices")
if not os.path.exists(path_indices):
    os.mkdir(path_indices)

path_cloudmasks= os.path.join(path_output, "cloudmasks")
if not os.path.exists(path_cloudmasks):
    os.mkdir(path_cloudmasks)

# Search for scene directories  --------------------------------------------------------------
if sensor == "Sentinel":
    sceneDirs = rsu.search_scene_directories(path_imagery, "S2[AB]*")
    sceneDirs_ID = []

    for i,sceD in enumerate(sceneDirs):
        metadatafiles = []

        try:

            if sceD.endswith(".zip"):
                zpfile = zipfile.ZipFile(sceD)
                metadatafiles = fnmatch.filter(zpfile.namelist(), "*/GRANULE/*/MTD_TL.xml")
                if len(metadatafiles) == 0:
                    metadatafiles = fnmatch.filter(zpfile.namelist(), "*/GRANULE/*/S2*.xml")
                    if len(metadatafiles) > 0:
                        sceneDirs_ID += [[sceD, x[-10:-4]] for x in metadatafiles]
                else:
                    sceneDirs_ID.append([sceD, None])
                del zpfile
            else:
                granuleDir = ""
                for (root, dirnames, filenames) in os.walk(sceD):
                    for d in fnmatch.filter(dirnames, "GRANULE"):
                        granuleDir = os.path.join(root, d)
                for d in os.listdir(granuleDir):
                    for f in fnmatch.filter(os.listdir(os.path.join(granuleDir,d)), "*.xml"):
                        metadatafiles.append(os.path.join(granuleDir, d,f))

                if len(metadatafiles) > 1:
                    sceneDirs_ID += [[sceD, os.path.basename(x)[-10:-4]] for x in metadatafiles]
                else:
                    sceneDirs_ID.append([sceD, None])

        except Exception,e:
            if not debug:
                raise GeoAlgorithmExecutionException("Invalid input data: Invalid file in 'Directory containing imagery': %s " % sceD)
            else:
                print("Invalid input data: Invalid file in 'Directory containing imagery': %s" % sceD)
                sys.exit(1)

    sceneDirs = sceneDirs_ID
else:
    sceneDirs = rsu.search_scene_directories(path_imagery, "L[COTE]*")
    sceneDirs = [[sD,None] for sD in sceneDirs]

# Search scenes --------------------------------------------------------------------------------
if not debug:
    progress.setText("Searching %s scenes ..." % sensor)

scenes = []
for sD in sceneDirs:

    if not debug:
        progress.setText("\n%s - %s" % (sD[0], sD[1]))
    else:
        print("\n%s - %s" % (sD[0], sD[1]))

    try:
        if sensor == "Landsat":
            newScene = rsu.LandsatScene(sD[0], ID=sD[1], tempDir=path_cloudmasks)
        elif sensor == "Sentinel":
            newScene = rsu.SentinelScene(sD[0], ID=sD[1], tempDir=path_cloudmasks)

        if tileID is None or (tileID in newScene.tileID):
            scenes.append(newScene)
        else:
            if not debug:
                progress.setText("Scene not included: Tile ID doesn't match %s." % tileID)

    except IOError as e:
        if not debug:
            progress.setText("WARNING: Scene is skipped. (%s)" % (e))
        else:
            print("WARNING: Scene skipped. (%s)" % (e))
    except Exception, e:
        if not debug:
            progress.setText("ERROR: Scene is skipped. (%s)" % (e))
        else:
            print "ERROR: Scene is skipped. (%s)" % e

if len(scenes) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No %s scenes found that match the selection criteria. Adjust the selection parameters." % sensor)
    else:
        print "No %s scenes found that match the selection criteria. Adjust the selection parameters." % sensor
        sys.exit(1)

# Check start and end dates ---------------------------------------------------------------------------------------------------------------------------------
if startDate == "":
    startDate = min([sce.date for sce in scenes])
else:
    try:
        startDate = dt.datetime.strptime(startDate, "%Y%m%d")
    except: 
        raise GeoAlgorithmExecutionException("Invalid input parameter: Format of 'Start date' is not valid.")

if endDate == "":
    endDate = max([sce.date for sce in scenes])
else:
    try:
        endDate = dt.datetime.strptime(endDate, "%Y%m%d")
    except: 
        raise GeoAlgorithmExecutionException("Invalid input parameter: Format of 'End date' is not valid.")

if endDate < startDate:
    raise GeoAlgorithmExecutionException("Invalid input parameters: 'Start date'  must be earlier than 'End date'.")

scenes = rsu.filter_scenes_by_date(scenes, startDate, endDate)

if len(scenes) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("No %s scenes found within the given time period. Adjust the start and end date." % sensor)
    else:
        print "No %s scenes found within the given time period. Adjust the start and end date." % sensor
        sys.exit(1)

# Get extent of AOI -----------------------------------------------------------------------------------------------------------------------------------------------------

if AOItype == 0:
    # use swir2 bands as reference for geotransformation. Must be provided
    exRaster = scenes[0].files[-1]
    extentAOI = rsu.convertShapeExtentToRasterExtent(path_AOI, exRaster)
    
    if extentAOI is None:
        if not debug:
            raise GeoAlgorithmExecutionException("Invalid input parameter: Shapefile contains more than one polygon or AOI does not intersect with input scenes.")
        sys.exit(1)
        
elif AOItype == 1:
    extentCoords = [float(coord) for coord in extentCoords.split(",")]
    if not debug:
        projWkt = qgis.utils.iface.mapCanvas().mapRenderer().destinationCrs().toWkt()

    projExtent = osr.SpatialReference()
    projExtent.ImportFromWkt(projWkt)
    extentAOI = rsu.extent(ulX=extentCoords[0], ulY=extentCoords[3], lrX=extentCoords[1], lrY=extentCoords[2], proj=projExtent)
    extentAOI.convertToRasterExtent(scenes[0])
    
    if extentAOI is None:
        if not debug:
            progress.setText("Invalid input parameter: User defined extent is invalid.")
        sys.exit(1)

    # Create folder for extents
    path_output_extents = os.path.join(path_output, "userDefinedExtents")
    if not os.path.exists(path_output_extents):
        os.mkdir(path_output_extents)
    extent_name = "userDefinedExtent_" + dt.datetime.today().strftime("%Y%m%d-h%Hm%M")
    extentAOI.save_as_shp(path_output_extents, extent_name)

    if not debug:
        layer = qgis.utils.iface.addVectorLayer(os.path.join(path_output_extents, extent_name+".shp"), extent_name, "ogr")

else:
    # use swir2 bands as reference for geotransformation. Must be provided.
    exBands = [sce.files[-1] for sce in scenes]
    extentAOI = rsu.getJointExtent(exBands)

    if extentAOI is None:
        if not debug:
            progress.setText("Invalid input data: Extent of AOI is invalid. Check if all scenes overlap each other.")
        sys.exit(1)


# Update extent of scenes -----------------------------------------------------------------------------------------
scenesWithinExtent = []
for sce in scenes:

    # Check whether scene is within AOI
    try:
        sce.updateExtent(extentAOI)
    except RuntimeWarning, e:
        print("%s not within AOI." % sce.ID)
        if not debug:
            progress.setText("%s not within AOI." % sce.ID)
        continue

    # Check whether scene has enough valid pixels within AOI
    #sce.calcCloudCoverage()
    if sce.cloudy > maxCloudCover:
        print "Cloud coverage too high. Scene is skipped."
        if not debug:
            progress.setText("Cloud coverage too high. Scene is skipped.")
            continue

    scenesWithinExtent.append(sce)
    sce.createVRT(path_vrt)

scenes = scenesWithinExtent
if len(scenes) == 0:
    if not debug:
        raise GeoAlgorithmExecutionException("Invalid input data: No %s scenes found that match the selection criteria." % sensor)
        sys.exit(1)
    else:
        print "Invalid input data: No %s scenes found that match the selection criteria." % sensor

if not debug:
    progress.setPercentage(10)
    progress.setText("\n%s %s scenes selected." % (len(scenes), sensor))

#  CALCULATION OF SPECTRAL INDICES -------------------------------------------------------------------------------------

if not debug:
    progress.setText("\nCalculating spectral indices ...")

for i, sce in enumerate(scenes):
    
    # Update progress bar
    if not debug:
        progress.setText("%s" % os.path.basename(sce.ID))

    if sensor == "Sentinel":
        calculate_indices_for_scene(sce, path_indices, [1, 2, 3, 7, 9, 10], extentAOI, WCRonly)
    elif sensor == "Landsat":
        calculate_indices_for_scene(sce, path_indices, [1, 2, 3, 4, 5, 6], extentAOI, WCRonly)

    # Update progress bar
    if not debug:
        state = ((i + 1.) / len(scenes)) * 90
        progress.setPercentage(state)

# Merge scenes from the same date to VRT file --------------------------------------------------------------------------

dates = [sce.date for sce in scenes]
uniqueDates = set(dates)
numDates = float(len(uniqueDates))
indexFolders = os.listdir(path_indices)

for date in uniqueDates:

    scenesAtdate = [sce for sce in scenes if (sce.date == date)]

    for index in indexFolders:
        inFiles = [os.path.join(path_indices, index, f) for f in fnmatch.filter(os.listdir(os.path.join(path_indices, index)), "*" + dt.datetime.strftime(date, "%Y%m%d") + "*.tif")]

        if len(inFiles) == 0:
            print("No indices found for %s" % dt.datetime.strftime(date, "%Y%m%d"))
            continue

        sensorCode = os.path.basename(inFiles[0])[:3]

        if sensor == "Landsat":
            mergedOutfile = os.path.join(path_indices, index, sensorCode + "_" + dt.datetime.strftime(date, "%Y%m%d") + "_" + index + ".vrt")
        else:
            mergedOutfile = os.path.join(path_indices, index, sensorCode + "_" + dt.datetime.strftime(date, "%Y%m%d") + "_" + index + ".vrt")

        try:
            rsu.createVRT(inFiles, mergedOutfile, separate=False)
        except Exception, e:
            if not debug:
                progress.setText("ERROR: %s" % e)
            else:
                print e

# Update progress bar
if not debug:
    progress.setPercentage(100)




