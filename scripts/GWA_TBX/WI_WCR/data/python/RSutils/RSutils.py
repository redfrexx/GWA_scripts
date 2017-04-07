# Copyright (c) 2016, GeoVille Information Systems GmbH
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.

"""
Support functions and classes for raster handling

Date created: 09/06/2016
Date last modified: 09/06/2016

"""

__author__ = "Christina Ludwig"
__version__ = "v1"


from osgeo import gdal, ogr
from gdalconst import GA_ReadOnly
import os, sys
import bottleneck as bn
import numpy as np
from osgeo import osr
from scipy import ndimage
from matplotlib import pyplot as plt
import logging
from datetime import datetime as dt
import math
import fnmatch
import operator
import datetime as dati
import subprocess

# CLASSES ===============================================================

def compareProjections(scenes):

    projs = [sce.proj for sce in scenes]
    #todo implement function

    # check if elements are all the same
    # if not choose projection of first scene as target projection
    # change projection of other scenes to target projection by transforming upperleft cooridnates of geotransform attribute of each scenes with reprojectPoint

class ImageComposite(object):
    """Object to create image composite."""

    def __init__(self, scenes, AOI=None, stepSize=None):

        if len(scenes) == 0:
            raise AttributeError ("No input scenes")
        elif len(scenes) > 255:
            raise AttributeError ("Too many input scenes for image composite.")
        else:
            self.scenes = scenes

        # Number of bands
        self.nbands = len(self.scenes[0].files)

        # Get extent
        if AOI == None:
            sceneSel = [sce.files[3] for sce in self.scenes]
            self.extent = getJointExtent(sceneSel)
            if self.extent is None:
                #logging.critical("Extent could not be calculated, because input files are not all overlapping.")
                raise ValueError ("Input files are not overlapping.")
                # todo how to return none instead of class
        else:
            # todo: Check if AOI intersects with other scenes
            self.extent = AOI

        # Get step size for index computation
        if stepSize == None or stepSize >= self.extent.ncol:
            self.stepSize = self.extent.ncol
        else:
            self.stepSize = stepSize

        # Get metadata
        metadata = gdal.Open(self.scenes[0].files[3], GA_ReadOnly)
        self.geotrans = list(metadata.GetGeoTransform())
        self.geotrans[0] = self.extent.ulX
        self.geotrans[3] = self.extent.ulY
        self.pixSize = self.geotrans[1]

        # todo:Check projections of input scenes
        self.proj = metadata.GetProjection()
        del metadata

        # initilize output arrays
        self.validObs_index = np.empty((0, self.extent.ncol), dtype=np.int8)
        self.masks = np.empty((len(self.scenes), 0, self.extent.ncol), dtype=np.int8)

        # Dictionary of layer index and sceneID
        self.sceneIDs = {}
        for j,f in enumerate(self.scenes):
            self.sceneIDs[j] = f.ID

        self.sceneLookup = {}
        for j,f in enumerate(self.scenes):
            self.sceneLookup[j] = f
            f.indexNr = j

        self.subextents = getSubExtents(self.extent, self.stepSize)

        # Transformation
        self.referenceScene = ""
        self.targetScenes = []
        self.transformedScenes = []
        self.indexScenes = []

    def calc_validPixels(self):
        """Calculate number of valid observations per pixel."""
        valid = []
        for sce in self.scenes:
            jointExtent = getJointExtent(sce.files[3], AOIextent=self.extent)
            mask = sce.getMask(extent=jointExtent)
            geotrans = [jointExtent.ulX, jointExtent.pixSize, 0.0, jointExtent.ulY, 0.0, -jointExtent.pixSize]
            if mask.shape[1] < self.extent.ncol or mask.shape[0] < self.extent.nrow:
                mask = padArray(mask, geotrans, sce.proj, self.extent)
            valid.append(np.where(mask == 0, 1, 0))
        valid = np.array(valid)
        self.validObs = bn.nansum(valid, axis=0)
        del valid

    def calc_index(self, inputScenes, method="HOT", bands=[1,2,3,4,5,6], transformed=True):
        """Calculate image composite index based on specified method."""
        print "\nCalculating composite index using method: %s" % method

        self.index = np.empty((0, self.extent.ncol), dtype=np.uint8)
        self.validObs_index = np.empty((0, self.extent.ncol), dtype=np.int8)

        # todo: remove negative values in normalized images by adding the minimum of all images
        # minimum = min([file.GetRasterBand(band).GetMinimum() for band in range(1, file.RasterCount+1)])

        for i, subext in enumerate(self.subextents):

            # Stack scenes with data and scene index number
            dims = (len(inputScenes), subext.nrow, subext.ncol)
            refStack = np.zeros(dims, dtype=[("data", "f4"),("ID", "u1"),("NDVI", "f4")])

            # Reference bands
            for j, sce in enumerate(inputScenes):

                # Get joint extent of file and subextent
                jointExtent = getJointExtent(sce.files[3], AOIextent=subext)

                # Read in bands depending on composite method
                if jointExtent == None:
                    ref = np.empty((subext.nrow, subext.ncol))
                    ref[:] = np.nan
                else:
                    if method == "blue":
                        ref = sce.getBand(bands[0], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    elif method == "maxNDVI":
                        nir = sce.getBand(bands[3], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        red = sce.getBand(bands[2], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        ref = (nir - red) / (nir + red)
                        del nir, red
                    elif method == "HOT":
                        blue = sce.getBand(bands[0], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        red = sce.getBand(bands[2], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        ref = blue - 0.5*red - 800.
                        del blue, red
                    elif method == "maxRatio":
                        blue = sce.getBand(bands[0], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        nir = sce.getBand(bands[3], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        swir1 = sce.getBand(bands[4], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        ref = np.nanmax(np.array([nir, swir1]), axis=0) / blue
                        del blue, nir, swir1
                    elif method == "maxNDWI":
                        nir = sce.getBand(bands[3], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        swir1 = sce.getBand(bands[4], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                        ref = (nir - swir1) / (nir + swir1)
                        del nir, swir1
                    else:
                        print "Method not valid."
                        return None

                    # Adjust size of array to size of subextent
                    geotrans = [jointExtent.ulX, self.pixSize, 0.0, jointExtent.ulY, 0.0, -self.pixSize]
                    if ref.shape[1] < subext.ncol or ref.shape[0] < subext.nrow:
                        ref = padArray(ref, geotrans, sce.proj, subext)

                # Write data to stack
                refStack["data"][j] = ref
                refStack["ID"][j] = sce.indexNr
                del ref

            # calculate composite index
            if method == "blue":
                index, valid_obs = indexByPercentile(refStack, 25)[1:3]
            elif method == "maxNDVI":
                median = bn.nanmedian(refStack["data"], axis=0)
                # todo check if results are the same
                index_max, valid_obs = indexByPercentile(refStack, 100)[1:3]
                index_min, valid_obs = indexByPercentile(refStack, 0)[1:3]
                index = np.where(median < -0.05, index_min, index_max)
            elif method == "HOT":
                # todo check if results are the same
                index, valid_obs = indexByPercentile(refStack, 25)[1:3]
            elif method == "maxRatio":
                index, valid_obs = indexByPercentile(refStack, 100)[1:3]
            elif method == "maxNDWI":
                index, valid_obs = indexByPercentile(refStack, 100)[1:3]

            # Append subextents to output arrays
            self.index = np.vstack([self.index, index])
            self.validObs_index = np.vstack([self.validObs_index, valid_obs])

            print "%.2f%% done ..." % ((float(i)+1)/float(len(self.subextents)) * 100)

            del refStack, index, valid_obs

        try:
            self.indexScenes = [self.sceneLookup[idx] for idx in np.unique(self.index) if idx != 255]
        except:
            print "ID not in sceneLookup."
            logging.critical("ID in index not in scene Lookup.")
            sys.exit(1)

    def calc_index_combi(self, inputScenes, bands=[1,2,3,4,5,6], transformed=True):
        """Calculate image composite index based a combined method."""

        print "\nCalculating composite index using method: Combi"

        self.index = np.empty((0, self.extent.ncol), dtype=np.int8)
        self.validObs_index = np.empty((0, self.extent.ncol), dtype=np.int8)

        # todo: remove negative values in normalized images by adding the minimum of all images
        # minimum = min([file.GetRasterBand(band).GetMinimum() for band in range(1, file.RasterCount+1)])

        # Scene ranking by cloud coverage
        inputScenes.sort(key=operator.attrgetter("cloudCoverage"), reverse=False)

        # Initial scene
        initScene = inputScenes[0]
        additionalScenes = [sce for sce in inputScenes if sce != initScene]

        for i, subext in enumerate(self.subextents):

            # Create empty stack for scene and features
            dims = (1, subext.nrow, subext.ncol)
            compStack = np.zeros(dims, dtype=[("ID", "u2"),("HOT", "f4"),("RATIO", "f4"),("NDVI", "f4"),("water", "i2"),("BLUE", "f4"), ("TCBI", "f4")])
            # Get joint extent of file and subextent
            jointExtent_init = getJointExtent(initScene.files[3], AOIextent=subext)

            # Read in bands depending on composite method
            if jointExtent_init == None:
                nan = np.empty((subext.nrow, subext.ncol))
                nan[:] = np.nan
                # Write data to stack
                compStack["ID"] = 255
                compStack["HOT"], compStack["NDVI"], compStack["water"], compStack["BLUE"], compStack["RATIO"], compStack["TCBI"] = nan, nan, nan, nan, nan, nan
            else:
                # Read bands
                blue = initScene.getBand(bands[0], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")
                green = initScene.getBand(bands[1], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")
                red = initScene.getBand(bands[2], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")
                nir = initScene.getBand(bands[3], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")
                swir1 = initScene.getBand(bands[4], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")
                swir2 = initScene.getBand(bands[5], extent=jointExtent_init, masked=True, transformed=transformed).astype("float32")

                # Read cfmask
                waterMask = initScene.getWatermask(extent=jointExtent_init)

                # Calculate features
                HOT = blue * math.sin(1.3) - red * math.cos(1.3)
                NDVI = (nir - red) / (nir + red)
                RATIO = swir1 / blue # np.nanmax(np.array([nir, swir1]), axis=0) / (blue + red)
                TCBI = (0.3029 * blue) + (0.2786 * green)+(0.4733 * red)+(0.5599 * nir)+(0.508 * swir1)+(0.1872 * swir2)

                # Adjust size of array to size of subextent
                geotrans = [jointExtent_init[0], self.pixSize, 0.0, jointExtent_init[1], 0.0, -self.pixSize]
                if HOT.shape[1] < subext.ncol or HOT.shape[0] < subext.nrow:
                    HOT = padArray(HOT, geotrans, initScene.proj, subext)
                    NDVI = padArray(NDVI, geotrans, initScene.proj,subext)
                    blue = padArray(blue, geotrans, initScene.proj,subext)
                    RATIO = padArray(RATIO, geotrans, initScene.proj,subext)
                    TCBI = padArray(TCBI, geotrans, initScene.proj,subext)
                    waterMask = padArray(waterMask, geotrans, initScene.proj,subext)

                # Write data to stack
                compStack["ID"] = np.where(~np.isnan(HOT), initScene.indexNr, 255)
                compStack["HOT"] = HOT
                compStack["NDVI"] = NDVI
                compStack["water"] = waterMask
                compStack["RATIO"] = RATIO
                compStack["BLUE"] = blue
                compStack["TCBI"] = TCBI

                # median of ND red Swir1
                # medianFile = os.path.join(self.scenes[0].tempDir, "NDRedSwir_median.tif")
                # if os.path.exists(medianFile):
                #     median = raster2array(medianFile, extent=subext)[0].astype("int16")
                # else:
                #     median = NDVI

                del HOT, RATIO, TCBI, blue, green, red, nir, swir1, swir2

            # Reference bands
            for j, sce in enumerate(additionalScenes):

                # Create empty stack for scene and features
                newStack = np.zeros(dims, dtype=[("ID", "u2"),("HOT", "f4"),("RATIO", "f4"),("water", "i2"),("NDVI", "f4"),("BLUE", "f4"), ("TCBI", "f4")])

                # Get joint extent of file and subextent
                jointExtent = getJointExtent(sce.files[3], AOIextent=subext)

                # Read in bands depending on composite method
                if jointExtent == None:
                    nan = np.empty((subext.nrow, subext.ncol))
                    nan[:] = np.nan
                    # Write data to stack
                    newStack["ID"] = sce.indexNr
                    newStack["HOT"], newStack["watermask"], newStack["NDVI"], newStack["BLUE"], newStack["RATIO"], newStack["TCBI"] = nan, nan, nan, nan, nan, nan
                else:
                    # Read bands
                    blue = sce.getBand(bands[0], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    green = sce.getBand(bands[1], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    red = sce.getBand(bands[2], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    nir = sce.getBand(bands[3], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    swir1 = sce.getBand(bands[4], extent=jointExtent, masked=True, transformed=transformed).astype("float32")
                    swir2 = sce.getBand(bands[5], extent=jointExtent, masked=True, transformed=transformed).astype("float32")

                    waterMask_sce = initScene.getWatermask(extent=jointExtent)
                    if waterMask_sce is None:
                        waterMask_sce = np.zeros((subext.nrow, subext.ncol))

                    # Calculate Features
                    HOT = blue * math.sin(1.3) - red * math.cos(1.3)
                    NDVI = (nir - red) / (nir + red)
                    RATIO = swir1 / blue #bn.nanmax(np.array([nir, swir1]), axis=0) / blue # (swir1 - red) / (swir1 + red) #bn.nanmax(np.array([nir, swir1]), axis=0) / blue
                    TCBI = (0.3029 * blue) + (0.2786 * green)+(0.4733 * red)+(0.5599 * nir)+(0.508 * swir1)+(0.1872 * swir2)

                    # Adjust size of array to size of subextent
                    geotrans = [jointExtent.ulX, self.pixSize, 0.0, jointExtent.ulY, 0.0, -self.pixSize]
                    if HOT.shape[1] < subext.ncol or HOT.shape[0] < subext.nrow:
                        HOT = padArray(HOT, geotrans, sce.proj, subext)
                        NDVI = padArray(NDVI, geotrans, sce.proj,subext)
                        blue = padArray(blue, geotrans, sce.proj,subext)
                        RATIO = padArray(RATIO, geotrans, sce.proj,subext)
                        TCBI = padArray(TCBI, geotrans, sce.proj,subext)
                        waterMask_sce = padArray(waterMask_sce, geotrans, sce.proj, subext)

                    # Write data to stack
                    newStack["ID"] = sce.indexNr
                    newStack["HOT"] = HOT
                    newStack["NDVI"] = NDVI
                    newStack["water"] = waterMask_sce
                    newStack["RATIO"] = RATIO
                    newStack["BLUE"] = blue
                    newStack["TCBI"] = TCBI

                    del HOT, RATIO, NDVI, TCBI, blue, green, red, nir, swir1, swir2, waterMask_sce

                # Choose best pixel out of compositeStack and newStack
                # add new pixels to composite where composite is nan
                compStack = np.where((compStack["ID"]==255) & (~np.isnan(newStack["HOT"])), newStack, compStack)

                # MASKS
                # Water mask
                water = np.where((compStack["water"] == 1) | (newStack["water"] == 1), 1, 0)
                # Shadow mask
                shadow = np.where((compStack["TCBI"] < 2000) | (newStack["TCBI"] < 2000), 1, 0)
                # HOT mask (HOT < 1000)
                HOTvalid = np.squeeze((compStack["HOT"] > 1000) | (newStack["HOT"] > 1000)).astype("bool")

                # REPLACE VALUES
                # On WATER use min(blue)
                compStack = np.where( water & (newStack["BLUE"] < compStack["BLUE"]), newStack, compStack)

                # On land and HOT is valid use min(HOT)
                compStack = np.where( ~water & HOTvalid & (newStack["HOT"] < compStack["HOT"]), newStack, compStack)
                #
                # # If HOT is not valid use ratio
                compStack = np.where( ~water & ~HOTvalid & shadow & (newStack["NDVI"] > compStack["NDVI"]), newStack, compStack)
                compStack = np.where( ~water & ~HOTvalid & ~shadow & (newStack["RATIO"] < compStack["RATIO"]), newStack, compStack)

                # Lueck 2016
                # compStack = np.where((compStack["ID"]==255) & ~np.isnan(newStack["HOT"]), newStack, compStack)
                # # ... using HOT
                # HOTvalid = (compStack["HOT"] > 1000) & (newStack["HOT"] > 1000)
                # compStack = np.where( HOTvalid & (newStack["HOT"] < compStack["HOT"]), newStack, compStack)
                # # ... using maxRatio
                # TCBIvalid = (compStack["TCBI"] < 5000) & (newStack["TCBI"] < 5000)
                # compStack = np.where( ~HOTvalid & TCBIvalid & (newStack["RATIO"] > compStack["RATIO"]), newStack, compStack)
                # # ... using maxNDVI
                # compStack = np.where( ~HOTvalid & ~TCBIvalid & (newStack["NDVI"] > compStack["NDVI"]), newStack, compStack)

            # Append subextents to output arrays
            self.index = np.vstack([self.index, np.squeeze(compStack["ID"])])
            #self.validObs_index = np.vstack([self.validObs_index, valid_obs])

            print "%.2f%% done ..." % ((float(i)+1)/float(len(self.subextents)) * 100)

            del compStack # , index, valid_obs

        self.indexScenes = [self.sceneLookup[idx] for idx in np.unique(self.index) if idx != 255]

    def export_composite(self, outfileDir, outfileName, transformed=False, validPixels=True):
        """Export composite image."""

        print "\nExporting composite:"

        # Create output files
        # ------------------
        outdriver = gdal.GetDriverByName("GTIFF")
        nodata = -32768

        # Output file for valid observations
        if validPixels:
            if not hasattr(self, 'validObs'):
                self.calc_validPixels()

            if not os.path.exists(outfileDir):
                raise IOError("Output folder does not exist.")

            dest = os.path.join(outfileDir, outfileName + "_validObs.tif")
            if os.path.exists(dest):
                os.unlink(dest)
            validFile = outdriver.Create(str(dest), self.extent.ncol, self.extent.nrow, 1, gdal.GDT_Byte)
            validFile.GetRasterBand(1).SetNoDataValue(255)
            validFile.SetGeoTransform(self.geotrans)
            validFile.SetProjection(self.proj)
            validFile.GetRasterBand(1).WriteArray(self.validObs)
            validFile.FlushCache()
            del validFile

        # output file for scene index
        dest = os.path.join(outfileDir, outfileName + "_sceneIndex.tif")
        if os.path.exists(dest):
            os.unlink(dest)
        indexFile = outdriver.Create(str(dest), self.extent.ncol, self.extent.nrow, 1, gdal.GDT_Byte)
        indexFile.GetRasterBand(1).SetNoDataValue(255)
        indexFile.SetGeoTransform(self.geotrans)
        indexFile.SetProjection(self.proj)
        indexFile.GetRasterBand(1).WriteArray(self.index)
        indexFile.FlushCache()
        del indexFile

        perc25 = np.empty(shape=(self.nbands, self.extent.nrow, self.extent.ncol), dtype="float32")
        perc25[:] = np.nan

        for j, sce in enumerate(self.indexScenes):

            print "%.2f%% done ..." % ((float(j)+1)/float(len(self.indexScenes)) * 100)

            jointExtent = getJointExtent(sce.files[3], AOIextent=self.extent)

            if jointExtent is not None:
                bands = []
                for b in range(1, self.nbands+1):
                    band = sce.getBand(bandNo=b, extent=jointExtent, masked=True, transformed=transformed).astype("Int16")
                    if band is None:
                        logging.error("Requested invalid band from image %s" % sce.ID)
                    else:
                        band = band.astype("int16")
                    if band.shape[1] < self.extent.ncol or band.shape[0] < self.extent.nrow:
                        band = padArray(band, sce.geotrans, sce.proj, self.extent)
                    bands.append(band)
                bands = np.array(bands)
                perc25 = np.where(self.index == sce.indexNr, bands, perc25)

        # Export band composite
        # Spectral bands
        dest = os.path.join(outfileDir, outfileName + ".tif")
        if os.path.exists(dest):
            os.unlink(dest)
        spectralFile = outdriver.Create(str(dest), self.extent.ncol, self.extent.nrow, self.nbands, gdal.GDT_Int16)
        spectralFile.SetGeoTransform(self.geotrans)
        spectralFile.SetProjection(self.proj)
        perc25 = np.where(np.isnan(perc25), nodata, perc25)

        for b in range(1, self.nbands+1):
            spectralFile.GetRasterBand(b).SetNoDataValue(nodata)
            spectralFile.GetRasterBand(b).WriteArray(perc25[b-1,:,:])
            spectralFile.FlushCache()

        del perc25

        del spectralFile

class Scene(object):
    """Basic Object for satellite scene."""
    def __init__(self, sceneDir, tempDir=None, extentAOI=None):

        # Scene description
        self.dir = sceneDir
        self.ID = os.path.basename(sceneDir)[:16]

        self.extent = extentAOI

        # Temp directory
        if tempDir == None:
            self.tempDir = self.dir
        else:
            self.tempDir = tempDir

        if os.path.exists(os.path.join(self.tempDir, "masks")):
            self.maskDir = os.path.join(self.tempDir, "masks")
        else:
            self.maskDir = self.tempDir

        if os.path.exists(os.path.join(self.tempDir, "histogramMatching")):
            self.histoDir = os.path.join(self.tempDir, "histogramMatching")
        else:
            self.histoDir = self.tempDir

    def getMetadata(self, band=0):

        # Get extent
        #todo raise ValueError if extent is None

        self.extent = getJointExtent(self.files[band], AOIextent=self.extent)
        if self.extent == None:
            raise Exception("Scene not in extent.")

        # Get metadata
        metadata = gdal.Open(self.files[band], GA_ReadOnly)
        self.geotrans = list(metadata.GetGeoTransform())
        self.geotrans[0] = self.extent.ulX
        self.geotrans[3] = self.extent.ulY
        self.pixSize = self.geotrans[1]
        self.proj = metadata.GetProjection()
        self.nodata = metadata.GetRasterBand(1).GetNoDataValue()
        del metadata

        # Initialize variables for later reference
        transformedScene = fnmatch.filter(os.listdir(self.histoDir), "*" + self.ID + "_matched.tif")
        if len(transformedScene) == 1:
            self.transformed = True
            self.transformedScene = os.path.join(self.histoDir, transformedScene[0])
        else:
            self.transformed = False
            self.transformedScene = ""

        self.VRTfile = ""
        self.indexNr = None
        self.cloudCoverage = None
        self.deltaDays = None
        self.bandMin = None

    def updateExtent(self, extent, band=0):
        # Get extent
        #todo raise ValueError if extent is None
        self.extent = getJointExtent(self.files[band], AOIextent=extent)
        if self.extent == None:
            raise Exception("Scene not in extent.")

    def getMinimumSpectralValue(self):

        if self.bandMin is None:
            # Band minimum value
            mins = []
            for f in self.files:
                openF = gdal.Open(f)
                mins.append(openF.GetRasterBand(1).GetStatistics(0,1)[0])
                del openF
            self.bandMin = min(mins)

        return self.bandMin

    def getBand(self, bandNo=1, masked=True, extent=None, transformed=False ):
        """Get specified band of scene as numpy array."""

        if extent == None:
            extent = self.extent

        if transformed and os.path.exists(self.transformedScene):
            band = raster2array(self.transformedScene, AOIExtent=extent, bandNo=bandNo)[0]
        else:
            band = raster2array(self.files[bandNo-1], AOIExtent=extent)[0]

        if band is None:
            raise ValueError("Requested band %s for scene %s could not be opened." % (bandNo, self.ID) )

        if masked:
            mask = self.getMask(extent=extent)
            if mask is not None:
                band = np.where((mask == 1) | (band == self.nodata), np.nan, band)

        return band

    def createVRT(self, outDir=None):

        if outDir is None:
            outDir = self.tempDir

        self.VRTfile = os.path.join(outDir, self.ID + ".vrt")

        if not os.path.exists(self.VRTfile) or (os.path.exists(self.VRTfile) and ((gdal.Open(self.VRTfile).RasterXSize != self.extent.ncol) or (gdal.Open(self.VRTfile).RasterYSize != self.extent.nrow))):
            cellSize = self.geotrans[1]
            extentVRT = [self.extent.ulX, self.extent.ulY - cellSize * self.extent.nrow, self.extent.ulX + cellSize * self.extent.ncol, self.extent.ulY]
            if os.path.exists(self.fmask):
                cmd = "gdalbuildvrt -resolution lowest -separate -te " + str(extentVRT[0]) + " " + str(extentVRT[1]) + " " + str(extentVRT[2]) + " " + str(extentVRT[3]) + " " + self.VRTfile + " " +  " ".join(self.files) + " " + self.fmask
            else:
                cmd = "gdalbuildvrt -resolution lowest -separate -te " + str(extentVRT[0]) + " " + str(extentVRT[1]) + " " + str(extentVRT[2]) + " " + str(extentVRT[3]) + " " + self.VRTfile + " " +  " ".join(self.files)
            #os.system(cmd)
            subprocess.call(cmd, shell=True) #universal_newlines=True, stderr=subprocess.STDOUT,

class SentinelScene(Scene):
    """Object for Sentinel scene """

    def __init__(self, sceneDir, tempDir=None, extentAOI=None):

        Scene.__init__(self, sceneDir, tempDir, extentAOI)

        # Scene description
        self.ID = os.path.basename(sceneDir) + "_" + os.path.basename(os.path.dirname(sceneDir))
        self.date = dt.strptime(os.path.basename(os.path.dirname(sceneDir)), "%Y%m%d")

        self.getFiles()
        self.getMetadata(band=4)

        self.nodata = 0.0

    def getFiles(self):

        # todo: error handling when if files/masks are found

        # Input files
        self.files = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*B0[2345678]*.jp2")]
        self.files += [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*B8A*.jp2")]
        self.files += [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*B1[12]*.jp2")]

        # Mask files
        fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_fmask.tif")]
        if len(fmask) != 0:
            self.fmask = fmask[0]
        else:
            fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_fmask.img")]
            if len(fmask) == 0:
                self.fmask = ""
            else:
                self.fmask = fmask[0]

        shadowMask = [os.path.join(self.maskDir, b) for b in fnmatch.filter(os.listdir(self.maskDir), self.ID + "_shadowMask.tif")]
        if len(shadowMask) != 0:
            self.shadowMask = shadowMask[0]
        else:
            self.shadowMask = ""

    def getMask(self, extent=None, withShadow=True):
        """Get cloud mask of scene as numpy array."""

        if extent is None:
            extent = self.extent

        if os.path.exists(self.fmask):
            fmask = raster2array(self.fmask, extent)[0]
        else:
            fmask = None

        if withShadow:
            if os.path.exists(self.shadowMask):
                shadowMask = raster2array(self.shadowMask, extent)[0]
            else:
                shadowMask = None

            # Mask out 2 - clouds, 3 - shadow and (4 - snow --> Not for Africa)
            # todo include snow in cloud mask for non african sites
            if fmask is not None and shadowMask is not None:
                mask = np.where((fmask == 2) | (fmask == 3) | (shadowMask == 1), 1, 0)
                mask = np.where((np.isnan(fmask) & np.isnan(shadowMask)), np.nan, mask)
            elif fmask is None and shadowMask is not None:
                mask = shadowMask
            elif shadowMask is None and fmask is not None:
                mask = np.where((fmask == 2) | (fmask == 3), 1, 0)
                mask = np.where(np.isnan(fmask), np.nan, mask)
            else:
                logging.warning("No masks found.")
                return None
            del fmask, shadowMask
        else:
            if fmask is not None:
                mask = np.where((fmask == 2) | (fmask == 3), 1, 0)
                mask = np.where(np.isnan(fmask), np.nan, mask)
            else:
                logging.warning("No masks found.")
                return None
            del fmask

        # Buffer masks
        mask = binaryBuffer(mask, size=1)
        mask = removeNoise_oneSide(mask)

        return mask

    def getWatermask(self, extent=None):
        """Get cloud mask of scene as numpy array."""

        if extent is None:
            extent = self.extent

        if os.path.exists(self.fmask):
            fmask = raster2array(self.fmask, extent)[0]
        else:
            fmask = None

        # Mask out 2 - clouds, 3 - shadow and 4 - snow
        if fmask is not None:
            mask = np.where((fmask == 5), 1, 0).astype("int8")

        else:
            logging.warning("No fmask found.")
            return None

        del fmask
        return mask

    def calcCloudCoverage(self, AOIextent=None):

        if self.cloudCoverage == None:

            if AOIextent is None:
                AOIextent = self.extent

            mask = self.getMask(extent=AOIextent)
            if mask is None:
                return

            # Compute cloud
            # todo: repace with bottleneck if running on qgis in version 1.0
            cols = AOIextent.ncol
            rows = AOIextent.nrow
            self.cloudCoverage = 100. - (float(np.nansum(mask == 0)) / (float(cols) * float(rows))) * 100.

class LandsatScene(Scene):
    """Object for Landsat scene (Inherits from Scene)"""

    def __init__(self, sceneDir, tempDir=None, extentAOI=None):

        Scene.__init__(self, sceneDir, tempDir, extentAOI)

        # Get files of spectral bands
        self.ID = os.path.basename(sceneDir)[:16]
        self.date = dt.strptime(self.ID[9:16], "%Y%j")

        self.getFiles()
        self.getMetadata(band=0)

    def getFiles(self):

        # Files
        dirname = os.path.basename(self.dir)
        self.files = []
        if fnmatch.fnmatch(dirname, "L[CO]*"):
            for filename in fnmatch.filter(os.listdir(self.dir), "*_toa_band[234567]*.tif"):
                self.files.append(os.path.join(self.dir, filename))
        elif fnmatch.fnmatch(dirname, "L[ET]*"):
            for filename in fnmatch.filter(os.listdir(self.dir), "*_toa_band[123457]*.tif"):
                self.files.append(os.path.join(self.dir, filename))
        self.files = sorted(self.files)

        # Mask files
        fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_cfmask_epsg3035.tif")]
        if len(fmask) != 0:
            self.fmask = fmask[0]
        else:
            self.fmask = ""

        shadowMask = [os.path.join(self.maskDir, b) for b in fnmatch.filter(os.listdir(self.maskDir), "*_shadowmask.tif")]
        if len(shadowMask) != 0:
            self.shadowMask = shadowMask[0]
        else:
            self.shadowMask = ""

    def getBand(self, bandNo=1, masked=True, extent=None, transformed=False ):
        """Get specified band of scene as numpy array."""

        if extent == None:
            extent = self.extent

        if transformed and os.path.exists(self.transformedScene):
            band = raster2array(self.transformedScene, extent, bandNo=bandNo)[0]
        else:
            band = raster2array(self.files[bandNo-1], extent)[0]

        if band is None:
            return None

        if masked:
            mask = self.getMask(extent=extent)
            band = np.where((mask == 1) | (band == self.nodata) | (band == 20000), np.nan, band)

        return band

    def getMask(self, extent=None):
        """Get cloud mask of scene as numpy array."""

        if extent is None:
            extent = self.extent

        if os.path.exists(self.fmask):
            fmask = raster2array(self.fmask, extent)[0]
        else:
            fmask = None

        if os.path.exists(self.shadowMask):
            shadowMask = raster2array(self.shadowMask, extent)[0]
        else:
            shadowMask = None

        # Mask out 2 - clouds, 3 - shadow and 4 - snow
        if fmask is not None and shadowMask is not None:
            mask = np.where((fmask == 2) | (fmask == 3) | (fmask == 4 ) | (shadowMask == 1), 1, 0)
            mask = np.where((np.isnan(fmask) | np.isnan(shadowMask)), np.nan, mask)
        elif fmask is None and shadowMask is not None:
            mask = shadowMask
        elif shadowMask is None and fmask is not None:
            mask = np.where((fmask == 2) | (fmask == 3) | (fmask == 4 ), 1, 0)
            mask = np.where(np.isnan(fmask), np.nan, mask)
        else:
            logging.warning("No masks found.")
            return None

        # Buffer masks
        mask = binaryBuffer(mask, size=1)
        mask = removeNoise_oneSide(mask)

        del fmask, shadowMask

        return mask

    def calcCloudCoverage(self):

        if self.cloudCoverage == None:

            fmask, shadowMask = None, None

            if os.path.exists(self.fmask):
                fmask = raster2array(self.fmask, self.extent)[0]

            if os.path.exists(self.shadowMask):
                shadowMask = raster2array(self.shadowMask, self.extent)[0]

            # Mask out 2 - clouds, 3 - shadow and 4 - snow
            if (shadowMask is not None) and (fmask is not None):
                mask = np.where((fmask == 2) | (fmask == 3) | (fmask == 4 ) | (shadowMask == 1) | np.isnan(shadowMask) | np.isnan(fmask), 1, 0)
            elif (fmask is None) and (shadowMask is not None):
                mask = np.where( np.isnan(shadowMask) | (shadowMask==1), 1, 0)
            elif (shadowMask is None) and (fmask is not None):
                mask = np.where((fmask == 2) | (fmask == 3) | (fmask == 4 ) | np.isnan(fmask), 1, 0)
            else:
                logging.warning("No masks found.")
                return None

            # Compute cloud coverate
            self.cloudCoverage = 100. - (float(bn.nansum(mask == 0)) / (float(self.extent.ncol) * float(self.extent.nrow))) * 100.

class extent(object):

    def __init__(self, ulX, ulY, ncol, nrow, pixSize, proj, lrX=None, lrY=None):
        self.ulX = ulX
        self.ulY = ulY
        self.ncol = ncol
        self.nrow = nrow
        self.pixSize = pixSize
        self.proj = proj
        self.lrX = lrX
        self.lrY = lrY

        # FUNCTIONS =============================================================

def improveFmask(sce):
    """ Remove commission errors over water from fmask"

    :param scene:
    :return:
    """
    try:
        fmask = raster2array(sce.fmask)[0]
        red = sce.getBand(3, masked=False, transformed=False).astype("float32")
        nir = sce.getBand(8, masked=False, transformed=False).astype("float32")
        swir1 = sce.getBand(10, masked=False, transformed=False).astype("float32")
    except ValueError,e:
        logging.warning("Scene %s requested band could not be opened.\n Error message:\n %s" % (sce.ID,e))

    # Calculate NDVI
    NDVI = (nir - red) / (red + nir)

    # Replace values in cloud mask
    betterFmask = np.where( (swir1 < 600) & ((fmask == 2) | (fmask == 4)), 5, fmask) #(mNDWI < -0.2) &
    betterFmask = np.where( (NDVI < -0) & (fmask == 3), 5, betterFmask)

    # Export file
    try:
        dest = os.path.join(sce.dir, sce.ID + "_fmask.tif")
        array2raster(betterFmask, sce.geotrans, sce.proj, dest, gdal.GDT_Byte, 255)
    except:
        print "Scene %s: could not write new fmask to file" % sce.ID
        logging.warning("Scene %s: could not write new fmask to file" % sce.ID)

    sce.fmask = dest

    return 0

# cloud masking

def createCloudWaterMask(redFile, nirFile, extent):

    # Mask clouded pixels over water
    red = raster2array(redFile, extent)[0]
    nir = raster2array(nirFile, extent)[0]

    NDVI = (nir.astype("float32") - red.astype("float32")) / (nir.astype("float32") + red.astype("float32"))
    invalidPixelMask = np.where((NDVI < 0.1) & (nir > 250), 1, 0)

    del red, nir, NDVI

    return invalidPixelMask

def createCloudMask(maskDir, extent):
    """ Create joint cloud mask based on fmask and BQA (LS8 only) mask."""

    maskFiles = [os.path.join(maskDir, f) for f in os.listdir(maskDir) if f.endswith('_BQA.TIF') or f.endswith('_cfmask_epsg3035.tif')]
    if len(maskFiles)==0:
        return None

    BQAfile = fnmatch.filter(maskFiles, "*_BQA*.TIF")
    cfmaskFile = fnmatch.filter(maskFiles, "*_cfmask*.tif")

    jointExtent = getJointExtent(maskFiles, AOIextent=extent)

    # Create binary mask from BQA file
    if BQAfile:
        BQAfile_ = os.path.join(maskDir, BQAfile[0])
        bqa, geotrans_bqa = raster2array(BQAfile_, jointExtent)
        if not np.all(np.isnan(bqa)):
            bqamask = createBQAMask(bqa)
        if bqamask.shape[1] < extent.ncol or bqamask.shape[0] < extent.nrow:
            bqamask = padArray(bqamask, geotrans_bqa, extent)
    else:
        bqamask = np.ones((extent.nrow, extent.ncol), dtype="float16")
        bqamask[:] = np.nan

    # Create binary mask form fmask file
    if cfmaskFile:
        cfmaskfile_ = os.path.join(maskDir, cfmaskFile[0])
        fmask, geotrans_cf = raster2array(cfmaskfile_, jointExtent)
        if not np.all(np.isnan(fmask)):
            fmask = createfmask(fmask)

        if fmask.shape[1] < extent.ncol or fmask.shape[0] < extent.nrow:
            fmask = padArray(fmask, geotrans_cf, extent)
    else:
        fmask = np.ones((extent.nrow, extent.ncol), dtype="float16")
        fmask[:] = np.nan

    return np.nanmax(np.array([bqamask, fmask]), axis=0)

def createfmask_sentinel(fmaskFile, extent):
        """Create Cloud mask from cfmask"""

        jointExtent = getJointExtent(fmaskFile, AOIextent=extent)

        # Read file
        fmask, geotrans_cf = raster2array(fmaskFile, jointExtent)[:2]
        if not np.all(np.isnan(fmask)):
            fmask = np.where((fmask == 2) | (fmask == 3) | (fmask == 255) | np.isnan(fmask), 1, 0).astype("int8")

        if fmask.shape[1] < extent.ncol or fmask.shape[0] < extent.nrow:
            fmask = padArray(fmask, geotrans_cf, jointExtent.proj, extent)

        return (fmask)

def createBQAMask(bqa):
    """Create Cloud Mask from Landsat 8 BQA band."""
    mask = np.where((bqa >= 24576) | (bqa == 1), 1, 0).astype("int8")
    mask = np.where(np.isnan(bqa), np.nan, mask)
    return mask

def createfmask(cfmask):
    """Create Cloud mask from cfmask"""
    mask = np.where((cfmask == 2) | (cfmask == 4), 1, 0).astype("int8")
    mask = np.where((cfmask == 255) | (np.isnan(cfmask)), np.nan, mask)
    return mask

# Binary image methods

def removeNoise(binary_img, iterations=2):
    """Remove noise in binary image."""
    eroded_img = ndimage.binary_erosion(binary_img, iterations=iterations)
    reconstruct_img = ndimage.binary_propagation(eroded_img, mask=binary_img)
    tmp = np.logical_not(reconstruct_img)
    eroded_tmp = ndimage.binary_erosion(tmp, iterations=iterations)
    reconstruct_final = np.logical_not(ndimage.binary_propagation(eroded_tmp, mask=tmp))
    return np.where(np.isnan(binary_img), np.nan, reconstruct_final )

def removeNoise_oneSide(binary_img):
    """Remove noise in binary image. Only pixels that are 0 are removed. Pixels of value 1 are all preserved. (used for cloud masks) """
    #eroded_img = ndimage.binary_erosion(binary_img, iterations=3)
    #reconstruct_img = ndimage.binary_propagation(eroded_img, mask=binary_img)
    tmp = np.logical_not(binary_img)
    eroded_tmp = ndimage.binary_erosion(tmp, iterations=5)
    reconstruct_final = np.logical_not(ndimage.binary_propagation(eroded_tmp, mask=tmp))
    return np.where(np.isnan(binary_img), np.nan, reconstruct_final )

def binaryBuffer(binary_img, size=1):
    """Buffer pixels with value 1 in binary image."""
    struct2 = ndimage.generate_binary_structure(2, 2)
    buffered_img = ndimage.morphology.binary_dilation(binary_img, iterations=size, structure=struct2)
    return np.where(np.isnan(binary_img), np.nan, buffered_img )

def binaryBuffer_negative(binary_img, size=1):
    """Buffer pixels with value 0 in binary image."""
    struct2 = ndimage.generate_binary_structure(2, 2)
    buffered_img = np.logical_not(ndimage.morphology.binary_dilation(np.logical_not(binary_img), iterations=size, structure=struct2))
    return np.where(np.isnan(binary_img), np.nan, buffered_img )

# Raster Input / Output

def raster2array(rasterPath, AOIExtent=None, bandNo=1):
    """ Load and clip band from raster
    Load a band from a raster image and clip it to extent of shapefile

    :param: rasterPaths, extent

    :rtype: object: ndarray, newGeoTrans
    """


    #todo check if raster and rasterExtent have same projection

    # Open file
    try:
        rasterF = gdal.Open(rasterPath, GA_ReadOnly )
    except:
        return (None, None, None)

    # Check if band number exists
    if bandNo > rasterF.RasterCount:
        logging.warning("Requested band %s of raster file %s does not exist." % (bandNo, rasterPath))
        raise ValueError("Requested band %s of raster file %s does not exist." % (bandNo, rasterPath))

    # Get geotransformation
    geoTrans = rasterF.GetGeoTransform()
    proj = rasterF.GetProjection()
    sourceProj = osr.SpatialReference()
    sourceProj.ImportFromWkt(proj)
    pixSize = geoTrans[1]

    if AOIExtent:

        if AOIExtent.proj != sourceProj:
            ext_ulX, ext_ulY = reprojectPoint(AOIExtent.proj, sourceProj, AOIExtent.ulX, AOIExtent.ulY)
            ext_lrX, ext_lrY = reprojectPoint(AOIExtent.proj, sourceProj, AOIExtent.lrX, AOIExtent.lrY)
        else:
            ext_ulX = AOIExtent.ulX
            ext_ulY = AOIExtent.ulY
            ext_lrX = AOIExtent.lrX
            ext_lrY = AOIExtent.lrY

        ext_ncol = AOIExtent.ncol
        ext_nrow = AOIExtent.nrow
        ext_pixSize = AOIExtent.pixSize

        # LR coordinates of raster
        lrX = geoTrans[0] + rasterF.RasterXSize * pixSize
        lrY = geoTrans[3] - rasterF.RasterYSize * pixSize

        # Convert geogr. coordinates to image coordinates
        ras_ulX, ras_ulY = world2Pixel(geoTrans, ext_ulX, ext_ulY)

        # Check if extent is valid
        if (ras_ulY < 0) or (ras_ulX < 0) or (lrX < ext_lrX) or (lrY > ext_lrY):
            logging.warning("Requested invalid extent for raster file %s." % rasterPath)
            raise ValueError("Requested invalid extent for raster file %s." % rasterPath)

        # Number of cols and rows
        ncol = int((ext_lrX - ext_ulX) / pixSize)
        nrow = int((ext_ulY - ext_lrY) / pixSize)

        # Import clip extent from raster file
        band = rasterF.GetRasterBand(bandNo).ReadAsArray(ras_ulX, ras_ulY, ncol, nrow)

        if ext_pixSize < pixSize:
            zoom = pixSize / ext_pixSize
            band = ndimage.zoom(band, zoom, order=0)
        elif ext_pixSize > pixSize:
            band = downsample(band, pixSize, ext_pixSize)

    else:
        ras_ulX, ras_ulY = 0, 0
        ncol = rasterF.RasterXSize
        nrow = rasterF.RasterYSize

        # Import clip extent from raster file
        band = rasterF.GetRasterBand(bandNo).ReadAsArray(ras_ulX, ras_ulY, ncol, nrow)

    # Mask nodata pixel
    nodataVal = rasterF.GetRasterBand(bandNo).GetNoDataValue()
    if nodataVal is not None:
        band = np.where((band == nodataVal) | (band == 20000.), np.nan, band)

    # Update Geotransformation to export clipped raster image
    newGeoTrans = list(geoTrans)
    newGeoTrans[0] = (ras_ulX * pixSize) + geoTrans[0]
    newGeoTrans[3] = (ras_ulY * (-1) * pixSize) + geoTrans[3]

    del rasterF

    return (band, newGeoTrans, proj)

def array2raster(array, trans, proj, dest, dataType, nodataVal):
    """Array > Raster
    Save a raster from a C order array.

    :param array: ndarray, trans, proj, pathOut, dataType, nan value
    """

    # Replace nodata values
    array = np.where(np.isfinite(array), array, nodataVal)

    # Check dimensions of array
    dims = len(array.shape)
    if dims == 2:
        bandCount = 1
        bandCols = array.shape[1]
        bandRows = array.shape[0]
        array = array[np.newaxis,...]
    elif dims == 3:
        bandCount = array.shape[0]
        bandCols = array.shape[2]
        bandRows = array.shape[1]
    else:
        print "Too many dimensions."
        return -1

    # Check if file exists, delete it if it does
    if os.path.exists(dest):
         os.remove(dest)

    outdriver = gdal.GetDriverByName("GTIFF")

    # create empty output file
    outFile = outdriver.Create(
                            str(dest),
                            bandCols,
                            bandRows,
                            bandCount,
                            dataType)

    outFile.SetGeoTransform(trans)
    outFile.SetProjection(proj)

    # replace np.nan by nodata value
    array = np.where(np.isnan(array), nodataVal, array)

    for band in range(bandCount):
        outFile.GetRasterBand(band+1).SetNoDataValue(nodataVal)
        outFile.GetRasterBand(band+1).WriteArray(array[band])
        outFile.FlushCache()  # Write to disk.

    del outFile

    return True

def scene2raster(scene, dest, dataType, nodataVal):
    """Array > Raster
    Save a raster from a C order array.

    :param array: ndarray, trans, proj, pathOut, dataType, nan value
    """

    bandCount = len(scene.files)
    array = []
    mask = scene.getMask()
    for b in range(1, bandCount+1):
        band = scene.getBand(b, masked=False, transformed=True).astype("float32")
        band = np.where(mask == 1, np.nan, band)
        array.append(band)

    array = np.array(array)

    # Replace nodata values
    array = np.where(np.isfinite(array), array, nodataVal)

    # Check dimensions of array
    dims = len(array.shape)
    if dims == 2:
        bandCount = 1
        bandCols = array.shape[1]
        bandRows = array.shape[0]
        array = array[np.newaxis,...]
    elif dims == 3:
        bandCount = array.shape[0]
        bandCols = array.shape[2]
        bandRows = array.shape[1]
    else:
        print "Too many dimensions."
        return -1

    # Check if file exists, delete it if it does
    if os.path.exists(dest):
         os.remove(dest)

    outdriver = gdal.GetDriverByName("GTIFF")

    # create empty output file
    outFile = outdriver.Create(
                            str(dest),
                            bandCols,
                            bandRows,
                            bandCount,
                            dataType)

    outFile.SetGeoTransform(scene.geotrans)
    outFile.SetProjection(scene.proj)

    for band in range(bandCount):
        outFile.GetRasterBand(band+1).SetNoDataValue(nodataVal)
        outFile.GetRasterBand(band+1).WriteArray(array[band])
        outFile.FlushCache()  # Write to disk.

    del outFile, band, mask

    return True

def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    pixel = int(np.floor((x - ulX) / xDist))
    line = int(np.floor((ulY - y) / xDist))
    return (pixel, line)

def reprojectPoint(sourceProj, targetProj, xCoord, yCoord):

    transform = osr.CoordinateTransformation(sourceProj, targetProj)

    point = ogr.CreateGeometryFromWkt("POINT (" + str(xCoord) + " " + str(yCoord) + ")")

    point.Transform(transform)

    return round(point.GetX()), round(point.GetY())

def padArray(array, geoTransOfArray, projArray, AOIextent):

    projRaster = osr.SpatialReference()
    projRaster.ImportFromWkt(projArray)

    if not projRaster.IsSame(AOIextent.proj):
        ext_ulX, ext_ulY = reprojectPoint(AOIextent.proj, projRaster, AOIextent.ulX, AOIextent.ulY)
    else:
        ext_ulX = AOIextent.ulX
        ext_ulY = AOIextent.ulY

    ext_cols = AOIextent.ncol
    ext_rows = AOIextent.nrow

    arr_ulX = geoTransOfArray[0]
    arr_ulY = geoTransOfArray[3]
    arr_rows = array.shape[0]
    arr_cols = array.shape[1]

    pixSize = geoTransOfArray[1]

    if ext_rows != arr_rows:
        if ext_ulY > arr_ulY:
            N_o = int((ext_ulY - arr_ulY) / pixSize)
            N_u = max(ext_rows - N_o - arr_rows, 0)
        else:
            N_o = 0
            N_u = max(ext_rows - arr_rows, 0)
    else:
        N_o = 0
        N_u = 0

    if ext_cols != arr_cols:
        if ext_ulX < arr_ulX:
            N_l = int((arr_ulX - ext_ulX) / pixSize)
            N_r = max(ext_cols - N_l - arr_cols, 0)
        else:
            N_l = 0
            N_r = max(ext_cols - arr_cols, 0)
    else:
        N_l = 0
        N_r = 0

    newArray = np.pad(array.astype(float),((N_o,N_u),(N_l,N_r)), mode='constant', constant_values=np.nan)

    return newArray

# File system search for scene directories

def searchSceneDirs_sentinel(inDir):
    sceneDirs = []
    for root, dirs, files in os.walk(inDir):
        for f in files:
            if f.endswith("metadata.xml"):
                sceneDirs.append(root)
    return sceneDirs

# todo adapt to longer periods
def getSceneDirsForSeason_sentinel(sceneDirs, year, month, period=0):

    #start_month = (month - period)
    #end_month = (month + period + 1)
    months = [12,1,2,3,4,5,6,7,8,9,10,11,12,1]

    # Images for season
    seasonFiles =[]
    for d in sceneDirs:
        date = dt.strptime(os.path.basename(os.path.dirname(d)), "%Y%m%d")
        if year == 0:
            if date.month in months[month-period:month+period+1]:
                seasonFiles.append(d)
        elif month == 0:
            if (date.year == year):
                seasonFiles.append(d)
        elif period != 0:
            if month == 1:
                if (date.year == year-1 and date.month == 12) or (date.year == year and date.month == month) or (date.year == year and date.month == month+1):
                    seasonFiles.append( d)
            elif month == 12:
                if (date.year == year+1 and date.month == 1) or (date.year == year and date.month == month) or (date.year == year and date.month == month-1):
                    seasonFiles.append( d)
            elif (date.year == year) and (date.month in months[month-period: month+period+1]):
                seasonFiles.append(d)
        else:
            if(date.year == year) and (date.month == month):
                seasonFiles.append(d)

    return sorted(seasonFiles)

def getSceneDirsForSeason(sceneDirs, year, month, period=1):

    months = [12,1,2,3,4,5,6,7,8,9,10,11,12,1]

    # Images for season
    seasonFiles =[]
    for dir in sceneDirs:
        date = dt.strptime(os.path.basename(dir)[9:16], "%Y%j")
        if year == 0:
            if (date.month in months[month-period:month+period+1]):
                seasonFiles.append( dir)
        elif month == 0:
            if (date.year == year):
                seasonFiles.append( dir)
        else:
            if month == 1:
                if (date.year == year-1 and date.month == 12) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month+1):
                    seasonFiles.append( dir)
            elif month == 12:
                if (date.year == year+1 and date.month == 1) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month-1):
                    seasonFiles.append( dir)
            else:
                if (date.year == year) and (date.month in range(month-period, month+period+1)):
                    seasonFiles.append( dir)

    return sorted(seasonFiles)

def getSceneDirs(rootDir, sceneList):
    """Get all scene directories of scenes in a list of scene IDs."""
    # Create full scene paths
    sceneDirs = []
    for sce in sceneList:
        scenePath = sce[3:6]
        sceneRow = sce[6:9]
        sensor = sce[0:3]
        sceneRoot = os.path.join(rootDir, scenePath + "_" + sceneRow, sensor)
        #sceneDir = [os.path.join(sceneRoot, dir) for dir in os.listdir(sceneRoot) if dir.startswith(sce[0:16]) and os.path.isdir(os.path.join(sceneRoot, dir))]
        sceneDir = os.path.join(sceneRoot, sce[0:16])
        if os.path.exists(sceneDir):
            sceneDirs.append(sceneDir)
        else:
            logging.warning("%s not found." % sce[0:16])

    return sceneDirs


# adapt to sentinel
def filterSceneDirs_Sentinel(directory, year, month, sceneDirs=None):

    # Images for season
    sceneDirs_season =[]

    if sceneDirs is None:
        sceneDirs = [dir for dir in fnmatch.filter(os.listdir(directory), "L*") if os.path.isdir(os.path.join(directory, dir))]

    for dir in sceneDirs:
        date = dt.strptime(os.path.basename(os.path.dirname(dir)), "%Y%m%d")
        if month == 1:
            if (date.year == year-1 and date.month == 12) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month+1):
                sceneDirs_season.append(os.path.join(directory, dir))
        elif month == 12:
            if (date.year == year+1 and date.month == 1) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month-1):
                sceneDirs_season.append(os.path.join(directory, dir))
        elif month == 0:
            if (date.year == year):
                sceneDirs_season.append(os.path.join(directory, dir))
        else:
            if (date.year == year) and (date.month in range(month-1, month+2)):
                sceneDirs_season.append(os.path.join(directory, dir))

    return sorted(sceneDirs_season)

def filterSceneDirs_Landsat(directory, year, month, sceneDirs=None):

    # Images for season
    sceneDirs_season =[]

    if sceneDirs is None:
        sceneDirs = [dir for dir in fnmatch.filter(os.listdir(directory), "L*") if os.path.isdir(os.path.join(directory, dir))]

    for dir in sceneDirs:
        date = dt.strptime(os.path.basename(dir)[9:16], "%Y%j")
        if month == 1:
            if (date.year == year-1 and date.month == 12) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month+1):
                sceneDirs_season.append(os.path.join(directory, dir))
        elif month == 12:
            if (date.year == year+1 and date.month == 1) or (date.year == year   and date.month == month) or (date.year == year   and date.month == month-1):
                sceneDirs_season.append(os.path.join(directory, dir))
        elif month == 0:
            if (date.year == year):
                sceneDirs_season.append(os.path.join(directory, dir))
        else:
            if (date.year == year) and (date.month in range(month-1, month+2)):
                sceneDirs_season.append(os.path.join(directory, dir))

    return sorted(sceneDirs_season)

def findSceneDirsByList_Landsat(rootDir, sceneList):

    # Create full scene paths
    sceneDirs = []
    for sce in sceneList:
        #scenePath = sce[3:6]
        #sceneRow = sce[6:9]
        # todo adapt if folder structure changes (subfolders for row and path)
        #sceneRoot = os.path.join(rootDir, scenePath + "_" + sceneRow)
        sceneDir = [os.path.join(rootDir, dir) for dir in os.listdir(rootDir) if dir.startswith(sce) and os.path.isdir(os.path.join(rootDir, dir))]
        if len(sceneDir) > 1:
            raise Warning ("Too many directories for scene %s found" % sce)
        elif len(sceneDir) < 1:
            raise Warning ("%s not found." % sceneDir)
        else:
            if os.path.exists(sceneDir[0]):
                sceneDirs.append(sceneDir[0])
            else:
                raise Warning ("%s not found." % sceneDir)

    return sceneDirs

def findSceneDirsByList_Sentinel(rootDir, sceneList):

    # Create full scene paths
    sceneDirs = []
    for sce in sceneList:
        sceneID = sce[:5]
        sceneDate = sce[6:14]
        sceneDir = [os.path.join(rootDir, dir) for dir in fnmatch.filter(os.listdir(rootDir), sceneID + "*" + sceneDate + "*")]
        if len(sceneDir) > 1:
            raise Warning ("Too many directories for scene %s found" % sce)
        elif len(sceneDir) < 1:
            raise Warning ("%s not found." % sceneDir)
        else:
            if os.path.exists(sceneDir[0]):
                sceneDirs.append(sceneDir[0])
            else:
                raise Warning ("%s not found." % sceneDir)

    return sceneDirs

# Extent calculation

def intersectRasterExtents(extent1, extent2):

    # todo check if projections of extents are the same
    ext1_ulX = extent1.ulX
    ext1_ulY = extent1.ulY
    ext1_ncol = extent1.ncol
    ext1_nrow = extent1.nrow

    if extent1.proj != extent2.proj:
        ext2_ulX, ext2_ulY = reprojectPoint(extent2.proj, extent1.proj, extent2.ulX, extent2.ulY)
    else:
        ext2_ulX = extent2.ulX
        ext2_ulY = extent2.ulY

    ext2_ncol = extent2.ncol
    ext2_nrow = extent2.nrow

    if extent2.pixSize != extent2.pixSize:
        raise ValueError ("Pixel size of extents must be equal. ")
    else:
        pixSize = extent2.pixSize

    # get lower right coordinates
    ext1_lrX = ext1_ulX + ext1_ncol * pixSize
    ext1_lrY = ext1_ulY - ext1_nrow * pixSize
    ext2_lrX = ext2_ulX + ext2_ncol * pixSize
    ext2_lrY = ext2_ulY - ext2_nrow * pixSize

    # Get minimum upper left coordinates
    ulX = max(ext1_ulX, ext2_ulX)
    ulY = min(ext1_ulY, ext2_ulY)

    # Get minimum lower right coordinates
    lrX = min(ext1_lrX, ext2_lrX)
    lrY = max(ext1_lrY, ext2_lrY)

    # Get minimum number of rows and columns min(shapefile, rasterfiles)
    ncol = (lrX - ulX) / pixSize
    nrow = (ulY - lrY) / pixSize

    # Check if extent is valid
    if ncol <= 0 or nrow <= 0:
        return None
    else:
        return extent( ulX, ulY, int(ncol), int(nrow), pixSize, extent1.proj, lrX, lrY )

def getJointExtent(rasterPaths, shapePath="", attributeName=None, attributeValue=None, AOIextent=None, buffer=None):
    """ Get overlapping extent of input files. Input files must all have the same resolution. If calculating the joint extent
    of images of different resolution, use the image with the lowest resolution.

    :param: rasterPaths, shapePath, attributeName (optional), attributeValue (optional), extent (optional)

    :rtype: object: extent
    """

    if isinstance(rasterPaths, basestring):
        rasterPaths = rasterPaths.split()

    # Open raster files
    rasterFiles = []
    geoTransforms = []
    if len(rasterPaths) > 1:
        rasterFiles = [gdal.Open(rasterPath, GA_ReadOnly) for rasterPath in rasterPaths]
        geoTransforms = [raster.GetGeoTransform() for raster in rasterFiles ]
    else:
        rasterFiles.append(gdal.Open(rasterPaths[0], GA_ReadOnly))
        geoTransforms.append(rasterFiles[0].GetGeoTransform())

    # Metadata of raster file
    pixSize = geoTransforms[0][1]
    projRaster = osr.SpatialReference()
    projRaster.ImportFromWkt(rasterFiles[0].GetProjectionRef())

    # Get UL coordinates and number of cols/rows of shapfile/feature if given
    if os.path.isfile(shapePath) and shapePath.endswith(".shp"):
        shp_extent = convertShapeExtentToRasterExtent(shapePath, rasterPaths[0], attributeName, attributeValue, buffer)
        shp_ulX = shp_extent.ulX
        shp_ulY = shp_extent.ulY
        shp_ncol = shp_extent.ncol
        shp_nrow = shp_extent.nrow
        shp_lrX = shp_extent.lrX
        shp_lrY = shp_extent.lrY
    else:
        shp_ulX, shp_ulY, shp_ncol, shp_nrow, shp_lrX, shp_lrY = None, None, None, None, None, None

    # Get UL coordinates and number of cols/rows of extent if given
    if AOIextent != None:

        if not projRaster.IsSame(AOIextent.proj):
            ext_ulX, ext_ulY = reprojectPoint(AOIextent.proj, projRaster, AOIextent.ulX, AOIextent.ulY)
            ext_lrX, ext_lrY = reprojectPoint(AOIextent.proj, projRaster, AOIextent.lrX, AOIextent.lrY)
        else:
            ext_ulX = AOIextent.ulX
            ext_ulY = AOIextent.ulY
            ext_ncol = AOIextent.ncol
            ext_nrow = AOIextent.nrow
            ext_lrX = AOIextent.lrX
            ext_lrY = AOIextent.lrY
    else:
        ext_ulX, ext_ulY, ext_ncol, ext_nrow, ext_lrX, ext_lrY = None, None, None, None, None, None

    # Get minimum upper left coordinates
    ulX = max(filter(lambda x: x is not None, [geoTrans[0] for geoTrans in geoTransforms] + [shp_ulX, ext_ulX]))
    ulY = min(filter(lambda x: x is not None, [geoTrans[3] for geoTrans in geoTransforms] + [shp_ulY, ext_ulY]))

    lrX = min(filter(lambda x: x is not None, [img.GetGeoTransform()[0] + img.RasterXSize*pixSize for img in rasterFiles] + [shp_lrX, ext_lrX]))
    lrY = max(filter(lambda x: x is not None, [img.GetGeoTransform()[3] - img.RasterYSize*pixSize for img in rasterFiles] + [shp_lrY, ext_lrY]))

    # Get minimum number of rows and columns
    ncol = int( (lrX - ulX) / pixSize)
    nrow = int( (ulY - lrY) / pixSize)

    lrX = ulX + ncol * pixSize
    lrY = ulY - nrow * pixSize

    # Close raster files
    del rasterFiles, geoTransforms

    # Check if extent is valid
    if ncol <= 0 or nrow <= 0:
        return None
    else:
        return extent(ulX, ulY, int(ncol), int(nrow), pixSize, projRaster, lrX, lrY)

def convertShapeExtentToRasterExtent(shapePath, rasterPath, attributeName=None, attributeValue=None, buffer=None):
    """ Convert extent taken from a shapefile (feature or shapefile extent) to coordinates of a reference raster file.
    """

    # Meta data of raster file
    rasterF = gdal.Open(rasterPath)
    if rasterF is None:
        print "File not found."
        logging.warning("File not found.")
        return 1

    geoTrans = rasterF.GetGeoTransform()
    pixSize = geoTrans[1]
    projRaster = osr.SpatialReference()
    projRaster.ImportFromWkt(rasterF.GetProjectionRef())

    if os.path.isfile(shapePath) and shapePath.endswith(".shp"):

        # Open shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        source = driver.Open(shapePath, 0)
        shapeLay = source.GetLayer()
        projShape = shapeLay.GetSpatialRef()

        if attributeValue == None or attributeName == None:

            # Extract coordinates of upper left (ul) and lower right (lr) edges of shapefile extent
            shp_ulX, shp_lrX, shp_lrY, shp_ulY = shapeLay.GetExtent()

            if projShape.IsSame(projRaster) == False:
                shp_ulX, shp_ulY = reprojectPoint(projShape, projRaster, shp_ulX, shp_ulY)
                shp_lrX, shp_lrY = reprojectPoint(projShape, projRaster, shp_lrX, shp_lrY)

        else:

            # Filter features
            query =  str(attributeName) + " = '" + str(attributeValue)+ "'" # "Location = 'E47N50' "
            shapeLay.SetAttributeFilter(query)

            # Get feature and geometry
            feature = shapeLay.GetNextFeature()
            if feature is None:
                logging.critical("AOI feature not found.")
                print("AOI feature not found.")
                sys.exit(1)

            geom = feature.GetGeometryRef()

            # reproject the geometry if necessary
            if projShape.IsSame(projRaster) == False:
                coordTrans = osr.CoordinateTransformation(projShape, projRaster)
                geom.Transform(coordTrans)

            # Extract coordinates of upper left (ul) and lower right (lr) edges of feature geometry extent
            shp_ulX, shp_lrX, shp_lrY, shp_ulY = geom.GetEnvelope()

        if buffer != None:
            shp_ulX -= pixSize * buffer
            shp_ulY += pixSize * buffer
            shp_lrX += pixSize * buffer
            shp_lrY -= pixSize * buffer

        # Get raster coordinates of ul and lr corners of feature
        ras_ulX, ras_ulY = world2Pixel(geoTrans, shp_ulX, shp_ulY)
        ras_lrX, ras_lrY = world2Pixel(geoTrans, shp_lrX, shp_lrY)

        # Get geographic coordinates of center of upper left and lower right raster cell
        ulX = (ras_ulX * pixSize) + geoTrans[0]
        ulY = (ras_ulY * (-1) * pixSize) + geoTrans[3]
        lrX = (ras_lrX * pixSize) + geoTrans[0]
        lrY = (ras_lrY * (-1) * pixSize) + geoTrans[3]

        # Get number of cols and rows
        ncol = int(ras_lrX - ras_ulX)
        nrow = int(ras_lrY - ras_ulY)

        # Check if extent is valid
        if ncol <= 0 or nrow <= 0:
            # todo: raise ValueError, delete return None
            return None
        else:
            return extent(ulX, ulY, int(ncol), int(nrow), pixSize, projRaster, lrX, lrY)

def getSubExtents(mainExtent, step):

    ulX = mainExtent.ulX
    startY = mainExtent.ulY
    ncols = mainExtent.ncol
    nrows = mainExtent.nrow
    lrX = mainExtent.lrX
    lrY = mainExtent.lrY
    pixSize = mainExtent.pixSize

    if step <= 0:
        print "Error: Step size must be positive."
        return

    subextents =[]
    if nrows <= step:
        subextents.append(mainExtent)
        return subextents

    endY = startY + (-1) * pixSize * nrows
    stepY = (-1) * pixSize * step

    # Define size of data chunks (number of lines to process within each iteration)
    lolim = np.arange(startY, endY, stepY)
    uplim = np.arange(startY+stepY, endY, stepY)
    nLoops = len(uplim)

    # Get coordinatsubextents
    for i in range(0, nLoops):
        subextent = extent(ulX, lolim[i], ncols, step,  pixSize, mainExtent.proj, lrX, uplim[i])
        subextents.append(subextent)

    # Last subextent
    if uplim[nLoops-1] != endY:
        lastStep = int((endY - uplim[nLoops-1]) / ((-1)*pixSize))
        subextent = extent(ulX, uplim[nLoops-1], ncols, lastStep, pixSize,mainExtent.proj, lrX, endY)
        subextents.append(subextent)

    return subextents

def transformProjections(scenes):
    # Check if raster files have same projection
    projs = [sce.proj for sce in scenes]
    uniqueProjs = set(projs)
    if len(uniqueProjs) == 1:
        print "All files share same projection. Nothing to do."
        return 0

    # Get one projection as reference
    refProj = osr.SpatialReference()
    refProj.ImportFromWkt(next(iter(uniqueProjs)))
    refProj_wkt = refProj.ExportToWkt()

    # Convert all rater files (change geotrans and proj attributes)
    for sce in scenes:
        targetProj = osr.SpatialReference()
        targetProj.ImportFromWkt(sce.proj)

        if not targetProj.IsSame(refProj):
            sce.proj = refProj_wkt
            x,y = sce.geotrans[0], sce.geotrans[3]

            newX, newY = reprojectPoint(targetProj, refProj, x, y)

            sce.geotrans[0] = newX
            sce.geotrans[3] = newY

    # Check if all projections are the same now
    projs = [sce.proj for sce in scenes]
    uniqueProjs = set(projs)
    if len(uniqueProjs) == 1:
        print "Transformation successful!"
        return 0

# Miscellaneous

def createVRT(inFiles, outfile, AOIextent=None):
    """Create a vrt file of several input files for specified extent. """
    if AOIextent != None:
        cmd = "gdalbuildvrt -separate -te " + str(AOIextent.ulX) + " " + str(AOIextent.lrY) + " " + str(AOIextent.lrX) + " " + str(AOIextent.ulY) + " " + outfile + " " + " ".join(inFiles)
    else:
        cmd = "gdalbuildvrt -separate " + outfile + " " + " ".join(inFiles)

    os.system(cmd)

    # Choose masked band
    # outfile2 = outfile[:-4] + "_masked.vrt"
    # idx_mask = len(inFiles)
    # cmd = "gdal_translate -of vrt " + outfile + " -mask " + str(idx_mask) + " " + outfile2 + " --config GDAL_TIFF_INTERNAL_MASK YES"
    # os.system(cmd)

    return (outfile)

def downsample(a, a_pixSize, out_pixSize):
    """Resample image to lower resolution

        :param a array
        :param a_pixSize Pixel size of source array
        :param out_pixSize Pixel size of output array

        :rtype ndarray
    """
    shape = (int(a.shape[0] * a_pixSize / out_pixSize), int(a.shape[1] * a_pixSize / out_pixSize))

    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def calcMetrics(stack):
    """Compute metrics of numpy stack

    :param: stack

    :return: metrics
    """
    metrics = []
    # Mean
    metrics.append(bn.nanmean(stack, axis=0))
    # Median
    metrics.append(bn.nanmedian(stack, axis=0))
    # Standard Deviation
    metrics.append(bn.nanstd(stack, axis=0))
    # Sum
    metrics.append(bn.nansum(stack, axis=0))
    # Maximum and point of time
    metrics.append(bn.nanmax(stack, axis=0))
    # TODO: fix argmin und argmax -> value error when all-nan-slice encountered
    #try:
    #    metrics.append(bn.nanargmax(stack, axis=0))
    #except:
    #    pass
    # Minimum and point of time
    metrics.append(bn.nanmin(stack, axis=0))
    #try:
    #    metrics.append(bn.nanargmin(stack, axis=0))
    #except:
    #    pass
    # Percentile
    allNaN = bn.allnan(stack, axis=0)
    percentiles = nan_percentile(np.array(stack, copy=True), [25,5,95])
    metrics.append(np.where(allNaN, np.nan, percentiles[0]))
    metrics.append(np.where(allNaN, np.nan, percentiles[1]))
    metrics.append(np.where(allNaN, np.nan, percentiles[2]))

    return metrics

# Percentile calculation

def nan_percentile(arr, q):
    """ Function to calcualte percentil of numpy stack faster than using numpy.percentile

        taken from http://krstn.eu/np.nanpercentile()-there-has-to-be-a-faster-way/ and
        http://stackoverflow.com/questions/32089973/numpy-index-3d-array-with-index-of-last-axis-stored-in-2d-array

    """

    # valid (non NaN) observations along the first axis
    no_obs = bn.allnan(arr, axis=0)
    valid_obs = np.sum(np.isfinite(arr), axis=0)
    # replace NaN with maximum
    max_val = np.nanmax(arr)
    arr[np.isnan(arr)] = max_val
    # sort - former NaNs will move to the end
    arr = np.sort(arr, axis=0)

    # loop over requested quantiles
    if type(q) is list:
        qs = []
        qs.extend(q)
    else:
        qs = [q]
    if len(qs) < 2:
        quant_arr = np.zeros(shape=(arr.shape[1], arr.shape[2]))
    else:
        quant_arr = np.zeros(shape=(len(qs), arr.shape[1], arr.shape[2]))

    result = []
    for i in range(len(qs)):
        quant = qs[i]
        # desired position as well as floor and ceiling of it
        k_arr = (valid_obs - 1) * (quant / 100.0)
        k_arr = np.where(k_arr < 0, 0, k_arr)
        f_arr = np.floor(k_arr).astype(np.int32)
        c_arr = np.ceil(k_arr).astype(np.int32)
        fc_equal_k_mask = (f_arr == c_arr)

        # linear interpolation (like numpy percentile) takes the fractional part of desired position
        floor_val = _zvalue_from_index_ogrid(arr=arr, ind=f_arr) * (c_arr - k_arr)
        ceil_val = _zvalue_from_index_ogrid(arr=arr, ind=c_arr) * (k_arr - f_arr)

        # If choose does not work, use _zvalue_from_index (not working however)
        #floor_val = f_arr.choose(arr) * (c_arr - k_arr)
        #ceil_val = c_arr.choose(arr) * (k_arr - f_arr)

        quant_arr = floor_val + ceil_val
        quant_arr[fc_equal_k_mask] = _zvalue_from_index_ogrid(arr=arr, ind=k_arr.astype(np.int32))[fc_equal_k_mask]  # if floor == ceiling take floor value
        # quant_arr[fc_equal_k_mask] = k_arr.astype(np.int32).choose(arr)[fc_equal_k_mask]
        quant_arr = np.where(no_obs, np.nan, quant_arr)

        result.append(quant_arr)

    return result

def _zvalue_from_index_ogrid(arr, ind):
    y = arr.shape[1]
    x = arr.shape[2]
    y,x =np.ogrid[0:y, 0:x]

    return arr[ind, y, x]

def indexByPercentile(arr, quant):
    """Calculate percentile of numpy stack and return the index of the chosen pixel. """
    #valid (non NaN) observations along the first axis
    arr_tmp = np.array(arr, copy=True)

    #no_obs = bn.allnan(arr_tmp["data"], axis=0)
    valid_obs = np.sum(np.isfinite(arr_tmp["data"]), axis=0)
    # replace NaN with maximum
    max_val = np.nanmax(arr_tmp["data"]) + 1
    arr_tmp["data"][np.isnan(arr_tmp["data"])] = max_val
    # sort - former NaNs will move to the end
    arr_tmp = np.sort(arr_tmp, kind="mergesort", order="data", axis=0)
    arr_tmp["data"] = np.where(arr_tmp["data"]==max_val, np.nan, arr_tmp["data"])

    # desired position as well as floor and ceiling of it
    k_arr = (valid_obs - 1) * (quant / 100.0)
    k_arr = np.where(k_arr < 0, 0, k_arr)
    f_arr = np.floor(k_arr + 0.5)
    f_arr = f_arr.astype(np.int)

    # get floor value of reference band and index band
    floor_val = _zvalue_from_index_ogrid(arr=arr_tmp, ind=f_arr.astype("int16"))
    idx = np.where(valid_obs==0, 255, floor_val["ID"])

    quant_arr = floor_val["data"]

    del arr_tmp

    return (quant_arr, idx, valid_obs)

def index_interpolation(arr, quant):
    #valid (non NaN) observations along the first axis
    arr_tmp = np.array(arr, copy=True)

    no_obs = bn.allnan(arr_tmp["data"], axis=0)
    valid_obs = np.sum(np.isfinite(arr_tmp["data"]), axis=0)
    # replace NaN with maximum
    max_val = np.nanmax(arr_tmp["data"]) + 1
    arr_tmp["data"][np.isnan(arr_tmp["data"])] = max_val
    # sort - former NaNs will move to the end
    arr_tmp = np.sort(arr_tmp, kind="mergesort", order="data", axis=0)
    arr_tmp["data"] = np.where(arr_tmp["data"]==max_val, np.nan, arr_tmp["data"])

    # desired position as well as floor and ceiling of it
    k_arr = (valid_obs - 1) * (quant / 100.0)
    k_arr = np.where(k_arr < 0, 0, k_arr)
    f_arr = np.floor(k_arr).astype(np.int32)
    c_arr = np.ceil(k_arr).astype(np.int32)

    #z_arr = np.where(k_arr%1 >= 0.5, c_arr, f_arr)

    # get floor value of reference band and index band
    floor_val = _zvalue_from_index_ogrid(arr=arr_tmp, ind=f_arr)
    f_idx = np.where(valid_obs==0, 255, floor_val["ID"])
    ceil_val = _zvalue_from_index_ogrid(arr=arr_tmp, ind=c_arr)
    c_idx = np.where(valid_obs==0, 255, ceil_val["ID"])
    k_idx = k_arr%1

    del arr_tmp

    return (f_idx, c_idx, k_idx, valid_obs)

def get_interpolatedIndexValue(arr, f_idx, c_idx, k_idx):

    # linear interpolation (like numpy percentile) takes the fractional part of desired position
    floor_val = np.nanmax(np.where(arr["ID"] == f_idx, arr["data"], np.nan), axis=0) * (1-k_idx)
    ceil_val = np.nanmax(np.where(arr["ID"] == c_idx, arr["data"], np.nan), axis=0) * k_idx

    quant_arr = floor_val + ceil_val

    return quant_arr

def downsample(a, a_pixSize, out_pixSize):
    """Resample image to lower resolution

        :param a array
        :param a_pixSize Pixel size of source array
        :param out_pixSize Pixel size of output array

        :rtype ndarray
    """
    shape = (int(a.shape[0] * a_pixSize / out_pixSize), int(a.shape[1] * a_pixSize / out_pixSize))

    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]

    return a.reshape(sh).mean(-1).mean(1)

def rescale(index):
    # Take 1st and 99th percentile for clipping the range
    minVal = np.nanpercentile(index, 1) #rsu.nan_percentile(index, 1)
    maxVal = np.nanpercentile(index, 99) # rsu.nan_percentile(index, 99)
    # Replace min/max values
    index = np.where((index > maxVal), maxVal, index)
    index = np.where((index < minVal), minVal, index)
    # Normalization
    index_scaled = (index - minVal) / (maxVal - minVal) * 2000. - 1000.

    return index_scaled
