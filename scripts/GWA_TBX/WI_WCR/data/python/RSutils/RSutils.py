# Copyright (c) 2016, GeoVille Information Systems GmbH
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.

"""
Support functions and classes for raster handling

Date created: 09/06/2016
Date last modified: 09/05/2017

"""

__author__ = "Christina Ludwig"
__version__ = "v1"

import glob
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
import zipfile
import xml.etree.ElementTree as eTree
import warnings

# CLASSES ===============================================================

class Scene(object):
    """Basic Object for satellite scene."""
    def __init__(self, sceneDir, ID=None, tempDir=None, extentAOI=None):

        # Check scene directory
        if os.path.isdir(sceneDir):
            self.zipped = False
            self.dir = sceneDir
        elif sceneDir.endswith(".zip"):
            self.zipped = True
            self.dir = sceneDir
        else:
            raise RuntimeError("Input directory is invalid.")

        # Temp directory
        if tempDir is None:
            self.tempDir = self.dir
        else:
            self.tempDir = tempDir

        self.extent = extentAOI

        self.cloudCoverage = None

    def getGeotrans(self):

        # Get extent
        self.extent = getJointExtent(self.files[-1], AOIextent=self.extent)
        if self.extent is None:
            raise RuntimeError("Scene not within AOI.")

        # Get metadata
        metadata = gdal.Open(self.files[-1], GA_ReadOnly)
        self.geotrans = list(metadata.GetGeoTransform())
        self.geotrans[0] = self.extent.ulX
        self.geotrans[3] = self.extent.ulY
        self.pixSize = self.geotrans[1]
        self.proj = metadata.GetProjection()
        self.nodata = metadata.GetRasterBand(1).GetNoDataValue()
        del metadata

    def getBand(self, bandNo=1, masked=True, extent=None):
        """Get specified band of scene as numpy array."""

        if extent == None:
            extent = self.extent

        band = raster2array(self.files[bandNo-1], self.extent)[0]

        if band is None:
            raise ValueError("Requested band %s for scene %s could not be opened." % (bandNo, self.ID))

        if masked:
            mask = self.getMask(extent=extent)
            if mask is not None:
                band = np.where((mask == 1) | (band == self.nodata), np.nan, band)

        return band

    def getMask(self, extent=None):
        """Get cloud mask of scene as numpy array."""

        if os.path.exists(self.fmask):
            fmask = raster2array(self.fmask, extent)[0]
        else:
            return None

        if extent is None:
            extent = self.extent

        # Mask out 2 - clouds, 3 - shadow and (4 - snow --> Not for Africa)

        mask = np.where((fmask == 2) | (fmask == 3), 1, 0)
        mask = np.where(np.isnan(fmask), np.nan, mask)

        del fmask

        # Buffer masks
        mask = binaryBuffer(mask, size=1)
        mask = removeNoise_oneSide(mask)

        return mask

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

    def updateExtent(self, extent):

        intersectedExtent = intersectRasterExtents(self.extent, extent)

        if intersectedExtent is None:
            raise RuntimeWarning("Given extent does not overlap with scene.")

        self.extent = intersectedExtent

        # Get metadata
        self.geotrans[0] = self.extent.ulX
        self.geotrans[3] = self.extent.ulY

    def calcCloudCoverage(self, AOIextent=None):

        if self.cloudCoverage is None:

            if not os.path.exists(self.fmask):
                return 1

            fmask = None
            fmask = raster2array(self.fmask, self.extent)[0]

            nonNANpixels = np.nansum(np.where(np.isnan(fmask), 0, 1))

            # Mask out 2 - clouds, 3 - shadow and 4 - snow
            mask = np.where((fmask == 2) | (fmask == 3) | (fmask == 4), 1, 0)

            # Compute cloud coverate
            # self.cloudCoverage = 100. - (float(bn.nansum(mask == 0)) / (float(self.extent.ncol) * float(self.extent.nrow))) * 100.
            self.cloudCoverage = 100. - (float(bn.nansum(mask == 0)) / nonNANpixels) * 100.

class SentinelScene(Scene):
    """Class to handle Sentinel-2 scenes 
    
    It does not include the 60 meter bands 1,9 and 10.
    
    """

    def __init__(self, sceneDir, ID=None, tempDir=None, extentAOI=None):

        Scene.__init__(self, sceneDir, tempDir, extentAOI)

        # Scene description
        self.ID = ID
        if ID is not None:
            self.tileID = self.ID[-6:]
        else:
            self.tileID = None
        self.getMetadata()
        self.getFiles()
        self.extent = extentAOI
        self.getGeotrans()

        self.nodata = 0.0

    def getMetadata(self):

        if self.dir.endswith(".zip"):
            zpfile = zipfile.ZipFile(self.dir)
            metadatafile = fnmatch.filter(zpfile.namelist(), "*/GRANULE/*/MTD_TL.xml")

            if len(metadatafile) != 1 and self.tileID is not None:
                metadatafile = fnmatch.filter(zpfile.namelist(), "*/GRANULE/*/*%s.xml" % self.tileID)

            if len(metadatafile) != 1:
                raise RuntimeError("No metadata file found!")

            # Open metadata file
            metaF = zpfile.open(metadatafile[0])
            metainfo = eTree.parse(metaF)
            rootElement = metainfo.getroot()

            # Get IDs
            self.ID = [child.text for child in rootElement.iter("TILE_ID")][0]
            tileIDindex = self.ID.rfind("_T")+1
            self.tileID = self.ID[tileIDindex:tileIDindex+6]

            # Get date
            datestring = [child.text for child in rootElement.iter("SENSING_TIME")][0]
            self.date = dt.strptime(datestring[:10], "%Y-%m-%d")

            metaF.close()
            zpfile.close()

        else:
            metadatafile = ""
            for root, dirs, files in os.walk(self.dir):
                for f in files:
                    if f == "MTD_TL.xml" or (self.tileID is not None and f.endswith("*"+self.tileID+".xml")):
                        metadatafile = os.path.join(root, f)

            if metadatafile == "":
                raise RuntimeError("No metadata file found!")

            metaF = open(metadatafile)
            metainfo = eTree.parse(metaF)
            rootElement = metainfo.getroot()

            datestring = [child.text for child in rootElement.iter("SENSING_TIME")][0]
            self.date = dt.strptime(datestring[:10], "%Y-%m-%d")

            # Get IDs
            if self.ID is None:
                self.ID = [child.text for child in rootElement.iter("TILE_ID")][0]
                tileIDindex = self.ID.rfind("_T")+1
                self.tileID = self.ID[tileIDindex:tileIDindex+6]

            metaF.close()

    def getFiles(self):
        """ 
        BandNo 1: BLUE (Sentinel Band number: 2)
        BandNo 2: GREEN (3)
        BandNo 3: RED (4)
        BandNo 4: Red Edge (5)
        BandNo 5: Red Edge (6)
        BandNo 6: Red Edge (7)
        BandNo 7: NIR (8)
        BandNo 8: Red Edge 8A (8A)
        BandNo 9: SWIR1 (11)
        BandNo 10: SWIR2 (12)
        """

        bandEndings = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
        optionalBands = ["B02", "B05", "B06", "B07", "B8A"]
        self.files = []

        if self.dir.endswith(".zip"):
            zpfile = zipfile.ZipFile(self.dir)

            for b in bandEndings:
                if self.tileID is None:
                    band = ["/".join(["/vsizip", self.dir, b]) for b in fnmatch.filter(zpfile.namelist(), "*%s*.[Jj][Pp]2" % b)]
                else:
                    band = ["/".join(["/vsizip", self.dir, b]) for b in fnmatch.filter(zpfile.namelist(), "*%s*%s*.[Jj][Pp]2" % (self.tileID, b))]
                if len(band) == 0:
                    band = None
                    if b in optionalBands:
                        warnings.warn("Band %s is missing" % b)
                    else:
                        raise IOError("Band %s is missing" % b)
                else:
                    band = band[0]
                self.files.append(band)

            self.fmask = ""

            zpfile.close()
        else:

            for b in bandEndings:
                band = glob.glob(self.dir + "*/GRANULE/*/IMG_DATA/*%s.[Jj][Pp]2"%b)
                if len(band) == 0:
                    band = glob.glob(self.dir + "/*/GRANULE/*/IMG_DATA/*%s.[Jj][Pp]2" % b)
                    if len(band) == 0:
                        band = None
                        if b in optionalBands:
                            warnings.warn("Band %s is missing" % b)
                        else:
                            raise IOError("Band %s is missing" % b)
                else:
                    band = band[0]
                self.files.append(band)

            #self.files = glob.glob(self.dir + "/*/GRANULE/*/IMG_DATA/*B0[2345678].[Jj][Pp]2")
            #self.files += glob.glob(self.dir + "/*/GRANULE/*/IMG_DATA/*B8A.[Jj][Pp]2")
            #self.files += glob.glob(self.dir + "/*/GRANULE/*/IMG_DATA/*B1[12].[Jj][Pp]2")
            #self.files.sort()

            # Mask files
            fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_fmask.[Tt][Ii][Ff]")]
            if len(fmask) != 0:
                self.fmask = fmask[0]
            else:
                fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_fmask.img")]
                if len(fmask) == 0:
                    self.fmask = ""
                else:
                    self.fmask = fmask[0]

    def getBand(self, bandNo=1, masked=True, extent=None):
        """Get specified band of scene as numpy array.
        
        BandNo 1: BLUE (Sentinel Band number: 2)
        BandNo 2: GREEN (3)
        BandNo 3: RED (4)
        BandNo 4: Red Edge (5)
        BandNo 5: Red Edge (6)
        BandNo 6: Red Edge (7)
        BandNo 7: NIR (8)
        BandNo 8: Red Edge 8A (8A)
        BandNo 9: SWIR1 (11)
        BandNo 10: SWIR2 (12)
        """

        if extent == None:
            extent = self.extent

        band = raster2array(self.files[bandNo-1], extent)[0]

        if band is None:
            raise ValueError("Requested band %s for scene %s could not be opened." % (bandNo, self.ID))

        if masked:
            mask = self.getMask(extent=extent)
            if mask is not None:
                band = np.where((mask == 1) | (band == self.nodata), np.nan, band)

        return band

class LandsatScene(Scene):
    """Class to handle Landsat scenes
    
    This class leaves out the thermal and coastal blue bands (Landsat 8) and
        
    """

    def __init__(self, sceneDir, ID=None, tempDir=None, extentAOI=None):

        Scene.__init__(self, sceneDir, tempDir, extentAOI)

        # Get files of spectral bands
        self.ID = os.path.basename(sceneDir)[:16]
        self.tileID = self.ID[3:9]
        self.date = dt.strptime(self.ID[9:16], "%Y%j")
        self.sensor = self.ID[2]

        self.getFiles()
        self.getGeotrans()

    def getFiles(self):
        """
        BandNo 1: BLUE (LS8: 2, LS5/7: 1)
        BandNo 2: GREEN (LS8: 3, LS5/7: 2)
        BandNo 3: RED (LS8: 4, LS5/7: 3)
        BandNo 4: NIR (LS8: 5, LS5/7: 4)
        BandNo 5: SWIR1 (LS8: 6, LS5/7: 5)
        BandNo 6: SWIR2 (LS8: 7, LS5/7: 7)
        """

        # Files
        self.files = []

        if self.sensor == "8":
            bandEndings = ["B2", "B3", "B4", "B5", "B6", "B7"]
            optionalBands = ["B2"]
            for b in bandEndings:
                band = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*%s*.[Tt][Ii][Ff]" % b)]
                if len(band) == 0:
                    band = None
                    if b in optionalBands:
                        warnings.warn("Band %s is missing" % b)
                    else:
                        raise IOError("Band %s is missing" % b)
                else:
                    band = band[0]
                self.files.append(band)

        else:
            bandEndings = ["B1", "B2", "B3", "B4", "B5", "B7"]
            optionalBands = ["B1"]
            for b in bandEndings:
                band = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*%s*.[Tt][Ii][Ff]" % b)]
                if len(band) == 0:
                    band = None
                    if b in optionalBands:
                        warnings.warn("Band %s is missing" % b)
                    else:
                        raise IOError("Band %s is missing" % b)
                else:
                    band = band[0]
                self.files.append(band)

        # Mask files
        fmask = [os.path.join(self.dir, b) for b in fnmatch.filter(os.listdir(self.dir), "*_fmask.[Tt][Ii][Ff]")]
        if len(fmask) != 0:
            self.fmask = fmask[0]
        else:
            self.fmask = ""

        # shadowMask = [os.path.join(self.maskDir, b) for b in fnmatch.filter(os.listdir(self.maskDir), "*_shadowmask.tif")]
        # if len(shadowMask) != 0:
        #     self.shadowMask = shadowMask[0]
        # else:
        #     self.shadowMask = ""

    def getBand(self, bandNo=1, masked=True, extent=None, transformed=False ):
        """Get specified band of scene as numpy array.
        
        BandNo 1: BLUE (LS8: 2, LS5/7: 1)
        BandNo 2: GREEN (LS8: 3, LS5/7: 2)
        BandNo 3: RED (LS8: 4, LS5/7: 3)
        BandNo 4: NIR (LS8: 5, LS5/7: 4)
        BandNo 5: SWIR1 (LS8: 6, LS5/7: 5)
        BandNo 6: SWIR2 (LS8: 7, LS5/7: 7)
        
        """

        if extent == None:
            extent = self.extent

        band = raster2array(self.files[bandNo-1], extent)[0]

        if band is None:
            return None

        if masked:
            mask = self.getMask(extent=extent)
            band = np.where((mask == 1) | (band == self.nodata) | (band == 20000), np.nan, band)

        return band

class extent(object):

    def __init__(self, ulX, ulY, ncol=None, nrow=None, pixSize=None, proj=None, lrX=None, lrY=None):
        self.ulX = ulX
        self.ulY = ulY
        self.ncol = ncol
        self.nrow = nrow
        self.pixSize = pixSize
        self.proj = proj
        self.lrX = lrX
        self.lrY = lrY

        # if ncol is None:
        #     self.ncol = int((self.lrX - self.ulX) / self.pixSize)
        #
        # if nrow is None:
        #     self.nrow = int((self.ulY - self.lrY) / self.pixSize)
        #
        # if lrX is None:
        #     self.lrX = self.ulX + self.pixSize * self.ncol
        #
        # if lrY is None:
        #     self.lrY = self.ulY - self.pixSize * self.nrow

    def getGeotrans(self):
        return [self.ulX, self.pixSize, 0.0, self.ulY, 0.0, -self.pixSize]

    def getProj(self):
        return self.proj.ExportToWkt()

    def convertToRasterExtent(self, refScene):

        rasterProj = osr.SpatialReference()
        rasterProj.ImportFromWkt(refScene.proj)

        if self.proj.IsSame(rasterProj) == False:
            self.ulX, self.ulY = reprojectPoint( self.proj, rasterProj, self.ulX, self.ulY)
            self.lrX, self.lrY = reprojectPoint( self.proj,rasterProj, self.lrX, self.lrY)

        geoTrans = refScene.geotrans

        # Get raster coordinates of ul and lr corners of feature
        ras_ulX, ras_ulY = world2Pixel(geoTrans,self.ulX, self.ulY)
        ras_lrX, ras_lrY = world2Pixel(geoTrans, self.lrX,self.lrY)

        # Get geographic coordinates of center of upper left and lower right raster cell
        self.ulX = (ras_ulX * refScene.pixSize) + geoTrans[0]
        self.ulY = (ras_ulY * (-1) * refScene.pixSize) + geoTrans[3]
        self.lrX = (ras_lrX * refScene.pixSize) + geoTrans[0]
        self.lrY = (ras_lrY * (-1) * refScene.pixSize) + geoTrans[3]

        # Get number of cols and rows
        self.ncol = int(ras_lrX - ras_ulX)
        self.nrow = int(ras_lrY - ras_ulY)

        self.proj = rasterProj
        self.pixSize = refScene.pixSize

        # Check if extent is valid
        if self.ncol <= 0 or self.nrow <= 0:
            print("Extent is invalid!")

# FUNCTIONS =============================================================

def filter_scenes_by_date(scenes, start_date, end_date):

    filteredScenes = []
    for sce in scenes:
        if start_date <= sce.date <= end_date:
            filteredScenes.append(sce)
    return filteredScenes

def search_scene_directories(inDir, searchPattern):
    """ Search for scene directores in input direcotry. Duplicate folders are removed. No recursive search."""

    sceneDirs = []
    for f in os.listdir(inDir):
        if fnmatch.fnmatch(f, searchPattern):
            sceneDirs.append(os.path.join(inDir, f))
    sceneIDs = set([os.path.basename(f).strip(".zip") for f in sceneDirs])

    # Filter out duplicates
    duplicates = []
    for id in sceneIDs:
        dirs = fnmatch.filter(sceneDirs, "*" + id + "*")
        if len(dirs) > 1:
            duplicates += dirs[1:]

    sceneDirs = [d for d in sceneDirs if d not in duplicates]

    return sceneDirs

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

def createVRT(inFiles, outfile, AOIextent=None, separate=True):
    """Create a vrt file of several input files for specified extent. """

    if separate:
        if AOIextent != None:
            cmd = "gdalbuildvrt -separate -te " + str(AOIextent.ulX) + " " + str(AOIextent.lrY) + " " + str(AOIextent.lrX) + " " + str(AOIextent.ulY) + " " + outfile + " " + " ".join(inFiles)
        else:
            cmd = "gdalbuildvrt -separate " + outfile + " " + " ".join(inFiles)
    else:
        if AOIextent != None:
            cmd = "gdalbuildvrt -te " + str(AOIextent.ulX) + " " + str(AOIextent.lrY) + " " + str(AOIextent.lrX) + " " + str(AOIextent.ulY) + " " + outfile + " " + " ".join(inFiles)
        else:
            cmd = "gdalbuildvrt " + outfile + " " + " ".join(inFiles)

    subprocess.call(cmd, shell=True)  # universal_newlines=True, stderr=subprocess.STDOUT,

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
