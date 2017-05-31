# Copyright (c) 2016, GeoVille Information Systems GmbH
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, is prohibited for all commercial applications without
# licensing by GeoVille GmbH.

"""
Calculation of topographic Wetness Index using SAGA

Keyword arguments:
    TWI: TWI file

Output:
    TWI binary mask

Dependencies:

Date created: 09/06/2016
Date last modified: 21/09/2016

"""

__author__ = "Christina Ludwig"
__version__ = "v1"

##Create TWI binary mask=name
##Wetland Inventory=group
##ParameterRaster|path_TWI|TWI
#OutputDirectory|out_dir|Output directory


import os, sys
import numpy as np
from scipy import ndimage
from osgeo import gdal

# Load additional library
here = r"C:\Users\ludwig\.qgis2\processing\scripts\GWA_TBX\WI_WCR"

here = os.path.dirname(scriptDescriptionFile)
pyDir = os.path.join(here, 'data', 'python')
if pyDir not in sys.path:
    sys.path.append(pyDir)


import RSutils.RSutils as rsu


#path_TWI = "I:/temp/TWI/TWI.tif"

out_dir = os.path.dirname(path_TWI)
print out_dir

attrName = None
attrVal = None
extentBuffer=0

# Create TWI binary mask ------------------------------------------------

outfile_binary = os.path.join(out_dir, "TWIbinary.tif")
if not os.path.exists(outfile_binary):

    print "Creating binary mask ..."

    # Median filter
    #TWI, geotrans, proj = rsu.raster2array(path_TWI)
    #TWI_smooth = ndimage.median_filter(TWI, size=7) * 1000
    #dest = os.path.join(out_dir, "TWI.tif")
    #rsu.array2raster(TWI_smooth, geotrans, proj, dest, gdal.GDT_Int16, -9999)

    # Read TWI file
    TWI_file = gdal.Open(path_TWI)
    TWI = TWI_file.GetRasterBand(1).ReadAsArray()
    nodata_TWI = TWI_file.GetRasterBand(1).GetNoDataValue()
    TWI = np.where(TWI == nodata_TWI, np.nan, TWI)

    # Find threshold
    TWI = (TWI-np.nanmin(TWI))/(np.nanmax(TWI)-np.nanmin(TWI))
    TWI_binary = np.where(TWI > 0.4, 0, 1)

    #Sieve mask
    TWI_binary = rsu.removeNoise(TWI_binary, 2)
    TWI_binary = rsu.binaryBuffer_negative(TWI_binary, 3)

    #Export to file
    x = rsu.array2raster(TWI_binary, TWI_file.GetGeoTransform(), TWI_file.GetProjection(), outfile_binary, gdal.GDT_Byte, 255)
    if x == False:
        print("Exporting TWI mask failed.")
        raise RuntimeError("Exporting TWI mask failed.")
    
    del TWI_file, TWI_binary, TWI

