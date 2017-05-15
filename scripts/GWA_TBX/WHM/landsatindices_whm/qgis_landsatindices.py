#Definition of inputs and outputs
#==================================
##Landsat 8 Indices=name
##Landsat Tools=group
##ParameterRaster|input|Input Reflectance Stack|False
##OutputDirectory|outputDirectory|Folder to save the stack of Indices
#OutputRaster|output|Name for Index Stack

# Call the function for Landsat 8 index calculation
#==================================
import os
import sys
here = os.path.dirname(scriptDescriptionFile)
if here not in sys.path:
    sys.path.append(here)
import landsatindices

print 'Starting index calculation...'
landsatindices.landsatindices(input, outputDirectory)
print 'Finished writing to disk...'
