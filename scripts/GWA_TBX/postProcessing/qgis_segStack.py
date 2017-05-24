#Definition of inputs and outputs
#==================================
##Generic Tools=group
##PreparePostProcessingStack=name
##ParameterRaster|inRst1|Input Stack|False
##ParameterBoolean|B1|Band 1|False
##ParameterBoolean|B2|Band 2|False
##ParameterBoolean|B3|Band 3|False
##ParameterBoolean|B4|Band 4|False
##ParameterBoolean|B5|Band 5|False
##ParameterBoolean|B6|Band 6|False
##ParameterBoolean|B7|Band 7|False
##ParameterBoolean|B8|Band 8|False
##ParameterBoolean|B9|Band 9|False
##ParameterBoolean|B10|Band 10|False

##ParameterRaster|inRst2|Classification Raster|False

##OutputRaster|outPath|Output Stack for Segmentation

#==================================
# import modules
import os
import sys
from processing.tools import dataobjects

here = os.path.dirname(scriptDescriptionFile)
if here not in sys.path:
    sys.path.append(here)
import segStack


outBands = []
if B1 == True:
    outBands.append(1)
if B2 == True:
    outBands.append(2)
if B3 == True:
    outBands.append(3)
if B4 == True:
    outBands.append(4)
if B5 == True:
    bandList.append(5)
if B6 == True:
    outBands.append(6)
if B7 == True:
    outBands.append(7)
if B8 == True:
    outBands.append(8)
if B9 == True:
    outBands.append(9)
if B10 == True:
    outBands.append(10)
if len(outBands) == 0:
    raise GeoAlgorithmExecutionException("ERROR: No bands selected")

progress.setConsoleInfo('Starting preparation of stack...')

# run script
segStack.prepSegStack(inRst1, outBands, inRst2, outPath)

progress.setConsoleInfo('Finished!')

#==================================
