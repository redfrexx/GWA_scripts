#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_03_OLCI_MphChl=name
##ParameterString|MPHvalidPixelExpression|Valid pixel expression|not quality_flags.invalid and (not pixel_classif_flags.IDEPIX_LAND or quality_flags.fresh_inland_water) and not (pixel_classif_flags.IDEPIX_CLOUD or pixel_classif_flags.IDEPIX_CLOUD_BUFFER)
##ParameterNumber|MPHcyanoMaxValue|Cyano maximum value|0|None|1000.0
##ParameterNumber|MPHchlThreshForFloatFlag|CHL threshold for float flag|0|None|350.0
##ParameterBoolean|MPHexportMph|Export MPH|false


import os
import glob
import tempfile

tempfolder = 'wq_scripts_'

def folder_check(tempfolder):
    try:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0]
        return False
    except IndexError:
        progress.setConsoleInfo('ERROR: Parameter folder could not be found. Please execute step 1 first!')
        return True
        
def convert(MPHcyanoMaxValue, MPHchlThreshForFloatFlag):
    MPHcyanoMaxValueS = '%.2f' % MPHcyanoMaxValue
    MPHchlThreshForFloatFlagS = '%.2f' % MPHchlThreshForFloatFlag
    return MPHcyanoMaxValueS, MPHchlThreshForFloatFlagS

def create_parameterfile(tempdir, MPHvalidPixelExpression, MPHcyanoMaxValueS, MPHchlThreshForFloatFlagS, MPHexportMph):
    with open(tempdir + "WaterQualityParametersOLCI03.txt", "w") as text_file:
        text_file.write('mphValidExpression='+ MPHvalidPixelExpression + '\n')
        text_file.write('mphCyanoMaxValue='+ MPHcyanoMaxValueS + '\n') 
        text_file.write('mphChlThreshForFloatFlag='+ MPHchlThreshForFloatFlagS + '\n')
        text_file.write('mphExportMph='+ str(MPHexportMph).lower() + '\n')

def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        MPHcyanoMaxValueS, MPHchlThreshForFloatFlagS = convert(MPHcyanoMaxValue, MPHchlThreshForFloatFlag)
        create_parameterfile(tempdir, MPHvalidPixelExpression, MPHcyanoMaxValueS, MPHchlThreshForFloatFlagS, MPHexportMph)

execution(tempfolder)