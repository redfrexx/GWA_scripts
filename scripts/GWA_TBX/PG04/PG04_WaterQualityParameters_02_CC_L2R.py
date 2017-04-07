#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_CoastColour_L2R=name
##ParameterBoolean|SnTMap|Use SnT map|false
##ParameterNumber|AverageSalinity|Average Salinity|0|None|1.0
##ParameterNumber|AverageTemperature|Average Temperature|0|None|15.0
##ParameterBoolean|ExtremeCaseMode|Use extreme case mode|true
##ParameterString|LandExpression|Land expression|l1p_flags.CC_LAND
##ParameterString|SnowIceExpression|Snow ice expression|l1p_flags.CC_CLOUD or l1p_flags.CC_SNOW_ICE
##ParameterBoolean|L2RToa|Calculate L2R TOA|false
##ParameterSelection|L2RReflectAs|Reflection as|RADIANCE_REFLECTANCES; IRRADIANCE_REFLECTANCES

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
        
def convert(L2RReflectAs, AverageSalinity, AverageTemperature):
    if L2RReflectAs == 0:
        L2RReflectAsS = 'RADIANCE_REFLECTANCES'
    else:
        L2RReflectAsS = 'IRRADIANCE_REFLECTANCES'
    AverageSalinityS = '%.2f' % AverageSalinity
    AverageTemperatureS = '%.2f' % AverageTemperature
    
    return L2RReflectAsS, AverageSalinityS, AverageTemperatureS

def create_parameterfile(tempdir, SnTMap, AverageSalinityS, AverageTemperatureS, ExtremeCaseMode, LandExpression, SnowIceExpression, L2RToa, L2RReflectAsS):
    with open(tempdir + "WaterQualityParameters02.txt", "w") as text_file:
        text_file.write('snTMap='+ str(SnTMap).lower() + '\n')
        text_file.write('averageSalinity='+ AverageSalinityS + '\n') 
        text_file.write('averageTemperature='+ AverageTemperatureS + '\n')
        text_file.write('extremeCaseMode='+ str(ExtremeCaseMode).lower() + '\n')
        text_file.write('landExpression='+ LandExpression + '\n')
        text_file.write('snowIceExpression='+ SnowIceExpression + '\n')
        text_file.write('l2RToab='+ str(L2RToa).lower() + '\n')
        text_file.write('l2RReflectAs='+ L2RReflectAsS + '\n')


def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        L2RReflectAsS, AverageSalinityS, AverageTemperatureS = convert(L2RReflectAs, AverageSalinity, AverageTemperature)
        create_parameterfile(tempdir, SnTMap, AverageSalinityS, AverageTemperatureS, ExtremeCaseMode, LandExpression, SnowIceExpression, L2RToa, L2RReflectAsS)

execution(tempfolder)