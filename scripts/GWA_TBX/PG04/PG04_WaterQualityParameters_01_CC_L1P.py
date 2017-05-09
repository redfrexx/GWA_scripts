#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_CoastColour_L1P=name
##ParameterBoolean|Icol|Icol|false
##ParameterBoolean|Calibration|Calibration|false
##ParameterBoolean|Smile|Smile|true
##ParameterBoolean|Equalization|Equalization|true
##ParameterBoolean|IgnoreSeaIceClim|Ignore sea ice climatology|false
##ParameterNumber|CloudBufferWidth|lnsert size of cloud buffer in pixel|0|None|2
##ParameterNumber|CloudScreeningAmbiguous|Cloud screening ambiguous threshold|0|None|1.4
##ParameterNumber|CloudScreeningSure|Cloud screening sure threshold|0|None|1.8
##ParameterNumber|GlintCloudThresholdAddition|Glint cloud screening addition|0|None|0.1
##ParameterBoolean|OutputCloudProbabilityFeatureValue|Write cloud probability feature value|false

import os
import glob
import tempfile

tempfolder = 'wq_scripts_'

def folder_create(tempfolder):
    try:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0]
        return tempdir
    except:
        progress.setConsoleInfo('Temporary folder:' + tempfolder + ' does not exist and will be created.')
        tempfile.mkdtemp(prefix=tempfolder)
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0]
        return tempdir

def folder_check(tempfolder):
    try:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0]
        return False
    except IndexError:
        progress.setConsoleInfo('ERROR: Temporary folder:' + tempfolder + ' cloud not be created. Check for administration rights to create folder.')
        return True

def create_parameterfile(tempdir, Icol, Calibration, Smile, Equalization, IgnoreSeaIceClim, CloudBufferWidth, CloudScreeningAmbiguous, CloudScreeningSure, GlintCloudThresholdAddition, OutputCloudProbabilityFeatureValue):
    with open(tempdir + "WaterQualityParameters01.txt", "w") as text_file:
        text_file.write('icol='+ str(Icol).lower() + '\n')
        text_file.write('calibration='+ str(Calibration).lower() + '\n')
        text_file.write('smile='+ str(Smile).lower() + '\n')
        text_file.write('equalization='+ str(Equalization).lower() + '\n')
        text_file.write('ignoreSeaIceClim='+ str(IgnoreSeaIceClim).lower() + '\n')
        text_file.write('cloudBufferWidth='+ str(CloudBufferWidth) + '\n')
        text_file.write('cloudScreeningAmbiguous='+ str(CloudScreeningAmbiguous) + '\n')
        text_file.write('cloudScreeningSure='+ str(CloudScreeningSure) + '\n')
        text_file.write('glintCloudThresholdAddition='+ str(GlintCloudThresholdAddition) + '\n')
        text_file.write('outputCloudProbabilityValue='+ str(OutputCloudProbabilityFeatureValue).lower() + '\n')

def execution(tempfolder):
    tempdir = folder_create(tempfolder) + '/'
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        create_parameterfile(tempdir, Icol, Calibration, Smile, Equalization, IgnoreSeaIceClim, CloudBufferWidth, CloudScreeningAmbiguous, CloudScreeningSure, GlintCloudThresholdAddition, OutputCloudProbabilityFeatureValue)

execution(tempfolder)