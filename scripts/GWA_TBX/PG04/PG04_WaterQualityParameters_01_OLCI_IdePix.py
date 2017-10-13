#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_01_OLCI_IdePix=name
##ParameterBoolean|CloudBuffer|Process cloud buffer|true
##ParameterNumber|CloudBufferWidth|lnsert size of cloud buffer in pixel|0|None|2
##ParameterBoolean|OutputCloudProbabilityFeatureValue|Write cloud probability feature value|false

import os
import glob
import tempfile

tempfolder = 'wq_scripts_'
OutputCloudProbabilityFeatureValue = 'false'

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

def create_parameterfile(tempdir, CloudBuffer, CloudBufferWidth,  OutputCloudProbabilityFeatureValue):
    with open(tempdir + "WaterQualityParametersOLCI01.txt", "w") as text_file:
        text_file.write('cloudBuffer='+ str(CloudBuffer) + '\n')
        text_file.write('cloudBufferWidth='+ str(CloudBufferWidth) + '\n')
        text_file.write('outputCloudProbabilityValue='+ str(OutputCloudProbabilityFeatureValue).lower() + '\n')

def execution(tempfolder):
    tempdir = folder_create(tempfolder) + '/'
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        create_parameterfile(tempdir, CloudBuffer, CloudBufferWidth, OutputCloudProbabilityFeatureValue)

execution(tempfolder)