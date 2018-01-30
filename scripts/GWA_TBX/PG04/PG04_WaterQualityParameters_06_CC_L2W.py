#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_CoastColour_L2W=name
##ParameterString|invalidPixelExpression|Invalid pixel expression|l2r_flags.INPUT_INVALID
#ParameterBoolean|outputAOT550|write AOT 550|false
##ParameterSelection|owtType|Select OWT type|INLAND; COASTAL; INLAND_NO_BLUE_BAND; GLASS_5C; GLASS_6C, GLASS_6C_NORMALIZED 

outputAOT550=False

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

def convert(owtType):
    if owtType == 0:
        owtTypeS = 'INLAND'
    if owtType == 1:
        owtTypeS = 'COASTAL'
    if owtType == 2:
        owtTypeS = 'INLAND_NO_BLUE_BAND'
    if owtType == 3:
        owtTypeS = 'GLASS_5C'
    if owtType == 4:
        owtTypeS = 'GLASS_6C'
    if owtType == 5:
        owtTypeS = 'GLASS_6C_NORMALIZED'
    return owtTypeS

def create_parameterfile(tempdir, invalidPixelExpression, outputAOT550, owtTypeS):
    with open(tempdir + "WaterQualityParameters05.txt", "w") as text_file:
        text_file.write('L2WinvalidPixelExpression='+ invalidPixelExpression + '\n')
        text_file.write('L2WoutputAOT550='+ str(outputAOT550).lower() + '\n') 
        text_file.write('L2WowtType='+ owtTypeS + '\n')

def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        owtTypeS = convert(owtType)
        create_parameterfile(tempdir, invalidPixelExpression, outputAOT550, owtTypeS)

execution(tempfolder)