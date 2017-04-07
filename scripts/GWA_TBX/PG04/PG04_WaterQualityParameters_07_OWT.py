#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_OWT=name
##ParameterString|OWTreflectancesPrefix|Reflectance prefix|reflec
##ParameterBoolean|writeInputReflectances|Write input reflectances|false

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

def create_parameterfile(tempdir, OWTreflectancesPrefix, writeInputReflectances):
    with open(tempdir + "WaterQualityParameters07.txt", "w") as text_file:
        text_file.write('OWTreflectancesPrefix='+ OWTreflectancesPrefix + '\n')
        text_file.write('OWTwriteInputReflectances='+ str(writeInputReflectances).lower() + '\n') 

def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        create_parameterfile(tempdir, OWTreflectancesPrefix, writeInputReflectances)

execution(tempfolder)