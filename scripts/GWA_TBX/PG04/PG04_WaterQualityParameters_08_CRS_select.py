#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_CRS_select=name
##Define_output_CRS=crs


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

def create_parameterfile(tempdir, Define_output_CRS):
    with open(tempdir + "WaterQualityParameters08.txt", "w") as text_file:
        text_file.write('crs='+ Define_output_CRS + '\n')

def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        create_parameterfile(tempdir, Define_output_CRS)

execution(tempfolder)