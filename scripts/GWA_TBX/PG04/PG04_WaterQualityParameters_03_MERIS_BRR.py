#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_MERIS_RayleighCorrection=name
##ParameterBoolean|TiePoints|Copy all tie points|true
##ParameterBoolean|L1Flags|Copy L1 flags|true
##ParameterBoolean|OutputToar|Output TOA|false
##ParameterSelection|CorrectionSurface|Correction Surface|ALL_SURFACES; LAND; WATER

TiePoints=True
L1Flags=True
OutputToar=False
CorrectionSurface='ALL_SURFACES'

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
        
def convert(CorrectionSurface):
    if CorrectionSurface == 0:
        CorrectionSurfaceS = 'ALL_SURFACES'
    if CorrectionSurface == 1:
        CorrectionSurfaceS = 'LAND'
    if CorrectionSurface == 2:
        CorrectionSurfaceS = 'WATER'
    return CorrectionSurfaceS

def create_parameterfile(tempdir, TiePoints, L1Flags, OutputToar, CorrectionSurfaceS):
    with open(tempdir + "WaterQualityParameters03.txt", "w") as text_file:
        text_file.write('tiePoints='+ str(TiePoints).lower() + '\n')
        text_file.write('l1Flags='+ str(L1Flags).lower() + '\n') 
        text_file.write('outputToar='+ str(OutputToar).lower() + '\n')
        text_file.write('correctionSurface='+ CorrectionSurfaceS + '\n')

def execution(tempfolder):
    if folder_check(tempfolder):
        return
    else:
        tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0] + '/'
        CorrectionSurfaceS = convert(CorrectionSurface)
        create_parameterfile(tempdir, TiePoints, L1Flags, OutputToar, CorrectionSurfaceS)

execution(tempfolder)