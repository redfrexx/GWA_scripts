#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityWorkflow_TempDataCleanUp=name

import shutil
import glob
import os
import tempfile

tempfolder = 'wq_scripts_'
tempdir = glob.glob(os.path.join(tempfile.gettempdir(), tempfolder + '*'))[0]
shutil.rmtree(tempdir)