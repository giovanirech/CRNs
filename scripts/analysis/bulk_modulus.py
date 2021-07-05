import subprocess
import math
import copy
import sys
import numpy as np
import glob
import os
import shutil

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_bulk_modulus():
    b = np.float(subprocess.check_output('grep \'Velocity S-wave\' gulp_full_opt.got | awk \'{print $7}\'', shell=True))
    return(b)
    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
list_crns = os.listdir('.')

for crn in list_crns:

    if os.path.isdir(crn):
        os.chdir(crn)
        b = get_bulk_modulus()
        os.chdir('..')
        print('%s %f' % (crn, b))
