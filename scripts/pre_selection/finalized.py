from ase.io import read, write
import ase.units as units
import numpy as np

from ase.calculators.gulp import GULP, Conditions
from multiprocessing import Process

import multiprocessing
import subprocess

import math
import sys
import os
import shutil


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def get_directories():

        directories = os.listdir('.')

        D = []

        for d in directories:
                if  os.path.isdir(d):
                        D.append(d)

        D_sort = sorted(D)

        return D_sort



#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------
list_crns = get_directories()

for crn in list_crns:

	os.chdir(crn)

	flag_exist = os.path.exists('opti.cif')
	
	os.chdir('..')

	if flag_exist == False:
		print('%s NOK: File Not found' % crn, flush=True)
		shutil.rmtree(crn)
	else:
		print('%s OK: The file exists' % crn, flush=True)

