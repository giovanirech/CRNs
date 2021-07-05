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



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def  negative_frequencies():

	fin = open('gulp.got', 'r')
	
	lines = fin.readlines()

	i = 0
	for l in lines:
		if l.find('Frequencies (cm-1) [NB: Negative implies an imaginary mode]:') != -1:
			if lines[i+2][0:8].find('*******') != -1:
				return True			
			else:			
				f  = round(float(lines[i+2][0:8]),2)
				if f < 0:
					return True 
		i = i + 1

	fin.close()

	return False

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def gulp_execution():

	atoms = read('opti.cif')
	c = Conditions(atoms)
	calc = GULP(keywords='conp prop phonon', library='rebo', conditions=c)
	atoms.set_calculator(calc)	
	print(atoms.get_potential_energy())
	

#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------
list_crns = get_directories()

for crn in list_crns:

	print(crn)

	os.chdir(crn)

	gulp_execution()

	flag_neg_freq = negative_frequencies()

	os.chdir('..')

	if flag_neg_freq == True:
		print('%s NOK: Negative Frequencies' % crn, flush=True)
		shutil.rmtree(crn)
	else:
		print('%s OK' % crn, flush=True)

