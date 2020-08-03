#!/usr/bin/python3
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import shutil
import os
import glob

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def  negative_frequencies(directories):
	
	for d in directories:
		f = glob.glob(d + '/*.out')[0]
		freq = np.loadtxt(f)

		for v in freq:
			if v < 0:
				return True
	return False

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def optimization_error(directories):

	for d in directories:
		try:
			ret = subprocess.check_output('grep \'**** Optimisation achieved ****\' ' + d + '/gulp_opti_*.got', shell=True)
		except:
			return True

	return False	

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_list_directories():
	D = []
	for l in os.listdir("."):
		if os.path.isdir(l) and l != 'amorph_arquivos':
			D.append(l)
	return sorted(D)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

if ( len(sys.argv) != 1 ):
	print (sys.argv[0], flush=True)
	sys.exit(0)

CRNS = get_list_directories()

for c in CRNS:

	os.chdir(c)

	D = get_list_directories()
	
	if optimization_error(D):
		print('%s: Optimization Error' % c )
		#shutil.rmtree(crn)
	elif negative_frequencies(D):
		print('%s: Negative Frequencies ' % c) 
		#shutil.rmtree(crn)
	else:
		print('%s: OK!' % c) 
			
	#	print(c)

	os.chdir('..')


	
		
		

