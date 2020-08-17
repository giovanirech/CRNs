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
def pre_optimization_error():
    try:
        ret = subprocess.check_output('grep \'**** Optimisation achieved ****\' gulp_pre_opt.got', shell=True)
    except:
        return True
    return False

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_list_directories():
	D = []
	for l in os.listdir("."):
		if os.path.isdir(l):
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

	flag_pre_opti_error = pre_optimization_error()

	os.chdir('..')
		
	if flag_pre_opti_error:
		print('%s: Pre optimization Error' % c )
		shutil.rmtree(c)
	else:
		print('%s: OK!' % c) 
		


