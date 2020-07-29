#!/usr/bin/env python

import sys, os, time, shutil, glob
import subprocess

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
if len(sys.argv) != 2:
	print ('%s  <initial value>' % sys.argv[0])
	sys.exit(0)

i = int(sys.argv[1])

list_crns = get_directories()

for crn in list_crns:
	os.system('mv ' + crn + ' crn_' + str(i).zfill(5) )
	i = i + 1





