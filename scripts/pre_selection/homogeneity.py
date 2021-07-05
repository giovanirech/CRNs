from operator import itemgetter
from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.data import covalent_radii
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
import statistics

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

DISTANCE = 1.85
STD = 7.0

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

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def percentage_sp2(a, octant):

	cutoffs =  covalent_radii[a.numbers]
	nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True) 
	nl.update(a)

	sp_number = [0,0,0,0,0,0,0,0,0,0]

	for i in octant: 

		indexes, offsets = nl.get_neighbors(i)
		
		hib = 0

		for j in indexes:

			d =  atoms.get_distance(i, j, True)

			if d <= DISTANCE:
				hib = hib + 1

			#hib = hib + hybridization_function(d)

		#ihib = int(round(hib))
		
		#sp_number[ihib] = sp_number[ihib] + 1
	
		sp_number[hib] = sp_number[hib] + 1

	return float(sp_number[3]) / (sp_number[2] + sp_number[3] + sp_number[4]) * 100

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def subcube(p, pi, pf):

	x, y, z = p
	xi, yi, zi = pi
	xf, yf, zf = pf

	if x>=xi and x <=xf and y>=yi and y<=yf and z>=zi and z<=zf:
		return True


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def homogeneity(atoms):

	parameter =  atoms.get_cell_lengths_and_angles()[0]

	middle = parameter / 2

	octant1 = []
	octant2 = []
	octant3 = []
	octant4 = []
	octant5 = []
	octant6 = []
	octant7 = []
	octant8 = []


	percentages = [] 

	for i in range(atoms.get_number_of_atoms()):

		x, y, z = atoms[i].position

		p = [ (x + parameter) % parameter, (y + parameter) % parameter, (z + parameter) % parameter]
					
		if subcube(p , [0, 0, 0], [middle, middle, middle]):
			octant1.append(i)

		if subcube(p, [0, 0, middle], [middle, middle, parameter]):
			octant2.append(i)

		if subcube(p, [middle, 0, 0], [parameter, middle, middle]):
			octant3.append(i)

		if subcube(p,  [middle, 0, middle], [parameter, middle, parameter]):
			octant4.append(i)

		if subcube(p,  [0, middle, 0], [middle, parameter, middle]):
			octant5.append(i)

		if subcube(p,  [0, middle, middle], [middle, parameter, parameter]):
			octant6.append(i)

		if subcube(p,  [middle, middle, 0], [parameter, parameter, middle]):
			octant7.append(i)

		if subcube(p,  [middle, middle, middle], [parameter, parameter, parameter]):
			octant8.append(i)
	
	perc1 = percentage_sp2(atoms, octant1)
	perc2 = percentage_sp2(atoms, octant2)
	perc3 = percentage_sp2(atoms, octant3)
	perc4 = percentage_sp2(atoms, octant4)
	perc5 = percentage_sp2(atoms, octant5)
	perc6 = percentage_sp2(atoms, octant6)
	perc7 = percentage_sp2(atoms, octant7)
	perc8 = percentage_sp2(atoms, octant8)


	
	average = (perc1 + perc2 + perc3 + perc4 + perc5 + perc6 + perc7 + perc8)/8

	#print('Perc1: %f  Average: %f' % (perc1, average))
	#print('Perc2: %f  Average: %f' % (perc2, average))
	#print('Perc3: %f  Average: %f' % (perc3, average))
	#print('Perc4: %f  Average: %f' % (perc4, average))
	#print('Perc5: %f  Average: %f' % (perc5, average))
	#print('Perc6: %f  Average: %f' % (perc6, average))
	#print('Perc7: %f  Average: %f' % (perc7, average))
	#print('Perc8: %f  Average: %f' % (perc8, average))

	percs = [perc1 ,perc2, perc3, perc4, perc5, perc6, perc7, perc8]

	#return statistics.stdev(percs)
	if  statistics.stdev(percs) <= STD:
		return True
	else:
		return False

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

list_crns = os.listdir('.')

for crn in list_crns:

	if os.path.isdir(crn):

		os.chdir(crn)
		
		atoms = read('opti.cif')
		
		
		flag_hom = homogeneity(atoms)

		os.chdir('..')
			
		if flag_hom == False:
			print('CRN %s: NOK - Incorrect homogeneity' % crn, flush=True)
			shutil.rmtree(crn)
		else:
			print('CRN %s!!: OK' % crn, flush=True)		

