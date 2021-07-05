from ase.calculators.lammpsrun import LAMMPS

from ase.io import read, write
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
from ase import Atoms

from ase.units import kJ
from ase.eos import EquationOfState

from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.io.trajectory import Trajectory
from ase.md import MDLogger

from ase.neighborlist import NeighborList
from ase.data import covalent_radii

import ase.units as units

import math
import copy
import sys
import numpy as np
import glob
import os
import shutil

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
DISTANCE_MIN = 1.7 
DISTANCE_MAX = 2.0

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def hybridization_function(distance):

	if ( distance <= DISTANCE_MIN):
		f = 1.0
	elif ( distance < DISTANCE_MAX): 
		c = (math.pi*(distance - DISTANCE_MIN)) / (DISTANCE_MAX - DISTANCE_MIN)
		f = ( 1.0 + math.cos(c))/2.0  
	else:
		f = 0.0
	return f

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def hybridization_calculation(atoms):
	
	cutoffs =  covalent_radii[atoms.numbers]
	nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True) 
	nl.update(atoms)

	sp_number = [0,0,0,0,0,0,0,0,0,0]

	for i in range(len(atoms)): 

		indexes, offsets = nl.get_neighbors(i)
		
		hib = 0

		for j in indexes:

			d =  atoms.get_distance(i, j, True)

			if d <= 1.85:
				hib = hib + 1

			#hib = hib + hybridization_function(d)

		#ihib = int(round(hib))
		
		#sp_number[ihib] = sp_number[ihib] + 1
	
		sp_number[hib] = sp_number[hib] + 1

	print("-------------------------------------",flush=True)
	print("------- HYBRIDIZATION ---------------",flush=True)
	print(sp_number,flush=True)
	print("-------------------------------------",flush=True)
	print("-------------------------------------",flush=True)

	return sp_number

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def hybridization( sp_number ):

	for i in range(len(sp_number)):
		if (i < 2 or i > 4) and sp_number[i]!=0:
			return False
	return True

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def percentage_sp ( sp_number, num_atoms ):

	MAX_SP = 3

	percentage_sp = float(sp_number[2])/num_atoms * 100

	print('Sp percentage: %f' % percentage_sp, flush=True)

	if percentage_sp <= MAX_SP:

		return True
	else:
		return False

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

list_crns = os.listdir('.')

for crn in list_crns:

	if os.path.isdir(crn):

		os.chdir(crn)
		
		atoms = read('opti.cif')

		sp_number = hybridization_calculation(atoms)

		flag_hib = hybridization(sp_number)
		
		flag_perc_sp = percentage_sp(sp_number, atoms.get_global_number_of_atoms())  

		os.chdir('..')
			
		if flag_hib == False:
			print('CRN %s: NOK - Incorrect hybridization' % crn, flush=True)
			shutil.rmtree(crn)
		elif flag_perc_sp == False:
			print('CRN %s: NOK - Incorrect SP Percentage' %crn, flush=True)
			shutil.rmtree(crn)
		else:
			print('CRN: OK', flush=True)		
		


