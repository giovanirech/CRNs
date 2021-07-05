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
def distance(atoms):
	
    cutoffs =  covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True) 
    nl.update(atoms)

    vol_tot = 0

    for i in range(len(atoms)): 

        indexes, offsets = nl.get_neighbors(i)
   
        dist = 0
        hib = 0
            
        for j in indexes:
            d =  atoms.get_distance(i, j, True)
            if d <= 1.85:
                dist = dist + d
                hib = hib + 1
   
        dist_ave = dist / hib
        r = dist_ave / 2
        va = (4 * math.pi * r**3) / 3 

        vol_tot = vol_tot + va

    fea = vol_tot / atoms.get_volume()
    
    return fea                  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
list_crns = os.listdir('.')

for crn in list_crns:

	if os.path.isdir(crn):

		os.chdir(crn)
		atoms = read('opti.cif')
		fea = distance(atoms)
		os.chdir('..')
		print('%s %f %f' % (crn, fea, atoms.get_volume()))
