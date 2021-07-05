from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.data import covalent_radii
import ase.units as units
import numpy as np

from ase.calculators.gulp import GULP, Conditions

import sys
import os
import shutil
import math
import subprocess
import statistics

#-------------------------------------------------------------------
#------------------------------------------------------------------------------

LIMIT = 30 

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def is_zero(constant):

	if math.fabs(constant) <= LIMIT:
		return True
	else:
		return False

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def born_criterion(C):

	C11 = ( C['C11'] + C['C22'] + C['C33'] ) / 3
	C12 = ( C['C12'] + C['C13'] + C['C23'] ) / 3
	C44 = ( C['C44'] + C['C55'] + C['C66'] ) / 3
                      
	if not(C11 - C12 > 0):
		return False

	if not(C11 + 2*C12 > 0):
		return False
	
	if not(C44 > 0):
		return False

	return True

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def isotropic_criterion(C):

	C11 = ( C['C11'] + C['C22'] + C['C33'] ) / 3
	C12 = ( C['C12'] + C['C13'] + C['C23'] ) / 3
	C44 = ( C['C44'] + C['C55'] + C['C66'] ) / 3

	if  math.fabs(C44-(C11-C12)/2) <= LIMIT:
		return True  		
	else:
		return False

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cubic_criterion(C):
	
	if np.std([C['C11'],C['C22'],C['C33']], ddof=1)>LIMIT:
		return False

	if np.std([C['C12'],C['C13'],C['C23']], ddof=1)>LIMIT:
		return False

	if np.std([C['C44'],C['C55'],C['C66']], ddof=1)>LIMIT:
		return False

	if is_zero(C['C14']) == False or is_zero(C['C15']) == False or \
	   is_zero(C['C16']) == False or is_zero(C['C24']) == False or \
	   is_zero(C['C25']) == False or is_zero(C['C26']) == False or \
	   is_zero(C['C34']) == False or is_zero(C['C35']) == False or \
	   is_zero(C['C36']) == False or is_zero(C['C45']) == False or \
	   is_zero(C['C46']) == False or is_zero(C['C56']) == False:
		return False

	return True

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def get_elastic_constants():

	constants = {}

	line = subprocess.check_output("grep -A 5 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)

	C11 = float(line[10:20])
	C12 = float(line[20:30])
	C13 = float(line[30:40])
	C14 = float(line[40:50])
	C15 = float(line[50:60])
	C16 = float(line[60:70])

	#print linha
	constants['C11'] = C11
	constants['C12'] = C12
	constants['C13'] = C13
	constants['C14'] = C14
	constants['C15'] = C15
	constants['C16'] = C16
	#print C11, C12, C13, C14, C15, C16

	line = subprocess.check_output("grep -A 6 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)

	C21 = float(line[10:20])
	C22 = float(line[20:30])
	C23 = float(line[30:40])
	C24 = float(line[40:50])
	C25 = float(line[50:60])
	C26 = float(line[60:70])
	
	#print linha
	#print C21, C22, C23, C24, C25, C26
	
	constants['C22'] = C22
	constants['C23'] = C23
	constants['C24'] = C24
	constants['C25'] = C25
	constants['C26'] = C26

	line = subprocess.check_output("grep -A 7 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)
	
	C31 = float(line[10:20])
	C32 = float(line[20:30])
	C33 = float(line[30:40])
	C34 = float(line[40:50])
	C35 = float(line[50:60])
	C36 = float(line[60:70])

	#print linha	
	#print C31, C32, C33, C34, C35, C36

	constants['C33'] = C33
	constants['C34'] = C34
	constants['C35'] = C35
	constants['C36'] = C36

	line = subprocess.check_output("grep -A 8 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)

	C41 = float(line[10:20])
	C42 = float(line[20:30])
	C43 = float(line[30:40])
	C44 = float(line[40:50])
	C45 = float(line[50:60])
	C46 = float(line[60:70])

	#print linha	
	#print C41, C42, C43, C44, C45, C46

	constants['C44'] = C44
	constants['C45'] = C45
	constants['C46'] = C46

	line = subprocess.check_output("grep -A 9 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)

	C51 = float(line[10:20])
	C52 = float(line[20:30])
	C53 = float(line[30:40])
	C54 = float(line[40:50])
	C55 = float(line[50:60])
	C56 = float(line[60:70])

	#print linha
	#print C51, C52, C53, C54, C55, C56

	constants['C55'] = C55
	constants['C56'] = C56

	line = subprocess.check_output("grep -A 10 'Elastic Constant Matrix: (Units=GPa)'  gulp_full_opt.got | tail -1", shell=True)
	
	C61 = float(line[10:20])
	C62 = float(line[20:30])
	C63 = float(line[30:40])
	C64 = float(line[40:50])
	C65 = float(line[50:60])
	C66 = float(line[60:70])

	#print linha	
	#print C61, C62, C63, C64, C65, C66
        
	constants['C66'] = C66

	return constants

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

	#gulp_execution()

	C = get_elastic_constants()

	#print(C)

	cubic =  cubic_criterion(C)

	born =  born_criterion(C)

	isotropic = isotropic_criterion(C)

	os.chdir('..')

	if not(cubic):
		print('NOK %s: cubic criterion' % crn)
		shutil.rmtree(crn)
	elif not(born):
		print('NOK %s: born criterion' % crn)
		shutil.rmtree(crn)
	elif not(isotropic):
		print('NOK %s: isotropic criterion' % crn)
		shutil.rmtree(crn)
	else:
		print('OK %s' % crn)

