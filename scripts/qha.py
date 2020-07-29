import subprocess
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import shutil
import os

from ase.spacegroup import crystal
from ase import Atoms
from ase.calculators.gulp import GULP, Conditions
from ase.io import read, write
from ase.units import kJ
from ase.eos import EquationOfState

from ase.neighborlist import NeighborList
from ase.data import covalent_radii

from scipy.optimize import minimize 
from scipy.optimize import curve_fit
from lmfit import fit_report, Minimizer, Parameters

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

DISTANCE_MIN = 1.7 
DISTANCE_MAX = 2.0

DISTANCE = 1.85

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def change_format_xyz(input_file, output_file):
	
	fin = open(input_file, 'r')
	fout = open(output_file, 'w')
	
	lines = fin.readlines()
	
	for i,l in enumerate(lines):

		if i == 1:
			a = float(l.split()[0])
			line_out = "Lattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" Properties=species:S:1:pos:R:3 spacegroup=\"P 1\" unit_cell=conventional pbc=\"T T T\"\n" % (a,a,a)
			fout.write(line_out)
		else:
			fout.write(l)
		
	fin.close()
	fout.close() 


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

			if d <= DISTANCE:
				hib = hib + 1

			#hib = hib + hybridization_function(d)

		#ihib = int(round(hib))
		
		#sp_number[ihib] = sp_number[ihib] + 1
	
		sp_number[hib] = sp_number[hib] + 1


	print("-------------------------------------", flush=True)
	print("------- HYBRIDIZATION ---------------", flush=True)
	print(sp_number, flush=True)
	print("-------------------------------------", flush=True)
	print("-------------------------------------", flush=True)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def molecular_dynamics(atoms):

	c = Conditions(atoms)
	calc = GULP(keywords='conp md isotropic', options=['integrator leapfrog verlet', 'ensemble npt 0.0005 0.0005', 'temperature 50','equilibration 5.0 ps', 'timestep 0.0001 ps', 'sample 0.0050', 'output cif dm.cif'], library='brenner', conditions=c)
	
	atoms.set_calculator(calc)
	print(atoms.get_potential_energy(), flush=True)

	shutil.copy2('gulp.gin', 'gulp_md.gin') 
	shutil.copy2('gulp.got', 'gulp_md.got') 

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def rfo_optimization(atoms, cif_opti_file):

	c = Conditions(atoms)
	calc = GULP(keywords='conp opti rfo isotropic prop phonon', options=['maxcyc 2000', 'output cif ' + cif_opti_file],  library='brenner ', conditions=c)
	atoms.set_calculator(calc)
	print(atoms.get_potential_energy(), flush=True)

	shutil.copy2('gulp.gin','gulp_rfo.gin')
	shutil.copy2('gulp.got','gulp_rfo.got')


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def full_optimization(atoms, cif_opti_file):

	c = Conditions(atoms)
	calc = GULP(keywords='conp opti lbfgs isotropic prop phonon', options=['maxcyc 2000', 'output cif ' + cif_opti_file],  library='brenner', conditions=c)
	atoms.set_calculator(calc)
	print(atoms.get_potential_energy(), flush=True)

	shutil.copy2('gulp.gin', 'gulp_full_opt.gin')
	shutil.copy2('gulp.got', 'gulp_full_opt.got')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def atomic_position_optimization(atoms, cif_opti_file):

	c = Conditions(atoms)
	calc = GULP(keywords='conv opti lbfgs isotropic', options=['maxcyc 2000', 'output cif ' + cif_opti_file], library='brenner', conditions=c)
	atoms.set_calculator(calc)
	pe = atoms.get_potential_energy()
	
	gulp_file_in = 'gulp_' + cif_opti_file[0:cif_opti_file.find(',')] + '.gin'
	gulp_file_out = 'gulp_' + cif_opti_file[0:cif_opti_file.find(',')] + '.got'
	shutil.copy2('gulp.gin', gulp_file_in)
	shutil.copy2('gulp.got', gulp_file_out)
	#return(pe)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def get_free_energy(atoms, factor_mul, temperature):
	
	c = Conditions(atoms)
	calc = GULP(keywords='phon noden prop', options=[temperature, 'shrink 8 8 8'], library='brenner', conditions=c)
	atoms.set_calculator(calc)
	atoms.get_potential_energy()
	
	gulp_file_in = 'gulp_free_' + str(factor_mul) + '.gin'
	gulp_file_out = 'gulp_free_' + str(factor_mul) + '.got'
		
	shutil.copy2('gulp.gin', gulp_file_in)
	shutil.copy2('gulp.got', gulp_file_out)
	
	fe_energy = np.float_(subprocess.check_output('grep \'Helmholtz free-energy\'  ' + gulp_file_out + '| awk \'{print $4}\'', shell=True).decode('utf-8').split())
	
	potential_energy =  float(subprocess.check_output('grep  \'Brenner potentials\' ' + gulp_file_out + '| awk \'{print $4}\'', shell=True).decode('utf-8'))
	
	zero_point_energy =  float(subprocess.check_output('grep  -m 1 \'Zero point energy\' ' + gulp_file_out + '| awk \'{print $5}\'', shell=True).decode('utf-8'))

	return potential_energy, zero_point_energy, fe_energy

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def  fit_potential_energy_vs_volume(V, P):
	eos = EquationOfState(V, P)
	v0, e0, B = eos.fit()
	print('Bulk modulus: %f GPa' % (B / kJ * 1.0e24), flush=True )
	eos.plot('potential_energy.pdf')
	plt.close('all')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def  fit_free_energy_vs_volume(V, P,t ):
	
	eos = EquationOfState(V, P)
	v0, e0, B = eos.fit()
	file_plot = 'free_energy_' + str(t) + '.pdf'
	eos.plot(file_plot)
	plt.close('all')	
	return(v0, e0, B)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def func(x, a0, a1, a2, a3):
	return a0 + (a1/x + a2*x) * np.exp(-a3/x)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def derived_function(x, a0, a1, a2, a3):
    return ((a2 - a1/x**2)*np.exp(-a3/x) + (a3*(a1/x + a2*x))/(np.exp(a3/x)*x**2)) 

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def func2min(params, T, V):
    a0, a1, a2, a3 = params['a0'], params['a1'], params['a2'], params['a3']
    model = func(T, a0, a1, a2, a3)
    return model - V

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def straight_equation(x, m, b):
	return m*x + b

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def fit_zero_point_energy_vs_volume(V,ZPE):

	popt, pcov = curve_fit(straight_equation, V, ZPE)
	m = popt[0]
	b = popt[1]
	
	print('m=%f b=%f' % (m,b))
	ZPE_FIT = m * np.array(V) + b
	plt.plot(V,ZPE,'bo', V, ZPE_FIT, 'b-')
	plt.xlabel('Volume $(\AA^3)$')
	plt.ylabel('Zero point energy (eV)')
	plt.savefig('volume_vs_zero_point_eneryg.pdf')
	plt.close('all')

	return ZPE_FIT 
	
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def fit_volume_vs_temperature(T, V):

	params = Parameters()
	params.add('a0', value=45, min=-100000, max=100000)
	params.add('a1', value=0, min=-1000, max=1000)
	params.add('a2', value=0, min=-1000, max=1000)
	params.add('a3', value=0, min=-1000, max=1000)

	min = Minimizer(func2min, params, fcn_args=(T, V))
	out = min.minimize(method='differential_evolution')
	fit = func2min(out.params, T, V)

	print(fit_report(out), flush=True)

	min = Minimizer(func2min, params=out.params, fcn_args=(T, V))
	out = min.minimize(method='leastsq')
	fit = func2min(out.params, T, V)

	print(fit_report(out), flush=True)

	a0 = out.params['a0'].value
	a1 = out.params['a1'].value
	a2 = out.params['a2'].value
	a3 = out.params['a3'].value

	return a0, a1, a2, a3
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def thermal_expansion_coefficient(T, a0, a1, a2, a3):

	return ( 1e6 * derived_function(T, a0, a1, a2, a3)/func(T, a0, a1, a2, a3) )
	

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_volume_vs_temperature(file_name):

	x, y = np.loadtxt(file_name, unpack=True)
	plt.plot(x,y,'ro-')
	plt.xlabel('Temperature (K)')
	plt.ylabel('Volume $(\AA^3)$')
	plt.savefig('volume_vs_temperature.pdf')
	plt.close('all')
	
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_thermal_expansion_vs_temperature(file_name):

	x, y = np.loadtxt(file_name, unpack=True)
	plt.plot(x,y,'ro-')
	plt.xlabel('Temperature (K)')
	plt.ylabel('Thermal expansion Coefficient ($\mu$K$^{-1}$)')
	plt.savefig('thermal_expansion_vs_temperature.pdf')
	plt.close('all')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def copy_intermediate_files():

	if os.path.exists('gulp_files'):
		shutil.rmtree('gulp_files')

	if os.path.exists('cif_files'):
		shutil.rmtree('cif_files')

	if os.path.exists('pdf_free_files'):
		shutil.rmtree('pdf_free_files')

	os.mkdir('gulp_files')
	os.system('mv *.gin  gulp_files')
	os.system('mv *.got  gulp_files')

	os.mkdir('cif_files')
	os.system('mv *.cif cif_files')

	os.mkdir('pdf_free_files')
	os.system('mv free_energy_*.pdf pdf_free_files')
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

if ( len(sys.argv) != 2 ):
	print (sys.argv[0] + ' <file.xyz>', flush=True)
	sys.exit(0)


change_format_xyz(sys.argv[1], 'amorph.xyz')

atoms = read('amorph.xyz')
print('Hibridizacao apos Amorph')

hybridization_calculation(atoms)

molecular_dynamics(atoms)

atoms = read('dm.cif')
print('Hibridizacao apos dinamica')
hybridization_calculation(atoms)

rfo_optimization(atoms, cif_opti_file = 'rfo.cif')	
atoms = read('rfo.cif')
print('Hibridizacao apos RFO')
hybridization_calculation(atoms)

full_optimization(atoms, cif_opti_file = 'opti.cif')

atoms = read('opti.cif')
print('Hibridizacao apos otimizacao completa')
hybridization_calculation(atoms)


displacement = 0
number_points = 9
start = - (number_points//2) + displacement
end =   number_points//2 + displacement

percentage = 3

V = []
PE = []
ZPE = []
FE = []

for p in list(range(start,end+1)):
		
	factor_mul = 1 + (p * (2.0*percentage)/(number_points-1))/100.0
	
	atoms = read('opti.cif')
	
	vol = factor_mul * atoms.get_volume() * 1.014
	a = vol ** (1./3)
	atoms.set_cell([a, a, a, 90, 90, 90], scale_atoms=True)
	
	v = atoms.get_volume()

	cif_opti_file = 'opti_' + str(factor_mul) + '.cif' 	
	atomic_position_optimization(atoms, cif_opti_file)
		
	atoms = read(cif_opti_file)
	print('Hibridizacao apos otimizacao a volume constante: %f' % factor_mul)
	hybridization_calculation(atoms)
	
	pe, zpe, fe_partial = get_free_energy(atoms, str(factor_mul), temperature='temperature 10 10 29')
	
	V.append(v)
	ZPE.append(zpe)	
	PE.append(pe)
	FE.append(fe_partial)


ZPE_FIT = fit_zero_point_energy_vs_volume(V,ZPE)

np.savetxt('potential_energy_vs_volume.dat',  np.transpose([V,PE]), fmt='%.6f')

fit_potential_energy_vs_volume(V,PE)

FE = np.array(FE)

L,C = np.shape(FE)

T = []
VA = []

t = 10
for k in range(C):
	F = FE[:,k] - ZPE + ZPE_FIT	
	v0, e0, B = fit_free_energy_vs_volume(V,F, t)
	T.append(t)
	t = t + 10
	VA.append(v0)

np.savetxt('volume_vs_temperature.dat', np.transpose([T,VA]), fmt='%.6f')

plot_volume_vs_temperature('volume_vs_temperature.dat')

T = np.array(T)
V = np.array(VA)

a0, a1, a2, a3 = fit_volume_vs_temperature(T,V)

A = thermal_expansion_coefficient(T, a0, a1, a2, a3)

np.savetxt('thermal_expansion_coefficient.dat', np.transpose([T,A]), fmt='%.6f')

plot_thermal_expansion_vs_temperature('thermal_expansion_coefficient.dat')

copy_intermediate_files()

