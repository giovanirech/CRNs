#!/usr/bin/python3
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import shutil
import os
import glob

from ase.spacegroup import crystal
from ase import Atoms
from ase.calculators.gulp import GULP, Conditions
from ase.io import read, write
from ase.units import kJ
from ase.eos import EquationOfState

from ase.neighborlist import NeighborList
from ase.data import covalent_radii

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

def atomic_position_optimization(atoms, factor_dir):

	c = Conditions(atoms)

	cif_opti_file = 'opti_' + factor_dir + '.cif'
	calc = GULP(keywords='conv opti lbfgs isotropic', options=['maxcyc 2000', 'output cif ' + cif_opti_file], library='brenner', conditions=c)
	atoms.set_calculator(calc)

	pe = atoms.get_potential_energy()

	gulp_file_in = 'gulp_opti_' + factor_dir + '.gin'
	gulp_file_out = 'gulp_opti_' + factor_dir + '.got'

	shutil.copy2('gulp.gin', gulp_file_in)
	shutil.copy2('gulp.got', gulp_file_out)

	#return(pe)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def phonon_calculation(atoms, factor_dir):

	c = Conditions(atoms)

	freq_file = 'freq_' + factor_dir + '.out'
	calc = GULP(keywords='phon noden prop', options=['shrink 8 8 8', 'output freq text ' + freq_file ], library='brenner', conditions=c)
	atoms.set_calculator(calc)

	atoms.get_potential_energy()

	gulp_file_in = 'gulp_phonon_' + factor_dir + '.gin'
	gulp_file_out = 'gulp_phonon_' + factor_dir + '.got'

	shutil.copy2('gulp.gin', gulp_file_in)
	shutil.copy2('gulp.got', gulp_file_out)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_volumes(directories):
	V = []
	for d in directories:
		vol = np.float(subprocess.check_output('grep \'Initial cell volume =\' ' +  d + '/gulp_opti_*.got | awk \'{print $5}\'', shell=True))
		V.append(vol)
	return V

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_potential_energy(directories):
	PE = []
	for d in directories:
		pe = np.float(subprocess.check_output('grep \'Final energy = \' ' +  d + '/gulp_opti_*.got | awk \'{print $4}\'', shell=True))
		PE.append(pe)
	return PE


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_Nkpoints(directories):
	N = []
	for d in directories:
		n = np.int(subprocess.check_output('grep \' Number of k points for this configuration =\' ' +  d + '/gulp_phonon_*.got | awk \'{print $9}\'', shell=True))
		N.append(n)
	return N

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_frequencies(directories):
	F = []
	for d in directories:
		f = glob.glob(d + '/*.out')[0]
		freq = np.loadtxt(f)
		F.append(freq)
	MF = np.transpose(np.array(F))
	return MF
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
def calc_free_energy(MF_FIT, N_kpoints, PE, T):
	h_eV_s = 4.135667696e-15
	c_cm_per_s = 2.99792458e+10
	Kb_eV_per_K = 8.617333262145e-5

	assert (np.shape(MF_FIT)[1]==len(N_kpoints))

	F = []
	for v in range(len(N_kpoints)):
		N = N_kpoints[v]
		freqs = MF_FIT[:,v]
		pe = PE[v]
		w = 1/N
		FT = []
		for t in T:
            		fvib = np.sum(w*0.5*h_eV_s*c_cm_per_s*freqs + w*Kb_eV_per_K*t*np.log(1-np.exp(-(h_eV_s*c_cm_per_s*freqs)/(Kb_eV_per_K*t))))
            		FT.append(fvib + pe)
		F.append(FT)

	return np.transpose(np.array(F))

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
def fit_zero_point_energy_vs_volume(V,ZPE):

	m,b = np.polyfit(V, ZPE, 1)

	print('m=%f b=%f' % (m,b))
	ZPE_FIT = m * np.array(V) + b
	plt.plot(V,ZPE,'ko', V, ZPE_FIT, 'k-')
	plt.xlabel('Volume $(\AA^3)$')
	plt.ylabel('Zero point energy (eV)')
	plt.savefig('volume_vs_zero_point_eneryg.pdf')
	plt.close('all')

	return ZPE_FIT

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def fit_frequency_vs_volume(V, MF):

	num_freq = np.shape(MF)[0]

	F_FIT = []

	for i in range(num_freq):
		F = MF[i,:]
		a, b, c = np.polyfit(V, F, 2)
		f_fit = a * np.power(np.array(V),2) + b * np.array(V) + c
		F_FIT.append(f_fit)

	MF_FIT = np.array(F_FIT)
	return MF_FIT
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def fit_volume_vs_temperature(T, V):

	params = Parameters()
	params.add('a0', min=0, max=100000)
	params.add('a1', min=-1000, max=1000)
	params.add('a2', min=-1000, max=1000)
	params.add('a3', min=-1000, max=1000)

	min = Minimizer(func2min, params, fcn_args=(T, V))
	kws={'popsize':30}
	out = min.minimize(method='differential_evolution', **kws)
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
def get_list_directories():
	D = []
	for l in os.listdir("."):
		if os.path.isdir(l):
			D.append(l)
	return sorted(D)

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

perc_list = list(range(start,end+1))
for p in perc_list:

	atoms = read('opti.cif')
	factor_mul = 1 + (p * (2.0*percentage)/(number_points-1))/100.0

	factor_dir = "{:.4f}".format(factor_mul)
	cur_dir = os.getcwd()
	os.mkdir(factor_dir)
	os.chdir(factor_dir)

	vol = factor_mul * atoms.get_volume() * 1.014
	a = vol ** (1./3)
	atoms.set_cell([a, a, a, 90, 90, 90], scale_atoms=True)
	atomic_position_optimization(atoms, factor_dir)

	cif_opti_file = 'opti_' + factor_dir + '.cif'

	atoms = read(cif_opti_file)
	print('Hibridizacao apos otimizacao a volume constante: %s' % factor_dir)
	hybridization_calculation(atoms)

	phonon_calculation(atoms, factor_dir)
	os.chdir(cur_dir)

D = get_list_directories()

V = get_volumes(D)
PE = get_potential_energy(D)
N_kpoints = get_Nkpoints(D)

np.savetxt('potential_energy_vs_volume.dat',  np.transpose([V,PE]), fmt='%.6f')
fit_potential_energy_vs_volume(V, PE)

MF = get_frequencies(D)
MF_FIT = fit_frequency_vs_volume(V, MF)

T = np.linspace(10,300,30)

MFREE_FIT = calc_free_energy(MF_FIT, N_kpoints, PE, T)

V_FIT = []
num_temp = np.shape(MFREE_FIT)[0]
for t in range(num_temp):
	F = MFREE_FIT[t,:]
	v0, e0, B = fit_free_energy_vs_volume(V,F, T[t])
	V_FIT.append(v0)

np.savetxt('volume_vs_temperature.dat', np.transpose([T,V_FIT]), fmt='%.6f')
plot_volume_vs_temperature('volume_vs_temperature.dat')

a0, a1, a2, a3 = fit_volume_vs_temperature(T,V_FIT)

A = thermal_expansion_coefficient(T, a0, a1, a2, a3)

np.savetxt('thermal_expansion_coefficient.dat', np.transpose([T,A]), fmt='%.6f')

plot_thermal_expansion_vs_temperature('thermal_expansion_coefficient.dat')



