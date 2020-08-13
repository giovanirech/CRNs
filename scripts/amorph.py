from multiprocessing import Process

import os
import shutil
import sys

#-------------------------------------------------------------------
#--------------------- AMORPH DEFINITION ---------------------------
#-------------------------------------------------------------------

file_crn_in = 'crn.in'
file_pre_in = 'pre.in'
file_import = 'crn.import'
file_run = 'run.sh'
amorph_directory_output = 'amorph_arquivos'
path_amorph = '/opt/amorph'
steps_pre_in = 1.28e7
steps_crn_in = 2.56e8
amorph_stdout = 'amorph.txt'
file_crn_xyz = 'crn.xyz'

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def get_params_script(sys):

	params = {}
	argv1, value1 = sys.argv[1].split('=')
	argv2, value2 = sys.argv[2].split('=')
	argv3, value3 = sys.argv[3].split('=')
	argv4, value4 = sys.argv[4].split('=')

	params[argv1]= float(value1)
	params[argv2]= float(value2)
	params[argv3]= float(value3)
	params[argv4]= float(value4)

	return params

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def execution_instances_amorph(name_crn_directory, ncrns, natoms, sp2, sp3):

	process = []

	for i in range(0,ncrns):
		try:	
			p = Process(target=amorph_execution, args=(name_crn_directory, i, natoms, sp2, sp3))
			p.start()
			process += [p]
		except:
			print ("-----------------------------------------------------", flush=True)
			print ("-------------- Erro: Amorph execution ---------------", flush=True )
			print ("-----------------------------------------------------", flush=True)
			sys.stdout.flush()
			os._exit(0)

	for p in process:
		p.join() 


#-------------------------------------------------------------------
#-------------------------------------------------------------------
def create_file_amorph_crn_in(file_crn_in):

	f = open(file_crn_in,'w')

	f.write('@import crn.import\n\n')
	f.write('[main]\n')
	f.write('output=1\n\n')

	f.write('[annealing]\n')
	f.write('source=pre.out\n')
	f.write('steps=' + str(steps_crn_in) + '\n')
#	f.write('steps= (' + str(passos_crn_in * 0.05 ) + ',' + str(passos_crn_in * 0.65 ) + ',' + str(passos_crn_in * 0.30 ) + ')' + '\n')
#	f.write('#(0.05,0.90,0.05)\n')
	f.write('acc=(0.95,0.75,0.40,1e-2)\n')

	f.close()

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def create_file_amorph_import(file_import, num_atoms, sp3, sp2):
	
	f = open(file_import,'w')

	str_num_atoms = str(num_atoms)
	str_sp2 = str(sp2/100.0)
	str_sp3 = str(sp3/100.0)

	title = str_num_atoms + '_sp3=' + str_sp3 +'_sp2=' + str_sp2

	f.write('[main]\n')
	f.write('title=' + title +'\n')
	f.write('atoms=' + str_num_atoms + '\n')
	f.write('log_surface=0\n\n')

	f.write('[params]\n')
	f.write('lambda_C=2.5\n') 
	f.write('lambda_E=1.0\n\n') 

	f.write('frac_sp3=' + str_sp3 + '\n')
	f.write('frac_sp2=' + str_sp2 + '\n\n')

	f.write('cost_nb=10\n')
	f.write('cost_sp0=5\n')
	f.write('cost_sp=2\n')
	f.write('cost_sp2=1.5\n')
	f.write('cost_sp3=1\n')
	f.write('cost_sp4=10\n\n')

	f.write('k_r=5.0\n')
	f.write('k_a=3.0\n') #corrigido conforme dissertacao FJ
	#f.write('k_a=2.0\n') valor usado ate 25/11/2016
	f.write('k_t=1.5\n')
	f.write('k_h=0.0\n')

	f.close()

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def create_file_amorph_pre_in(file_pre_in):

	f = open(file_pre_in,'w')

	f.write('@import crn.import\n\n')
	f.write('[main]\n')
	f.write('output=0\n\n')

	f.write('[annealing]\n')
	f.write('temp=1.0(1e3,1e-5)\n')
	f.write('steps=' + str(steps_pre_in) + '\n')

	f.close()

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def create_file_amorph_run(file_run):

	f = open(file_run,'w')

	f.write('#!/bin/bash\n\n')
	f.write(path_amorph + '/amorph pre.in\n')
	f.write(path_amorph + '/amorph crn.in\n')

	f.close()

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def amorph_run(file_run, destination_directory):

	#print os.getcwd()

	os.system('sh ' + file_run + ' > out.txt' )

	shutil.copy(file_crn_xyz, destination_directory ) 

#-------------------------------------------------------------------
#-------------------------------------------------------------------
def amorph_execution(name_crn_directory, i, natoms, sp2, sp3):

	current_directory = os.getcwd()

	crn_directory_output = current_directory + '/' + name_crn_directory + '_' + str(i)

	if not os.path.exists(crn_directory_output):
		
		os.makedirs(crn_directory_output)

	os.chdir(crn_directory_output)
	
	sys.stdout = open(amorph_stdout, 'w')
	
	if not os.path.exists(amorph_directory_output):
	
		os.makedirs(amorph_directory_output)

		os.chdir(amorph_directory_output)

		create_file_amorph_crn_in(file_crn_in)

		create_file_amorph_import(file_import, natoms, sp3, sp2)
	
		create_file_amorph_pre_in(file_pre_in)

		create_file_amorph_run(file_run)

		amorph_run(file_run, crn_directory_output)

		sys.stdout = sys.__stdout__
		

	os.chdir(current_directory)
	
#-------------------------------------------------------------------
#-------------------------------------------------------------------
if len(sys.argv) != 5 :
	print ('------------------------------------------------------------------------------------------------', flush=True)
	print (sys.argv[0] + ' ncrns=<number_of_crns> natoms=<number_of_atoms> sp2=<percentage> sp3=<percentage>', flush=True)
	print ('------------------------------------------------------------------------------------------------', flush=True)
	exit(0); 

params = get_params_script(sys)

ncrns = int(params['ncrns'])
natoms = int(params['natoms'])
sp3 = params['sp3']
sp2 = params['sp2']

execution_instances_amorph('crn', ncrns, natoms, sp2, sp3)




