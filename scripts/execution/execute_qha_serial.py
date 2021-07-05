import os, sys
import subprocess

directories = os.listdir('.')

np = os.cpu_count()

process_list  = []

ncrns = 0

for d in directories:

	if os.path.isdir(d):
		ncrns = ncrns + 1
		os.chdir(d)
		print(os.getcwd(), flush=True)
		p = subprocess.Popen('python3 /home/almartin/CRNs/scripts/qha.py crn.xyz >> output.txt', shell=True)
		process_list.append(p)

		if len(process_list) == np or ncrns == len(directories):
			for p in process_list:
				p.wait()
			process_list = []

		os.chdir('..')

