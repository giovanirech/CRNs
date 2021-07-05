import os, sys
import subprocess


directories = os.listdir('.')

for d in directories:
	if os.path.isdir(d):
		os.chdir(d)
		print(os.getcwd(), flush=True)
		os.system('python3 /home/almartin/CRNs/scripts/qha.py crn.xyz >> output.txt')
		os.chdir('..')
