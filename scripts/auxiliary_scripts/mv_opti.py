import os, sys
import subprocess

if len(sys.argv)!=2:
	print ('%s <destination directory>' % sys.argv[0])
	sys.exit(0)

dest = sys.argv[1]
	
directories = os.listdir('.')

for d in directories:
	if os.path.isdir(d):
		os.chdir(d)
		exists = os.path.exists('opti.cif')
		os.chdir('..')
		if exists:
			cmd = 'mv ' + d + ' ' + dest
			os.system(cmd)
