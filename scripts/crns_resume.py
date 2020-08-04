#!/usr/bin/python3
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import shutil
import os
import glob
import xlsxwriter

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def get_volume(crn):
	v = np.float(subprocess.check_output('grep \'Initial cell volume\' ' + crn +'/1.0000/gulp_opti_1.0000.got | awk \'{print $5}\'', shell=True))
	return v

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_hybridization(crn):
	string = subprocess.check_output('grep \'\[0\' ' +  crn + '/output.txt | tail -n 1', shell=True).decode('utf-8')
	res = string.replace('[', ' ').replace(']', ' ')
	sp = int(res.split(',')[2])
	sp2 = int(res.split(',')[3])
	sp3 = int(res.split(',')[4])

	return sp, sp2, sp3
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def sp2_hybridization_percentage(sp, sp2, sp3):
	sp2_perc = sp2 / ( sp + sp2 + sp3) * 100
	return sp2_perc

	
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_bulk_modulus(crn):
	b = np.float(subprocess.check_output('grep \'Bulk modulus:\' ' + crn +'/output.txt | awk \'{print $3}\'', shell=True))
	return b

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_cet(c, temp):
	mcet  = np.loadtxt(c + '/thermal_expansion_coefficient.dat')

	for t,cet in mcet:
		if t == temp:
			return cet  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_list_directories():
	D = []
	for l in os.listdir("."):
		if os.path.isdir(l) and l != 'amorph_arquivos':
			D.append(l)
	return sorted(D)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

if ( len(sys.argv) != 1 ):
	print (sys.argv[0], flush=True)
	sys.exit(0)

CRNS = get_list_directories()

workbook = xlsxwriter.Workbook('crn.xlsx')

worksheet = workbook.add_worksheet() 

headers = ['CRN', 'VOL', 'B0', 'CET_100', 'CET_300', 'SP2_PERC', 'SP' , 'SP2', 'SP3' ]

headers_format = workbook.add_format({'bold': True, 'bg_color': 'yellow', 'valign': 'center', 'border': 1})

column = 0
for h in headers:
	worksheet.write(0, column, h, headers_format)
	column += 1

cell_format = workbook.add_format({'border': 1})
row = 1
for c in CRNS:
	column = 0

	name = c
	worksheet.write(row, column, name, cell_format)
	column += 1
	
	v = get_volume(c)
	worksheet.write(row, column, v, cell_format)
	column += 1
	
	b = get_bulk_modulus(c)
	worksheet.write(row, column, b, cell_format)
	column += 1
	
	cet100 = get_cet(c, 100)
	worksheet.write(row, column, cet100, cell_format)
	column += 1

	cet300 = get_cet(c, 300)
	worksheet.write(row, column, cet300, cell_format)
	column += 1

	sp, sp2, sp3 = get_hybridization(c)
	
	sp2_perc = sp2_hybridization_percentage(sp, sp2, sp3)
	worksheet.write(row, column, sp2_perc, cell_format)
	column += 1

	worksheet.write(row, column, sp, cell_format)
	column += 1

	worksheet.write(row, column, sp2, cell_format)
	column += 1

	worksheet.write(row, column, sp3, cell_format)
	column += 1
	
	row += 1

workbook.close() 	
