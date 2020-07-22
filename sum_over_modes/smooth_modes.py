import glob
import sys
import numpy as np
import argparse
import copy

#ARGPARSER STUFF
parser = argparse.ArgumentParser(description='Smooths frequency dependency with volume using a 2nd degree polynomial')

parser.add_argument('file_list', metavar='file_list', type=str, nargs='*', help='List of files to be analysed. Assumes that all files are outputs of GULP phonon calculations containing the same number of atoms and the same Monkhorst-Pack mesh. Wild cards are accepted, so instead of passing \n out_v1.gout out_v2.gout out_v3.gout out_v4.gout\n you can pass just \nout_v*.gout')

parser.add_argument('--zpe', help='Calcualtes the zero point energy for each volume',action='store_true' )

parser.add_argument('--nosave', help='Don\'t create text files with frequency modes. The default is to save all modes, the read and the smoothed in separate files, for each volume given.',action='store_false' )

parser.add_argument('-F', metavar='Temperature', nargs='*', help='Calcualtes the Helmholtz free-energy for the given temperature (in Kelvin). When this option is requested, the zero point energy is automatically calculated.', action='store' )

args = parser.parse_args()
####################################################################################
#ASSUMPTIONS:
#    ALL FILES HAVE THE SAME MONKHORST-PACK MESH
#    ALL FILES HAVE THE SAME NUMBER OF ATOMS PER UNIT CELL
#    FREQUENCIES ARE READ IN CM-1
####################################################################################

def get_kpoints_weights_and_coords(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        
    N = 0
    N_begin_klist = len(lines)
    kpoints_weights = []
    kpoints_coords = []
    
    for n, line in enumerate(lines):
        if 'Brillouin zone sampling points :' in line:
            N_begin_klist = n+5
        
        if '--------------------------------' in line and n>N_begin_klist:
            N_end_klist = n-1
            N, x, y, z, w = [float(a) for a in lines[N_end_klist].split()]
            break
        
    for line in lines[N_begin_klist:N_end_klist+1]:
        n, x, y, z, w = [float(a) for a in line.split()]
        kpoints_weights.append(w)
        kpoints_coords.append((x, y, z))
    
    if N==0:
        print('Could not find K points list while reading file {0}.'.format(file))
        return kpoints_weights, kpoints_coords
    
    return kpoints_weights, kpoints_coords



def get_frequencies(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    n_begin_freq = len(lines)
    getting_frequencies = False
    frequency_list = []
    
    Kpoints_freqs_list = []
    
    
    for n,line in enumerate(lines):
        
        if getting_frequencies and len(line.split())<1 and n>n_begin_freq:
            getting_frequencies = False
            frequency_list = np.array(frequency_list)
            Kpoints_freqs_list.append(frequency_list)
            frequency_list = []
            
        if 'Frequencies (cm-1) [NB: Negative implies an imaginary mode]:' in line:
            n_begin_freq = n+2
            getting_frequencies = True
                                  
        if getting_frequencies and n>=n_begin_freq:
            f = line.split()
            for i in f:
                frequency_list.append(float(i))
            
            
    return Kpoints_freqs_list

def get_volume(file):
    with open(file, 'r') as f:
        for line in f:
            if 'Initial cell volume =' in line:
                volume = float(line.split()[4])
                break    
    return volume

def get_zpe_gulp(file):
    with open(file, 'r') as f:
        for line in f:
            if 'Zero point energy            =' in line:
                zpe_gulp = float(line.split()[4])
                break 
    return zpe_gulp

def get_F_gulp(file):
    with open(file, 'r') as f:
        for line in f:
            if 'Helmholtz free-energy        =' in line:
                F_gulp = float(line.split()[3])
                break 
    return F_gulp

def get_temp_gulp(file):
    with open(file, 'r') as f:
        for line in f:
            if 'Phonon properties (per mole of unit cells): Temperature =' in line:
                temp = float(line.split()[-2])
                break 
    return temp

def get_Epot_gulp(file):
    with open(file, 'r') as f:
        for line in f:
            if ('Total lattice energy       =' in line) and ('eV' in line):
                Epot = float(line.split()[-2])
                break 
    return Epot

def get_data_gulp(file):
    '''
    Reads GULP output file and gather several informations. Returns several values:
    volume, zpe, F, temp, Epot = get_data_gulp(file)
    
    where:
       volume = unit cell volume (A^3)
       zpe = Zero point energy (eV)
       F = Helmholtz free-energy (eV)
       temp = Temperature (K)
       Epot = Potential energy (eV)
    '''
    with open(file, 'r') as f:
        Epot_gulp = False
        zpe_gulp = False
        F_gulp = False
        temp_gulp = False
        volume_gulp = False
        for line in f:
            if ('Total lattice energy       =' in line) and ('eV' in line):
                Epot_gulp = float(line.split()[-2])
            if 'Phonon properties (per mole of unit cells): Temperature =' in line:
                temp_gulp = float(line.split()[-2])
            if 'Helmholtz free-energy        =' in line:
                F_gulp = float(line.split()[3])
            if 'Zero point energy            =' in line:
                zpe_gulp = float(line.split()[4])
            if 'Initial cell volume =' in line:
                volume_gulp = float(line.split()[4])
                
        if not (Epot_gulp and zpe_gulp and F_gulp and temp_gulp and volume_gulp):
            print('!!WARNING: One or more properties were not found in the output file {0}. Please check the file!'.format(file))
            answer = input('Would you like to continue anyway (Values not found will be treated as zero)? (y/[n])')
            if answer == 'y' or answer == 'Y':
                pass
            else:
                exit()
    return volume_gulp, zpe_gulp, F_gulp, temp_gulp, Epot_gulp
##################################################################################

##################################################################################
#reads files and construct data structure (a list of lists of lists)

if len(sys.argv) < 2:
    parser.print_help()
    exit()


    
if __name__ == '__main__':
    #Test if arguments of -F are numbers and automatically sets --zpe option to True
    if args.F != None:
        for i in range(len(args.F)):
            try:
                float(args.F[i])
            except:
                print('All the temperature values for the -F option must be numerical!')
                exit()
        args.zpe = True
        
        
    file_list = args.file_list
    
    if len(file_list) < 3:
        print('You need at least 3 files!!!')
        exit()
    
    print('List of files to be used:')
    for file in file_list:
        print(file)
    
    answer = input('Are these the intended files? (y/[n])')
    
    if answer!='y' and answer!='Y':
        exit()
    else:
        pass
    
    data = [] #data[volume][kpoint][mode_number]
    v_list = []
    zpe_gulp_list = []
    F_gulp_list = []
    Epot_gulp_list=[]
    temperature_gulp_list=[]
    for file in file_list:
        volume, zpe_gulp, F_gulp, temperature_gulp, Epot_gulp = get_data_gulp(file)
        Kpoints_freqs_list = get_frequencies(file)
        kpoints_weights, kpoints_coords = get_kpoints_weights_and_coords(file)
        
        v_list.append(volume)
        zpe_gulp_list.append(zpe_gulp)
        data.append(Kpoints_freqs_list)
        F_gulp_list.append(F_gulp)
        Epot_gulp_list.append(Epot_gulp)
        temperature_gulp_list.append(temperature_gulp)
    
    #smooths frequencies using a second degree polynomial
    
    data_smooth = copy.deepcopy(data)
    total_number_of_kpoints = len(data[0])
    total_number_of_modes = len(data[0][0])
    
    print('Total number of k-points: {0}\nTotal number of modes per k-point: {1}'.format(
        total_number_of_kpoints,total_number_of_modes))
    
    
    for kpoint in range(total_number_of_kpoints):
        for mode_number in range(total_number_of_modes):
    
            x = v_list #eixo x = volumes
            y = [data[i][kpoint][mode_number] for i in range(len(v_list))] #eixo y = um modo, de um determinado k-point, variando com o volume
            
            parameters = np.polyfit(x, y, 2)
            polynomial = np.poly1d(parameters)
            y_smooth = polynomial(x)
            
            for i in range(len(v_list)):
                data_smooth[i][kpoint][mode_number]=y_smooth[i]
    
    
    
    #writes all data structures to text files, unless --nosave is given
    if args.nosave:
        for i in range(len(v_list)):
            volume = '{0:.2f}'.format(v_list[i])
            
            #original modes calculated by gulp
            with open('modes_original_v{0}.txt'.format(volume), 'w') as f:
                f.write('#K-point \tx-coord \ty-coord \tz-coord \tweight \tfrequencies(cm-1)*3Ncolumns \n')
                
                for j in range(total_number_of_kpoints):
                    modes_string = '\t'.join(['{:9.2f}'.format(s) for s in data[i][j]])
                    initial_string = '{0}\t{1}\t{2}\t{3}\t{4}'.format(j+1, *kpoints_coords[j], kpoints_weights[j])
                    f.write(initial_string + '\t' + modes_string + '\n')
            
            #smoothed modes 
            with open('modes_smooth_v{0}.txt'.format(volume), 'w') as f:
                f.write('#K-pt \tx-coord \ty-coord \tz-coord \trel. weight \tfrequencies(cm-1)*3Ncolumns \n')
                
                for j in range(total_number_of_kpoints):
                    modes_string = '\t'.join(['{:9.2f}'.format(s) for s in data_smooth[i][j]])
                    initial_string = '{0:4}\t{1:8.6f}\t{2:8.6f}\t{3:8.6f}\t{4:10.8f}'.format(j+1, *kpoints_coords[j], kpoints_weights[j])
                    f.write(initial_string + '\t' + modes_string + '\n')

    
    
    ## CALCULATION OF ZERO POINT ENERGY, IF REQUIRED
    if args.zpe:
        h_eV_s = 4.135667696e-15 #planks constant in eV*s
        c_cm_per_s = 2.99792458e+10 # speed of light in cm/s
        print('\n{:-^65}'.format('-'))       
        print('Zero Point Energy calculation (in eV):')
        print('{0:^15} {1:^15} {2:^15} {3:^15}'.format('Volume (A^3)','GULP ZPE', 'Calc ZPE', 'Smooth ZPE'))
        print('{:-^65}'.format('-'))
        ZPE_list = []
        ZPE_smooth_list = []
        for i in range(len(v_list)):
            ZPE = 0.0
            ZPE_smooth = 0.0
            for k in range(total_number_of_kpoints):
                for f in range(total_number_of_modes):
                    ZPE += kpoints_weights[k]*0.5*h_eV_s*c_cm_per_s*data[i][k][f]
                    ZPE_smooth += kpoints_weights[k]*0.5*h_eV_s*c_cm_per_s*data_smooth[i][k][f]
            ZPE_list.append(ZPE)
            ZPE_smooth_list.append(ZPE_smooth)
        
            print('{0:^15f} {1:^15f} {2:^15f} {3:^15f}'.format(v_list[i],zpe_gulp_list[i], ZPE, ZPE_smooth))
        print('{:-^65}'.format('-'))
        
        
    ## CALCULATION OF HELMHOLTZ FREE-ENERGY, IF REQUIRED
    if args.F != None:
        Kb_eV_per_K = 8.617333262145e-5 #Boltzmann constant in eV per kelvin
        
        print('\n{:=^65}'.format('='))
        print('Helmholtz free-energy calculation:\n')
        
        for i in range(len(args.F)):
            temp = float(args.F[i])
            
            different_temperature_found = False
            for t in temperature_gulp_list:
                if t != temp:
                    different_temperature_found=True
            if different_temperature_found:
                print('''WARNING: The temperature of the output GULP files is different for the requested temperature\nfor the Helmholtz free-energy calculation. The values in the GULP column below refer to the temperatures\nfound in the output files, while the other columns refer to a temperature   of {0} K'''.format(temp))
            print('Free-energy Temperature: {0} K'.format(temp))        
            print('|{0:^14}|{1:^29}|{2:^29}|{3:^29}|'.format(' ','GULP', 'Sum over modes', 'Sum over smoothed modes'))
            print('|{0:^14}|{1:^14}|{2:^14}|{3:^14}|{4:^14}|{5:^14}|{6:^14}|'.format('Volume (A^3)','F (eV)','Epot (eV)', 'Fvib (eV)','F (eV)', 'Fvib (eV)', 'F (eV)'))
            print('{:-^106}'.format('-'))
            
            for j in range(len(v_list)):
                F_thermal = 0.0
                F_thermal_smooth = 0.0
                for k in range(total_number_of_kpoints):
                    for f in range(total_number_of_modes):
                        F_thermal += kpoints_weights[k] * Kb_eV_per_K*temp * np.exp(-(h_eV_s*c_cm_per_s*data[j][k][f])/(Kb_eV_per_K*temp))
                        
                        F_thermal_smooth += kpoints_weights[k]*Kb_eV_per_K*temp*np.exp(-(h_eV_s*c_cm_per_s*data_smooth[j][k][f])/(Kb_eV_per_K*temp))
                
            
                print('{0:^14f} {1:^14f} {2:^14f} {3:^14f} {4:^14f} {5:^14f} {6:^14f}'.format(v_list[j], F_gulp_list[j], Epot_gulp_list[j], ZPE_list[j]+F_thermal,ZPE_list[j]+F_thermal+Epot_gulp_list[j], ZPE_smooth_list[j]+F_thermal_smooth, ZPE_smooth_list[j]+F_thermal_smooth+Epot_gulp_list[j]))
            print('{:-^106}\n'.format('-'))
            
            
            
            