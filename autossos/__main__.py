#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 00:48:46 2020

@author: john
"""


import os
import sys
import shutil

from astropy.io import fits


if len(sys.argv) == 1:
    print("IT WORKS!")  
    sys.exit()


if sys.argv[1] in ['-d', '--default']:
    
    # Estructure example
    print('Creating example directory in: \n ' + os.getcwd())
    
    os.makedirs(os.path.join(os.getcwd(), 'multiple_nights', 'night1'), exist_ok=True)
    os.makedirs(os.path.join(os.getcwd(), 'multiple_nights', 'night2'), exist_ok=True)
    os.makedirs(os.path.join(os.getcwd(), 'SSOS'), exist_ok=True)
    
    
    config_files = ['config.ssos_calib', 'config.sex', 'config.scamp', 'default.param']
    path_to_module = os.path.dirname(__file__)
    for config in config_files:
        shutil.copy(os.path.join(path_to_module, config), os.path.join(os.getcwd(), config))
    
    os.chdir(os.path.join(os.getcwd(), 'SSOS'))
    os.system('ssos -d')
    
    conf = os.path.join(os.getcwd(), 'semp', 'ssos.scamp')
    with open(conf, 'r') as file:
        lines = file.readlines()
        lines = [i.replace('Y', 'N', 1) if "SOLVE_ASTROM  " == i[0:14] or "MATCH " == i[0:6] else i for i in lines]
    
    with open(conf, 'w') as file:
        file.writelines(lines)
    
    sys.exit()


if sys.argv[1] in ['-plot_sc', '--plot_sigmaclipping']:
    from ssos_calib.magcalib import sigmclip_plot
    sigmclip_plot(int(sys.argv[2])-1, int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]))
    sys.exit()


if sys.argv[1] in ['-check', '--check_calibrations']:
    from ssos_calib.filab import check
    check()
    sys.exit()


if sys.argv[1] in ['-period', '--periodogram']:
    from ssos_calib.magcalib import getperiod
    getperiod(sys.argv[2], sys.argv[3], sys.argv[4])
    sys.exit()


if sys.argv[1] in ['-re_sc', '--redo_sigmaclipping']:
    from ssos_calib import magcalib
    magcalib.sgmclip(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    magcalib.mag_calib()
    magcalib.ast_calib(sys.argv[6])
    sys.exit()



settings = {}
with open(os.path.join(os.getcwd(), 'config.ssos_calib'), 'r') as set_up:
    for line in set_up:
        if not line[0].isalpha():
            continue
        param, val, *_ = line.split()
        settings[param] = val
    
    locals().update(settings)

if sys.argv[1] in ['-SSOS']:
    from ssos_calib import magcalib
    magcalib.imgstats()
    magcalib.exec_ssos()
    magcalib.matchcats()
    magcalib.sgmclip(float(SIGM), float(MIN_MAG), float(MAX_MAG), float(MAX_MAGERR))
    magcalib.mag_calib()
    magcalib.ast_calib(R2)
    sys.exit()


from ssos_calib import filab
from ssos_calib import magcalib

main_dir = os.path.abspath(sys.argv[1]) 
os.chdir(main_dir)

def main():
    
    input_arguments = set(sys.argv)
    
    # Check nights for duplicates
    filab.duplicates()
    
    
    # Filabres execution
    skip_mainfil_arg = set(['-skip_filabres'])
    skip_mainfil = any(skip_mainfil_arg.intersection(input_arguments))
    
    if not skip_mainfil:
            filab.reduction(FLIPSTAT, INSTRUMENT)
    
    
    # Filabres calibration check
    outliers, notsolved, medianseps = filab.check()


    # Filabres recalibration of outliers
    skip_rec_arg = set(['-skip_recalibration'])
    skip_rec = any(skip_rec_arg.intersection(input_arguments))
    
    if not skip_rec:
        outliers, medianseps = filab.recalib(outliers, RECALIB)
    
    
    
    # Ask for SSOS recalibration
    filab.night_sepmed(medianseps)
    
    
    # Copy of reduced images without outliers
    
    ignored = outliers.copy()
    bad_nights = set(list(outliers.keys()) + list(notsolved.keys()))
    for n in bad_nights:
        if n in outliers.keys() and n in notsolved.keys():
            ignored[n] = outliers[n] + notsolved[n]
        elif n not in d1.keys():
            ignored[n] = notsolved[n]
    
    print("Copying the final images to SSOS folder...")
    filab.fil2ssos(ignored)
    
    
    # Check FILTER keyword in header
    print("Checking FILTER keyword in header...")
    sci_list = sorted([f for f in os.listdir(os.path.join(os.getcwd(), 'SSOS')) if ".fits" in f.lower()] )
    for i in sci_list:
        if 'FILTER' not in fits.getheader(os.path.join(os.getcwd(), 'SSOS', i)):
            fits.setval(os.path.join(os.getcwd(), 'SSOS', i), 'FILTER', value = 'NONE')    
    
    #####
    
    # Images statistics
    magcalib.imgstats()
    
    
    # SSOS execution
    magcalib.exec_ssos()

    # Matched catalogues    
    magcalib.matchcats()
    

    # Sigma clipping    
    magcalib.sgmclip(float(SIGM), float(MIN_MAG), float(MAX_MAG), float(MAX_MAGERR))
    
    
    # Magnitude calibration
    magcalib.mag_calib()
    
    
    # Asteroid light curve
    magcalib.ast_calib(R2)
    
    
    
if __name__ == '__main__':
    main()
