#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    Filabres related functions
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import warnings

from astroML.crossmatch import crossmatch_angular
from tqdm import tqdm
from pprint import pprint
from distutils.util import strtobool


def sys_exec(command):
    test = os.system(command)
    
    # Stop pipeline if finished with any error
    if test != 0:
        sys.exit()

        
def duplicates():
    
    ''' Find if there is any duplicated file inside nights folders '''
    
    print("Searching for duplicated files...")
    
    # Nights directories
    fil_folders = ['bias', 'flat-imaging', 'lists', 'science-imaging', 'SSOS']
    raw_folder = [f for f in sorted(os.listdir(os.getcwd())) if os.path.isdir(os.path.join(os.getcwd(), f)) and f not in fil_folders][0]
    raw_dir = os.path.join(os.getcwd(), raw_folder)                                        
    nights = [n for n in sorted(os.listdir(raw_dir)) if os.path.isdir(os.path.join(raw_dir, n))]          
    
    # Images names inside each night
    imgs  = []
    for n in nights:
        night_imgs = [i for i in sorted(os.listdir(os.path.join(raw_dir,n))) if 'bias' not in i.lower() and 'flat' not in i.lower() and (i.endswith('.fts') or i.endswith('.fits'))]
        imgs.append(night_imgs)
    
    # Find duplicates
    dupl = []
    for i in range(len(imgs)):
        for j in range(i+1, len(imgs)):
            dp = any([True for i in imgs[i] if i in imgs[j]])
            dupl.append(dp)
            if dp == True:
                print('There are duplicated files in: ' + nights[i] +'\t and \t'+ nights[j])
    
    # User choice if there are duplicates. Exit or append night name.
    if any(dupl):
        while True:
            yes = {'y','ye', 'yes'}
            no = {'n', 'no'}
            
            choice = input("\n PAUSED \nCan't proceed with duplicated files, please change the names. This can be done automatically. \nDo you want to append the folder name to the images files? 'No' will stop the pipeline. (yes/no): ").lower()
            if choice in yes:
                for root, dirs, files in os.walk(raw_dir):
                    
                    files = [i for i in files if 'bias' not in i.lower() and 'flat' not in i.lower() and (i.endswith('.fts') or i.endswith('.fits'))]
                    if not files:
                        continue
                    prefix = os.path.basename(root)
                    for f in files:
                        os.rename(os.path.join(root, f), os.path.join(root, "{}_{}".format(prefix,f)))
                        
                break
            elif choice in no:
                print('Stopping the pipeline...')
                sys.exit()
            else:
               sys.stdout.write("Please respond with 'yes' or 'no.")
    else:
        print('\n No duplicates found, continuing...\n')


def reduction(flipstat, instrument):
    
    ''' Filabres reduction. Reduce the images and add WCS header '''
    
    # True if there is a meridian flip
    flipstat = bool(strtobool(flipstat))
    
    # Nights directories    
    fil_folders = ['bias', 'flat-imaging', 'lists', 'science-imaging', 'SSOS']
    raw_folder = [f for f in sorted(os.listdir(os.getcwd())) if os.path.isdir(os.path.join(os.getcwd(), f)) and f not in fil_folders][0]
    raw_dir = os.path.join(os.getcwd(), raw_folder)                                        
    nights = [n for n in sorted(os.listdir(raw_dir)) if os.path.isdir(os.path.join(raw_dir, n))]          

    # Create flipped BIAS and FLAT if flipstat = True
    for n in nights:
        
        # Change .fts to .fits
        if any([True for i in os.listdir(os.path.join(raw_dir, n)) if i.endswith('.fts')]):
            cmd = "for f in " + os.path.join(raw_dir,n,'*.fts') + "; do mv $f ${f%.*}.fits; done"
            sys_exec(cmd)
        
        masterbias = [f for f in os.listdir(os.path.join(raw_dir, n)) if 'bias' in f.lower() and not f.endswith('r.fits')]
        masterflat = [f for f in os.listdir(os.path.join(raw_dir, n)) if 'flat' in f.lower() and not f.endswith('r.fits')]
                    
        if masterbias == [] or masterflat == []:
            print('Need reduction images in the night: ' + str(n))
            sys.exit()
        
        if (flipstat):
            if not os.path.isfile(masterbias[0].replace('.fits', 'r.fits')):
                cmd_bias_flip = 'filabres-rotate_flipstat ' + os.path.join(raw_dir, n, masterbias[0])
                sys_exec(cmd_bias_flip)
            if not os.path.isfile(masterflat[0].replace('.fits', 'r.fits')):
                cmd_flat_flip = 'filabres-rotate_flipstat ' + os.path.join(raw_dir, n, masterflat[0])
                sys_exec(cmd_flat_flip)
    
    
    # Remove old *.yaml
    os.system('rm ' + os.path.join(os.getcwd(),'*.yaml'))
    
    # Filabres setup
    cmd_1 = "filabres --setup " + instrument + " "+ raw_folder
    sys_exec(cmd_1)
    
    # Configurable files for filabres 
    setup_filabres = os.path.join(os.getcwd(), 'setup_filabres.yaml')
    with open(setup_filabres, 'r') as file:
        lines = file.readlines()
        lines.append('\n# sextractor and scamp files\n' +
                     'default_param: default.param\n' +
                     'config_sex: config.sex\n' +
                     'config_scamp: config.scamp\n')
    with open(setup_filabres, 'w') as file:
        file.writelines(lines)
    
    
    # Filabres initialize
    cmd_init = "filabres -rs initialize -v"
    sys_exec(cmd_init)

    # Filabres BIAS and FLAT reduction
    sys_exec("filabres -rs bias -v")
    sys_exec("filabres -rs flat-imaging -v")
    
    # Retry if Filabres runs out of memory
    retry = 35072
    while True:
        if retry == 35072:
            cmd = "filabres -rs science-imaging"
            retry = os.system(cmd)
        elif retry == 0:
            break
        else:
            sys.exit()
    
    
    
def check():
    
    ''' Check if any median separations with the reference stars is an outlier '''
    
    science_path = os.path.join(os.getcwd(), 'science-imaging')
    nights = [n for n in sorted(os.listdir(science_path)) if os.path.isdir(os.path.join(science_path, n))]          
    
    outliers = {}       # Median outliers
    notsolved = {}      # Not solved by Astrometry
    medianseps = {}     # Quartiles of separations
    
    # Check outliers inside each night using a boxplot
    for n in tqdm(nights, desc='Checking calibrations'):
        img_path = os.path.join(science_path, n)
        
        name = []; mediansep = []; notsol = []
        sci_path = [sci for sci in sorted(os.listdir(img_path)) if  os.path.isdir(os.path.join(img_path,sci)) and sci.startswith('science-imaging')]
        for science in sci_path:
            folder = os.path.join(science_path, n, science)

            if not os.path.isfile(os.path.join(folder, 'xxx.new')):
                notsol.append(science)
                continue
            else:
                dat_df = pd.read_csv(os.path.join(folder, 'full_1.cat'), delim_whitespace=True, header=None, comment = '#')
                
                ref_data = dat_df[dat_df[1] == 0]
                ref = np.empty((len(ref_data), 2), dtype=np.float64)
                ref[:, 0] = ref_data[10]
                ref[:, 1] = ref_data[11]
                
                sci_data = dat_df[dat_df[1] == 1]
                sci = np.empty((len(sci_data), 2), dtype=np.float64)
                sci[:, 0] = sci_data[10]
                sci[:, 1] = sci_data[11]
                
                
                max_radius = 1. / 3600  # 1 arcsec
                dist, ind = crossmatch_angular(sci, ref, max_radius)
                match = ~np.isinf(dist)
                
                dist_match = dist[match]
                dist_match *= 3600
                
                name.append(science)
                
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    mediansep.append(np.median(dist_match))


        sep = pd.DataFrame([name,mediansep]).T
        
        
        out =  ['NaN '+ str(i) for i in list(sep[0][np.isnan(list(sep[1]))])]
        
        sep_aux = sep.dropna().reset_index(drop=True)
        sep_aux[1] = [float("%.50f" % elem) for elem in list(sep_aux[1])]

        fl = plt.boxplot(sep_aux[1], showfliers=True); plt.close()
        fl_data = [i.get_data() for i in fl['fliers']]

        for i in [float("%.50f" % elem) for elem in fl_data[0][1]]:
            if i > fl['medians'][0].get_ydata()[0]:
                out.append(str(i) + ' '+ str(sep_aux[0][list(sep_aux[1]).index(i)]))
        
        
        outliers[n] = out
        notsolved[n] = notsol
        medianseps[n] = 'Median_of_seps: ' + str(fl['medians'][0].get_ydata()[0]) + ' | Q1: ' + str(fl['whiskers'][0].get_ydata()[0]) + ' | Q3: ' + str(fl['whiskers'][1].get_ydata()[0])
    
    if any(notsolved.values()):
        print("\n Astrometry couldn't find astrometric solution for:")
        NS = pd.DataFrame(notsolved.items(), columns=['Night', 'Images'])
        NS.to_csv(os.path.join(os.getcwd(), "list_notsolved.csv"), index=False, header=True)    
        pprint(notsolved, width=1)
    else:
        print(" - Astrometry solved all the images.")
    

    if any(outliers.values()):
        print('\n Bad astrometric solutions found in:')
        OL = pd.DataFrame(outliers.items(), columns=['Night', 'Images'])
        OL.to_csv(os.path.join(os.getcwd(), "list_outliers.csv"), index=False, header=True)
        pprint(outliers, width=1)
    else:
        print(" - There are no outliers.")
    
    pprint(medianseps, width=1)
    
    return [outliers, notsolved, medianseps]


def recalib(outliers, medianseps, RECALIB):

    ''' Recalibrate outliers if RECALIB = True. Runs Filabres individualy for each outlier '''    

    retry = True
    while retry:

        if bool(strtobool(RECALIB)) and any(outliers.values()):
            
            out = [item.replace('science-imaging_','').replace('_red','') +'.fits' for sublist in list(outliers.values()) for item in sublist]
            
            for o in tqdm(out, desc='\n Recalibrating'):
                cmd = "filabres -rs science-imaging --filename " + o.split()[1] + " --force -ng"
                sys_exec(cmd)
            
            outliers, *_, medianseps = check()
    
            while True:
                len_out = len([j for i in list(outliers.values()) for j in i])

                if len_out != 0:    
                    yes = {'yes','y', 'ye'}
                    no = {'no','n', ''}
                    ignore = {'ignore','ign'}
                    
                    choice = input('PAUSED.\n Detected outliers in mean separations between reference and extracted stars after individual recalibration of those images.\n You can edit now "config.scamp" parameters and try again the recalibration. \n Please respond:\n - "EXIT" to stop the pipeline.\n - "ignore" to use SSOS with the outliers\n - "yes" to try again the calibration.\n - "no" or <ENTER> to continue the pipeline removing the outliers from SSOS.\n').lower()
                    if choice in yes:
                        break
                    elif choice in no:
                        retry = False
                        break
                    elif choice in ignore:
                        outliers = {}
                        retry = False
                        break
                    elif choice == 'exit':
                        sys.exit()
                    else:
                       sys.stdout.write("Please respond with 'yes', 'no', 'ignore' or 'EXIT'")
                
                else:
                    retry = False
                    break
        else:
            retry = False
    
    return [outliers, medianseps]


def night_sepmed(medianseps):
    
    ''' Prints separation quartiles for each night and ask for a SSOS recalibration with scamp '''
    
    print('\n Median of separations for each night:')
    for n in medianseps.keys():
        print('\n' + n + '-> '+ medianseps[n])

    while True:
        yes = {'y','ye', 'yes'}
        no = {'n', 'no'}
        cont = {''}
        
        conf = os.path.join(os.getcwd(), 'SSOS', 'semp', 'ssos.scamp')
        
        choice = input("\n If you want to try a recalibration with SSOS respond (yes), if not respond (no) or press <ENTER> to continue with current SSOS parameters. \n ").lower()
        if choice in yes:
            with open(conf, 'r') as file:
                lines = file.readlines()
                lines = [i.replace('N', 'Y', 1) if (("SOLVE_ASTROM  " == i[0:14] or "MATCH " == i[0:6]) and i.split()[1] == 'N') else i for i in lines]
            
            with open(conf, 'w') as file:
                file.writelines(lines)
            
            break
        
        elif choice in no:        
            with open(conf, 'r') as file:
                lines = file.readlines()
                lines = [i.replace('Y', 'N', 1) if (("SOLVE_ASTROM  " == i[0:14] or "MATCH " == i[0:6]) and i.split()[1] == 'Y') else i for i in lines]
            
            with open(conf, 'w') as file:
                file.writelines(lines)
            
            break
        
        elif choice in cont:
            break
        
        else:
           sys.stdout.write("Please respond with a valid input.")
    

def fil2ssos(outliers):
    
    ''' Copy all calibrated images to SSOS folder except for the ones marked as outliers '''
    
    out = [item.split(' ')[-1] +'.fits' for sublist in list(outliers.values()) for item in sublist]
    science_path = os.path.join(os.getcwd(), 'science-imaging')
    nights = [n for n in sorted(os.listdir(science_path)) if os.path.isdir(os.path.join(science_path, n))]
    
    imgs_dir = []
    for n in nights:
        img_path = os.path.join(science_path, n)
        imgs_night = [i for i in sorted(os.listdir(img_path)) if i.endswith('.fits') and i not in out]
        imgs_dir.extend([os.path.join(img_path, i) for i in imgs_night])
    
    cp_path = os.path.join(os.getcwd(), 'SSOS')
    
    # Remove old images inside SSOS folder
    if any([True for i in os.listdir(cp_path) if '.fits' in i]):
        print('Removing older images...\n')
        os.system('rm ' + os.path.join(cp_path, '*.fits'))
    
    for i in tqdm(imgs_dir, desc = 'Copying images'):
        cmd = "cp " + i + " " + cp_path
        sys_exec(cmd)
