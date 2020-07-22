#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    SSOS related functions
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
import autossos.GLS as GLS

from tqdm import tqdm
from astropy import units as u
from astropy.nddata import CCDData

from astropy.table import Table
from astroML.crossmatch import crossmatch_angular
from astropy.io import fits

from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats


def sys_exec(command):
    test = os.system(command)
    
    # Stop pipeline if finished with any error
    if test != 0:
        sys.exit()

def imgstats():
    
    ''' Get statistics for the backgrounds of the images to remove the ones with bad photometry '''
    
    # Nights directories
    fil_folders = ['bias', 'flat-imaging', 'lists', 'science-imaging', 'SSOS']
    raw_folder = [f for f in sorted(os.listdir(os.getcwd())) if os.path.isdir(os.path.join(os.getcwd(), f)) and f not in fil_folders][0]
    raw_dir = os.path.join(os.getcwd(), raw_folder)                                        
    nights = [n for n in sorted(os.listdir(raw_dir)) if os.path.isdir(os.path.join(raw_dir, n))] 
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS')    
    
    # List of all the images
    lista = sorted([f for f in os.listdir(ssos_path) if ".fits" in f.lower()])
    
    # Check if the statistics have already been computed for those images
    if os.path.isfile(os.path.join(ssos_path, 'Stats.dat')):
        Stats = pd.read_csv(os.path.join(ssos_path, 'Stats.dat'))
        if list(set(Stats['Name']).symmetric_difference(set(lista))) == []:
            return
    

    img_mean =[]; img_sd = []; img_median = []; img_mad = []; names = []
    for i in tqdm(lista, desc = "Creating image statistics"):
        names.append(i)
        
        data = CCDData.read(os.path.join(ssos_path,i), unit = u.adu)
        
        mask = make_source_mask(data.data, nsigma=2, npixels=5, dilate_size=11)
        mean, median, std = sigma_clipped_stats(data.data, sigma=3.0, mask=mask)

        img_mean.append(mean)
        img_sd.append(std)
        img_median.append(median)
        img_mad.append(np.median(abs(data.data-median)))

    Stats = pd.DataFrame(data = {'Name': names , 'Mean': img_mean, 'SD': img_sd, 'Median': img_median, 'MAD': img_mad})
    Stats['Outlier'] = False
    Stats['Night'] = 'nan'
   
    # Mark as True bad photometry images using the outliers of a boxplot
    global_indx = []
    for n in nights:
        img_in_night = sorted(os.listdir(os.path.join(raw_dir, n)))
        Stats_night_indx = [list(Stats['Name']).index(i) for i in list(Stats['Name']) if any(os.path.splitext(s)[0] in i for s in img_in_night)]
        
        Stats.loc[Stats_night_indx, 'Night'] = n
        
        
        
        for st in range(1,len(Stats.columns)-1):
            dat = [float("%.50f" % elem) for elem in list(Stats.iloc[Stats_night_indx,st])]
            
            fl_mean = plt.boxplot(dat, showfliers=True)['fliers']; plt.close()
            fl_mean_data = [i.get_data() for i in fl_mean]
            
            indx = []
            for i in [float("%.50f" % elem) for elem in fl_mean_data[0][1]]:
                indx.append(list(dat).index(i))
            
            global_indx.extend(list(map(Stats_night_indx.__getitem__, indx)))
            
    # Save statistics data
    Stats.loc[list(set(global_indx)), 'Outlier'] = True
    Stats.to_csv(os.path.join(ssos_path, "Stats.dat"), index=False, header=True)
    

def exec_ssos():
    
    ''' Main SSOS execution '''
    
    os.makedirs(os.path.join(os.getcwd(), 'SSOS', 'cats','old'), exist_ok=True)
    old_results = [old for old in os.listdir(os.path.join(os.getcwd(), 'SSOS', 'cats')) if 'ssos_' in old]
    for o in old_results:
        shutil.move(os.path.join(os.getcwd(), 'SSOS', 'cats', o), os.path.join(os.getcwd(), 'SSOS', 'cats', 'old', o))
    
    cmd = "(cd " + os.path.join(os.getcwd(), 'SSOS') + " && exec ssos .)"
    sys_exec(cmd)
    
    # SSOS inspection    
    inpt = {'yes','y', 'ye', ''}
    while True:    
        choice = input('Inspection of the candidates. \n - LEFT ARROW: Artifact \n - UP ARROW: Unknown/unclear \n - RIGHT ARROW: Asteroid \n Press <ENTER> to continue. \n').lower()
        if choice in inpt:
            cmd = "ssos --inspect " + os.path.join(os.getcwd(),'SSOS')
            sys_exec(cmd)
            break

def matchcats():
    
    ''' Create matched catalogues with the reference stars of each catalogue'''
    
    print("Matching catalogues with reference stars...")
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS') 
    catspath = os.path.join(ssos_path, 'cats')
    
    # Remove old matched catalogues    
    os.system('rm ' + os.path.join(catspath, '*match_*'))
    
    # Match full catalogues
    full_cats = [f for f in sorted(os.listdir(catspath)) if f.startswith("full_")]
    for full in full_cats:
        
        full_match = pd.DataFrame()
        
        
        hdul = fits.open(os.path.join(catspath, full))
        dat_df = pd.DataFrame(np.array(hdul[2].data).byteswap().newbyteorder())
        
        
        reff_data = dat_df[dat_df['CATALOG_NUMBER'] == 0]
        
        reff = np.empty((len(reff_data['DELTA_J2000']), 2), dtype=np.float64)
        reff[:, 0] = reff_data['ALPHA_J2000']
        reff[:, 1] = reff_data['DELTA_J2000']
    
         
        num_cats = max(dat_df['CATALOG_NUMBER'])
        meansep = [];mediansep=[]
        for i in range(1,num_cats+1):
            
            catt_data = dat_df[dat_df['CATALOG_NUMBER'] == i]
            
            catt = np.empty((len(catt_data['ALPHA_J2000']), 2), dtype=np.float64)
            catt[:, 0] = catt_data['ALPHA_J2000']
            catt[:, 1] = catt_data['DELTA_J2000']
            
            
            max_radius = 1. / 3600  # 1 arcsec
            dist, ind = crossmatch_angular(catt, reff, max_radius)
            match = ~np.isinf(dist)
            
            dist_match = dist[match]
            dist_match *= 3600
            
            meansep.append(np.mean(dist_match))
            mediansep.append(np.median(dist_match))
            
            ind_match = ind[match]   
            
            
            sep = pd.DataFrame(data=dist_match, columns=['Sep'])
            
            reff_match = reff_data.iloc[ind_match,:].reset_index(drop=True)  
            reff_match.columns = [col.lower() for col in reff_match.columns]        
            cattg_match = catt_data[match].reset_index(drop=True)
            
            result = pd.concat([reff_match,cattg_match, sep], axis=1)
            full_match = full_match.append(result).reset_index(drop=True)
        
        fullmatch_tb = Table.from_pandas(full_match)
        
        name = 'match_' + full.split('.')[0]
        fullmatch_tb.write(os.path.join(catspath, name + '.fits'))
        os.rename(os.path.join(catspath, name + '.fits'), os.path.join(catspath, name + '.cat'))


def sigm_clipping(cat_MAG, sigm, min_mag, max_mag, max_magerr):
    
    ''' Sigma clipping of a single catalogue '''
    
    # FILTER MASKS
    ## Remove 99 values
    mask_99 = cat_MAG['maginst']!=99
    
    ## Remove magnitudes with error greater than threshold
    mask_magerr_thres = cat_MAG['mag_err'] < max_magerr
    mask_maginsterr_thres = cat_MAG['maginst_err'] < max_magerr
    
    ## Remove detection limits
    mask_maglim_max = cat_MAG['mag'] < max_mag
    mask_maglim_min = cat_MAG['mag'] > min_mag
    
    # Final mask
    mask_lin = mask_99 & mask_magerr_thres & mask_maginsterr_thres & mask_maglim_max & mask_maglim_min
    
    # Data for linear region
    x=np.array(cat_MAG['maginst'][mask_lin])
    y=np.array(cat_MAG['mag'][mask_lin])
    
    iters = 0
    len_valid = 0
    len_x = len(x)
    
    # If there are less that 4 values it will return "0" values
    if len_x < 4:
        return(len(cat_MAG), 0, 0, 0, 0, 0, cat_MAG['epoch'][0], len(x), mask_lin, x, y, 'nan', 0)

    # Iterates until all remaining points are inside the threshold
    while len_valid < len_x:
        coef, V = np.polyfit(x, y, 1, cov = True)
        poly1d_fn = np.poly1d(coef) 
        
        diff = abs(y - poly1d_fn(x))
        
        threshold = np.mean(diff) + sigm * np.std(diff)
        
        valid = diff < threshold
        
        len_x = len(x)
        
        x = x[valid]
        y = y[valid]
        len_valid = len(valid[valid==True])
        iters = iters +1

        
        
    SSreg = np.sum( (y - poly1d_fn(x))**2 )
    SStot = np.sum( (y - np.mean(y))**2 )
    
    R2 = 1-SSreg/SStot
    
    
    return(len(cat_MAG), coef[1], np.sqrt(V[1][1]), coef[0], np.sqrt(V[0][0]), R2, cat_MAG['epoch'][0], len(x), mask_lin, x, y, poly1d_fn, iters)


def sgmclip(sigm, min_mag, max_mag, max_magerr):
    
    ''' Sigma clipping of all matched catalogues '''
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS') 
    catspath = os.path.join(ssos_path, 'cats')
    
    # Matched catalogues
    matchs = [f for f in sorted(os.listdir(catspath)) if f.startswith("match_")]
    
    intercept = []; slope = []; CoefReg = []; Epoch = []; catalognumber = []; Nsourc = []
    intercept_err = []; slope_err = []; Nfit = []
    for match_cat in matchs:
        
        print('Sigma clipping catalogues in ' + match_cat + ' ...')
        
        current_epoch = []; current_R2 = []
        
        with fits.open(os.path.join(catspath, match_cat)) as hdul:
            
            merg_cat = hdul[1].data

            mag = list(merg_cat['mag'])
            mag_err = list(merg_cat['magerr'])
            maginst = list(merg_cat['MAG'])
            maginst_err = list(merg_cat['MAGERR'])
            epoch = list(merg_cat['EPOCH'])
            
            MAGS = pd.DataFrame([mag, mag_err, maginst, maginst_err, epoch]).T
            MAGS.columns = ['mag', 'mag_err', 'maginst', 'maginst_err', 'epoch']

            num_cats = list(set(merg_cat['CATALOG_NUMBER']))
            
            for i in tqdm(num_cats):
                
                cat_mask = merg_cat['CATALOG_NUMBER'] == i
                cat_MAG = MAGS[cat_mask].reset_index(drop=True)
                
                
                Coeff = sigm_clipping(cat_MAG, sigm, min_mag, max_mag, max_magerr)
                
                
                catalognumber.append(i)
                
                Nsourc.append(Coeff[0])
                intercept.append(Coeff[1])
                intercept_err.append(Coeff[2])
                slope.append(Coeff[3])
                slope_err.append(Coeff[4])
                CoefReg.append(Coeff[5])
                Epoch.append(Coeff[6])
                Nfit.append(Coeff[7])
                
                current_epoch.append(Coeff[6])
                current_R2.append(Coeff[5])
            
            # R2 evolution plot
            plt.figure()
            plt.plot(current_epoch, current_R2, label = str(sigm) + ' sigma')
            plt.title(match_cat)
            plt.xlabel('Epoch')
            plt.ylabel('R2')
            plt.grid(linewidth = 0.5)
            plt.legend()
            plt.savefig(os.path.join(ssos_path, 'checkplots', 'R2_ev_'+ match_cat +'.png'))
    
    # Store coefficients data
    Coeff_df = pd.DataFrame(data = {'Catalog Number': catalognumber, 'Intercept': intercept, 'Int_err': intercept_err, 'Slope': slope, 'Slope_err': slope_err, 'R2': CoefReg, 'N tot': Nsourc,'N fit': Nfit, 'Epoch': Epoch})
    Coeff_df.to_csv(os.path.join(catspath, 'coeff_fit.csv'), sep = ',', index = False)

        
            
def sigmclip_plot(full, cat, sigm, min_mag, max_mag, max_magerr):
    
    ''' Sigma clipping plot of a single catalogue '''
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS') 
    catspath = os.path.join(ssos_path, 'cats')    
    
    match_cat = [f for f in os.listdir(catspath) if f.startswith("match_")][full]            

    with fits.open(os.path.join(catspath, match_cat)) as hdul:

        merg_cat = hdul[1].data
        mag = list(merg_cat['mag'])
        mag_err = list(merg_cat['magerr'])
        maginst = list(merg_cat['MAG'])
        maginst_err = list(merg_cat['MAGERR'])
        epoch = list(merg_cat['EPOCH'])
        
        MAGS = pd.DataFrame([mag, mag_err, maginst, maginst_err, epoch]).T
        MAGS.columns = ['mag', 'mag_err', 'maginst', 'maginst_err', 'epoch']
        
        numb = cat
        
        cat_mask_ind = merg_cat['CATALOG_NUMBER'] == numb

        cat_MAG_ind = MAGS[cat_mask_ind].reset_index(drop=True)

        Coeff = sigm_clipping(cat_MAG_ind, sigm, min_mag, max_mag, max_magerr)
        
        if Coeff[12] == 0:
            print("Sigma clipping returns 0 values. Please check sigma clipping parameters")
            sys.exit()
            
        plt.figure()
        plt.scatter(cat_MAG_ind['maginst'][-Coeff[8]], cat_MAG_ind['mag'][-Coeff[8]], s = 1, color = 'red')
        plt.scatter(cat_MAG_ind['maginst'][Coeff[8]], cat_MAG_ind['mag'][Coeff[8]], s = 1, color = 'orange')
        plt.scatter(Coeff[9], Coeff[10], s = 1)
         
        plt.plot(np.linspace(min(cat_MAG_ind['maginst']),max(cat_MAG_ind['maginst']),1000), Coeff[11](np.linspace(min(cat_MAG_ind['maginst']),max(cat_MAG_ind['maginst']),1000)), '-', color='black',
                  label= f'y = {Coeff[1]} + {Coeff[3]} * x \n R2 = {Coeff[5]} \n sigm = {sigm} \n nº iters = {Coeff[12]} \n nº fit = {len(Coeff[9])} \n nº tot_sourc = {Coeff[0]}')
        
        plt.title(f"Catalog nº {numb}, with {sigm} sigm")
        plt.xlabel('Mag Instr')
        plt.ylabel('Mag Ref')
        plt.legend()
        plt.grid()
        plt.show()


def mag_calib():
    
    ''' Obtain calibrated magnitudes for the SSOS results '''
    
    # Avoid warnings when working on slices of dataframes
    pd.options.mode.chained_assignment = None
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS') 
    catspath = os.path.join(ssos_path, 'cats')
    
    # Coefficients data and SSOS results
    coeff_csv = [f for f in os.listdir(catspath) if f.startswith("coeff_")][0]
    ssos_csv = [f for f in os.listdir(catspath) if f.startswith("ssos_")][0]
    
    Coeff = pd.read_csv(os.path.join(catspath, coeff_csv))
    ssos_data = pd.read_csv(os.path.join(catspath, ssos_csv))
    
    ssos_mag = ssos_data['MAG']
    
    ssos_mag_cal = np.zeros(len(ssos_mag)); ssos_R2 = np.zeros(len(ssos_mag)); ssos_mag_cal_err = np.zeros(len(ssos_mag))
    for n in ssos_data['CATALOG_NUMBER']:
        Cpos = Coeff['Catalog Number'] == n
        Spos = ssos_data['CATALOG_NUMBER'] == n
        
        ssos_mag_cal[Spos] = np.array(Coeff['Intercept'][Cpos]) + np.array(Coeff['Slope'][Cpos]) * np.array(ssos_data['MAG'][Spos])
        ssos_R2[Spos] = Coeff['R2'][Cpos]
        ssos_mag_cal_err[Spos] = np.sqrt( (np.array(ssos_data['MAG'][Spos])*np.array(Coeff['Slope_err'][Cpos]))**2 + (np.array(Coeff['Int_err'][Cpos]))**2 + (np.array(Coeff['Slope'][Cpos])*np.array(ssos_data['MAGERR'][Spos]))**2 )
    
    
    ssos_data['MAG_CALIB'] = ssos_mag_cal
    ssos_data['MAG_CALIB_ERR'] = ssos_mag_cal_err
    ssos_data['R2'] = ssos_R2
   
    # Creates a column with True for outliers
    ssos_data['PHOT_OUTLIER'] = False
    stats_data = pd.read_csv(os.path.join(ssos_path, 'Stats.dat'))
    bad_phot = list(stats_data['Name'][stats_data['Outlier'] == True])
    for bad in bad_phot:
        ssos_data['PHOT_OUTLIER'][ssos_data['IMAGE_FILENAME'] == bad] = True
    
    # Assign "names" to possible asteroides without name or unknowns
    ssos_data_aux = ssos_data[ssos_data['MATCHED'] == False]
    
    ## For unknowns
    ssos_data_unk = ssos_data_aux[ssos_data_aux['UNKNOWN'] == True]
    ssos_unk = ssos_data_unk
    catnumb_unk = ssos_unk['CATALOG_NUMBER']
    catnumb_unk = catnumb_unk.append(pd.DataFrame([-1]))
    
    j = 0; n=1
    for i in range(len(catnumb_unk)-1):
        if any(catnumb_unk.iloc[i+1] > catnumb_unk.iloc[i]):
            continue
        else:
            indx = list(ssos_unk.iloc[j:i+1].index)
            ssos_data['SKYBOT_NAME'][indx] = 'UNK_'+str(n)
            j = i+1; n += 1
    
    ## For asteroids
    ssos_data_ast = ssos_data_aux[ssos_data_aux['ASTEROID'] == True]
    ssos_ast = ssos_data_ast
    catnumb_ast = ssos_ast['CATALOG_NUMBER']
    catnumb_ast = catnumb_ast.append(pd.DataFrame([-1]))    
    
    j = 0; n=1
    for i in range(len(catnumb_ast)-1):
        if any(catnumb_ast.iloc[i+1] > catnumb_ast.iloc[i]):
            continue
        else:
            indx = list(ssos_ast.iloc[j:i+1].index)
            ssos_data['SKYBOT_NAME'][indx] = 'AST_'+str(n)
            j = i+1; n += 1

    
    ssos_data.to_csv(os.path.join(catspath, 'calib_' + ssos_csv), sep = ',', index = False)

def ast_calib(R2):
    
    ''' Create data only with valid data without photometry outliers and R2 within threshold '''
    
    ssos_path = os.path.join(os.getcwd(), 'SSOS') 
    catspath = os.path.join(ssos_path, 'cats')
    
    img_stats = sorted([f for f in os.listdir(ssos_path) if "Stats" in f])[0]
    
    Stats = pd.read_csv(os.path.join(ssos_path, img_stats))
    ssos_calib_name = sorted([f for f in os.listdir(catspath) if f.startswith("calib_")])[0]
    ssos_data = pd.read_csv(os.path.join(catspath, ssos_calib_name))
    
    # Bad photometry images
    img_notvalid_phot = sorted(list(set(Stats['Name'][Stats['Outlier'] == True])))
    
    # Bad R2 images
    img_notvalid_R2 = sorted(list(set(ssos_data['IMAGE_FILENAME'][ssos_data['R2'] < float(R2)])))
    
    # Not valid images
    img_notvalid = sorted(list(set(img_notvalid_R2 + img_notvalid_phot)))
    
    
    ssos_data_index = [i for i in range(len(ssos_data['IMAGE_FILENAME'])) if ssos_data['IMAGE_FILENAME'][i] not in img_notvalid]
    ssos_data_valid = ssos_data.iloc[ssos_data_index,:]
    
    # Creates light curves with valid data
    ast = [i for i in list(set(ssos_data['SKYBOT_NAME'])) if str(i) != 'nan']
    ast_data = []
    for i in range(len(ast)):
        ast_data.append(ssos_data_valid[ssos_data_valid['SKYBOT_NAME'] == ast[i]])
                
            
        plt.figure()
        plt.title(ast[i])
        plt.plot(ast_data[i]['EPOCH'], ast_data[i]['MAG_CALIB'], 'o', markersize = 1.5)
        plt.xlabel('Epoch')
        plt.ylabel('Mag Calib')
        plt.grid()
        plt.savefig(os.path.join(ssos_path, 'checkplots', 'mag_' +str(ast[i]) + '.png'))
    
    # Save valid data
    ssos_data_valid.to_csv(os.path.join(catspath, 'valid' + ssos_calib_name), sep = ',', index = False)

    # Create plot with image statistics and not valid data
    Stats_nights = list(set(Stats['Night']))
    for n in Stats_nights:
        
        dat = Stats[Stats['Night'] == n]
        
        plt.figure()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        indx_dat_badR2 = [list(dat.Name).index(j) for j in img_notvalid_R2 if j in list(dat.Name)]
        indx_dat_nodetect = [list(dat.Name).index(j) for j in list(dat.Name) if j not in list(ssos_data['IMAGE_FILENAME'])]

        for i in range(1,5):
            
            plt.subplot(2,2,i)
            plt.plot(dat.Name, dat.iloc[:,i])
            plt.scatter(dat.Name[dat.Outlier == True], dat[dat.Outlier == True].iloc[:,i], color = 'red', marker = "x", label = 'Photometry outlier')
            plt.scatter(dat.iloc[indx_dat_badR2,0], dat.iloc[indx_dat_badR2,i], color = 'orange', marker = "o", s = 10, label = 'R2 < '+ R2)
            plt.scatter(dat.iloc[indx_dat_nodetect,0], dat.iloc[indx_dat_nodetect,i], color = 'green', marker = "o", s = 10, label = 'No detections')
            plt.title(dat.columns[i]); plt.xlabel('Images'); plt.ylabel('ADUs');plt.grid(axis='y')
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        
        plt.legend()
        plt.savefig(os.path.join(ssos_path, 'checkplots', 'bad_phot_' + str(n) + '.png'))
        



def getperiod(file, name, iterations):
    
    ''' Use GLS to obtain periods '''
    
    file_dir = os.path.abspath(file)
    
    data_all = pd.read_csv(file_dir)
    data_ast = data_all[data_all['SKYBOT_NAME'] == name]
    
    time = list(data_ast['EPOCH']*365.25*24)
    mag = list(data_ast['MAG_CALIB'])
    
    lines = []
    for i in range(1, int(iterations)+1):
        if i == 1:
            yres = mag
        else:
            yfit = periodograma.sinmod()
            yres = periodograma.y-yfit
        
        periodograma = GLS.Gls((time,yres))
        periodograma.plot(period=True,fap=(0.1,0.01,0.001), block=True)
        
        best_period_res = periodograma.best['P']
        period_error_res = periodograma.best['e_P']
        fap_res = periodograma.FAP()
        
        line1 = 'Best period: {0} +/- {1} days'.format(best_period_res,period_error_res)
        line2 = 'False alarm probability: {0}%'.format(fap_res*100)
       
        lines.extend(['\n Iteration:', str(i),'\n', line1, '\n', line2])
        print('\n', line1, '\n', line2)
        
    # Creates a file to store data
    f = open(os.path.join(os.getcwd(), name + "_period.txt"), "w")
    for l in lines:
        f.write(l)
