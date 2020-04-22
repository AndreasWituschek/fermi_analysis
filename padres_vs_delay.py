# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: Andi

THis program takes the delay, i0, and the average over certain peak areas from 
the ion and electron time of flight traces and saves them to another hd5 file, 
so that they can be processed later on in order to check for saturation.
"""
import numpy as np
import os
import errno
import h5py

'''Which beamtime?'''
beamtime = 2

"""Experimental Parameters"""
delay_zero_pos = 11025.66

runs = [31, 693] #[893,1209]#[4651,5020]#[1164,1424] #first and last run of a delay scan
run_remove = [226,227,292,293,294,308] #[1079,1080,1081,1078,1153,1161,1178,1205,1206]#[4964]#[1250]
run_missing_i0 = [3044] # runs that are missing iom_sh_a

run_list = list(range(runs[0], runs[1] + 1))
for rr in run_remove:
    if rr in run_list:
        run_list.remove(rr)

""" Analysis parameter """
#scan31
window = [30,160]


""" File paths """
# path where it saves the preanalysed data
analyse_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime1/combined/'
if beamtime==2:
    analyse_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
analyse_path += 'scan_{:03}/'.format(runs[0])
# path where the raw data is: 
root_path = '//10.5.71.28/fermiserver/beamtime1/Day_3/'
if beamtime==2:
    root_path = '//10.5.71.28/fermiserver/beamtime2/'


try:
    os.makedirs(analyse_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

analyse_path += 'scan_{:03}_padres_spectra.h5'.format(runs[0])

try:
    analysis_file = h5py.File(analyse_path, 'r')
    run_numbers = np.array(analysis_file.get('run_numbers'))
    print(run_numbers)
    analysis_file.close()
    file_exists = True
except:
    run_numbers = np.array([])
    file_exists = False

print(file_exists)


# LDM data variables
delay = np.array([])
padres = np.array([])


# Run numbers
new_run_numbers = np.array([])


for run in run_list:
    if run in run_numbers:
        print('Skip run {}'.format(run))
        continue
    else:
        new_run_numbers = np.append(new_run_numbers, int(run))
    print('Analysing run {}'.format(int(run)))
# Generate run path
    run_path = root_path + 'Run_{:03}/'.format(run)

# Handle LDM data
    ldm_file_path = run_path + 'rawdata/'
    ldm_delay = np.array([])
    ldm_padres = np.array([])
    for run_file in os.listdir(ldm_file_path):
        ldm_file = h5py.File(ldm_file_path + run_file, 'r')
        #reading delay from files
        ldm_delay = np.append(ldm_delay, np.array(ldm_file['/photon_source/SeedLaser/trls_sl_03_pos']))

        ldm_padres = np.append(ldm_padres, np.mean(np.array(ldm_file['/photon_diagnostics/Spectrometer/hor_spectrum'])[window[0]:window[1],::],axis=0))
        
        
        ldm_file.close()
#    ldm_padres = np.reshape(ldm_padres,(2,-1))
    delay = np.append(delay, ldm_delay - delay_zero_pos)
    padres = np.append(padres, ldm_padres)

padres = padres.reshape((len(ldm_delay)*len(run_list),-1))


if file_exists:
    if not new_run_numbers.size == 0:
#        print('run_numbers: {}'.format(new_run_numbers))
        analysis_file = h5py.File(analyse_path, 'a')
        # Append run_numbers
        analysis_file['run_numbers'].resize((analysis_file['run_numbers'].shape[0] + new_run_numbers.shape[0]), axis = 0)
        analysis_file['run_numbers'][-new_run_numbers.shape[0]:] = new_run_numbers

        # Append LDM data
        ldm_group = analysis_file.get('LDM')
        # Delay
        ldm_group['delay'].resize((ldm_group['delay'].shape[0] + delay.shape[0]), axis = 0)
        ldm_group['delay'][-delay.shape[0]:] = delay


else:
    analysis_file = h5py.File(analyse_path, 'w')
    # Create LDM dataset
    ldm_group = analysis_file.create_group('LDM')
    # Delay
    ldm_group.create_dataset('delay', data=delay, maxshape=(None,))
    # I0
    ldm_group.create_dataset('padres', data=padres)


    # Run numbers
    analysis_file.create_dataset('run_numbers', data=new_run_numbers, maxshape=(None,))
    analysis_file.create_dataset('removed_run_numbers', data=np.asarray(run_remove), maxshape=(None,))


if not new_run_numbers.size == 0:
    analysis_file.close()
