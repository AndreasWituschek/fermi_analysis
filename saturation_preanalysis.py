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
import matplotlib.pyplot as plt

'''Which beamtime?'''
beamtime = 2

"""Experimental Parameters"""
delay_zero_pos = 11025.66
mfli_name = ['dev3265', 'dev3269']  # number of MFLI demodulator to be analyzed

runs = [4651, 5019]#693] #first and last run of a delay scan
run_remove = [4964]#[5251,5272,5328]#[2981,1079,1080,1081,1078,1153,1161,1178,1205,1206]#[5251,5272,5328] # #[1248] [4964] #[226, 227, 292, 293, 294, 308] #[1385, 1387, 1413, 1420, 6052,6053,6054,6055,6056]# #[3991] #3874] #[3367, 3368,3471,3472,3486,3485,3605,3606,3607,3608,3609,3610,3611,3612,3613] #[2386,2504,2505,2506,2507,2508,2509,2510,2511,2529] #[1575,1576,1577] #[1153, 1161, 1178, 1205, 1206]
run_missing_i0 = [3044] # runs that are missing iom_sh_a

run_list = range(runs[0], runs[1] + 1)
for rr in run_remove:
    if rr in run_list:
        run_list.remove(rr)

""" Analysis parameter """
cut_frac = 0.0   # constant fraction of data points cut from MFLI data vector, starting at 0 -> [0:cut_frac*length] is cut
transfer_ref = beamtime==2 # decides wether reference wavelenght is written into hd5 files. 

#windows for TOF traces
##scan893
#ion_window = [8520,8630]
#elc_window = [7900,8450]

#scan4651
ion_window = [15600,16400]
elc_window = [5900,6200]

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

analyse_path += 'scan_{:03}_saturation_preanalysis.h5'.format(runs[0])

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
i0 = np.array([])
ion = np.array([])
elc = np.array([])


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
    ldm_i0 = np.array([])
    ldm_ion = np.array([])
    ldm_elc = np.array([])
    for run_file in os.listdir(ldm_file_path):
        ldm_file = h5py.File(ldm_file_path + run_file, 'r')
        #reading delay from files
        ldm_delay = np.append(ldm_delay, np.array(ldm_file['/photon_source/SeedLaser/trls_sl_03_pos']))
        
        #reading i0 from files
        if run in run_missing_i0:
            ldm_i0 = np.zeros(ldm_delay.shape)
        else:
            ldm_i0 = np.append(ldm_i0, np.array(ldm_file['/photon_diagnostics/FEL01/I0_monitor/iom_uh_a']))
        
#        reading ion and electron spectrum from files:
#        raw_ion = np.mean(np.array(ldm_file['/digitizer/channel1']),axis=0)
#        raw_elc = np.mean(np.array(ldm_file['/digitizer/channel3']),axis=0)
        ldm_ion = np.append(ldm_ion, np.mean(np.array(ldm_file['/digitizer/channel1'])[::,ion_window[0]:ion_window[1]],axis=1))
        ldm_elc = np.append(ldm_elc, np.mean(np.array(ldm_file['/digitizer/channel3'])[::,elc_window[0]:elc_window[1]],axis=1))                
        
        ldm_file.close()
    delay = np.append(delay, ldm_delay[1::2] - delay_zero_pos)
    i0 = np.append(i0, ldm_i0)
    ion = np.append(ion, ldm_ion)
    elc = np.append(elc, ldm_elc)


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

        # I0
        ldm_group['I0'].resize((ldm_group['I0'].shape[0] + i0.shape[0]), axis = 0)
        ldm_group['I0'][-i0.shape[0]:] = i0
        # ion average
        ldm_group['s_I0'].resize((ldm_group['ion_avg'].shape[0] + ion.shape[0]), axis = 0)
        ldm_group['s_I0'][-ion.shape[0]:] = ion
        # electron average
        ldm_group['l_seed'].resize((ldm_group['elc_avg'].shape[0] + elc.shape[0]), axis = 0)
        ldm_group['l_seed'][-elc.shape[0]:] = elc

else:
    analysis_file = h5py.File(analyse_path, 'w')
    # Create LDM dataset
    ldm_group = analysis_file.create_group('LDM')
    # Delay
    ldm_group.create_dataset('delay', data=delay, maxshape=(None,))
    # I0
    ldm_group.create_dataset('I0', data=i0, maxshape=(None,))
    # ion avgerage
    ldm_group.create_dataset('ion_avg', data=ion, maxshape=(None,))
    # electron average
    ldm_group.create_dataset('elc_avg', data=elc, maxshape=(None,))
    # l_ref
#    if transfer_ref:
#        ldm_group.create_dataset('l_ref', data=ldm_l_ref_m, maxshape=(None,))
#    # i0 fft
#    ldm_group.create_dataset('i0_fft', data=ldm_i0_fft, maxshape=(None,None))

    # Run numbers
    analysis_file.create_dataset('run_numbers', data=new_run_numbers, maxshape=(None,))
    analysis_file.create_dataset('removed_run_numbers', data=np.asarray(run_remove), maxshape=(None,))


if not new_run_numbers.size == 0:
    analysis_file.close()

#execfile("c_file.py")