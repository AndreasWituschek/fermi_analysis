# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: ldm

02/04/19: 
- Changed i0 from iom_sh_a to iom_uh_a
- ref laser wavelength goes to hd5 files
- removed runs go to hd5 files

10/04/19:
- added FFT of I0 to hd5 files
"""
import numpy as np
import os
import errno
import h5py

'''Which beamtime?'''
beamtime = 2

"""Experimental Parameters"""
delay_zero_pos = 11025.66
mfli_name = ['dev3265', 'dev3269']  # number of MFLI demodulator to be analyzed

runs = [1164,1424] #first and last run of a delay scan
run_remove = [1250]#[5251,5272,5328]#[2981,1079,1080,1081,1078,1153,1161,1178,1205,1206]#[5251,5272,5328] #[226, 227, 292, 293, 294, 308] #[1248] [4964] #[226, 227, 292, 293, 294, 308] #[1385, 1387, 1413, 1420, 6052,6053,6054,6055,6056]# #[3991] #3874] #[3367, 3368,3471,3472,3486,3485,3605,3606,3607,3608,3609,3610,3611,3612,3613] #[2386,2504,2505,2506,2507,2508,2509,2510,2511,2529] #[1575,1576,1577] #[1153, 1161, 1178, 1205, 1206]
run_missing_i0 = [3044] # runs that are missing iom_sh_a

run_list = range(runs[0], runs[1] + 1)
for rr in run_remove:
    if rr in run_list:
        run_list.remove(rr)

""" Analysis parameter """
cut_frac = 0.2   # constant fraction of data points cut from MFLI data vector, starting at 0 -> [0:cut_frac*length] is cut
transfer_ref = beamtime==2 # decides wether reference wavelenght is written into hd5 files. 

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
#root_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/'
#root_path = '/home/ldm/ExperimentalData/Online4LDM/RBT-UOF_4/Data/'
#root_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_4/'
#/home/ldm/ExperimentalData/Online4LDM/20149020/results/PM_Data/PM_Data_Day5/Late/scan_1834
#/home/ldm/ExperimentalData/Online4LDM/20149020/results/PM_Data/PM_Data_Day5/Late
#smb://online4ldm.esce.elettra.trieste.it/store/20149020/Day_5/Run_2181/rawdata
#smb://online4ldm.esce.elettra.trieste.it/store/20149020/Day_5/Run_2181/rawdata/Run_2181_223236826.h5


try:
    os.makedirs(analyse_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

analyse_path += 'scan_{:03}.h5'.format(runs[0])

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
ldm_delay_m = np.array([])
ldm_delay_s = np.array([])
ldm_i0_m = np.array([])
ldm_i0_s = np.array([])
ldm_i0_fft = np.array([])
ldm_l_fel_m = np.array([])
ldm_l_ref_m = np.array([])


# Run numbers
new_run_numbers = np.array([])

# MFLI data variables
mfli_data = {
        'dev3265':{
                'x0': np.array([]),
                'y0': np.array([]),
                'x1': np.array([]),
                'y1': np.array([]),
                'x2': np.array([]),
                'y2': np.array([]),
                's_x0': np.array([]),
                's_y0': np.array([]),
                's_x1': np.array([]),
                's_y1': np.array([]),
                's_x2': np.array([]),
                's_y2': np.array([]),
                'harmonic': np.array([]),},
        'dev3269':{
                'x0': np.array([]),
                'y0': np.array([]),
                'x1': np.array([]),
                'y1': np.array([]),
                'x2': np.array([]),
                'y2': np.array([]),
                's_x0': np.array([]),
                's_y0': np.array([]),
                's_x1': np.array([]),
                's_y1': np.array([]),
                's_x2': np.array([]),
                's_y2': np.array([]),
                'harmonic': np.array([]),}}

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
    ldm_l_ref = np.array([])
    for run_file in os.listdir(ldm_file_path):
        ldm_file = h5py.File(ldm_file_path + run_file, 'r')
        ldm_delay = np.append(ldm_delay, np.array(ldm_file['/photon_source/SeedLaser/trls_sl_03_pos']))
        
        if run in run_missing_i0:
            ldm_i0 = np.zeros(ldm_delay.shape)
        else:
            ldm_i0 = np.append(ldm_i0, np.array(ldm_file['/photon_diagnostics/FEL01/I0_monitor/iom_uh_a']))
        ldm_l_fel = np.array(ldm_file['/photon_source/SeedLaser/Wavelength'])
        if transfer_ref:
            ldm_l_ref = np.array(ldm_file['/photon_source/Reference_Laser_UV_Spectrum/GaussFitCenter'])
        ldm_file.close()
    ldm_delay_m = np.append(ldm_delay_m, np.mean(ldm_delay) - delay_zero_pos)
    ldm_delay_s = np.append(ldm_delay_s, np.std(ldm_delay))
    ldm_i0_m = np.append(ldm_i0_m, np.mean(ldm_i0))
    ldm_i0_s = np.append(ldm_i0_s, np.std(ldm_i0))
    ldm_l_fel_m = np.append(ldm_l_fel_m, ldm_l_fel) 
    #i0 harmonics:
    ldm_i0_fft_raw = abs(np.fft.fftshift(np.fft.fft(ldm_i0)))
    ldm_i0_fft = np.append(ldm_i0_fft, ldm_i0_fft_raw)
    ldm_i0_fft = ldm_i0_fft.reshape((-1,int(ldm_i0_fft_raw.size)))
#    ldm_i0_fft_freq = np.fft.fftshift(np.fft.fftfreq(len(ldm_i0),1./50.))

    
            
    
    if transfer_ref:
        ldm_l_ref_m = np.append(ldm_l_ref_m, ldm_l_ref)
    
# Handle MFLI data
    for mfli in mfli_name:
        mfli_file_path = run_path +'work/' + mfli + '/'
        dir_list = os.listdir(mfli_file_path)
        run_file = dir_list[0]
        if run == runs[0]:
            mfli_file = h5py.File(mfli_file_path + run_file, 'r')
            mfli_data[mfli]['harmonic'] = np.array(mfli_file['harmonic'])
            mfli_file.close()

        # loop through demodulators
        for demod in range(3):
            # loop through files
            mfli_x = np.array([])
            mfli_y = np.array([])
            for run_file in dir_list:
                mfli_file = h5py.File(mfli_file_path + run_file, 'r')
                mfli_x = np.append(mfli_x, np.array(mfli_file['x' + str(demod)]))
                mfli_y = np.append(mfli_y, np.array(mfli_file['y' + str(demod)]))
                mfli_file.close()
            data_size = np.min([mfli_x.size, mfli_y.size])
            cut = int(cut_frac * data_size)
            mfli_data[mfli]['x' + str(demod)] = np.append(mfli_data[mfli]['x' + str(demod)], np.mean(mfli_x[cut:]))
            mfli_data[mfli]['s_x' + str(demod)] = np.append(mfli_data[mfli]['s_x' + str(demod)], np.std(mfli_x[cut:]))
            mfli_data[mfli]['y' + str(demod)] = np.append(mfli_data[mfli]['y' + str(demod)], np.mean(mfli_y[cut:]))
            mfli_data[mfli]['s_y' + str(demod)] = np.append(mfli_data[mfli]['s_y' + str(demod)], np.std(mfli_y[cut:]))


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
        ldm_group['delay'].resize((ldm_group['delay'].shape[0] + ldm_delay_m.shape[0]), axis = 0)
        ldm_group['delay'][-ldm_delay_m.shape[0]:] = ldm_delay_m
        # Std Delay
        ldm_group['s_delay'].resize((ldm_group['s_delay'].shape[0] + ldm_delay_s.shape[0]), axis = 0)
        ldm_group['s_delay'][-ldm_delay_s.shape[0]:] = ldm_delay_s
        # I0
        ldm_group['I0'].resize((ldm_group['I0'].shape[0] + ldm_i0_m.shape[0]), axis = 0)
        ldm_group['I0'][-ldm_i0_m.shape[0]:] = ldm_i0_m
        # Std I0
        ldm_group['s_I0'].resize((ldm_group['s_I0'].shape[0] + ldm_i0_s.shape[0]), axis = 0)
        ldm_group['s_I0'][-ldm_i0_s.shape[0]:] = ldm_i0_s
        # l_seed
        ldm_group['l_seed'].resize((ldm_group['l_seed'].shape[0] + ldm_l_fel_m.shape[0]), axis = 0)
        ldm_group['l_seed'][-ldm_l_fel_m.shape[0]:] = ldm_l_fel_m
        # l_ref
        if transfer_ref:
            ldm_group['l_ref'].resize((ldm_group['l_ref'].shape[0] + ldm_l_ref_m.shape[0]), axis = 0)
            ldm_group['l_ref'][-ldm_l_ref_m.shape[0]:] = ldm_l_ref_m

        # Append MFLI data
        for mfli in mfli_name:
            mfli_group = analysis_file.get(mfli)
            mfli_device = mfli_data[mfli]

            for demod in range(3):
                # X values
                mfli_group['x' + str(demod)].resize((mfli_group['x' + str(demod)].shape[0] + mfli_device['x' + str(demod)].shape[0]), axis = 0)
                mfli_group['x' + str(demod)][-mfli_device['x' + str(demod)].shape[0]:] = mfli_device['x' + str(demod)]
                mfli_group['s_x' + str(demod)].resize((mfli_group['s_x' + str(demod)].shape[0] + mfli_device['s_x' + str(demod)].shape[0]), axis = 0)
                mfli_group['s_x' + str(demod)][-mfli_device['s_x' + str(demod)].shape[0]:] = mfli_device['s_x' + str(demod)]
                # Y values
                mfli_group['y' + str(demod)].resize((mfli_group['y' + str(demod)].shape[0] + mfli_device['y' + str(demod)].shape[0]), axis = 0)
                mfli_group['y' + str(demod)][-mfli_device['y' + str(demod)].shape[0]:] = mfli_device['y' + str(demod)]
                mfli_group['s_y' + str(demod)].resize((mfli_group['s_y' + str(demod)].shape[0] + mfli_device['s_y' + str(demod)].shape[0]), axis = 0)
                mfli_group['s_y' + str(demod)][-mfli_device['s_y' + str(demod)].shape[0]:] = mfli_device['s_y' + str(demod)]
else:
    analysis_file = h5py.File(analyse_path, 'w')
    # Create LDM dataset
    ldm_group = analysis_file.create_group('LDM')
    # Delay
    ldm_group.create_dataset('delay', data=ldm_delay_m, maxshape=(None,))
    # Std Delay
    ldm_group.create_dataset('s_delay', data=ldm_delay_s, maxshape=(None,))
    # I0
    ldm_group.create_dataset('I0', data=ldm_i0_m, maxshape=(None,))
    # Std I0
    ldm_group.create_dataset('s_I0', data=ldm_i0_s, maxshape=(None,))
    # l_seed
    ldm_group.create_dataset('l_seed', data=ldm_l_fel_m, maxshape=(None,))
    # l_ref
    if transfer_ref:
        ldm_group.create_dataset('l_ref', data=ldm_l_ref_m, maxshape=(None,))
    # i0 fft
    ldm_group.create_dataset('i0_fft', data=ldm_i0_fft, maxshape=(None,None))

    # Run numbers
    analysis_file.create_dataset('run_numbers', data=new_run_numbers, maxshape=(None,))
    analysis_file.create_dataset('removed_run_numbers', data=np.asarray(run_remove), maxshape=(None,))

    # Append MFLI data
    for mfli in mfli_name:
        mfli_group = analysis_file.create_group(mfli)
        mfli_device = mfli_data[mfli]

        # Append harmonic
        mfli_group.create_dataset('harmonic', data=mfli_device['harmonic'], maxshape=(None,))

        for demod in range(3):
            # X values
            mfli_group.create_dataset('x' + str(demod), data=mfli_device['x' + str(demod)], maxshape=(None,))
            mfli_group.create_dataset('s_x' + str(demod), data=mfli_device['s_x' + str(demod)], maxshape=(None,))
            # Y values
            mfli_group.create_dataset('y' + str(demod), data=mfli_device['y' + str(demod)], maxshape=(None,))
            mfli_group.create_dataset('s_y' + str(demod), data=mfli_device['s_y' + str(demod)], maxshape=(None,))

if not new_run_numbers.size == 0:
    analysis_file.close()

execfile("c_file.py")