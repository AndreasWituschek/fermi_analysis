# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: ldm
"""
import numpy as np
import os
import errno
import h5py


"""Experimental Parameters"""
# parameters for theoretical curve
delay_zero_pos = 11025.66
mfli_name = ['DEV3265', 'DEV3269']  # number of MFLI demodulator to be analyzed

runs = [1367, 1434]

run_list = range(runs[0], runs[1] + 1)

run_remove = [1385, 1387, 1413, 1420]
for rr in run_remove:
    run_list.remove(rr)

""" Analysis parameter """
cut_frac = 0.2   # constant fraction of data points cut from MFLI data vector, starting at 0 -> [0:cut_frac*length] is cut

""" File paths """
analyse_path = 'C:/Users/FemtoMeasure/Desktop/HeDroplets/'
analyse_path += 'scan_{}/'.format(runs[0])


root_path = '//online4ldm.esce.elettra.trieste.it/store/20149020/Day_4/'
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

analyse_path += 'scan_{}.h5'.format(runs[0])

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

# Run numbers
new_run_numbers = np.array([])

# MFLI data variables
mfli_data = {
        'DEV3265':{
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
                's_y2': np.array([]),},
        'DEV3269':{
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
                's_y2': np.array([]),}}

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
    for run_file in os.listdir(ldm_file_path):
        ldm_file = h5py.File(ldm_file_path + run_file, 'r')
        ldm_delay = np.append(ldm_delay, np.array(ldm_file['photon_source']['SeedLaser']['trls_sl_03_pos']))
        ldm_i0 = np.append(ldm_i0, np.array(ldm_file['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a']))
        ldm_file.close()
    ldm_delay_m = np.append(ldm_delay_m, np.mean(ldm_delay) - delay_zero_pos)
    ldm_delay_s = np.append(ldm_delay_s, np.std(ldm_delay))
    ldm_i0_m = np.append(ldm_i0_m, np.mean(ldm_i0))
    ldm_i0_s = np.append(ldm_i0_s, np.std(ldm_i0))

# Handle MFLI data
    for mfli in mfli_name:
        mfli_file_path = run_path +'work/' + mfli + '/'

        # loop through demodulators
        for demod in range(3):
            # loop through files
            mfli_x = np.array([])
            mfli_y = np.array([])
            for run_file in os.listdir(mfli_file_path):
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
        print('run_numbers: {}'.format(new_run_numbers))
        analysis_file = h5py.File(analyse_path, 'a')
#        analysis_file.create_dataset('run_numbers', data=np.arange(10))
        # Append LDM data
        # Delay
        analysis_file['delay'].resize((analysis_file['delay'].shape[0] + ldm_delay_m.shape[0]), axis = 0)
        analysis_file['delay'][-ldm_delay_m.shape[0]:] = ldm_delay_m
        # Std Delay
        analysis_file['s_delay'].resize((analysis_file['s_delay'].shape[0] + ldm_delay_s.shape[0]), axis = 0)
        analysis_file['s_delay'][-ldm_delay_s.shape[0]:] = ldm_delay_s
        # I0
        analysis_file['I0'].resize((analysis_file['I0'].shape[0] + ldm_i0_m.shape[0]), axis = 0)
        analysis_file['I0'][-ldm_i0_m.shape[0]:] = ldm_i0_m
        # Std I0
        analysis_file['s_I0'].resize((analysis_file['s_I0'].shape[0] + ldm_i0_s.shape[0]), axis = 0)
        analysis_file['s_I0'][-ldm_i0_s.shape[0]:] = ldm_i0_s

        # Append run_numbers
        analysis_file['run_numbers'].resize((analysis_file['run_numbers'].shape[0] + new_run_numbers.shape[0]), axis = 0)
        analysis_file['run_numbers'][-new_run_numbers.shape[0]:] = new_run_numbers


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
    # Delay
    analysis_file.create_dataset('delay', data=ldm_delay_m, maxshape=(None,))
    # Std Delay
    analysis_file.create_dataset('s_delay', data=ldm_delay_s, maxshape=(None,))
    # I0
    analysis_file.create_dataset('I0', data=ldm_i0_m, maxshape=(None,))
    # Std I0
    analysis_file.create_dataset('s_I0', data=ldm_i0_s, maxshape=(None,))

    # Run numbers
    analysis_file.create_dataset('run_numbers', data=new_run_numbers, maxshape=(None,))

    # Append MFLI data
    for mfli in mfli_name:
        mfli_group = analysis_file.create_group(mfli)
        mfli_device = mfli_data[mfli]
        for demod in range(3):
            # X values
            mfli_group.create_dataset('x' + str(demod), data=mfli_device['x' + str(demod)], maxshape=(None,))
            mfli_group.create_dataset('s_x' + str(demod), data=mfli_device['s_x' + str(demod)], maxshape=(None,))
            # Y values
            mfli_group.create_dataset('y' + str(demod), data=mfli_device['y' + str(demod)], maxshape=(None,))
            mfli_group.create_dataset('s_y' + str(demod), data=mfli_device['s_y' + str(demod)], maxshape=(None,))

analysis_file.close()

plt.plot(ldm_delay_m)
