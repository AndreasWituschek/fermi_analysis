# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: ldm
"""
import numpy as np
import os
import fermi_analysis.functions as fk
import h5py
import matplotlib.pyplot as plt
import scipy.constants as spc
import pickle
import os


"""Experimental Parameters"""
# parameters for theoretical curve
l_He = 52.22  # helium 1s-4p transition in nm
l_ref = 265.98  # reference laser wavelenght in nm
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
mfli_name = ['DEV3265', 'DEV3269']  # MFLI device number
demod_harmonic = [['5_0', '5_1', '5_2'],['9_0', '9_1', '9_2']]  # demodulated harmonic corresponding to MFLI channels 0-2 (1-3 MFLI surface), first vector for DEV3265, second for DEV3269
demod_channel = [0, 1, 2]  # number of MFLI demodulator to be analyzed

Runs = [2397, 2404]
delays = [200, 400]


""" Parameters for theoretical curve """
phi = 0.  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
offset = 0.  # offset

""" Analysis parameter """
cut_frac = 0.2   # constant fraction of data points cut from MFLI data vector, starting at 0 -> [0:cut_frac*length] is cut

pickle_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/results/PM_Data/PM_Data_Day5/Night/'
#pickle_path += 'run_{0}-{1}_ion_{2}-{3}fs/'.format(Runs[0], Runs[1] - 1, delays[0], delays[1])
pickle_path += 'ion_33K_{0}-{1}fs/'.format(delays[0], delays[1])

# Create path if not exists
if not os.path.exists(pickle_path):
    os.mkdir(pickle_path)

for mfli in mfli_name:
    if mfli == 'DEV3265':
        demod_harmonic_dev = demod_harmonic[0]
    else:
        demod_harmonic_dev = demod_harmonic[1]
    for demod in demod_channel:
        mfli_x_m = np.array([])
        mfli_y_m = np.array([])
        mfli_x_s = np.array([])
        mfli_y_s = np.array([])
        ldm_delay_m = np.array([])
        ldm_delay_s = np.array([])
        ldm_i0_m = np.array([])
        ldm_i0_s = np.array([])
        demod_str_x = 'x' + str(demod)
        demod_str_y = 'y' + str(demod)


        run_list = range(Runs[0], Runs[1])
#        run_remove = [1786, 1768]
#        for rr in run_remove:
#            run_list.remove(rr)
        for run in run_list:  # Attention: range(11,13) executes run 11 and 12 (but not 13)
            ldm_file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_5/Run_' + str(run) +'/rawdata/'.format(int(run))
            mfli_file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_5/Run_' + str(run) +'/work/' + mfli + '/'.format(int(run))

            ldm_files = os.listdir(ldm_file_path)
            mfli_files = os.listdir(mfli_file_path)

            mfli_x = np.array([])
            mfli_y = np.array([])
            ldm_delay = np.array([])
            ldm_i0 = np.array([])

            for mfli_file, ldm_file in zip(mfli_files, ldm_files):
                mfli_data = fk.ReadMfliData(mfli_file_path + mfli_file)
                mfli_x = np.append(mfli_x, mfli_data[demod_str_x])
                mfli_y = np.append(mfli_y, mfli_data[demod_str_y])

                ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
                ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
                ldm_i0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a'])
            mfli_data_cutindex = int(round(cut_frac * np.shape(mfli_x)[0]))
            mfli_x_m = np.append(mfli_x_m, np.mean(mfli_x[mfli_data_cutindex :]))
            mfli_y_m = np.append(mfli_y_m, np.mean(mfli_y[mfli_data_cutindex :]))
            mfli_x_s = np.append(mfli_x_s, np.std(mfli_x[mfli_data_cutindex :]))
            mfli_y_s = np.append(mfli_y_s, np.std(mfli_y[mfli_data_cutindex :]))
            ldm_delay_m = np.append(ldm_delay_m, ldm_delay - delay_zero_pos)
            ldm_i0_m = np.append(ldm_i0_m, np.mean(ldm_i0))
            ldm_i0_s = np.append(ldm_i0_s, np.std(ldm_i0))
        print demod_harmonic_dev[demod] , mfli

        #length = ldm_delay_m.size
        #
        #Xtd, Ytd, X, Y = fk.Curve(l_He, l_ref, harmonic, phi, A,
        #                          offset, start, stop, length)
        #
        #print(ldm_delay_m)
        #plt.plot(Xtd, Xtd, 'b-')

        delay =[]
        X = []
        X_s = []
        Y = []
        Y_s = []
        i0 = []
        i0_s = []

#        if os.path.exists('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Night/delayscan_' + demod_harmonic_dev[demod]):
#            with open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Night/delayscan_' + demod_harmonic_dev[demod],'rb') as fp:
#        if os.path.exists(pickle_path + 'delayscan_run_{0}-{1}_ion_{2}-{3}fs_{4}'.format(Runs[0], Runs[1] - 1, delays[0], delays[1], demod_harmonic_dev[demod])):
#            with open(pickle_path + 'delayscan_run_{0}-{1}_ion_{2}-{3}fs_{4}'.format(Runs[0], Runs[1] - 1, delays[0], delays[1], demod_harmonic_dev[demod]),'rb') as fp:
        if os.path.exists(pickle_path + 'delayscan_ion_{0}-{1}fs_{2}'.format(delays[0], delays[1], demod_harmonic_dev[demod])):
            with open(pickle_path + 'delayscan_ion_{0}-{1}fs_{2}'.format(delays[0], delays[1], demod_harmonic_dev[demod]),'rb') as fp:
                delay = pickle.load(fp)
                X = pickle.load(fp)
                X_s = pickle.load(fp)
                Y = pickle.load(fp)
                Y_s = pickle.load(fp)
                i0 = pickle.load(fp)
                i0_s = pickle.load(fp)

        delay = np.concatenate((delay,ldm_delay_m),axis=None)
        X = np.concatenate((X,mfli_x_m),axis=None)
        X_s = np.concatenate((X_s,mfli_x_s),axis=None)
        Y = np.concatenate((Y,mfli_y_m),axis=None)
        Y_s = np.concatenate((Y_s,mfli_y_s),axis=None)
        i0 = np.concatenate((i0,ldm_i0_m),axis=None)
        i0_s = np.concatenate((i0_s,ldm_i0_s),axis=None)

#        with open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Night/delayscan_' + demod_harmonic_dev[demod],'wb') as fp:
#        with open(pickle_path + 'delayscan_run_{0}-{1}_ion_{2}-{3}fs_{4}'.format(Runs[0], Runs[1] - 1, delays[0], delays[1], demod_harmonic_dev[demod]),'wb') as fp:
        with open(pickle_path + 'delayscan_ion_{0}-{1}fs_{2}'.format(delays[0], delays[1], demod_harmonic_dev[demod]),'wb') as fp:
            pickle.dump(delay,fp)
            pickle.dump(X,fp)
            pickle.dump(X_s,fp)
            pickle.dump(Y,fp)
            pickle.dump(Y_s,fp)
            pickle.dump(i0,fp)
            pickle.dump(i0_s,fp)


#        fig, axs = plt.subplots(2, 1)
#
#        axs[0].errorbar(delay, X, yerr=X_s, color='k', linestyle='')
#        axs[0].plot(delay, X, 'ko-')
#        axs[0].errorbar(delay, Y, yerr=Y_s, color='r', linestyle='')
#        axs[0].plot(delay, Y, 'ro-')
#        axs[0].set_ylabel('Amplitude in a.u.')
#
#        axs[1].errorbar(delay, i0, yerr=i0_s, color='k', linestyle='')
#        axs[1].plot(delay, i0, 'ko-')
#        axs[1].set_xlabel('Delay in fs')
#        axs[1].set_ylabel('i0')
#
#
#        #axs[0].errorbar(ldm_delay_m, mfli_x_m, yerr=mfli_x_s, color='k', linestyle='')
#        #axs[0].plot(ldm_delay_m, mfli_x_m, 'ko-')
#        #axs[0].errorbar(ldm_delay_m, mfli_y_m, yerr=mfli_x_s, color='r', linestyle='')
#        #axs[0].plot(ldm_delay_m, mfli_y_m, 'ro-')
#        #axs[0].set_ylabel('Amplitude in a.u.')
#        #
#        #axs[1].errorbar(ldm_delay_m, ldm_i0_m, yerr=ldm_i0_s, color='k', linestyle='')
#        #axs[1].plot(ldm_delay_m, ldm_i0_m, 'ko-')
#        #axs[1].set_xlabel('Delay in fs')
#        #axs[1].set_ylabel('i0')
#plt.grid()
#plt.show()

print 'Runs: ' + str(Runs)