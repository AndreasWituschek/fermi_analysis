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


"""Experimental Parameters"""
# parameters for theoretical curve
l_He = 52.22  # helium 1s-4p transition in nm
l_ref = 266.  # reference laser wavelenght in nm
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66

""" Parameters for theoretical curve """
phi = 0.  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
offset = 0.  # offset

""" Analzsis paramter """
cut_frac = 0.25 # constant fraction of data points cut from MFLI data vector, starting at 0 -> [0:cut_frac*length] is cut


mfli_x_m = np.array([])
mfli_y_m = np.array([])
mfli_x_s = np.array([])
mfli_y_s = np.array([])
ldm_delay_m = np.array([])
ldm_delay_s = np.array([])
ldm_i0_m = np.array([])
ldm_i0_s = np.array([])

for run in range(59, 79):
#for run in range(11, 12):
    ldm_file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/Run_{:03}/rawdata/'.format(int(run))
    mfli_file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/Run_{:03}/work/mfli/'.format(int(run))

    ldm_files = os.listdir(ldm_file_path)
    mfli_files = os.listdir(mfli_file_path)

    mfli_x = np.array([])
    mfli_y = np.array([])
    ldm_delay = np.array([])
    ldm_i0 = np.array([])

    for mfli_file, ldm_file in zip(mfli_files, ldm_files):
        mfli_data = fk.ReadMfliData(mfli_file_path + mfli_file)
        mfli_x = np.append(mfli_x, mfli_data['x2'])
        mfli_y = np.append(mfli_y, mfli_data['y2'])

        ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
        ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
        ldm_i0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a'])
#    mfli_data_cutindex = int(round(cut_frac * np.shape(mfli_x)[0]))
#    mfli_x_m = np.append(mfli_x_m, np.mean(mfli_x[mfli_data_cutindex :]))
#    mfli_y_m = np.append(mfli_y_m, np.mean(mfli_y[mfli_data_cutindex :]))
#    mfli_x_s = np.append(mfli_x_s, np.std(mfli_x[mfli_data_cutindex :]))
#    mfli_y_s = np.append(mfli_y_s, np.std(mfli_y[mfli_data_cutindex :]))
    mfli_x_m = np.append(mfli_x_m, np.mean(mfli_x))
    mfli_y_m = np.append(mfli_y_m, np.mean(mfli_y))
    mfli_x_s = np.append(mfli_x_s, np.std(mfli_x))
    mfli_y_s = np.append(mfli_y_s, np.std(mfli_y))
    ldm_delay_m = np.append(ldm_delay_m, ldm_delay - delay_zero_pos)
    ldm_i0_m = np.append(ldm_i0_m, np.mean(ldm_i0))
    ldm_i0_s = np.append(ldm_i0_s, np.std(ldm_i0))

#length = ldm_delay_m.size
#
#Xtd, Ytd, X, Y = fk.Curve(l_He, l_ref, harmonic, phi, A,
#                          offset, start, stop, length)
#
#print(ldm_delay_m)
#plt.plot(Xtd, Xtd, 'b-')
fig, axs = plt.subplots(2, 1)
axs[0].errorbar(ldm_delay_m, mfli_x_m, yerr=mfli_x_s, color='k', linestyle='')
axs[0].plot(ldm_delay_m, mfli_x_m, 'ko-')
axs[0].errorbar(ldm_delay_m, mfli_y_m, yerr=mfli_x_s, color='r', linestyle='')
axs[0].plot(ldm_delay_m, mfli_y_m, 'ro-')
axs[0].set_ylabel('Amplitude in a.u.')

axs[1].errorbar(ldm_delay_m, ldm_i0_m, yerr=ldm_i0_s, color='k', linestyle='')
axs[1].plot(ldm_delay_m, ldm_i0_m, 'ko-')
axs[1].set_xlabel('Delay in fs')
axs[1].set_ylabel('i0')
plt.grid()
plt.show()