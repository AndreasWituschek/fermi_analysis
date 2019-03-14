# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:22:20 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

i = 150
run = 217
delay_zero_pos = 11025.66

ldm_file_path = '//online4ldm.esce.elettra.trieste.it/store/20149020/Day_2/Run_' + str(run) +'/rawdata/'.format(int(run))

ldm_file = os.listdir(ldm_file_path)[0]

def gaus(x,amp,x0,sigma):
    return amp*np.exp(-(x-x0)**2/(2*sigma**2))

def padres_fit_function(x,amp,x0,sigma,freq, phase, c):
    return gaus(x,amp,x0,sigma) * (1-np.sin(2*np.pi+(freq * x + phase))) + c


ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')

padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
print('padres_wavelength: {}'.format(padres_wavelength))

padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
print('padres_span: {}'.format(padres_span))

ldm_padres = np.array(ldm_data['photon_diagnostics']['Spectrometer']['hor_spectrum'])
print('ldm_padres shape: {}'.format(ldm_padres.shape))


ldm_data.close()

ldm_w_padres = ldm_padres[i][400:600]

#est_amp = np.max(ldm_w_padres) * 3
est_amp = 199829.
#est_x0 = np.where(ldm_w_padres == np.max(ldm_w_padres))[0][0]
#print(np.where(ldm_w_padres == np.max(ldm_w_padres))[0][0])
est_x0 = 136.5
est_sigma = 12.
est_freq = 0.2
est_phase = 1.1
est_c = 8000

p0=[est_amp, est_x0, est_sigma, est_freq, est_phase, est_c]

popt,pcov = curve_fit(padres_fit_function,np.arange(ldm_w_padres.size),ldm_w_padres,p0=p0)

fig, ax = plt.subplots(1,1)
ax.plot(np.arange(ldm_w_padres.size),ldm_w_padres)
ax.plot(np.arange(ldm_w_padres.size), padres_fit_function(np.arange(ldm_w_padres.size), *popt))