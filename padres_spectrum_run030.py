# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.signal as sig


run= 30
delay_zero_pos = 11025.66


padres_spec = [] #padres spectrum of all harmonics(5,7,10)
padres_roi = []  #padres spectrum for all harmonics(5,7,10)

''' Loop over all harmonics '''
ldm_file_path ='//10.5.71.28/FermiServer/Beamtime2/Run_030/rawdata/'
ldm_files = os.listdir(ldm_file_path)

ldm_padres = np.array([]) #raw padres spectrum from ldm files
padres = np.zeros((0, 1000))#padres spectrum

  
for ldm_file in ldm_files:
    ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
    ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
#    ldm_i0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a'])
#    i0 = np.append(i0, ldm_i0) # I0 of FEL
    
    #padres spectrum
    ldm_padres = np.array(ldm_data['photon_diagnostics']['Spectrometer']['hor_spectrum'])
#    ldm_padres = ldm_padres/np.max(ldm_padres,axis=1).reshape(400,1).astype(float)
    padres = np.append(padres, ldm_padres,axis=0) 
    #padres spectrometer metadata
    padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
    padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
    padres_pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])

padres_lambda = padres_wavelength + (np.arange(1,1001)-500)*padres_pixel2micron*1E-3*padres_span #padres spectrometer calibration
padres_spec.append(padres)
padres_roi.append(padres_lambda)



''' Proper Data windows for 5H '''
padres_low=100
padres_high=310
padres_spec_5H = np.transpose(padres_spec[0])[padres_low:padres_high,::]
padres_roi_5H = padres_roi[0][padres_low:padres_high]