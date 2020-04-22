# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

run = 6774 #run with 30degree phase steps in 7th harmonic
delay_zero_pos = 11025.66


padres_spec = [] #padres spectrum of all harmonics(5,7,10)
padres_roi = []  #padres spectrum for all harmonics(5,7,10)

sl_spec = [] #seed spectrum of all runs
sl_roi = [] #seed ROI of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files
sl = np.zeros((0, 2068)) #seed laser spectum


ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_' + str(run) +'/rawdata/'.format(int(run))
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
    ldm_padres = ldm_padres/np.max(ldm_padres,axis=1).reshape(400,1).astype(float)
    padres = np.append(padres, ldm_padres,axis=0) 
    #padres spectrometer metadata
    padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
    padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
    padres_pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])
    
    #seed laser spectrum
    ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)
    ldm_sl = ldm_sl/np.max(ldm_sl,axis=1).reshape(400,1) #normation of seed laser specturm
    print(sl.shape,ldm_sl.shape)
#            sl = np.append(sl, ldm_sl,axis=0) 
    sl = np.append(sl, ldm_sl,axis=0) 
    #Seed laser metadata
    ldm_roi = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['RegionOfInterest']) # ROI for seed laser spectrum
    ldm_data.close

padres_lambda = padres_wavelength + (np.arange(1,1001)-500)*padres_pixel2micron*1E-3*padres_span #padres spectrometer calibration
padres_spec.append(padres)
padres_roi.append(padres_lambda)
    
#sl_lambda= np.linspace(ldm_roi[0],ldm_roi[1],ldm_sl.shape[1]) # this seed laser spectrometer calibration gives wrong values, apparently something wrong with DAQ
sl_lambda = np.arange(ldm_roi[0],395,0.06818072)-80 # this calibration is a post mortem correction of the upper one with the help of the spectra Paolo sent me from the Control room



''' Proper Data windows for all harmonics '''

shots= 8000

#seed laser 
sl_low=1220
sl_high=1330
sl_spec = np.transpose(sl)[sl_low:sl_high,0:shots]
sl_roi = sl_lambda[sl_low:sl_high]


#data window 7H run 6774
padres_low=200
padres_high=310
padres_spec_7H = np.transpose(padres_spec[0])[padres_low:padres_high,0:shots]
padres_roi_7H = padres_roi[0][padres_low:padres_high]


''' Plotting '''
ticksize= 3.0
ticklength = 8.
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['axes.linewidth'] = ticksize

fig = plt.figure(figsize=(10,8))

#seed laser spectrum
ax = fig.add_subplot(211)
ax.pcolormesh(np.arange(sl_spec.shape[1]),sl_roi, sl_spec)
#ax.set_ylabel(r'Seed $\lambda$ [nm]')
ax.ticklabel_format(useOffset=False)
ax.tick_params(length=ticklength, width = ticksize)
ax.set_xlim([0,shots-1])
ax.set_ylim([259.,263.])
ax.set_xticklabels([])
ax.set_yticks([260.,261.,262.])
ax.set_yticklabels(['260.0','261.0','262.0'])


#7H spectrum
ax = fig.add_subplot(212)
ax.pcolormesh(np.arange(padres_spec_7H.shape[1]), padres_roi_7H, padres_spec_7H)
ax.tick_params(length=ticklength, width = ticksize)
#ax.set_ylabel(r'7H $\lambda$ [nm]')
ax.set_xlabel('laser shots')
ax.set_xlim([0,shots-1])
ax.set_ylim([37.2,37.4])
#ax.set_xticklabels([])
ax.set_yticks([37.25,37.30,37.35])
ax.set_yticklabels(['37.25','37.30','37.35'])


