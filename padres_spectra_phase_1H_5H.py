# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

runs = [6777,6773,6768] #[5H, 7H, 10H]
#runs = [6777,6774,6768] #[5H, 7H, 10H]
#run=3213
#run 6774 is with 30degree phase steps in 7th harmonic
delay_zero_pos = 11025.66


padres_spec = [] #padres spectrum of all harmonics(5,7,10)
padres_roi = []  #padres spectrum for all harmonics(5,7,10)

sl_spec = [] #seed spectrum of all runs
sl_roi = [] #seed ROI of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files
sl = np.zeros((0, 2068)) #seed laser spectum


''' Loop over all harmonics '''
for run in runs:
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
        ldm_padres = ldm_padres/np.mean(ldm_padres,axis=1).reshape(400,1).astype(float) #normation of padres spectrum
        padres = np.append(padres, ldm_padres,axis=0)  
        #padres spectrometer metadata
        padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
        padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
        padres_pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])
        
        #seed laser spectrum
        if run == 6777: #only analyze first run
            ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)
#            ldm_sl = ldm_sl/np.max(ldm_sl,axis=1).reshape(400,1) #normation of seed laser specturm
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
#seed laser 
sl_low=1220
sl_high=1330
sl_spec = np.transpose(sl)[sl_low:sl_high,0:600]
sl_roi = sl_lambda[sl_low:sl_high]


#data window 5H run 6777
padres_low=100
padres_high=310
padres_spec_5H = np.transpose(padres_spec[0])[padres_low:padres_high,0:600]
padres_roi_5H = padres_roi[0][padres_low:padres_high]


#data window 7H run 6773
padres_low=200
padres_high=310
padres_spec_7H = np.transpose(padres_spec[1])[padres_low:padres_high,0:600]
padres_roi_7H = padres_roi[1][padres_low:padres_high]

#data window 10H run 6769
padres_low=130
padres_high=280
padres_spec_10H = np.transpose(padres_spec[2])[padres_low:padres_high,1900:2500] #3400:4000
padres_roi_10H = padres_roi[2][padres_low:padres_high]



''' Plotting '''
#ticksize= 3.0
#ticklength = 8.
#fontsize=16.
#plt.rcParams['xtick.labelsize'] = fontsize
#plt.rcParams['ytick.labelsize'] = fontsize
#plt.rcParams['axes.labelsize'] = fontsize
#plt.rcParams['xtick.labelsize'] = 25
#plt.rcParams['ytick.labelsize'] = 25
#plt.rcParams['axes.linewidth'] = ticksize
#plt.rcParams['lines.linewidth'] = 2.



#plot params
ticksize= 2.
ticklength = 5.
fontsize=16.
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.major.width'] = ticksize
plt.rcParams['xtick.major.size'] = ticklength
plt.rcParams['ytick.major.width'] = ticksize
plt.rcParams['ytick.major.size'] = ticklength
plt.rcParams['axes.linewidth'] = ticksize
plt.rcParams['lines.linewidth'] = 2.



fig = plt.figure(figsize=(9,8))

#seed laser spectrum
ax = fig.add_subplot(411)
ax.pcolormesh(np.arange(sl_spec.shape[1]),sl_roi, sl_spec)
#ax.set_ylabel(r'Seed $\lambda$ [nm]')
ax.ticklabel_format(useOffset=False)
ax.tick_params(length=ticklength, width = ticksize)
ax.set_xlim([0,599])
ax.set_ylim([259.,263.])
ax.set_xticklabels([])
ax.set_yticks([260.,261.,262.])
ax.set_yticklabels(['260.0','261.0','262.0'])

#5H spectrum
ax = fig.add_subplot(412)
ax.pcolormesh(np.arange(padres_spec_5H.shape[1]), padres_roi_5H, padres_spec_5H) 
ax.tick_params(length=ticklength, width = ticksize)
#ax.set_ylabel(r'5H $\lambda$ [nm]')
ax.set_xlim([0,599])
ax.set_ylim([52.1,52.4])
ax.set_xticklabels([])
ax.set_yticks([52.15,52.25,52.35])
ax.set_yticklabels(['52.15','52.25','52.35'])


#7H spectrum
ax = fig.add_subplot(413)
ax.pcolormesh(np.arange(padres_spec_7H.shape[1]), padres_roi_7H, padres_spec_7H)
ax.tick_params(length=ticklength, width = ticksize)
#ax.set_ylabel(r'7H $\lambda$ [nm]')
#ax.set_xlabel('# Shots')
ax.set_xlim([0,599])
ax.set_ylim([37.2,37.4])
ax.set_xticklabels([])
ax.set_yticks([37.25,37.30,37.35])
ax.set_yticklabels(['37.25','37.30','37.35'])

#10H spectrum
ax = fig.add_subplot(414)
ax.pcolormesh(np.arange(padres_spec_10H.shape[1]), padres_roi_10H, padres_spec_10H)
ax.tick_params(length=ticklength, width = ticksize)
#ax.set_ylabel(r'10H $\lambda$ [nm]')
ax.set_xlabel('laser shot')
ax.set_xlim([0,599])
ax.set_ylim([26.06,26.16])
ax.set_yticks([26.08,26.10,26.12,26.14])
ax.set_yticklabels(['26.08','26.10','26.12','26.14'])

#plt.tight_layout()






#fig, ax = plt.subplots(4, 1)
#
#fig.suptitle('Modulated Seed and FEL spectra')
##seed laser spectrum
#ax[0].pcolormesh(np.arange(sl_spec.shape[1]),sl_roi, sl_spec)
#ax[0].set_ylabel('Seed')
#ax[0].ticklabel_format(useOffset=False)
#ax[0].set_ylim([264.2,264.8])
#ax[0].set_xticklabels([])
#
##5H spectrum
#ax[1].pcolormesh(np.arange(padres_spec_5H.shape[1]), padres_roi_5H, padres_spec_5H)
#ax[1].set_ylabel('5H')
#ax[1].set_ylim([52.1,52.4])
#ax[1].set_xticklabels([])
#
##7H spectrum
#ax[2].pcolormesh(np.arange(padres_spec_7H.shape[1]), padres_roi_7H, padres_spec_7H)
#ax[2].set_ylabel('7H')
##ax[2].set_xlabel('# Shots')
#ax[2].set_ylim([37.2,37.4])
#ax[2].set_xticklabels([])
#
##10H spectrum
#ax[3].pcolormesh(np.arange(padres_spec_10H.shape[1]), padres_roi_10H, padres_spec_10H)
#ax[3].set_ylabel('10H')
#ax[3].set_xlabel('# Shots')
#ax[3].set_ylim([26.05,26.18])


