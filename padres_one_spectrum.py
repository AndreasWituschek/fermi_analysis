# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import fermi_analysis.functions as fk

runs = [4750] #[5H, 7H, 10H]
#run=3213
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
        ldm_padres = ldm_padres/np.max(ldm_padres,axis=1).reshape(400,1).astype(float)
        padres = np.append(padres, ldm_padres,axis=0) 
        #padres spectrometer metadata
        padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
        padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
        padres_pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])
        
        #seed laser spectrum
        
        ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)
        ldm_sl = ldm_sl/np.max(ldm_sl,axis=1).reshape(400,1) #normation of seed laser specturm
        print ldm_sl.shape
        sl = np.append(sl, ldm_sl,axis=0) 
#        Seed laser metadata
        ldm_roi = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['RegionOfInterest']) # ROI for seed laser spectrum
        ldm_data.close
    
    padres_lambda = padres_wavelength + (np.arange(1,1001)-500)*padres_pixel2micron*1E-3*padres_span #padres spectrometer calibration
    padres_spec.append(padres)
    padres_roi.append(padres_lambda)
    
#sl_lambda= np.linspace(ldm_roi[0],ldm_roi[1],ldm_sl.shape[1]) #seed laser spectrometer calibration


''' Proper Data windows for all harmonics '''
#seed laser 
sl_low=1220
sl_high=1330
sl_spec = np.transpose(sl)[sl_low:sl_high]
sl_roi = sl_lambda[sl_low:sl_high]


#data window 5H run 6777
padres_low=20
padres_high=300
padres_spec_5H = np.transpose(padres_spec[0])[padres_low:padres_high]
padres_roi_5H = padres_roi[0][padres_low:padres_high]



padres_single = np.mean(padres_spec_5H,axis=1)

#save_file = np.array([padres_single,padres_roi_5H])
#np.savetxt("run4644_padres.csv", save_file, delimiter="\t")
#pd.DataFrame(save_file.transpose()).to_csv("run4644_padres.csv",delimiter="\t",index=False)

#==============================================================================
# Fitting gaussion to spectrum average
#==============================================================================

def Gauss(x,amp,mean,fwhm,offs):
    return amp*np.exp(-4*np.log(2)*(x-mean)**2/fwhm**2)+offs

popt_gauss, pcov_gauss = curve_fit(Gauss, padres_roi_5H, padres_single, p0=[1.,padres_roi_5H[fk.find_index(padres_single,max(padres_single))],0.07,0.05], absolute_sigma=False)
perr_gauss = np.sqrt(np.diag(pcov_gauss))
peak_pos_sl = popt_gauss[1]

print popt_gauss[2]

#==============================================================================
# ''' Plotting '''
#==============================================================================
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(padres_roi_5H,padres_single)
ax.plot(padres_roi_5H,Gauss(padres_roi_5H,popt_gauss[0],popt_gauss[1],popt_gauss[2],popt_gauss[3]))



fig = plt.figure(figsize=(10,5))

#seed laser spectrum
ax = fig.add_subplot(211)
ax.pcolormesh(np.arange(sl_spec.shape[1]),sl_roi, sl_spec)
ax.set_ylabel(r'Seed $\lambda$ [nm]')
ax.ticklabel_format(useOffset=False)
ax.set_ylim([264.2,264.8])
ax.set_xticklabels([])

#5H spectrum
ax = fig.add_subplot(212)
ax.pcolormesh(np.arange(padres_spec_5H.shape[1]), padres_roi_5H, padres_spec_5H)
ax.set_ylabel(r'5H $\lambda$ [nm]')
#ax.set_ylim([52.1,52.4])
ax.set_xticklabels([])

##7H spectrum
#ax = fig.add_subplot(223)
#ax.pcolormesh(np.arange(padres_spec_7H.shape[1]), padres_roi_7H, padres_spec_7H)
#ax.set_ylabel(r'7H $\lambda$ [nm]')
#ax.set_xlabel('# Shots')
#ax.set_ylim([37.2,37.4])
##ax.set_xticklabels([])
#
##10H spectrum
#ax = fig.add_subplot(224)
#ax.pcolormesh(np.arange(padres_spec_10H.shape[1]), padres_roi_10H, padres_spec_10H)
#ax.set_ylabel(r'10H $\lambda$ [nm]')
#ax.set_xlabel('# Shots')
#ax.set_ylim([26.05,26.18])
#
#plt.tight_layout()
#





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


