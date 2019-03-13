# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import fermi_analysis.functions as fk
import h5py
import matplotlib.pyplot as plt
import scipy.constants as spc
import pickle
import os


run = 217
delay_zero_pos = 11025.66

ldm_file_path = '//online4ldm.esce.elettra.trieste.it/store/20149020/Day_2/Run_' + str(run) +'/rawdata/'.format(int(run))

ldm_files = os.listdir(ldm_file_path)
ldm_i0 = np.array([])
i0 = np.array([])


ldm_padres = np.array([])#padres spectrum
padres = np.zeros((0, 1000))#padres spectrum

ldm_slspec = np.array([]) #seed laser spectrum
slspec = np.zeros((0, 2068)) #seed laser spectum


for ldm_file in ldm_files:
    ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
    ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
    ldm_i0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a'])
    i0 = np.append(i0, ldm_i0) # I0 of FEL

    ldm_padres = np.array(ldm_data['photon_diagnostics']['Spectrometer']['hor_spectrum']) 
    padres = np.append(padres, ldm_padres,axis=0) 
    padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
    padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
    
    ldm_slspec = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum']) 
    slspec = np.append(slspec, ldm_slspec,axis=0) 


#plot whole run
padres = np.fliplr(np.transpose(padres)[440:570,:])
slspec = np.transpose(slspec)[1230:1320,:]

##only 400 shots:
#padres = np.fliplr(np.transpose(padres)[440:570,1000:1400])
#slspec = np.transpose(slspec)[1230:1320,1000:1400]

fig, axs = plt.subplots(2, 1)

axs[0].imshow(slspec)
axs[0].set_ylabel('Seed')
fig.suptitle('Modulated Seed and FEL spectra')

axs[1].imshow(padres)
axs[1].set_ylabel('10H')
axs[1].set_xlabel('# Shots')






#plt.imshow(slspec)



#plt.plot(slspec[:,1])

