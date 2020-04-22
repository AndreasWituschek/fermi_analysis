# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: Andi

Plots Padres spektrum from preanalyzed files (skript for that: padres_vs_delay.py) vs. the delay. so one can see shifts in the laser frequency.
"""
import numpy as np
#import os
#import errno
import h5py
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import scipy.constants as spc

#plt.close('all')

file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined\scan_031/scan_031_padres_spectra.h5'
data = h5py.File(file_path, 'r')


#indices where below and above only zeros are written into padres spectrum respectively
low = 8
high = 310

E_r = 23.74207019 #eV Helium 1s-4p transition

padres = np.array(data['/LDM/padres']).transpose()[low:high,::]
delay = np.array(data['/LDM/delay'])

padres -= np.min(padres)
padres /= np.max(padres)

data.close()

#spectrometer calibration  as obtained from LDM files
padres_lambda = 51.49174461 + (np.arange(1,1001)-500)*15.4639*1E-3*-0.15670998176276718 #scan 31 
#padres_lambda = 51.49957925234294 + (np.arange(1,1001)-500)*15.4639*1E-3*-0.15672136832796663 #scan 1164
#padres_lambda = 42.649328216442186 + (np.arange(1,1001)-500)*15.4639*1E-3*-0.141875698427862 #scan 4651
#padres_lambda = 52.220296899060386 + (np.arange(1,1001)-500)*15.4639*1E-3*-0.15787277243579323 #scan 893

padres_lambda = padres_lambda[low:high]

#laser and transition frequency
wl_laser = 261.07/5. #scan 31 and scan893
#wl_laser = 262.07/5. #scan1164
#wl_laser = 260.87/6. #scan4651
wl_trans = 52.221309 #Helium transition wavelength
#wl_trans = 43.48763 #argon transition wavelength

wn_lim = [23.64,23.9] #Helium in eV
#wn_lim = [43.3,43.65] #Argon in nm

factor = spc.h/spc.e*spc.c*1e9 #conversion from nm to eV
padres_lambda = factor/padres_lambda #calibrated spectrometer energy axis
delay_plot = np.linspace(min(delay),max(delay),padres.shape[1]) #delay axis 

maxvls = np.max(padres,axis=0)
max_idx = np.array([])
for i in np.arange(padres.shape[1]):
    a = padres_lambda[np.argwhere(padres[::,i]==maxvls[i])]
    max_idx = np.append(max_idx,a)


#cut = np.max(np.argwhere(delay_plot>-200))

plt.figure('padres spectrum vs delay',figsize=(7,4))
plt.pcolormesh(delay_plot,padres_lambda,padres*1.1,vmax=1 , rasterized = True)
plt.colorbar()
plt.ylim([E_r-0.13,E_r+0.11])
plt.xlim([-300,300])
#plt.axhline(wl_laser,color='grey',linewidth=2, label = 'wl_laser')
plt.axhline(E_r,color='grey',linestyle='--') #He 1s 4p transition energy
plt.xlabel('delay [fs]')
plt.ylabel(r'h$\nu$ [eV]')
plt.plot(delay_plot,savgol_filter(max_idx,19,3),color='k')    
plt.tight_layout()


