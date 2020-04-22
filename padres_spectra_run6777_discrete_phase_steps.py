# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as spc
from scipy.ndimage import convolve1d

run=6777

delay_zero_pos = 11025.66


padres_spec = [] #padres spectrum of all harmonics(5,7,10)
padres_roi = []  #padres spectrum for all harmonics(5,7,10)

sl_spec = [] #seed spectrum of all runs
sl_roi = [] #seed ROI of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files
sl = np.zeros((0, 2068)) #seed laser spectum


''' Loop over all harmonics '''
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
#    ldm_padres = ldm_padres/np.max(ldm_padres,axis=1).reshape(400,1).astype(float)
    padres = np.append(padres, ldm_padres,axis=0) 
    #padres spectrometer metadata
    padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
    padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
    padres_pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])
    
    #seed laser spectrum
    ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)
#    ldm_sl = ldm_sl/np.max(ldm_sl,axis=1).reshape(400,1) #normation of seed laser specturm
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
padres_high=311
padres_spec_5H = np.transpose(padres_spec[0])[padres_low:padres_high,0:600]
padres_roi_5H = padres_roi[0][padres_low:padres_high]



## averaged padres and sl spectrum, normalized
padres_spec_avg = np.mean(padres_spec_5H,axis=1)/max(abs(np.mean(padres_spec_5H,axis=1)))




#==============================================================================
# Theoretical waveform of fringe potters
#==============================================================================

l_fel = 52.243840860758795 # in nm
delay = 250. #delay in fs
phase = 1.2
spec_reso = 2.5*1e-4 #padres spectrometer resolution in nm

def Envelope(l,amp,l0,fwhm,offs):
    return amp*np.exp(-4*np.log(2)*(l-l0)**2/fwhm**2)+offs

def Fringes(l,amp,tau,phase,offs):
    return amp*(1.+np.cos(2.*np.pi*spc.c/(l*1e-9)*tau*1e-15+phase))
    
popt, pcov = curve_fit(Envelope, padres_roi_5H, padres_spec_avg, p0=[1,52.2,0.09,0.01],absolute_sigma=False)
perr = np.sqrt(np.diag(pcov))

#l_theo = np.linspace(51.9,52.4,100000)
l_theo = padres_roi_5H
fringes= Fringes(l_theo,1,250,phase,0)
envelope= Envelope(l_theo,popt[0],popt[1],popt[2],popt[3])-popt[3]
pattern = 0.5*(fringes*envelope+popt[3])
resolution =  Envelope(l_theo,1,l_theo[int(len(l_theo)/2)],popt[1]*spec_reso,0)
resolution /= resolution.max()
pattern2= np.convolve(pattern,resolution,mode='same')
#pattern2 = convolve1d(pattern,resolution)


plt.plot(padres_roi_5H, padres_spec_avg)
#plt.plot(padres_roi_5H,padres_spec_5H[::,555]/max(padres_spec_5H[::,555]))
plt.plot(l_theo,envelope+popt[3])
plt.plot(l_theo,pattern)
plt.plot(l_theo,resolution)
plt.plot(l_theo,pattern2/max(pattern2))

#==============================================================================
# plotting data
#==============================================================================

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




#==============================================================================
# ''' Plotting moving fringes '''
#==============================================================================
fig = plt.figure(figsize=(10,8))

#seed laser spectrum
ax = fig.add_subplot(211)
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
ax = fig.add_subplot(212)
ax.pcolormesh(np.arange(padres_spec_5H.shape[1]), padres_roi_5H, padres_spec_5H) 
ax.tick_params(length=ticklength, width = ticksize)
ax.set_ylabel(r'5H $\lambda$ [nm]')
ax.set_xlim([0,599])
ax.set_ylim([52.1,52.4])
#ax.set_xticklabels([])
ax.set_yticks([52.15,52.25,52.35])
ax.set_yticklabels(['52.15','52.25','52.35'])

#==============================================================================
# ''' plotting  4 discrete phase steps '''
#==============================================================================
steps = [0,21,42,63] #corresponds phase of 5H: 0, pi/2, pi, 3pi/2
start_idx = 199
idx = [i + start_idx for i in steps]

fig2 = plt.figure(figsize=(7,8))

ax= fig2.add_subplot(211)
for i in idx:
    ax.plot(sl_roi,sl_spec[::,i]/59220.+((i-start_idx)/21.)*0.44,color='blue')
ax.legend([r'$0$',r'$\pi/10$','$\pi/5$','$3\pi/10$'])
ax.set_xlim([258.8,263.2])
ax.set_ylim([-0.1,2.4])
ax.set_xlabel(r'seed wavelength [nm]')
ax.set_ylabel('intensity [a.u.]')

ax= fig2.add_subplot(212)
for i in idx:
    ax.plot(padres_roi_5H,padres_spec_5H[::,i]/367800.+((i-start_idx)/21.)*0.54,color='darkviolet')
ax.legend([r'$0$',r'$\pi/2$','$\pi$','$3\pi/2$'])
ax.set_xlim([52.1,52.4])
ax.set_ylim([-0.1,2.4])
ax.set_xlabel(r'5H wavelength [nm]')
ax.set_ylabel('intensity [a.u.]')
plt.tight_layout()

#==============================================================================
# plotting  2 discrete phase steps '''
#==============================================================================
steps = [42,0] #corresponds phase of 5H: 0, pi/2
start_idx = 230
idx = [i + start_idx for i in steps]

fig2 = plt.figure(figsize=(7,8))

ax= fig2.add_subplot(211)
for i in idx:
    ax.plot(sl_roi,sl_spec[::,i]/59220.+((i-start_idx)/21.)*0.54,color='blue')
ax.legend([r'$0$',r'$\pi/10$','$\pi/5$','$3\pi/10$'])
ax.set_xlim([258.8,263.2])
ax.set_ylim([-0.3,2.4])
ax.set_xlabel(r'seed wavelength [nm]')
ax.set_ylabel('intensity [a.u.]')

ax= fig2.add_subplot(212)
for i in idx:
    ax.plot(padres_roi_5H,padres_spec_5H[::,i]/367800.+((i-start_idx)/21.)*0.54,color='darkviolet')
ax.legend([r'$0$',r'$\pi/2$','$\pi$','$3\pi/2$'])
ax.set_xlim([52.08,52.38])
ax.set_ylim([-0.3,2.4])
ax.set_xlabel(r'5H wavelength [nm]')
ax.set_ylabel('intensity [a.u.]')
plt.tight_layout()

