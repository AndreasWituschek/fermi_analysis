# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy import constants as spc
import scipy.signal as sig
from matplotlib.gridspec import GridSpec

run = 6776 #6756#1431  #5H 0.0Hz modulation
#run = 6777 #5H 0.12Hz modulation

#run = 215 #6H 0.12Hz modulation, CAREFUL: First beamtime!

#run= 6757 # 7H 0Hz Modulation
#run= 6773 # 7H 0.12Hz modulation

#run = 6760 #8H 0.0Hz modulation
#run = 6761 #8H 0.12Hz modulation

#run = 6765 #9H 0.0Hz modulation
#run = 6764 #8H 0.12Hz modulation

#run= 6766 #10H 0Hz modulation, tau = 350
#run= 6768 #10H 0.12Hz modulation

padres_spec = [] #padres spectrum of all harmonics(5,7,10)
padres_roi = []  #padres spectrum for all harmonics(5,7,10)

sl_spec = [] #seed spectrum of all runs
sl_roi = [] #seed ROI of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files
sl = np.zeros((0, 2068)) #seed laser spectum


''' Loop over all harmonics '''
#ldm_file_path = 'D:/FermiServer/Beamtime2/Run_' + str(run) +'/rawdata/'.format(int(run))
ldm_file_path = 'C:/Users/andreas/Desktop/Run_' + str(run) +'/rawdata/'.format(int(run))

ldm_files = os.listdir(ldm_file_path)

ldm_padres = np.array([]) #raw padres spectrum from ldm files
padres = np.zeros((0, 1000))#padres spectrum
bam = np.array([]) #beam arrival monitor [fs]

  
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
    
    #bunch arrival monitor
    ldm_bam = np.array(ldm_data['photon_source']['bc01_diagnostics']['bam_bc01.01']['Arrival'])
    bam = np.append(bam,ldm_bam)
    ldm_data.close

padres_lambda = padres_wavelength + (np.arange(1,1001)-500)*padres_pixel2micron*1E-3*padres_span #padres spectrometer calibration
padres_spec.append(padres)
padres_roi.append(padres_lambda)

padres_spec = np.transpose(padres_spec[0])
#filtering the zero FEL power or ROI problem shots out:
filt = np.nonzero(np.mean(padres_spec,axis=0)>1000.)
padres_spec = padres_spec[::,filt[0]]
padres_roi = padres_roi[0]
padres_spec = padres_spec[::,]

sl = np.transpose(sl)

#==============================================================================
 #''' Determining phase stability via determination of maximum position of fringe in spectrometer data '''
#==============================================================================
#boundaries for max_idx
low_idx=0
high_idx=195
#max_idx gives the indices where the max of the padres spectrum is located (in pixels of the spectrometer)
max_idx = np.asarray([np.min(np.argwhere(padres_spec[low_idx:high_idx,i]==np.max(padres_spec[low_idx:high_idx,i]))).flatten() for i in range(padres_spec.shape[1])]).flatten()+low_idx
plt.figure()
plt.plot(max_idx,'bo')
plt.xlabel('shots')
plt.ylabel('index of maxiumum')

l_fel = 261.1/5. # in nm
delay = 250. #delay in fs

fringe_spacing =  (l_fel*1e-9)**2/(delay*1e-15*spc.c)/1e-9 #in nm
print(1e9*fringe_spacing*spc.h*spc.c/l_fel**2/spc.e)
padres_reso =  np.mean(-np.unique(np.diff(padres_roi))) # "stepsize" of padres spectrometer in nm
stab = (np.std(max_idx)*padres_reso)/fringe_spacing # converting from pixels to nm
stab *= 360.
print('The CEP stability by the max method is {}'.format(stab)+' degree')

plt.figure()
plt.plot(padres_spec[::,15:18])
plt.xlim([0,400])
plt.xlabel('index')

#
#==============================================================================
# ''' Determininge fringe stability by looking at phase of Fringe peak in DFT '''
#==============================================================================
#multiplying padres_spec with a arbitrary phase to avoid that phase always wraps around 2pi
padres_spec = padres_spec*np.exp(1j*np.pi/2)

#fourier transforming one spectrum to see where in dft the peak-maximum for the fringe spacing is at
dft =  np.fft.fftshift(np.fft.fft(padres_spec[::,10]))
plt.figure()
plt.plot(np.angle(dft)/np.pi)
plt.plot(abs(dft)/max(abs(dft)))
dft_max_idx = np.argwhere(abs(dft)==max(abs(dft[550:700]))) #index where abs(dft) has its maximum
print(dft_max_idx)

#phase is mean value of phase in the peak (located at dft_max_idx) of the DFT 
phase = np.array([np.mean(np.angle(np.fft.fftshift(np.fft.fft(padres_spec[::,i])))[dft_max_idx]) for i in range(padres_spec.shape[1])])
stab = np.std(phase)
#print('The CEP stability by the FT method is {}'.format(stab)+' rad')
stab = np.std(phase)/(2.*np.pi)*360. #phase stability in degree
print('The CEP stability by the FT method is {}'.format(stab)+' degree')
phase_filtered =  sig.savgol_filter(phase,window_length=31,polyorder=3)


# =============================================================================
# doing the same for the seed laser
# =============================================================================
#multiplying seed laser spectrum with a arbitrary phase to avoid that phase always wraps around 2pi
#sl = sl*np.exp(1j*np.pi/2)

#fourier transforming one spectrum to see where in dft the peak-maximum for the fringe spacing is at
dft_sl =  np.fft.fftshift(np.fft.fft(sl[::,10]))
plt.figure()
plt.plot(np.angle(dft_sl)/np.pi)
plt.plot(abs(dft_sl)/max(abs(dft_sl)))
dft_sl_max_idx = np.argwhere(abs(dft_sl[1130:1300])==max(abs(dft_sl[1130:1300])))+1130 #index where abs(dft) has its maximum
print(dft_sl_max_idx)

#phase is mean value of phase in the peak (located at dft_max_idx) of the DFT 
phase_sl = np.array([np.mean(np.angle(np.fft.fftshift(np.fft.fft(sl[::,i])))[dft_sl_max_idx]) for i in range(sl.shape[1])])
stab_sl = np.std(phase_sl)/(2.*np.pi)*360. #phase stability in degree
print('The seed CEP stability by the FT method is {}'.format(stab_sl)+' degree')


# =============================================================================
# ''' Plots for Paper '''
# =============================================================================
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

fig_dft = plt.figure('CEP_stability_Run_{}'.format(run),figsize=(8, 5))

gs = GridSpec(1, 11)
ax = fig_dft.add_subplot(gs[:,0:7])

#ax.set_title('CEP_stability_Run_{}'.format(run))
ax.plot(phase)
ax.plot(phase_filtered)
ax.set_xlabel('laser shot')
ax.set_ylabel('phase [rad]',alpha=0.7)
ax.legend([r'$n\phi_{21}$','filt. phase'])
ax.set_ylim([-3.5,3.5])
#props = dict(boxstyle='round', facecolor='white', alpha=1)
#ax.text(3000,-2.5, 'CEP stability: {}'.format(round(stab,1))+' degree',bbox=props)
plt.tight_layout()
#plt.plot(np.angle(np.fft.fftshift(np.fft.fft(np.sin(np.linspace(10,100,1000)+1.)))))

ax = fig_dft.add_subplot(gs[:,7:])
ax.hist(phase,60,orientation='horizontal',alpha=0.7)
ax.set_ylim([-3.5,3.5])
ax.get_yaxis().set_visible(False)
ax.set_xticklabels([])
ax.set_xlabel('counts')

# =============================================================================
# plotting BAM vs. fringe phase
# =============================================================================

start = 4000
stop = start + 2000

#some statistics:
print('The mean bam is {} fs. The bam std dev is {} fs'.format(np.mean(bam[start:stop]),np.std(bam[start:stop])))
print('The mean seed CEP is {} rad. The std dev is {} rad'.format(np.mean(phase_sl[start:stop]),np.std(phase_sl[start:stop])))
print('The mean FEL CEP is {} rad. The std dev is {} rad'.format(np.mean(phase[start:stop]),np.std(phase[start:stop])))



fig_stability = plt.figure('BAM vs. fringe phase',figsize=(8, 8))

ax = fig_stability.add_subplot(211)
ax.plot(-phase[start:stop])
ax.plot(-phase_sl[start:stop])
ax.set_ylabel('FEL fringe phase [rad]')
ax.set_ylim([-1,3])
ax.set_xlim([0,2000])
ax.legend(['FEL fringe phase', 'seed fringe phase'])
ax.set_xlabel('laser shot')

#ax2 = fig_stability.add_subplot(312)
#ax2.plot(-phase_sl[start:stop])
#ax2.set_ylim([-1,3])
#ax2.set_xlabel('laser shot')
#ax2.set_ylabel('FEL fringe phase [rad]',alpha=0.7)
#


ax3 = fig_stability.add_subplot(212)
ax3.scatter(-phase[start:stop],bam[start:stop],s=3)
#plt.plot(max_idx[start:stop],bam[start:stop],'o')
ax3.set_xlabel('FEL fringe phase [rad]')
ax3.set_ylabel('e-bunch jitter [fs]')
ax3.set_xlim([-0.5,3])
plt.tight_layout()

plt.figure()
plt.scatter(-phase[start:stop],bam[start:stop],s=3)
#plt.plot(max_idx[start:stop],bam[start:stop],'o')
plt.xlabel('Ramsey-fringe phase [rad]')
plt.ylabel('e-bunch jitter [fs]')
plt.xlim([-0.5,3])
plt.tight_layout()

plt.savefig('figS2_jitter',dpi=500)