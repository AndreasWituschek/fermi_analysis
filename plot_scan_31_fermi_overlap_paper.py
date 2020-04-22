# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:46:56 2019

@author: FemtoMeasure
"""

import numpy as np
import fermi_analysis.functions as fk 
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import scipy.constants as spc


#plt.close('all')
""" Data parameters """
run = 31 #first run of delay scan
demod = ['0','1','2']
device = ['dev3265','dev3269']
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'


''' plot params '''
plot_R = 1
plot_harmonics = 1
plot_paper = 1
export_Z = 0

"""Experimental Parameters"""
# parameters for theoretical curve
E_r = 23.74207019 # helium 1s-4p transition in eV 
l_trans =  spc.h*spc.c/(E_r*spc.e)*1e9  # helium 1s-4p transition in nm 
fwhm_FEL = 0.075 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Scan_{0:03}/'.format(int(run))
color = 'r'
data_window = [-300,300] #[-200,370]
gauss = 1 # gauss window on TD (true) or not (False)
modfreq = 3.96


'''Phase correction'''
phas_corr = harmonic*8. #phase of reference transmitter [degree] (gamma in Lukas phasing routine)
phas_corr -= 0.01*(harmonic*modfreq)*360 #phase of boxcar integrator [degree] (beta in Lukas phasing routine)
print('A phase correction of {} degrees due to electronics was applied'.format(phas_corr))
phas_res = 30. #more or less arbitrary phasing factor that makes real part absorptive and imaginary part dispersive
print('A phase correction of {} degrees was applied to yield a perfect absorptive line shape'.format(phas_res))
phas_corr += phas_res
print('Total phase correction: {} degree'.format(phas_corr))
phas_corr -= 180. #to take into account that we measure negative signals by ion detection
phas_corr += (180.*harmonic) % 360. #accounts for phase shift due to the fact that ref signal is other output of interferometer.

""" Parameters for theoretical curve """
phi = 0. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
zeroPaddingFactor = 1
suscept = False
peakMargin = 0.05
wn_lim = np.asarray([185000.,196000.])


# Load demodulated data
data = fk.ImportPreanalysedData(root_file_path,run)

#l_fel = data['LDM']['l_seed'][0]
#ldm_l_ref = data['LDM']['l_ref']
l_ref = 266.003 # reference laser frequency
E_undersamp = harmonic*spc.h*spc.c/(l_ref*spc.e)*1e9 # undersampling energy
T = data['LDM']['delay']-5.38




Z = np.empty((5,len(T)),dtype = complex)
for dev in device:
    # loop over the demodulators
    for d in demod:
        if not (dev == 'dev3269' and d == '2'): #excluding the "droplet" signal.
            mfli_harmonic = data[dev]['harmonic'][int(d)]
            i = str(mfli_harmonic) + 'H'
            Z[mfli_harmonic - 1] = data[dev]['x' + d] + 1j * data[dev]['y' + d]

#plotting all harmonics
if plot_harmonics:
    colorcycle = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
    nrow = len(Z)
    ncol = 1
    norm = np.max(np.abs(Z[0])) #normation on 1H
    fig, axs = plt.subplots(nrow, ncol,figsize=(6,9),sharex = 'col')
    for i, ax in enumerate(fig.axes):
        ax.plot(T,Z[i]/norm,label='{}H'.format(i+1),color=colorcycle[i])
        ax.plot(T,np.abs(Z[i]/norm),color='grey')
        ax.get_xaxis().set_visible(False)
        if i==nrow-1:
            ax.get_xaxis().set_visible(True)
        ax.legend()
        ax.set_xlim([-300,300])
        ax.axvline(0,color='black')
        ax.set_xlabel('delay [fs]')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

if export_Z:
    for i in [0,1,2,3,4]:
        np.savetxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/X_vs_delay_' +str(i+1) + '.txt',Z[i].real)
        np.savetxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/Y_vs_delay_' +str(i+1) + '.txt',Z[i].imag)



#plotting the absolute values
if plot_R:
    R1to5 = np.sum(np.abs(Z),axis = 0)
    R1to4 = np.sum(np.abs(Z[0:4]),axis = 0)
    plt.figure('R-values')
    plt.plot(T,R1to5)
    plt.plot(T,R1to4)

#seed laser AC    
sl_fwhm = 99. # seed intensity fwhm in fs
t_theo_seed = np.linspace(-250,250,5000)
sl_AC = np.exp(-4.0*np.log(2)*(t_theo_seed/(sl_fwhm*np.sqrt(2)))**2) # seed laser AC
sl_AC /= np.max(sl_AC)

# =============================================================================
# spectrogram
# =============================================================================

#analysing 5H
Z5 = Z[4]
Z5 *= np.exp(1j*np.pi/180.*phas_corr) #phasing data

T_d,X_d,Y_d = fk.CutDataSet(T, Z5.real, Z5.imag, data_window)
Z5_cut = X_d + 1j*Y_d
Z5_cut = fk.Phasing_TD(Z5_cut,T_d,E_r-E_undersamp)
scaling = 1E-4 #max(Z.real)/max(Ttheo)

# Fourier traffo
Td = np.mean(np.unique(np.diff(T)))
if gauss:
    Z5g = fk.GaussWindow(T_d, Z5_cut, suscept)
    wn, dft = fk.DFT(T_d, Z5g, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
else:
    wn, dft = fk.DFT(T_d, Z5_cut, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
FWHM , peakFullWidth = fk.PeakWidth(T_d, peakMargin, suscept, zeroPaddingFactor)


#conversion factor from wn to eV
factor = spc.h*spc.c/spc.e*100.

#spectrogram
slide_step = Td  #in fs
FWHM_slide = 45.0  #in fs
slide_positions = np.arange(T[0],T[-1], slide_step)
S = np.zeros((np.size(T), np.size(wn)))
i = 0
for slide_pos in slide_positions:
    Z5_slide = fk.slide_window(T, Z5, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
    wn,  DFT_slide = fk.DFT(T_d, Z5_slide, Td, l_ref, harmonic, zeroPaddingFactor = zeroPaddingFactor)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
    i += 1
S[:,:] /= np.size(slide_positions)
S /= np.max(abs(S))
S = S.transpose()
wn *= factor

maxvls = np.max(S,axis=0) #energy value where the spectrum peaks
max_idx = np.array([])
for i in np.arange(S.shape[1]):
    a = wn[np.argwhere(S[::,i]==maxvls[i])]
    max_idx = np.append(max_idx,a)
    
    

''' Plots for Fermi_pulse_overlap paper '''

if plot_paper:
    ticksize= 2.
    ticklength = 5.
    fontsize=13.
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.major.width'] = ticksize
    plt.rcParams['xtick.major.size'] = ticklength
    plt.rcParams['ytick.major.width'] = ticksize
    plt.rcParams['ytick.major.size'] = ticklength
    plt.rcParams['axes.linewidth'] = ticksize
    plt.rcParams['lines.linewidth'] = 2.
    
    
    #plotting spectrogram
    plt.figure('spectrogram_scan31',figsize=(7,4))
    plt.pcolormesh(T,wn,S, vmin=0, rasterized = True)
    plt.colorbar()
    plt.xlim([-300,300])
    plt.ylim([E_r-0.13,E_r+0.11])
    plt.axhline(E_r,color='grey',linestyle='--') #He 1s 4p transition energy
    plt.xlabel('delay [fs]')
    plt.ylabel(r'energy [eV]')
    plt.axhline(E_r,color='grey',linestyle='--',label ='He 1s4p') #He 1s 4p transition energy
    plt.axhline(spc.h*spc.c/(spc.e*52.23*1e-9),color='grey',linestyle='-',label = 'FEL')
#    plt.plot(T,savgol_filter(max_idx,19,3),color='k')    
    plt.legend()
    plt.tight_layout()

#version 01 of helium plot    
    #plotting time domain data
    fig, axs = plt.subplots(2, 1,figsize=(7,6), sharex = 'col')
    
    ax = axs[0]
    ax.plot(T, Z5.real)
    ax.plot(5,5,'silver')
    ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z5.real))*0.8,color='silver')
    ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z5.real))*0.8,color='silver')
    ax.set_ylabel('Re(S) [arb.units]')
#    ax.set_xlabel(r'$\tau$ [fs]')
    ax.set_xlim(-310,310)
    ax.set_ylim(-0.25,0.25)
    
    ax = axs[1]
    ax.plot(T, data['LDM']['I0'])
    ax.set_ylabel(r'FEL pulse energy [$\mu$J]')
    ax.set_xlabel(r'$\tau$ [fs]')
    ax.set_xlim(-310,310)
    ax.set_ylim(20,80)
    plt.tight_layout()
    
#version 02 ov helium plot
    #plotting time domain data
    fig, axs = plt.subplots(2, 1,figsize=(7,6))
    
    ax = axs[0]
    ax.plot(T, Z5.real)
    ax.plot(5,5,'silver')
    ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z5.real))*0.6,color='silver')
    ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z5.real))*0.6,color='silver')
    ax.set_ylabel('Re(S) [arb.units]')
    ax.set_xlabel(r'$\tau$ [fs]')
    ax.set_xlim(-310,310)
    ax.set_ylim(-0.2,0.25)
    
    axi = ax.twinx()
    axi.plot(T, data['LDM']['I0'],color='#d62728')
    axi.set_ylabel(r'FEL pulse energy [$\mu$J]')
    axi.set_xlabel(r'$\tau$ [fs]')
    axi.set_xlim(-310,310)
    axi.set_ylim(-15,80)
    
#    axi.spines['right'].set_color('#d62728')
#    ax.spines['left'].set_color('red')
    
    ax.yaxis.label.set_color('#1f77b4')                         
    ax.tick_params(axis='y', colors='#1f77b4')
    axi.yaxis.label.set_color('#d62728')
    axi.tick_params(axis='y', colors='#d62728')

    
    ax = axs[1]
    ax.plot(T, Z5.real)
    ax.set_ylabel('Re(S) [arb.units]')
    ax.set_xlabel(r'$\tau$ [fs]')
    ax.set_xlim(-100,100)
    ax.set_ylim(-0.025,0.025)
    plt.tight_layout()