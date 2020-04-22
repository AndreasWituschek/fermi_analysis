# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:29:08 2020

@author: andreas

This code takes simulation from Matlab script "num_sim_1D". It demodulates the
phase cycled population at each delay step and gives the demodulated amplitude
at the harmonics of the modulation frequency.
"""

import scipy.io as spio
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc

# =============================================================================
# functions
# =============================================================================
def slide_window(T, Z, t_center, FWHM):
    # multiplies a gaussion onto the dataset and returns the data set of the whole range of T
    # T, t_center, FWHM in ps
    # Z may be complex
    return Z*np.exp(-4*np.log(2)*((T-t_center)/FWHM)**2)


# =============================================================================
# parameters
# =============================================================================
harm = 5 #FEL harmonic
tile = 100 # how many cycles of 2pi are tiled 
plot_pop =1 #plot populations?
plot_paper = 0
E_r = 23.74207019 # helium 1s-4p transition in eV 

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

# =============================================================================
# #importing multiple population scans
# =============================================================================
#
#base_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/'
#base_name = '2020-02-07_amp125_'

base_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/full_Sa/'
base_name = '2020-02-24_full_Sa_amp130_'
file_pulse = '2020-02-24_full_Sa_pulse_amp130'

#base_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/full_Sa/'
#base_name = '2020-02-17_amp199_'
#file_pulse = base_name + 'pulse'

tau_steps_nr = 7001 #number of tau steps used in matlab script (=max_tau12 + 1 )
scan_nr = np.arange(0,1) #number pf phase cycling steps (defined in matlab script)
scans = np.zeros((len(scan_nr),tau_steps_nr)) #contains all the delay scans for all phase cycling steps
for n in scan_nr:
    mat = spio.loadmat(base_path + base_name + str(n) + '.mat', squeeze_me = True) # data population
    if n == 0: tau = mat['ttau12']
    scans[n] = mat['Popi']

mat2 = spio.loadmat(base_path + file_pulse + '.mat', squeeze_me = True) # simulated laser pulse
lambda_0 = mat2['WDp']# carrier wavelength in nm

amp = int(base_name[-4:-1])/100
#amp = 1/100

# =============================================================================
# importing data from Scan_31
# =============================================================================
X = np.zeros((harm,657))
Y = np.zeros((harm,657))
for i in [0,1,2,3,4]:
    XX = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/X_vs_delay_' +str(i+1) + '.txt').transpose()
    YY = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/Y_vs_delay_' +str(i+1) + '.txt').transpose()
    X[i] = XX
    Y[i] = YY
Z= X + 1j*Y
tau_data = np.linspace(-310.0393465165588,1001.9601165780342,657)

# =============================================================================
# plotting
# =============================================================================
plt.close('all')
if plot_pop:
    plt.figure('populations')
    for i in scan_nr:
        plt.plot(tau,scans[i],label=i)
#    plt.plot(tau,np.mean(scans,axis=0),label = 'mean',linewidth = 3,color = 'k')
    plt.legend()
    plt.xlabel('delay [fs]')
    plt.ylabel('population')
    plt.title('amp = {}'.format(amp))
    plt.tight_layout()


# =============================================================================
# Spectrogram analysis
# =============================================================================

X = np.append(np.flip(scans[0,1::]),scans[0])
T = np.append(-np.flip(tau[1::]),tau)
Td = np.mean(np.unique(np.diff(T)))

dft = np.fft.rfft(X)
wn = np.fft.rfftfreq(len(T),T[1]-T[0])/100.0/spc.c*1e15 # frequency axis in wn

#wn, dft = DFT(T, X, Td, np.inf , harm, zeroPaddingFactor = 2)

#spectrogram
slide_step = 2  #in fs
FWHM_slide = 45.0  #in fs
range_5H = np.arange(3800,4200) #range that gives 5H in dft
#range_5H = np.arange(0,len(dft))

slide_positions = np.arange(T[0],T[-1], slide_step)
S = np.zeros((len(slide_positions), np.size(wn)))
i = 0
for slide_pos in slide_positions:
    X_slide = slide_window(T, X, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
    DFT_slide = np.fft.rfft(X_slide)  # FFT of interferogram which has been truncated by mulitplying wiht gaussian
    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide[range_5H]))
    i += 1
S[:,:] /= np.size(slide_positions)
S /= np.max(abs(S))
S = S.transpose()

factor = spc.h*spc.c/spc.e*100.
wn *= factor #conversion to eV

maxvls = np.max(S,axis=0) #energy value where the spectrum peaks
max_idx = np.array([])
for i in np.arange(S.shape[1]):
    a = wn[np.argwhere(S[::,i]==maxvls[i])]
    max_idx = np.append(max_idx,a)
 
#plotting spectrogram

plt.figure('spectrogram_simu',figsize=(7,4))
plt.pcolormesh(slide_positions,wn[range_5H],S[range_5H]/np.max(abs(S[range_5H])), rasterized = True)
plt.colorbar()
plt.xlim([-300,300])
plt.ylim([E_r-0.13,E_r+0.11])
plt.axhline(E_r,color='grey',linestyle='--',label ='He 1s4p') #He 1s 4p transition energy
plt.axhline(spc.h*spc.c/(spc.e*52.23*1e-9),color='grey',linestyle='-',label = 'FEL')
plt.xlabel('delay [fs]')
plt.ylabel(r'energy [eV]')
plt.legend()
#plt.plot(T,savgol_filter(max_idx,19,3),color='k')    
plt.tight_layout()



# =============================================================================
# tests
## =============================================================================
##
#range_5H = np.arange(2800,2900)
#T = np.append(-np.flip(tau[1::]),tau)
#X = interp*np.cos(100*T)
#Td = np.mean(np.unique(np.diff(T)))
#
#dft = np.fft.rfft(X)
#wn = np.fft.rfftfreq(len(T),T[1]-T[0])/100.0/spc.c*1e15 # frequency axis in wn
#
##wn, dft = DFT(T, X, Td, np.inf , harm, zeroPaddingFactor = 2)
#
##spectrogram
#slide_step = 2  #in fs
#FWHM_slide = 25.0  #in fs
#
#
#slide_positions = np.arange(T[0],T[-1], slide_step)
#S = np.zeros((len(slide_positions), np.size(wn)))
#i = 0
#for slide_pos in slide_positions:
#    X_slide = slide_window(T, X, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
#    DFT_slide = np.fft.rfft(X_slide)  # FFT of interferogram which has been truncated by mulitplying wiht gaussian
#    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
#    i += 1
#S[:,:] /= np.size(slide_positions)
#S /= np.max(abs(S))
#S = S.transpose()
#
#factor = spc.h*spc.c/spc.e*100.
#wn *= factor #conversion to eV
#
#maxvls = np.max(S,axis=0) #energy value where the spectrum peaks
#max_idx = np.array([])
#for i in np.arange(S.shape[1]):
#    a = wn[np.argwhere(S[::,i]==maxvls[i])]
#    max_idx = np.append(max_idx,a)
# 
##plotting spectrogram
#
#plt.figure('spectrogram_scan31',figsize=(7,4))
#plt.pcolormesh(slide_positions,wn,S, rasterized = True)
#plt.colorbar()
#plt.xlim([-300,300])
##plt.ylim([E_r-0.13,E_r+0.11])
#plt.axhline(E_r,color='grey',linestyle='--') #He 1s 4p transition energy
#plt.xlabel('delay [fs]')
#plt.ylabel(r'energy [eV]')
##plt.plot(T,savgol_filter(max_idx,19,3),color='k')    
#plt.tight_layout()