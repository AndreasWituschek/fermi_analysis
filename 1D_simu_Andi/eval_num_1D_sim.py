# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:45:24 2020

@author: MarcelB
"""

import scipy.io as spio
from scipy import fftpack
import scipy as sp
import numpy as np
from scipy import constants as spc
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
from scipy import optimize 
from math import log10, floor
from scipy.signal import butter, lfilter, freqz, savgol_filter


# user input ==================================================================

file_name1 = '2020-02-05_test' # name of file to be loaded with simulated population
file_name2 = '2020-02-05_test_pulse' # name of file to be loaded with simulated laser pulses (in rotating frame)

plot_pulse = False

drop = 0.05 # drop Gauss to x
zero_pad_val = 2 # value for zeropadding
wlam = 1239.8424121      # nm --- eV conversion factor for nm to eV

fs = 20 # fontsize

#reference lines
freq_trans = 191492.711#transition
harm = 5

# defining some functions =====================================================


def Gaussian_1D(size, fwhm, center=0):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
# skalieren auf Länge!
    x = np.arange(0, size, 1, float)
    if center is None:
        x0  = size / 2
# warum? sollte immer bei 0,0 fs liegen
# Weil ich die Funktion aus einer Machine Learning Übung kopiert habe :P
    else:
        x0 = center

#    return 1.0
    return np.exp(-4*np.log(2) * ((x-x0)**2) / (fwhm)**2)


def absmax(data):
    return np.max(np.abs(data))

def zero_padding(T, factor):
    return 2**(int(round(np.log2(T.size)+0.5))+factor)  # Anzahl Datenpunkte auf die zur Not mit Zero-Padding aufgefuellt wird (für 1D FFT)

def redistribute(array, index):
    return np.concatenate((array[index:],array[0:index]))

def find_index(array, value):
    return np.argmin(abs(array-value))

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def idft_window(dft_ac,cent,sigma,harm,points,order):
    points  = int(points)
    b, a = butter_bandpass(harm*cent-sigma, harm*cent+sigma, points, order=order)
    w, h = freqz(b, a, worN= int(points/2))
    window = np.abs(np.append(np.roll(np.flip(h),1),h))
    return np.fft.ifft(np.fft.ifftshift(dft_ac*window)), window

# main script =================================================================
    
# import data from matlab simulation

mat1 = spio.loadmat(file_name1 + '.mat', squeeze_me = True) # data population
mat2 = spio.loadmat(file_name2 + '.mat', squeeze_me = True) # simulated laser pulse


#importing multiple population scans

base_name = '2020-02-07'
tau_steps_nr = 201 #number of tau steps used in matlab script (=max_tau12 + 1 )
scan_nr = np.arange(0,20) #number if phase cycling points
scans = np.zeros((len(scan_nr),tau_steps_nr))
for n in scan_nr:
    mat = spio.loadmat(base_name + str(n) + '.mat', squeeze_me = True) # data population
    if n == 0: tau = mat['ttau12']
    scans[n] = mat['Popi']

# assigning content to variables

# 1D scan data
tau = mat1['ttau12'] # tau vector
step_size = mat1['step_tau'] # step size tau scan
wp = mat1['wp'] # pulse carrier frequency in eV (Spectrum has to be shifted by this value later, as simulation is in rotating frame)
time_dat = mat1['Popi'] # time data, vector with population values dependent on pulse delay tau

laser_freq_nm = wlam/wp
laser_freq_wn = 1.0/(laser_freq_nm*10**(-7))# freq axis in wavenumbers

# simulated pulses
delta_t_limit = mat2['TDp'] # transform limited pulse duration (intensity)
lambda_0 = mat2['WDp']# carrier wavelength in nm
ch2 = mat2['ch2'] # second order dispersion in fs^2
ch3 = mat2['ch3'] # third order dispersion in fs^3
t = mat2['TRK'] # values for time axis
pulse_TD = mat2['Pulse1'] # temporal amplitude (normalized)
freq = mat2['WWW']# array with values for frequency axis in 2Pi*10^3*THz, in rotating frame of carrier freq
pulse_FD = mat2['PulseW'] # frequency domain data of laser pulse (spectral amplitude)

# shifting the frequency axis (eV) and then transforming it to nm, THz and cm^-1

freq_THz = ((freq*10**15/(2*np.pi))+(spc.c/(lambda_0*10**(-9))))*10**(-12)
freq_nm = (spc.c/(freq_THz*10**12))*10**(9)
freq_eV = wlam/freq_nm 
freq_wn = 1.0/(freq_nm*10**(-7))# freq axis in wavenumbers

temporal_intens = abs(pulse_TD)**2
spectral_intens = abs(pulse_FD)**2


# plotting E-field data (rotating frame, normalized)

if plot_pulse:
    #TD data: temporal amplitued
    fig1 = plt.figure(file_name2 + '_TD', figsize=(10,6))
    ax1 = fig1.add_subplot(111)
    ax1.plot(t, pulse_TD.real, '-o', label='real', color ='r', linewidth = 2)
    ax1.plot(t, pulse_TD.imag, '-o', label='imag', color ='b', linewidth = 2)
    ax1.plot(t, abs(pulse_TD), '--o', label = 'abs' , color = 'k', linewidth = 2) 
    ax1.set_ylabel('temporal amplitude [a.u.]', fontsize = fs)
    ax1.set_xlabel('time [fs]', fontsize = fs)
    ax1.tick_params(labelsize=fs)
    ax1.legend(loc='best', fontsize = fs)


    #FD data: spectral amplitude
    fig2 = plt.figure(file_name2 + '_FD', figsize=(10,10))
    ax2 = fig2.add_subplot(221)
    ax2.plot(freq_eV, pulse_FD.real, '-o', label='real', color ='r', linewidth = 2)
    ax2.plot(freq_eV, pulse_FD.imag, '-o', label='imag', color ='b', linewidth = 2)
    ax2.plot(freq_eV, abs(pulse_FD), '--o', label = 'abs' , color = 'k', linewidth = 2) 
    ax2.set_ylabel('spectral amplitude [a.u.]', fontsize = fs)
    ax2.set_xlabel('energy [eV]', fontsize = fs)
    ax2.tick_params(labelsize=fs)
    ax2.legend(loc='best', fontsize = fs)

    
    ax3 = fig2.add_subplot(222)
    ax3.plot(freq_nm, pulse_FD.real, '-o', label='real', color ='r', linewidth = 2)
    ax3.plot(freq_nm, pulse_FD.imag, '-o', label='imag', color ='b', linewidth = 2)
    ax3.plot(freq_nm, abs(pulse_FD), '--o', label = 'abs' , color = 'k', linewidth = 2) 
    #ax3.set_ylabel('amplitude [a.u.]', fontsize = fs)
    ax3.set_xlabel('wavelength [nm]', fontsize = fs)
    ax3.tick_params(labelsize=fs)
    ax3.legend(loc='best', fontsize = fs)
    
    
    ax4 = fig2.add_subplot(223)
    ax4.plot(freq_THz, pulse_FD.real, '-o', label='real', color ='r', linewidth = 2)
    ax4.plot(freq_THz, pulse_FD.imag, '-o', label='imag', color ='b', linewidth = 2)
    ax4.plot(freq_THz, abs(pulse_FD), '--o', label = 'abs' , color = 'k', linewidth = 2) 
    ax4.set_ylabel('spectral amplitude [a.u.]', fontsize = fs)
    ax4.set_xlabel('frequency [THz]', fontsize = fs)
    ax4.tick_params(labelsize=fs)
    ax4.legend(loc='best', fontsize = fs)

    
    ax5 = fig2.add_subplot(224)
    ax5.plot(freq_wn/10**3, pulse_FD.real, '-o', label='real', color ='r', linewidth = 2)
    ax5.plot(freq_wn/10**3, pulse_FD.imag, '-o', label='imag', color ='b', linewidth = 2)
    ax5.plot(freq_wn/10**3, abs(pulse_FD), '--o', label = 'abs' , color = 'k', linewidth = 2) 
    #ax5.set_ylabel('amplitude [a.u.]', fontsize = fs)
    ax5.set_xlabel('wavenumber [$10^3$cm${}^{-1}$]', fontsize = fs)
    ax5.tick_params(labelsize=fs)
    ax5.legend(loc='best', fontsize = fs)

    # temporal intensity

    fig3 = plt.figure(file_name2 + '_TD_intensity', figsize=(10,6))
    ax6 = fig3.add_subplot(111)
    ax6.plot(t, temporal_intens, '-o', label = 'intensity' , color = 'k', linewidth = 2) 
    ax6.set_ylabel('temporal intensity [a.u.]', fontsize = fs)
    ax6.set_xlabel('time [fs]', fontsize = fs)
    ax6.tick_params(labelsize=fs)
    ax6.legend(loc= 'best', fontsize = fs/1.5)
    
    # spectral intensity
    fig4 = plt.figure(file_name2 + '_FD_intensity', figsize=(10,10))
    ax7 = fig4.add_subplot(221)
    ax7.plot(freq_eV, spectral_intens, '-o', label='intensity', color ='k', linewidth = 2)
    ax7.set_ylabel('spectral intensity [a.u.]', fontsize = fs)
    ax7.set_xlabel('energy [eV]', fontsize = fs)
    ax7.tick_params(labelsize=fs)
  
    ax8 = fig4.add_subplot(222)
    ax8.plot(freq_nm, spectral_intens, '-o', label='intensity', color ='k', linewidth = 2)
    ax8.set_xlabel('wavelength [nm]', fontsize = fs)
    ax8.tick_params(labelsize=fs)
  
    ax9 = fig4.add_subplot(223)
    ax9.plot(freq_THz, spectral_intens, '-o', label='intensity', color ='k', linewidth = 2)
    ax9.set_ylabel('spectral intensity [a.u.]', fontsize = fs)
    ax9.set_xlabel('frequency [THz]', fontsize = fs)
    ax9.tick_params(labelsize=fs)
    
    
    ax10 = fig4.add_subplot(224)
    ax10.plot(freq_wn/10**3, spectral_intens, '-o', label='intensity', color ='k', linewidth = 2)
    ax10.set_xlabel('wavenumber [$10^3$cm${}^{-1}$]', fontsize = fs)
    ax10.tick_params(labelsize=fs)


# evaluation 1D scan data

# calculate FWHM of Gaussian window for Fourier Trafo
fwhmGauss = np.sqrt(-4*np.log(2)*tau.size**2/np.log(drop))
# multiply time dat with Gaussian
time_dat_g = time_dat*Gaussian_1D(time_dat.shape[0], fwhm = fwhmGauss/1.5, center= len(tau)/2)

#Fourier Trafo (V1)

N = zero_padding(tau,1)

freq_raw = fftpack.fftfreq(N,d = step_size)
if step_size<0:
    freq_raw = freq_raw[::-1] # reverse array
cut_index = np.argmin(freq_raw)
freq_raw = redistribute(freq_raw, cut_index)
freq_raw = (freq_raw*10**15*spc.h)/spc.e # in eV, not shifted yet by pulse carrier frequency
freq_raw_shifted = (freq_raw + wp) # shift by wp, still in eV
freq = (freq_raw_shifted*spc.e/(spc.h*spc.c))/100.0 # in cm^-1
FFT_dat = redistribute(fftpack.fft(time_dat, n=N), cut_index) # sortiere Array so, dass passend zu freq array
FFT_dat_g = redistribute(fftpack.fft(time_dat_g, n=N), cut_index) # sortiere Array so, dass passend zu freq array

# plotting of simulated 1D scan

#time domain data

fig = plt.figure('TD_1Dscan', figsize=(10,6))
plt.plot(tau,time_dat, 'k', lw=2)
plt.ylabel('population', fontsize = fs)
plt.xlabel('tau [fs]', fontsize = fs)
plt.tick_params(labelsize =fs)

fig = plt.figure('TDg_1Dscan', figsize=(10,6))
plt.plot(tau,time_dat_g, 'k', lw=2)
plt.ylabel('population', fontsize = fs)
plt.xlabel('tau [fs]', fontsize = fs)
plt.tick_params(labelsize =fs)

fig = plt.figure('FFTg_1Dscan_window', figsize=(10,6))
plt.plot(freq,abs(FFT_dat_g), 'r', lw=2)
plt.ylabel('amplitude [a.u.]', fontsize = fs)
plt.xlabel('wavenumber [cm^-1]', fontsize = fs)
for i in np.arange(1,6,1):
    plt.axvline(x=freq_trans-(freq_trans-laser_freq_wn)*(5-i)/5, color = 'r', linestyle = '--')
plt.axvline(x=laser_freq_wn, color = 'k', linestyle = '--')
plt.tick_params(labelsize =fs)
plt.xlim([190000,192000])


# ======== IDFT of AC to get harmonic components (DÒESNT REALLY WORK!)
#cent = 5.5
#sigma =3
#
#Z_theo = np.array([]) #contains the individual UNsaturated harmonic demodulated AC Z components
#windows = np.array([]) #windows used for IDFT in each harmoinic
#harmonics = [1,2,3,4,5]
#
#for n in harmonics:
#    idft, window = idft_window(FFT_dat_g,cent,sigma,n,len(FFT_dat_g),order=4)
#    Z_theo = np.append(Z_theo,idft)
#    windows = np.append(windows,window)
#Z_theo = Z_theo.reshape((-1,int(len(Z_theo)/len(harmonics)))).transpose()
#windows = windows.reshape((-1,int(len(windows)/len(harmonics)))).transpose()  
#
#plt.figure()
#plt.plot(np.abs(FFT_dat_g)/max(abs(FFT_dat_g)))
#for n in harmonics:
#    plt.plot(windows[::,n-1])
#    
#plt.figure()
#for n in harmonics:
#    plt.plot(Z_theo[::,n-1],label = n)
#    plt.legend()
##

# =============================================================================
# Demodulation of phase cycling
# =============================================================================

plt.figure('populations')
for i in scan_nr:
    plt.plot(tau,scans[i],label=i)
    plt.legend()

plt.figure('population at delay vs. laboratory time')
plt.plot(np.tile(scans[::,80],4),'o-')

RvsTau = np.zeros((harm,len(tau)))

tile = 100
for i in np.arange(len(tau)):
    dft_t = np.abs(np.fft.rfft(np.tile(scans[::,i],tile)))
    for h in np.arange(0,harm):
        RvsTau[h,i] = dft_t[h*100]
        

#plt.plot(dft_t,'o-')
for i in np.roll(np.arange(0,harm),-1):
    plt.plot(RvsTau[i],label='{}H'.format(i))
    plt.legend()












