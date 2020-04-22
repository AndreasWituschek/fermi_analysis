# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:48:24 2019

@author: Primoz
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import jv
import fermi_analysis.functions as fk

#File in out
i = 1
run = 250
delay_zero_pos = 11025.66

#define constants
c=299792458
planck=6.62607004e-34
eV=1.60217662E-19

plt.close('all')
E0 = 1.298E9*eV #electron beam nominal energy
sigmaE = 100E3*eV #electron beam energy spread
R56 = 55E-6 #dispersive strength
ebeamquadchirp = 1.9E6*eV/(1E-12)**2 #electron beam quadratic chirp
ebeamcubechirp = 11.8E6*eV/(1E-12)**3 #electron beam cubic chirp

#ebeamquadchirp = 0*eV/(1E-12)**2 #electron beam quadratic chirp
#ebeamcubechirp = 0*eV/(1E-12)**3 #electron beam cubic chirp


lambdaseed = 261.09E-9 #seed laser wavelength
k1 = 2*np.pi/lambdaseed #seed laser wave number
omegaseed = 2.*np.pi*c/lambdaseed #seed laser angular frequency

n = 5 #harmonic number
lambdaFEL = lambdaseed/n #FEL central wavelength


tau10 = 80E-15 #first seed transform-limited pulse duration
dlambda1 = (c/(c/lambdaseed)**2)*(0.44/tau10) #transform-limited bandwidth
GDD1 = 0E-27 #first seed linear frequency (quadratic phase) chirp
tau1 = np.sqrt(1+(4*np.log(2)*GDD1/tau10**2)**2)*tau10 #first seed pulse duration


tau20 = tau10 # second seed transform-limited pulse duration
dlambda2 = (c/(c/lambdaseed)**2)*(0.44/tau20) #transform-limited bandwidth
GDD2 = 0E-27 #second seed linear frequency (quadratic phase) chirp
tau2 = np.sqrt(1+(4*np.log(2)*GDD2/tau20**2)**2)*tau20 #second seed pulse duration


deltaphi = 0. #relative phase between the seeds
ebeamtiming = -300e-15 #relative timing between the electron beam and the two seeds


start = -250e-15
stop = 250e-15
points = 12000

#tmin = -40*max(tau1,tau2) + stop/2 #time window for calculations
#tmax = 40*max(tau1,tau2) + stop/2
tmin = -800e-15 #time window for calculations
tmax = 800e-15
npoints = 100000 #number of points in the time window
step = (tmax-tmin)/npoints #step
time = np.linspace(tmin,tmax,npoints)
z = time*c #longitudinal coordinate in the electron beam


t = np.array([])
cc_t = np.array([])
for deltat in np.linspace(start,stop,points):
#    deltat = 150.4e-15 #separation between the seeds
    print(deltat*1e15)
    #seed phases accounting for the seed chirps
    Psi1= (1./(2.*GDD1+(tau10**4.)/(8.*(np.log(2.)**2.)*GDD1)))*time**2.
    Psi2 = (1./(2.*GDD2+(tau20**4.)/(8.*(np.log(2.)**2.)*GDD2)))*(time-deltat)**2.

    '''with fast phase'''
    #seed phases accounting for the seed chirps without fast oscillating omegaseed frequency
    Psi1_fast = (1./(2.*GDD1+(tau10**4.)/(8.*(np.log(2.)**2.)*GDD1)))*time**2.+omegaseed*time
    Psi2_fast = (1./(2.*GDD2+(tau20**4.)/(8.*(np.log(2.)**2.)*GDD2)))*(time-deltat)**2.+omegaseed*(time-deltat)


    C1 = 1. #relative seed amplitudes
    C2 = 1.


    seedfield_fast = (C1*np.exp(-2.*np.log(2)*time**2/tau1**2)*np.exp(1j*Psi1_fast) + C2*np.exp(-2*np.log(2)*(time-deltat)**2/tau2**2)*np.exp(1j*Psi2_fast)*np.exp(1j*deltaphi)) #seed electric field; first seed centered at time=0 fs
    seedfield= (C1*np.exp(-2.*np.log(2)*time**2/tau1**2)*np.exp(1j*Psi1) + C2*np.exp(-2*np.log(2)*(time-deltat)**2/tau2**2)*np.exp(1j*Psi2)*np.exp(1j*deltaphi)) #seed electric field; first seed centered at time=0 fs

    seedenvelope = abs(seedfield)**2 #seed envelope
    seedenvelope_fast = abs(seedfield_fast)**2
    #seedenvelope = seedfield.real #seed envelope
    seedphase = np.unwrap(np.angle(seedfield)) #seed phase
    seedphase_fast = np.unwrap(np.angle(seedfield_fast))


    A0 = 4. #amplitude of the energy modulation of the electron beam induced by the seeds
    A = A0*np.sqrt(seedenvelope_fast) #temporal profile of the energy modulation of the electron beam
    B = R56*k1*sigmaE/E0 #normalized dispersive strength

    #electron beam energy profile defined as the nominal energy E0 plus quadratic and cubic terms
    ebeamenergyprofile = (E0+ebeamquadchirp*(time-ebeamtiming)**2+ebeamcubechirp*(time-ebeamtiming)**3)
    ebeamphase = (B/sigmaE)*ebeamenergyprofile #electorn beam energy profile induces a phase onto the FEL pulse

    #bunching (proportional to the FEL electric field) in the time domain
    btime = np.exp(-(n*B)**2/2.)*jv(n, -n*B*A)*np.exp(1j*n*seedphase)*np.exp(1j*n*ebeamphase)
    btime_fast = np.exp(-(n*B)**2/2.)*jv(n, -n*B*A)*np.exp(1j*n*seedphase_fast)*np.exp(1j*n*ebeamphase)
    #plt.plot(btime_fast)
    crosscorr = np.conjugate(btime_fast)*btime_fast
    crosscorr = sum(crosscorr)
    cc_t = np.append(cc_t,crosscorr)
    t = np.append(t,deltat*1e15)
#    plt.plot(crosscorr)

plt.close('all')
plt.plot(t,cc_t)
plt.xlabel('delay [fs]')
plt.close('all')
cc_dft = np.fft.fftshift(np.fft.fft(np.fft.fftshift(cc_t)))
freqs = np.fft.fftshift(np.fft.fftfreq(points,d=(stop-start)/points))
plt.plot(freqs/1.14785e15,abs(cc_dft))
plt.xlabel('frequency [1.1478 PHz]')


#maxbunching = max(abs(btime)) #calculates maximum bunching, this depends on the energy modulation amplitude A and R56, should be several percent to correspond to realistic experimental conditions
#
#FELint = sum(abs(btime)**2) #total FEL intensity in a.u.
#
#FELtime = abs(btime)**2/max(abs(btime)**2) #normalized FEL intensity profile in the time domain
#FELphase = np.unwrap(np.angle(btime)) # FEL phase
#
#
#
##transformation into the spectral domain
#seedspectralenvelope = np.fft.fftshift(np.fft.fft(np.fft.fftshift(seedfield))) # seed spectrum
#seedspectralenvelope = abs(seedspectralenvelope)**2/max(abs(seedspectralenvelope)**2) # normalize the seed spectrum
#
#seedspectralenvelope_fast = np.fft.fftshift(np.fft.fft(np.fft.fftshift(seedfield_fast))) # seed spectrum
#seedspectralenvelope_fast = abs(seedspectralenvelope_fast)**2/max(abs(seedspectralenvelope_fast)**2) # normalize the seed spectrum
#plt.plot(seedspectralenvelope_fast[259270:259450])
#
#FELfreq = abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(btime))))**2 #FEL spectral envelope
#FELfreq = FELfreq/max(FELfreq) #normalize the FEL spectral enevelope
#timerange = max(time)-min(time) # time range for FFT
##sizetime =size(time)
##freqrange = c/(lambdaseed/n ) + (-(0.5)*(1/timerange)*(sizetime(2)-1):(1/timerange):(0.5)*(1/timerange)*(sizetime(2)-1)) # calculated frequency range for FFT
#sizetime = len(time)
#freqrange = c/(lambdaseed/n)+np.linspace(-(0.5)*1/timerange*(sizetime-1),(0.5)*(1/timerange)*(sizetime-1)+1/timerange,npoints)
#
#lambdarange = c/freqrange #calculated wavelength range for FFT
