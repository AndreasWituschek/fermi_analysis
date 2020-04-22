# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:48:24 2019

@author: Primoz
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.special import jv

#File in/out
i = 5
run =  236
delay_zero_pos = 11025.66

#define constants
c=299792458
planck=6.62607004e-34
eV=1.60217662E-19

'''plotting parameters'''
plot_seed=1
plot_FEL=1
plot_bunch=False
plot_FELspec=True

tminplot = -200e-15 #time range for plotting
tmaxplot = 200e-15

#plt.close('all')

'''File I/O'''
#ldm_file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/Run_' + str(run) +'/rawdata/'.format(int(run))
ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_' + str(run) +'/rawdata/'.format(int(run))

ldm_file = os.listdir(ldm_file_path)[0]
ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')

dat_Wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
print('padres_wavelength: {}'.format(dat_Wavelength))

dat_WavelengthSpan = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
print('padres_span: {}'.format(dat_WavelengthSpan))

dat_Pixel2micron = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Pixel2micron'])

spectrum = np.array(ldm_data['photon_diagnostics']['Spectrometer']['hor_spectrum'])

ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])

ldm_data.close()

spectrum = spectrum[i]
spectrum = spectrum/max(spectrum).astype(float) #normalize padres spectrum
lambda_use = dat_Wavelength + (np.arange(1,1001)-500)*dat_Pixel2micron*1E-3*dat_WavelengthSpan #calculated wavelength


E0 = 1.298E9*eV #electron beam nominal energy
sigmaE = 100E3*eV #electron beam energy spread
R56 = 75E-6 #dispersive strength
ebeamquadchirp = 1.9E6*eV/(1E-12)**2 #electron beam quadratic chirp
ebeamcubechirp = 11.8E6*eV/(1E-12)**3 #electron beam cubic chirp

#ebeamquadchirp = 0*eV/(1E-12)**2 #electron beam quadratic chirp
#ebeamcubechirp = 0*eV/(1E-12)**3 #electron beam cubic chirp


lambdaseed = 261.3E-9 #seed laser wavelength
k1 = 2*np.pi/lambdaseed #seed laser wave number
omegaseed = 2.*np.pi*c/lambdaseed #seed laser angular frequency

n = 5 #harmonic number
lambdaFEL = lambdaseed/n #FEL central wavelength


tau10 = 141E-15 #first seed transform-limited pulse duration (annotation: E-field!!)
dlambda1 = (c/(c/lambdaseed)**2)*(0.44/tau10) #transform-limited bandwidth
GDD1 = 0E-27 #first seed linear frequency (quadratic phase) chirp
tau1 = np.sqrt(1+(4*np.log(2)*GDD1/tau10**2)**2)*tau10 #first seed pulse duration


tau20 = tau10 # second seed transform-limited pulse duration
dlambda2 = (c/(c/lambdaseed)**2)*(0.44/tau20) #transform-limited bandwidth
GDD2 = 0E-27 #second seed linear frequency (quadratic phase) chirp
tau2 = np.sqrt(1+(4*np.log(2)*GDD2/tau20**2)**2)*tau20 #second seed pulse duration

'''SEED TIMING and Phase'''
#deltat = 250e-15 #separation between the seeds
deltat = (ldm_delay-delay_zero_pos)*1e-15
print(deltat*1e15)
deltaphi = 1.5*np.pi #relative phase between the seeds
ebeamtiming = -200e-15 #relative timing between the electron beam and the two seeds


tmin = -40*max(tau1,tau2) + deltat/2 #time window for calculations
tmax = 40*max(tau1,tau2) + deltat/2
npoints = 500000 #number of points in the time window
step = (tmax-tmin)/npoints #step
time = np.linspace(tmin,tmax,npoints)
z = time*c #longitudinal coordinate in the electron beam

tminplot = -200e-15 #time range for plotting
tmaxplot = 500e-15

#
#Psi1 = (1./(2.*GDD1+(tau10**4.)/(8.*(np.log(2.)**2.)*GDD1)))*time**2. #seed phases accounting for the seed chirps
#Psi2 = (1./(2.*GDD2+(tau20**4.)/(8.*(np.log(2.)**2.)*GDD2)))*(time-deltat)**2.


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
#seedenvelope = seedfield.real #seed envelope
seedphase = np.unwrap(np.angle(seedfield)) #seed phase
seedphase_fast = np.unwrap(np.angle(seedfield_fast))



A0 = 4 #amplitude of the energy modulation of the electron beam induced by the seeds
A = A0*np.sqrt(seedenvelope) #temporal profile of the energy modulation of the electron beam
B = R56*k1*sigmaE/E0 #normalized dispersive strength
print 'dispersive strength = {}'.format(B)

#electron beam energy profile defined as the nominal energy E0 plus quadratic and cubic terms
ebeamenergyprofile = (E0+ebeamquadchirp*(time-ebeamtiming)**2+ebeamcubechirp*(time-ebeamtiming)**3)
ebeamphase = (B/sigmaE)*ebeamenergyprofile #electorn beam energy profile induces a phase onto the FEL pulse

#bunching (proportional to the FEL electric field) in the time domain
btime = np.exp(-(n*B)**2/2.)*jv(n, -n*B*A)*np.exp(1j*n*seedphase)*np.exp(1j*n*ebeamphase)
btime_fast = np.exp(-(n*B)**2/2.)*jv(n, -n*B*A)*np.exp(1j*n*seedphase_fast)*np.exp(1j*n*ebeamphase)
#plt.plot(btime_fast)
maxbunching = max(abs(btime)) #calculates maximum bunching, this depends on the energy modulation amplitude A and R56, should be several percent to correspond to realistic experimental conditions
print 'maxbunching = {}'.format(maxbunching)
FELint = sum(abs(btime)**2) #total FEL intensity in a.u.

#FELtime = abs(btime)**2/max(abs(btime)**2) #normalized FEL intensity profile in the time domain
FELtime = abs(btime)**2 #un-normalized FEL intensity profile in the time domain
FELphase = np.unwrap(np.angle(btime)) # FEL phase



#transformation into the spectral domain
seedspectralenvelope = np.fft.fftshift(np.fft.fft(np.fft.fftshift(seedfield))) # seed spectrum
seedspectralenvelope = abs(seedspectralenvelope)**2/max(abs(seedspectralenvelope)**2) # normalize the seed spectrum

seedspectralenvelope_fast = np.fft.fftshift(np.fft.fft(np.fft.fftshift(seedfield_fast))) # seed spectrum
seedspectralenvelope_fast = abs(seedspectralenvelope_fast)**2/max(abs(seedspectralenvelope_fast)**2) # normalize the seed spectrum
#plt.plot(seedspectralenvelope_fast[259270:259450])

FELfreq = abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(btime))))**2 #FEL spectral envelope
#FELfreq = FELfreq/max(FELfreq) #normalize the FEL spectral enevelope
FELfreq = FELfreq/1E6 #normalize the FEL spectral enevelope
timerange = max(time)-min(time) # time range for FFT
#sizetime =size(time)
#freqrange = c/(lambdaseed/n ) + (-(0.5)*(1/timerange)*(sizetime(2)-1):(1/timerange):(0.5)*(1/timerange)*(sizetime(2)-1)) # calculated frequency range for FFT
sizetime = len(time)
freqrange = c/(lambdaseed/n)+np.linspace(-(0.5)*1/timerange*(sizetime-1),(0.5)*(1/timerange)*(sizetime-1)+1/timerange,npoints)

lambdarange = c/freqrange #calculated wavelength range for FFT

'''Plotting Results'''
'''Seed'''
if plot_seed:
    fsimages=16 # font size in images

    fig, ax = plt.subplots(3)
    fig.suptitle('SEED')

    ax[0].plot(time*1e15, seedenvelope, linewidth=2) #plots the seed envelope in the time domain
    #set(gca,'LineWidth',2)
    #set(gca,'FontSize',fsimages)
    ax[0].set_xlabel('TIME (fs)')
    ax[0].set_ylabel('SEED INTENSITY (a.u.)')
    ax[0].set_xlim(tminplot*1e15, tmaxplot*1e15)

    ax[1].plot(time*1e15,seedphase,linewidth=2)#plots the seed phase in the time domain
    #set(gca,'LineWidth',2)
    #set(gca,'FontSize',fsimages)
    ax[1].set_xlabel('TIME (fs)')
    ax[1].set_ylabel('SEED PHASE')
    ax[1].set_xlim(tminplot*1e15, tmaxplot*1e15)

    ax[2].plot(lambdarange*1e9*n,seedspectralenvelope,linewidth=2) #plots the seed envelope in the wavelength domain
    #ax[2].set_xlim((0.997)*lambdaseed/n*1e9, (1.003)*lambdaseed/n*1e9)
    ax[2].set_xlim(260.5, 262)
    #set(gca,'LineWidth',2)
    #set(gca,'FontSize',fsimages)
    ax[2].set_xlabel('WAVELENGTH (nm)')
    ax[2].set_ylabel('SEED SPECTRUM (a.u.)')
    plt.tight_layout()
    #
'''FEL'''
if plot_FEL:
    fig, ax = plt.subplots(3)
    fig.suptitle('FEL')

    ax[0].plot(time*1e15,FELtime,linewidth=2) #plots the FEL envelope in the time domain
    #set(gca,'LineWidth',3)
    #set(gca,'FontSize',fsimages)
    ax[0].set_xlabel('TIME (fs)')
    ax[0].set_ylabel('FEL intensity (a.u.)')
    ax[0].set_xlim(tminplot*1e15, tmaxplot*1e15)

    ax[1].plot(time*1e15,FELphase,linewidth=2) #plots the FEL phase in the time domain
    #set(gca,'LineWidth',2)
    #set(gca,'FontSize',fsimages)
    ax[1].set_xlabel('TIME (fs)')
    ax[1].set_ylabel('FEL PHASE')
    ax[1].set_ylim(-5,5*6.)
    ax[1].set_xlim(tminplot*1e15, tmaxplot*1e15)

    ax[2].plot(lambdarange*1e9,FELfreq,linewidth=2) #%plots the FEL envelope in the time domain
    #ax[2].set_xlim([(0.997)*lambdaseed/n*1e9 (1.003)*lambdaseed/n*1e9])
    #ax[2].set_xlim([-50/1000+c/(lambdaseed/n)*planck/eV 50/1000+c/(lambdaseed/n)*planck/eV])
    ax[2].set_xlim(52, 52.4)
    #set(gca,'LineWidth',2)
    #set(gca,'FontSize',fsimages)
    ax[2].set_xlabel('WAVELENGTH (nm)')
    ax[2].set_ylabel('FEL spectrum (a.u.)')
    plt.tight_layout()

if plot_FELspec:
    plt.figure('FEL spectrum')
    plt.plot(lambda_use,spectrum,lambdarange*1e9,FELfreq,linewidth=2) #plot the experimental and calculated FEL spectra onto the same plot
    plt.legend(['experiment','calculation'])
    plt.xlim(52., 52.4) #wavelength range in plots
    plt.xlabel('WAVELENGTH (nm)')
    plt.ylabel('FEL spectrum (a.u.)')

if plot_bunch:
    fig, ax = plt.subplots(1)
    ax.plot(time*1e12,ebeamenergyprofile/eV/1e9,linewidth=2) #shows the electron beam energy profile
    ax.set_xlim(-0.8+ebeamtiming, 0.8+ebeamtiming)
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('electron beam energy (GeV)')






#%
#% for kl=1:200
#% ebeamtiming=-1000e-15+kl*10e-15
#%
#% ebeamenergyprofile=(E0+ebeamquadchirp*(time+ebeamtiming).^2+ebeamcubechirp*(time+ebeamtiming).^3);
#% ebeamphase=(B/sigmaE)*ebeamenergyprofile;
#% %ebeamenergyprofile=0;
#%
#% btime=exp(-(n*B)^2/2).*besselj(n, -n*B*A).*exp(1i*n*seedphase).*exp(1i*n*ebeamphase);
#% maxbunching=max(abs(btime));
#%
#% FELint=sum(abs(btime).^2);
#%
#% FELtime=abs(btime).^2/max(abs(btime).^2);
#%
#% FELfreq=abs(fftshift(fft(fftshift(btime)))).^2;
#% FELfreq=FELfreq/max(FELfreq);
#%
#% timerange=max(time)-min(time);
#% sizetime=size(time);
#% freqrange=c/(lambdaseed/n)+(-(0.5)*(1/timerange)*(sizetime(2)-1):(1/timerange):(0.5)*(1/timerange)*(sizetime(2)-1));
#% lambdarange=c./freqrange;%-lambdaseed/n;
#%
#% spectraumdelayscan(kl,:)=FELfreq(:);
#% ebeamtimingscan(kl)=ebeamtiming;
#%
#% end
#%
#% figure(433)
#% surfc(ebeamtimingscan*1e12,lambdarange*1e9,spectraumdelayscan')
#% view(2)
#% shading interp
#% ylim([50.4 51.2])
#% xlim([-0.75 0.75])
#
#
#deltaphiscan=0;
#doscanphase=0;
#doscanduration=0;
#
#
#npointsphase=100;
#npointstime=100;
#scandeltaphi=(0:2*pi/npointsphase:4*pi);
#scanduration=(10e-15:180e-15/npointstime:400e-15);
#sizescandeltaphi=size(scandeltaphi);
#sizescanduration=size(scanduration);
#
#% if doscanphase==1
#%
#%
#% for ijk=1:sizescandeltaphi(2)
#%
#%     deltaphi=scandeltaphi(ijk);
#% seedfield=(C1*exp(-2*log(2)*time.^2/tau1^2).*exp(1i*Psi1)+C2*exp(-2*log(2)*(time-deltat).^2/tau2^2).*exp(1i*Psi2)*exp(1i*deltaphi));
#% seedenvelope=abs(seedfield).^2;
#% seedphase=unwrap(angle(seedfield));
#%
#%
#% seedspectralenvelope=fftshift(fft(fftshift(seedfield)));
#% seedspectralenvelope=abs(seedspectralenvelope).^2/max(abs(seedspectralenvelope).^2);
#%
#% scandeltaphiseedspectrum(ijk,:)=seedspectralenvelope;
#%
#% A=A0*sqrt(seedenvelope);
#%
#%
#% btime=besselj(n, -n*B*A).*exp(1i*n*seedphase);
#% %maxbunchingscan=max(abs(btime))
#% FELintensity(ijk)=sum(abs(btime).^2);
#% FELtime=abs(btime).^2/max(abs(btime).^2);
#% FELphase=unwrap(angle(btime));
#%
#% FELfreq=abs(fftshift(fft(fftshift(btime)))).^2;
#% FELfreq=FELfreq/max(FELfreq);
#%
#% timerange=max(time)-min(time);
#% sizetime=size(time);
#% freqrange=c/(lambdaseed/n)+(-(0.5)*(1/timerange)*(sizetime(2)-1):(1/timerange):(0.5)*(1/timerange)*(sizetime(2)-1));
#% lambdarange=c./freqrange-lambdaseed/n;
#%
#% fsimages=16;
#%
#% scandeltaphispectrum(ijk,:)=FELfreq;
#%
#% end
#%
#% figure(333)
#% imagesc(scandeltaphi,(freqrange-c/(lambdaseed/n))*planck/eV*1e3,scandeltaphispectrum')
#% ylim([-50 50])
#% set(gca,'LineWidth',2);
#% set(gca,'FontSize',fsimages)
#% xlabel('FEL PHASE DIFFERENCE');
#% ylabel('ENERGY DIFFERENCE (meV)');
#%
#% figure(888)
#% imagesc(scandeltaphi,(freqrange-c/(lambdaseed/n))*planck/eV*1e3,scandeltaphiseedspectrum')
#% ylim([-50 50])
#% set(gca,'LineWidth',2);
#% set(gca,'FontSize',fsimages)
#% xlabel('PHASE DIFFERENCE');
#% ylabel('ENERGY DIFFERENCE (meV)');
#%
#% figure(333111)
#% plot(scandeltaphi,FELintensity,'LineWidth',2)
#% set(gca,'LineWidth',2);
#% set(gca,'FontSize',fsimages)
#% xlabel('PHASE DIFFERENCE');
#% ylabel('FEL INTENSITY (a.u.)');
#%
#% end
#%
#%
#%
#%
#%
#%
#% if doscanduration==1
#%   deltaphi=deltaphiscan;
#%   clear FELfreq;
#% for ijk=1:sizescanduration(2)
#%
#%     tau1=scanduration(ijk);
#%     tau2=tau1;
#%
#% seedfield=(C1*exp(-2*log(2)*time.^2/tau1^2).*exp(1i*Psi1)+C2*exp(-2*log(2)*(time-deltat).^2/tau2^2).*exp(1i*Psi2)*exp(1i*deltaphi));
#% seedenvelope=abs(seedfield).^2;
#% seedphase=unwrap(angle(seedfield));
#%
#%
#% A=A0*sqrt(seedenvelope);
#%
#%
#% btime=besselj(n, -n*B*A).*exp(1i*n*seedphase);
#% FELintensity(ijk)=sum(abs(btime).^2);
#% FELtime=abs(btime).^2/max(abs(btime).^2);
#% FELphase=unwrap(angle(btime));
#%
#% FELfreq=abs(fftshift(fft(fftshift(btime)))).^2;
#% FELfreq=FELfreq/max(FELfreq);
#%
#% timerange=max(time)-min(time);
#% sizetime=size(time);
#% freqrange=c/(lambdaseed/n)+(-(0.5)*(1/timerange)*(sizetime(2)-1):(1/timerange):(0.5)*(1/timerange)*(sizetime(2)-1));
#% lambdarange=c./freqrange-lambdaseed/n;
#%
#% fsimages=16;
#%
#%
#% scandurationspectrum(ijk,:)=FELfreq;
#%
#% end
#%
#% figure(444)
#% imagesc(scanduration*1e15,(freqrange-c/(lambdaseed/n))*planck/eV*1e3,scandurationspectrum')
#% ylim([-50 50])
#% set(gca,'LineWidth',2);
#% set(gca,'FontSize',fsimages)
#% xlabel('SEED DURATION (fs)');
#% ylabel('ENERGY DIFFERENCE (meV)');
#%
#%
#% figure(444111)
#% plot(scanduration*1e15,FELintensity,'LineWidth',2)
#% set(gca,'LineWidth',2);
#% set(gca,'FontSize',fsimages)
#% xlabel('SEED DURATION (fs)');
#% ylabel('FEL INTENSITY (a.u.)');
#% end
#
#
