# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:46:56 2019

@author: FemtoMeasure
"""

import numpy as np
import fermi_analysis.functions as fk 
import h5py
from scipy import signal
import matplotlib.pyplot as plt
import scipy.constants as spc
from matplotlib.gridspec import GridSpec
import pandas as pd

plt.close('all')
""" Data parameters """
run = 5234 #first run of delay scan
demod = ['0','1','2']
device = ['dev3265','dev3269']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
#root_file_path = '//10.5.71.1/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
unwrap = False
export_6H = False
export_Z = True
phasing = False #decides wether phasing according to electronics is performed
# Load preanalysed data
data = fk.ImportPreanalysedData(root_file_path,run)

T = data['LDM']['delay']

#==============================================================================
# harmonic content of i0
#==============================================================================
i0_fft = data['LDM']['i0_fft'] #fft of i0
i0_fft_mean =  np.mean(i0_fft,axis=0)
i0_fftfreq = np.fft.fftshift(np.fft.fftfreq(n=800,d=1/50.))
#plotting the i0 fft:
figFFT = plt.figure('i0 FFT')
ax = figFFT.add_subplot(111)
ax.plot(i0_fftfreq,i0_fft_mean)
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('intensity [a.u.]')
ax.set_xlim([-2,27]) 

i0_0H = np.sum(i0_fft[::,395:405],axis=1)
i0_1H = np.sum(i0_fft[::,443:453],axis=1)
i0_2H = np.sum(i0_fft[::,492:502],axis=1)
i0_3H = np.sum(i0_fft[::,540:550],axis=1)
i0_4H = np.sum(i0_fft[::,589:599],axis=1)
i0_5H = np.sum(i0_fft[::,637:647],axis=1)
i0_6H = np.sum(i0_fft[::,686:696],axis=1)
i0_sum16 = i0_1H+i0_2H+i0_3H+i0_4H+i0_5H+i0_6H
i0_sum06 = i0_0H+i0_1H+i0_2H+i0_3H+i0_4H+i0_5H+i0_6H

figI0 = plt.figure('i0 Fourier transform components')
ax = figI0.add_subplot(111)
ax.plot(T,i0_1H)
ax.plot(T,i0_2H)
ax.plot(T,i0_3H)
ax.plot(T,i0_4H)
ax.plot(T,i0_5H)
ax.plot(T,i0_6H)
ax.plot(T,i0_sum16)
ax.plot(T,i0_sum06,color = 'brown')
ax.plot(T,i0_0H,color = 'grey')
ax.legend(['1H','2H','3H','4H','5H','6H','sum 1-6','sum 0-6','0H'])
ax.set_title('i0 Fourier transform components')
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.]')


# loop over the mfli devices

Zh = np.array([])
legend = np.array([])
Phih = np.array([])

for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = (data[dev]['x' + d] + 1j * data[dev]['y' + d])# / data['LDM']['I0']

        Zh = np.append(Zh,Z,axis=0)        
        legend =np.append(legend,i)

order = [0,1,2,3,4,5] #orders harmonics in ascending order        
Zh = Zh.reshape((len(device)*len(demod),-1)).transpose()[::,order]

#phasing:
if phasing:
    phi_corr = np.array([180,0,180,0,180,0]) # corrects for other output of interferometer
    gamma = 10.5*(np.arange(6)+1) #corrects for phase shift of reference transmitter
    beta = 0.01*3.03*360*(np.arange(6)+1) #phase of boxcar integrator [degree] (beta in Lukas phasing routine)
    
    Zh = Zh*np.exp(1j*(gamma-beta)/180.*np.pi)
    Zh = Zh*np.exp(1j*phi_corr/180.*np.pi)

Phih = np.angle(Zh)
t0 = 0 # indicates where actual time zero is [fs].

if unwrap:
    Phih =np.unwrap(Phih, axis =0)
    phi_0 = Phih[fk.find_index(T,t0),::] #phase at 0 delay of demodulate signal at each harmonic 
    phi_0_mod = Phih[fk.find_index(T,t0),::] % (2.*np.pi)#phase at 0 delay of demodulate signal at each harmonic 
    phi_0_mod = [a if a < np.pi else a-2*np.pi for a in phi_0_mod]
    print phi_0_mod
    Phih = Phih-phi_0+phi_0_mod


figphase = plt.figure('Phase of individual harmonic')
ax.set_title('Phase of individual harmonic')
ax = figphase.add_subplot(111)
ax.plot(T,Phih,'-')
ax.legend(legend,ncol=3)
ax.set_xlabel('delay [fs]')
ax.set_ylabel('phase [rad]')
ax.set_xlim([-10,10])


Rh = np.absolute(Zh)
Rh_sum =np.sum(Rh,axis=1) # sum of all harmonic contributions
legend = legend[order]
        
figTD = plt.figure('R of demodulated signal')
ax = figTD.add_subplot(111)
ax.set_title('R of demodulated signal')
ax.plot(T, Rh)
ax.plot(T,Rh_sum)
ax.legend(np.append(legend,'sum 1-6'))
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.]')


#plt.polar((phi_0)/180.*np.pi,np.ones_like(phi_0))
#==============================================================================
# export data into csv
#==============================================================================
if export_6H:
    df=pd.DataFrame({'delay': T, 'X_6H': Zh[::,5].real})
    df.to_csv('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/X_6H_vs_delay.csv',index=False)

if export_Z:
    np.savetxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/X_vs_delay.txt',Zh.real)
    np.savetxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/Y_vs_delay.txt',Zh.imag)
