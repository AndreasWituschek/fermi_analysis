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
run = 2829 #first run of delay scan
demod = ['0','1','2']
device = ['dev3269','dev3265']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'

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
[ax.axvline(_x,color='black') for _x in [13.5,23,9.5,4,17.5,19]]
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('intensity [a.u.]')
ax.set_xlim([-2,27])


i0_0H = np.sum(i0_fft[::,395:405],axis=1)
i0_1H = np.sum(i0_fft[::,458:468],axis=1)
i0_2H = np.sum(i0_fft[::,522:532],axis=1)
i0_3H = np.sum(i0_fft[::,585:595],axis=1)
i0_4H = np.sum(i0_fft[::,648:658],axis=1)
i0_5H = np.sum(i0_fft[::,711:721],axis=1)
i0_sum15 = i0_1H+i0_2H+i0_3H+i0_4H+i0_5H
i0_sum05 = i0_0H+i0_1H+i0_2H+i0_3H+i0_4H+i0_5H


figI0 = plt.figure('i0 Fourier transform components')
ax = figI0.add_subplot(111)
ax.plot(T,i0_1H)
ax.plot(T,i0_2H)
ax.plot(T,i0_3H)
ax.plot(T,i0_4H)
ax.plot(T,i0_5H)
ax.plot(T,i0_sum15)
ax.plot(T,i0_sum05,color = 'brown')
ax.plot(T,i0_0H)
ax.legend(['1H','2H','3H','4H','5H','sum 1-5','sum 0-5','0H'],loc='lower right')
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

        Phi = np.angle(Z, deg=True)
        
        Zh = np.append(Zh,Z,axis=0)
        legend =np.append(legend,i)
        Phih =  np.append(Phih,Phi)
        

order = [3,4,5,0,1,2] #orders harmonics in ascending order        
Zh = Zh.reshape((len(device)*len(demod),-1)).transpose()[::,order]
Phih = Phih.reshape((len(device)*len(demod),-1)).transpose()[::,order]
legend = legend[order]
        
figTD = plt.figure('R of demodulated signal')
ax = figTD.add_subplot(111)
ax.set_title('R of demodulated signal')
ax.plot(T, Zh[::,4].real)
ax.legend(np.append(legend,'sum 1-6'))
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.]')


#==============================================================================
# export data into csv
#==============================================================================

df=pd.DataFrame({'delay': T, 'X': Zh[::,4].real})
df.to_csv('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_2828/Z_5H_vs_delay.csv',index=False)
