# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 16:19:21 2019

@author: Andreas
"""

'''Program simulates Fano profile in TD and FD '''

import numpy as np
import fermi_analysis.functions as fk 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as spc
import pandas as pd

#plt.close('all')
#Parameters:
undersamp = True #undersampling or not
harmonic = 6
l_ref = harmonic*4.5607 #eV

#fano parameters
q=-0.2
E_r = 28.506 #eV
gamma = 0.016 #eV
E_fano = np.linspace(0,E_r+800,200000) # frequency domain energy values

if undersamp:
    E_r = E_r - l_ref #eV
    E_fano = np.linspace(0,E_r+70,800000)

fano_resonance = fk.fano(E_fano,E_r,q, gamma, 1.,0.)


fano_td_theo = np.fft.fftshift(np.fft.fft(fano_resonance))
fano_td_time = 1e15*spc.h/spc.e* np.fft.fftshift(np.fft.fftfreq(len(fano_resonance),E_fano[1]-E_fano[0]))

''' Plotting results '''
fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(311)
ax.set_title('undersampled fano q=' + str(q))
ax.plot(E_fano,fano_resonance)
if undersamp: 
    ax.set_xlabel('undersampled energy [eV]')
    ax.set_xlim(0,1)
else: 
    ax.set_xlabel('energy [eV]')
    ax.set_xlim(27.8,29)


ax = fig.add_subplot(312)
#ax.plot(fano_td_time,np.abs(fano_td_theo)/np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+20:])))
ax.set_title('undersampled pump probe transient')
ax.plot(fano_td_time,fano_td_theo.real)
ax.plot(fano_td_time,fano_td_theo.imag)
ax.plot(fano_td_time,abs(fano_td_theo))
ax.legend(['real','imaginary','abs'])
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.')
ax.set_xlim(-10,500)
#ax.set_ylim(-100,100)
#ax.set_ylim(np.min(np.abs(fano_td_theo)),np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+1:])))

ax = fig.add_subplot(313)
#ax.plot(fano_td_time,np.abs(fano_td_theo)/np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+20:])))
ax.set_title('undersampled pump probe transient')
ax.plot(fano_td_time,fano_td_theo.real)
ax.plot(fano_td_time,fano_td_theo.imag)
ax.plot(fano_td_time,abs(fano_td_theo))
ax.legend(['real','imaginary','abs'])
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.')
ax.set_xlim(-30,30)


plt.tight_layout()

#fig = plt.figure()
#ax = fig.add_subplot(111)
##ax.plot(abs(np.fft.fftshift(np.fft.fft(fano_td_theo.imag))))
#ax.plot(abs(np.fft.fftshift(np.fft.fft(fano_td_theo[51000:]))))


##writing the time domain profile of the resonance into a file:
#t=fano_td_time[50000:58000]
#Ir=fano_td_theo[50000:58000].real
#Ii=fano_td_theo[50000:58000].imag
#save_file = np.array([t,Ir,Ii])
##np.savetxt("run4644_padres.csv", save_file, delimiter="\t")
#pd.DataFrame(save_file.transpose()).to_csv("fano_theo_TD.csv",delimiter="\t",index=False)
