# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 16:19:21 2019

@author: Andreas
"""

'''Program simulates Fano profile in TD and FD '''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as spc
import pandas as pd

#plt.close('all')
#Parameters:
undersamp = False #undersampling or not
harmonic = 6
l_ref = harmonic*4.5607 #eV

#fano parameters
q=-0.17
E_r = 28.506 #eV
gamma = 0.012 #eV
n = 400
#E_fano = np.linspace(E_r-n,E_r+n,160000) # frequency domain 1energy values
E_fano = np.linspace(0,+n,60000) # frequency domain energy values
epsilon = (E_fano-E_r)/gamma

#creating the fano profile and preparing it for inverse fourier trafo
fano_resonance = ((q+1j)**2 + 1j*(epsilon-1j))/(epsilon - 1j) # formula for h(epsilon) as found in equation S21 of the Sup. Mat. for the 2D Fano Paper.
fano_resonance_s = np.append(fano_resonance,fano_resonance[::-1])
fano_resonance_s =  fano_resonance_s[0:len(fano_resonance_s)-1]


fano_td_theo = np.fft.fftshift(np.fft.ifft(fano_resonance_s))
fano_td_time = 1e15*spc.h/spc.e* np.fft.fftshift(np.fft.fftfreq(len(fano_resonance_s),E_fano[1]-E_fano[0]))


''' Plotting results '''
fig = plt.figure(figsize=(10, 13))

ax = fig.add_subplot(411)
ax.set_title('spectrum: fano q=' + str(q))
ax.plot(E_fano,fano_resonance.real)
ax.plot(E_fano,fano_resonance.imag)
ax.plot(E_fano,abs(fano_resonance))
ax.set_xlabel('Energy [eV]')
#if undersamp: 
#    ax.set_xlabel('undersampled energy [eV]')
#ax.set_xlim(0,1)
ax.set_xlim(27,30)
ax.legend(['real=dispersion','imaginary=absorption','abs'])

ax = fig.add_subplot(412)
#ax.plot(fano_td_time,np.abs(fano_td_theo)/np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+20:])))
ax.set_title('inverse FT of spectrum')
ax.plot(fano_td_time,fano_td_theo.real)
ax.plot(fano_td_time,fano_td_theo.imag)
ax.plot(fano_td_time,abs(fano_td_theo))
ax.legend(['real','imaginary','abs'])
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.')
ax.set_xlim(-10,300)
ylim =  1.4*max(abs(fano_td_theo[len(fano_td_theo)/2.+1::]))
#ax.set_ylim(-ylim,ylim)
#ax.set_ylim(np.min(np.abs(fano_td_theo)),np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+1:])))

ax = fig.add_subplot(413)
#ax.plot(fano_td_time,np.abs(fano_td_theo)/np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+20:])))
ax.set_title('inverse FT of spectrum (zoom 1)')
ax.plot(fano_td_time,fano_td_theo.real)
ax.plot(fano_td_time,fano_td_theo.imag)
ax.plot(fano_td_time,abs(fano_td_theo))
ax.legend(['real','imaginary','abs'])
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.')
ax.set_xlim(-10,300)
ylim =  1.4*max(abs(fano_td_theo[len(fano_td_theo)/2.+1::]))
ax.set_ylim(-ylim,ylim)

ax = fig.add_subplot(414)
#ax.plot(fano_td_time,np.abs(fano_td_theo)/np.max(np.abs(fano_td_theo[len(fano_td_theo)/2+20:])))
ax.set_title('inverse FT of spectrum (zoom 2)')
ax.plot(fano_td_time,fano_td_theo.real)
ax.plot(fano_td_time,fano_td_theo.imag)
ax.plot(fano_td_time,abs(fano_td_theo))
ax.legend(['real','imaginary','abs'])
ax.set_xlabel('delay [fs]')
ax.set_ylabel('intensity [a.u.')
ax.set_xlim(-1,1)
ax.set_ylim(-ylim,ylim)


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
