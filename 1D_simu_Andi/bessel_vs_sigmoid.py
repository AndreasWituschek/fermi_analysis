# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 08:58:45 2020

@author: andreas
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz, savgol_filter
from scipy.special import jv    
import time


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
plt.rcParams['lines.linewidth'] = 1.5

colors = ['#1f77b4','#d62728','#2ca02c','#ff7f0e','#9467bd','#8c564b']
          
#E field after harmonic generation with and without saturation.    
def sigmoid(pulse,harmonic):
    xmax = (1+np.sqrt(2/3)*harmonic**(-2/3))
    jmax = jv(harmonic,harmonic*xmax)
    f= np.sqrt(2)/xmax
    print(f)
    return np.real((f*pulse)**harmonic/(1+np.abs((f*pulse)**harmonic)))*jmax, jmax 
#E field after harmonic generation with and without saturation.
def bessel(pulse,harmonic):
    seedenvelope = np.abs(pulse)
    seedphase = np.unwrap(np.angle(pulse))
    return jv(harmonic,harmonic* seedenvelope)*np.exp(1j*seedphase*harmonic)

# =============================================================================
# parameters
# =============================================================================

#parameters for 6H AC
h = 6
ampb = 0.83 #amp for bessel function.
amps = 0.979 #amp for sigmoid function. amps =1.11 in old versions where sigmoind function was not adjusted to bessel function

##parameters for 5H WPI
#h = 5
#ampb = 1.05 #amp for bessel function.
#amps = 1.147 #amp for sigmoid function. amps =1.3 in old versions where sigmoind function was not adjusted to bessel function

# =============================================================================
# plotting
# =============================================================================

x = np.linspace(-0.5,1.8,1000)
sig, jmax = sigmoid(x,h)
bes= bessel(x,h)

plt.plot(x,bes,color=colors[0],label = 'Bessel')
plt.plot(x,sig,color=colors[1],label='Sigmoid')
plt.plot(x,x**6*jmax,color=colors[2],label = r'E$_{tot}^6$')
#plt.axvline(amps,color=colors[0],linestyle='dashed')
#plt.axvline(ampb,color=colors[1],linestyle='dashed')
plt.legend()
plt.xlim([0,1.6])
plt.ylim([-0.05,0.4])
plt.xlabel(r'E$_{tot}$ [arb. units]')
plt.ylabel(r'E$_{FEL}$ [arb. units]')
plt.tight_layout()