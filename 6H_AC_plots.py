# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:38:40 2019

@author: andreas
"""

import numpy as np
import matplotlib.pyplot as plt


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
plt.rcParams['lines.linewidth'] = 2.


def pulse(t,cen,fwhm):
    return np.exp(-4*np.log(2)*(t-cen)**2/fwhm**2)

tau = 1.
fwhm = 0.6
t = np.linspace(tau/2-1,tau/2+1,1000)
harmonic = 6.

 
fig = plt.figure('6H AC plots',figsize=(12,6))
ax = fig.add_subplot(111)
amps = np.array([1.,1.1,1.3,1.5,1.3,1.1,1.])
for i in np.arange(harmonic+1):
#    pulse1 = pulse(t,-tau,fwhm)+pulse(t,tau,fwhm)
    pulse_i = pulse(t,tau*(harmonic-i)/harmonic,fwhm)*amps[i]
    ax.plot(t,pulse_i,color=(i/harmonic, 0.3, (harmonic-i)/harmonic))
    ax.set_xlabel('delay')
    ax.set_ylim([-0.1,1.7])
    plt.gca().axes.get_yaxis().set_visible(False)
    ax.grid(b=True,axis='x')
    ax.set_xticks(np.arange(harmonic+1)/harmonic)