# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 23:15:04 2019

Calculates excitation probability based on Einstein coefficient,
first order perturbation theory and Bloch-equations.
For the latter, a rectangular pulse is assumed.

@author: LukasB
"""
import scipy as sp
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt

# input parameters ------------------------------------------------------------
wl_ex = 52.22E-9   # excitation wavelength of transition (m)
A = 0.24E9      # Einstein A coefficient (1/s)
E_L= 25.0E-6      # pulse energy (J)
dt_L = 43.0E-15   # pulse duration of intensity envelope (s)
d_L = 100.0E-6    # laser focus diameter (m)
att = 1.2E-3    # intensity attenuation of laser

I0 = att*E_L/(0.25*np.pi*d_L**2*dt_L)    # avg laser intensity (J/m^2*s)

# Bloch equations:

Omega = np.sqrt(3*I0*wl_ex**3*A/(2*np.pi*const.c*const.h))    # Rabi frequency
pp = np.sin(Omega/np.sqrt(2)*dt_L)**2                       # population probability
pc = np.sqrt(pp)                                            # probability excitation of coherence
print pp

# plot bloch
# define max intensity for plot
E_laser_max = 30.0E-6 # pulse energy (J)
d_L_min = 70.0E-6  # laser focus diameter (m)

Imax = att*E_laser_max/(0.25*np.pi*d_L_min**2*dt_L)
Imax *=1E-4
I = np.linspace(0, Imax, 200)# laser intensity in W/cm^2
Omega = [np.sqrt(3*i*1E4*wl_ex**3*A/(2*np.pi*const.c*const.h)) for i in I]
pp = [np.sin(o/np.sqrt(2)*dt_L)**2 for o in Omega]
pc = np.sqrt(pp)
plt.semilogx(I, pp)
plt.xlabel('Laser intensity (W/cm^2)')
plt.ylabel('Population probability')
plt.grid()