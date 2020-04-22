# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 17:43:00 2019

@author: andreas
"""

import numpy as np
import matplotlib.pyplot as plt

eps = np.linspace(-10,10,1000)
q = -0.17
phase_res = 0
phase_offres = -2
amp_offres = 0.5

res = q+eps/(eps+1j)*np.exp(1j*phase_res)
offres = amp_offres*np.exp(1j*phase_offres)
total = res + offres

fig1 = plt.figure('q={} phi_res{} phi_off{} amp_off{}'.format(q,phase_res,phase_offres,amp_offres),figsize=(5,10))
fig1.suptitle('q={} phi_res{} phi_off{} amp_off{}'.format(q,phase_res,phase_offres,amp_offres))

ax1 = fig1.add_subplot(211)
ax1.plot(res.real,res.imag,label='resonant')
ax1.plot(total.real,total.imag,label='total')
ax1.plot(0,0,'o',color='black')
ax1.legend()

ax2 = fig1.add_subplot(212)
ax2.plot(eps,np.angle(res),label='resonant')
ax2.plot(eps,np.angle(total),label='total')
ax2.legend()

wpi = np.conjugate(((q+1j)**2+1j*(eps-1j))/(eps-1j))
wpi = 1j*wpi
gruson = (eps+q)/(eps+1j)


fig2 = plt.figure()
fig2.suptitle('q={}'.format(q))
ax3 = fig2.add_subplot(211)
ax3.plot(eps,wpi.real, label = 'wpi.real')
#ax3.plot(eps,-wpi.real-1, label = '-wpi.real-1')
ax3.plot(eps,wpi.imag, label = 'wpi.imag')
ax3.plot(eps,gruson.real,label='gruson.real')
ax3.plot(eps,gruson.imag,label='gruson.imag')
ax3.legend()

ax4 = fig2.add_subplot(212)
ax4.plot(eps,np.angle(wpi), label = 'wpi')
ax4.plot(eps,np.angle(wpi.real-1+1j*wpi.imag), label = '-wpi-1')
ax4.plot(eps,np.angle(gruson),label='gruson')
ax4.legend()

#npormal lorenzian
gamma = 1

lorenz = (1j*eps-gamma)**-1

fig3 = plt.figure()
ax4 = fig3.add_subplot(111)
ax4.plot(eps,lorenz.real,label = 'lorenz.real')
ax4.plot(eps,lorenz.imag,label = 'lorenz.imag')
ax4.legend()

energy = np.linspace(0,20,100)
weg = 10
wkg = np.linspace(0,20,1000)

def Sw(w,wkg,weg):
    e = wkg-weg
    u = (e+q)/(e+1j)
    return np.abs(u)**2/(w-wkg-1j*0.1)

vector = np.array([])
for w in energy:
    test = np.array([Sw(w,a,weg) for a in wkg])
    vector = np.append(vector,np.sum(test))

plt.plot(energy,vector.real)
plt.plot(energy,vector.imag)
