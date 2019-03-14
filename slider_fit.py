# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:23:48 2019

@author: FemtoMeasure
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import os
import h5py


def gaus(x,amp,x0,sigma):
    return amp*np.exp(-(x-x0)**2/(2*sigma**2))

def padres_fit_function(x,amp,x0,sigma,freq, phase, c):
    return gaus(x,amp,x0,sigma) * (1-np.sin(2*np.pi+(freq * x + phase))) + c


i = 150
run = 218
delay_zero_pos = 11025.66

ldm_file_path = '//online4ldm.esce.elettra.trieste.it/store/20149020/Day_2/Run_' + str(run) +'/rawdata/'.format(int(run))

ldm_file = os.listdir(ldm_file_path)[0]

ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')

padres_wavelength = np.array(ldm_data['photon_diagnostics']['Spectrometer']['Wavelength'])
print('padres_wavelength: {}'.format(padres_wavelength))

padres_span = np.array(ldm_data['photon_diagnostics']['Spectrometer']['WavelengthSpan'])
print('padres_span: {}'.format(padres_span))

ldm_padres = np.array(ldm_data['photon_diagnostics']['Spectrometer']['hor_spectrum'])
print('ldm_padres shape: {}'.format(ldm_padres.shape))


ldm_data.close()

ldm_w_padres = ldm_padres[i][400:600]


fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.45)
t = np.linspace(0, 200, 500, endpoint=True)
a0 = 195588
f0 = 0.2
s0= 12.11
p0 = 1.14
m0 = 136.54
c0 = 8000
s = a0*np.sin(2*np.pi*f0*t)
s = padres_fit_function(t,a0,m0,s0,f0, p0, c0)
l, = plt.plot(t, s, lw=2, color='red')
plt.plot(ldm_w_padres)
#plt.axis([0, 200, 8000, 10000])


axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
axamp = plt.axes([0.25, 0.15, 0.65, 0.03])
axphase = plt.axes([0.25, 0.2, 0.65, 0.03])
axsigma = plt.axes([0.25, 0.25, 0.65, 0.03])
axm = plt.axes([0.25, 0.3, 0.65, 0.03])
axc = plt.axes([0.25, 0.35, 0.65, 0.03])

sfreq = Slider(axfreq, 'Freq', 0.01, 1.0, valinit=f0)
samp = Slider(axamp, 'Amp', 1.0e3, 3.5e5, valinit=a0)
sphase = Slider(axphase, 'Phase', 0.0, 2*np.pi, valinit=p0)
ssigma = Slider(axsigma, 'Sigma', 0.1, 100.0, valinit=s0)
sm = Slider(axm, 'm', 0.0, 200.0, valinit=m0)
sc = Slider(axc, 'c', 0.0, 1.0e4, valinit=c0)


def update(val):
    amp = samp.val
    freq = sfreq.val
    m = sm.val
    sigma = ssigma.val
    phase = sphase.val
    c = sc.val
    l.set_ydata(padres_fit_function(t,amp,m,sigma,freq, phase, c))
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)
sphase.on_changed(update)
ssigma.on_changed(update)
sm.on_changed(update)
sc.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    sfreq.reset()
    samp.reset()
    sphase.reset()
    ssigma.reset()
    sm.reset()
    sc.reset()
button.on_clicked(reset)

plt.show()