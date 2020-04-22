# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:46:56 2019

@author: FemtoMeasure
"""

import numpy as np
import fermi_analysis.functions as fk 
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as spc
import pandas as pd


#plt.close('all')
""" Data parameters """
run = 4651 #first run of delay scan
demod = '1'
device = 'dev3265'
#demod = ['0','1', '2']
#device = ['dev3265', 'dev3269']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
interactive_plots = 1
save_fig = 0
plotTheo = False
spectrogram = 0
plot_eV = 1

"""Experimental Parameters"""
l_trans =  28.510223  # [eV] fanoresonance in argon
harmonic = 6.  # harmonic of FEL
title = 'Scan_{0:03}/'.format(int(run))
color = 'b'
data_window = [-600,100]
gauss = 0 # gauss window on TD (true) or not (False)
modfreq = 3.028



""" Parameters for theoretical curve """
phi = 0. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = .5  # amplitude (=R from MFLI)
offset = 0.  # offset

""" analysis parameters """
i0_correction = False
i0_6H_correction = False
unwrap_phase = True
slide_step = 32
FWHM_slide = slide_step / np.sqrt(48)
zeroPaddingFactor = 1
suscept = True
wn_lim = np.asarray([228400.,231580.])

# Load preanalysed data
file_path = root_file_path + 'scan_{0:03}/'.format(int(run))
h5f = h5py.File(file_path + 'scan_{0:03}.h5'.format(int(run)), 'r')

# sort data by delay in ascending order
sort_inds = np.negative(np.array(h5f.get('LDM/delay'))).argsort()
#sort_inds = np.array(h5f.get('LDM/delay')).argsort()

data = {
    'run_numbers': np.array(h5f.get('run_numbers'))[sort_inds],
    'LDM': {
        'I0': np.array(h5f.get('LDM/I0'))[sort_inds],
        's_I0': np.array(h5f.get('LDM/s_I0'))[sort_inds],
        'delay': np.array(h5f.get('LDM/delay'))[sort_inds],
        's_delay': np.array(h5f.get('LDM/s_delay'))[sort_inds],
        'l_seed': np.array(h5f.get('LDM/l_seed'))[sort_inds],
        'l_ref': np.array(h5f.get('LDM/l_ref'))[sort_inds],
        'i0_fft': np.array(h5f.get('LDM/i0_fft'))[sort_inds],
            },
    'dev3265': {
        'harmonic': np.array(h5f.get('dev3265/harmonic')),
        'x0': np.array(h5f.get('dev3265/x0'))[sort_inds],
        'y0': np.array(h5f.get('dev3265/y0'))[sort_inds],
        'x1': np.array(h5f.get('dev3265/x1'))[sort_inds],
        'y1': np.array(h5f.get('dev3265/y1'))[sort_inds],
        'x2': np.array(h5f.get('dev3265/x2'))[sort_inds],
        'y2': np.array(h5f.get('dev3265/y2'))[sort_inds],
        's_x0': np.array(h5f.get('dev3265/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('dev3265/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('dev3265/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('dev3265/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('dev3265/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('dev3265/s_y2'))[sort_inds], },
    'dev3269': {
        'harmonic': np.array(h5f.get('dev3269/harmonic')),
        'x0': np.array(h5f.get('dev3269/x0'))[sort_inds],
        'y0': np.array(h5f.get('dev3269/y0'))[sort_inds],
        'x1': np.array(h5f.get('dev3269/x1'))[sort_inds],
        'y1': np.array(h5f.get('dev3269/y1'))[sort_inds],
        'x2': np.array(h5f.get('dev3269/x2'))[sort_inds],
        'y2': np.array(h5f.get('dev3269/y2'))[sort_inds],
        's_x0': np.array(h5f.get('dev3269/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('dev3269/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('dev3269/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('dev3269/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('dev3269/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('dev3269/s_y2'))[sort_inds], }, }
h5f.close()
l_fel = data['LDM']['l_seed'][0]
ldm_l_ref = data['LDM']['l_ref']
l_ref = np.mean(ldm_l_ref)
print l_ref
l_ref = 266.003
delay = (data['LDM']['delay']-0.59)*1.0001848 #calibrated delay. Factor arises from the difference in refractive index between the seed CWL used for He and for Ar.
Ttheo = np.linspace(data_window[0],data_window[1],1000)


#==============================================================================
# import UHLIA data
#==============================================================================

#h5f_uhlia = h5py.File('//nanoserver/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_4651/uhlia_scan_4651_phase155.h5', 'r')
#Z_uhlia = np.array(h5f_uhlia['6H']['1s_4o'])
#Z_uhlia /= max(abs(Z_uhlia))
#
#plt.plot(Z_uhlia.imag)

#==============================================================================
# '''' Fano Resonance profile from Paper '''
#==============================================================================
#fano parameters
q = -0.17
#q=-1000
E_r = 28.510223 #28.509 #eV 
gamma = 0.0126 #eV

E_undersamp = harmonic*spc.h*spc.c/(l_ref*spc.e)*1e9 # undersampling energy
E_r = E_r - E_undersamp  #eV. undersampled fano resonance
E_fano = np.linspace(0,E_r+30,500000)

fano_resonance = fk.fano(E_fano,E_r,q,gamma,1.,-1.)
I_fano_theo = np.conjugate(np.fft.fftshift(np.fft.fft(fano_resonance)))
t_fano_theo = 1e15*spc.h/spc.e* np.fft.fftshift(np.fft.fftfreq(len(fano_resonance),E_fano[1]-E_fano[0]))
#adjusting data window
if data_window[0]==0.:
    low = len(t_fano_theo)/2
else:
    low = max(np.argwhere(t_fano_theo<data_window[0]))
high = min(np.argwhere(t_fano_theo>data_window[1]))
t_fano_theo = t_fano_theo[low:high] 
t_fano_theo = t_fano_theo[::-1]
I_fano_theo = I_fano_theo[low:high]
I_fano_theo = I_fano_theo[::-1]
I_fano_theo = fk.Phasing_TD(I_fano_theo,data_window,E_r) #time domain phasing of theoretical dipole response so it starts with Im(d(t))=-1 at beginning of data_window in case of a lorenzian line shape

   


#==============================================================================
# Analysis
#==============================================================================
mfli_harmonic = data[device]['harmonic'][int(demod)]
i = str(mfli_harmonic) + 'H'
Z = (data[device]['x' + demod] + 1j * data[device]['y' + demod])# / data['LDM']['I0']


if i0_correction:
    Z = Z*(np.nanmean(data['LDM']['I0'])/data['LDM']['I0'])

R = np.abs(Z)
Phi = np.angle(Z, deg=False)
#        plt.plot(np.unwrap(Phi))


# Create theoretical curve
Ttheo = np.linspace(delay[0],delay[-1], 30000)
Xtd,Ytd,Xt,Yt = fk.Curve(spc.h*spc.c/(spc.e*l_trans)*1e9, l_ref, harmonic, phi, A, offset, delay[0], delay[-1], 30000) # theo curve for transition
#Xtd,Ytd,Xt,Yt = fk.Curve(l_fel/harmonic, l_ref, harmonic, phi, a*max(abs(Z)), offset, delay[0], delay[-1], 30000) # theo curve for laser CWL

T_d,X_d,Y_d = fk.CutDataSet(delay, Z.real, Z.imag, data_window)
Z_cut = X_d + 1j*Y_d
Z_plot = Z_cut
Z_cut = fk.Phasing_TD(Z_cut,T_d,E_r) #time domain phasing of data set to take into account position of data window


# Fourier traffo

Td = np.mean(np.unique(np.diff(delay)))
if gauss:
    Zg = fk.GaussWindow(T_d, Z_cut, suscept)
    wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
else:
    wn, dft = fk.DFT(T_d, Z_cut, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
FWHM , peakFullWidth = fk.PeakWidth(T_d, 0.05, False, zeroPaddingFactor)

       
#Fouriert Trafo of the theoretical dataset
wn_theo, dft_theo = fk.DFT(t_fano_theo, I_fano_theo, t_fano_theo[1]-t_fano_theo[0], l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)

''' Fitting the dephasing for 6th harmonic R'''

def decay(t,amp,tau):
    return amp*np.exp(-t/tau) + 0.001828
def decay_theo(t,amp,tau):
    return amp*np.exp(-t/tau)
    
t_start =-200 #time from where the fitting of the decay starts
idx_start = fk.find_index(T_d,t_start)

R_fit = abs(Z_cut) #region of R that will be fitted
popt, pcov = curve_fit(decay, T_d[idx_start:], R_fit[idx_start:], p0=[0.8,-110.],absolute_sigma=False)
perr = np.sqrt(np.diag(pcov))

#fitting theoretical curve:
popt_theo, pcov_theo = curve_fit(decay_theo, t_fano_theo, abs(I_fano_theo), p0=[0.8,-110.],absolute_sigma=False)
perr_theo = np.sqrt(np.diag(pcov_theo))

print('The time constant of the data fit is {} +- {} fs'.format(popt[1],perr[1]))
print('The time constant of the data fit is {} +- {} fs'.format(popt_theo[1],perr_theo[1]))

  

#importing 6H Signal from AC scan with 6th harmonic and Helium (scan5234)
df = pd.read_csv('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/X_6H_vs_delay.csv')
T5234 = np.array(df['delay'])
X_6H5234 = np.array(df['X_6H'])
X_6H5234 /= np.max(abs(X_6H5234))


#==============================================================================
# ''' Plotting '''
#==============================================================================
figTD = plt.figure('comparison Fano from data and theo',figsize=(14, 10))
# plot X
ax = figTD.add_subplot(111)
ax.set_title('comparison Fano from data and theo')
ax.plot(delay, Z.real/max(abs(Z.real)), '-', color='b',alpha=0.7)
ax.plot(Ttheo, Xt, 'r', alpha=0.3)
ax.plot(T5234, X_6H5234, color='g',alpha=0.7)
ax.set_ylabel('a.u.')
ax.set_xlabel('delay [fs]')
ax.legend(['X data','X theo','AC X_6H'])


#%%


#==============================================================================
# ''' plots for paper '''
#==============================================================================

#reversing the time domain arrays so everything goes in positive direction
Z_plot =  Z_plot[::-1]
Z_plot = np.conjugate(Z_plot)
t_fano_theo = np.negative(t_fano_theo)[::-1]
T_d = np.negative(T_d)[::-1]
I_fano_theo= I_fano_theo[::-1]
phase = np.unwrap(np.angle(Z_plot))

#fitting phase to get value at 0 delay:
def Lin(x,a,b):
    return a*x+b

start =  fk.find_index(T_d,368)
stop = fk.find_index(T_d,150)
popt, pcov = curve_fit(Lin, T_d[start:stop], phase[start:stop], absolute_sigma=False)
perr = np.sqrt(np.diag(pcov))

a= popt[0]
b = popt[1]
sb = perr[1]

print('Phase at zero delay is {} +- {} rad'.format(b % 2.*np.pi, sb/2.*np.pi))

#spectral phase at peak value:
phase0 = np.angle(dft)[fk.find_index(wn,28.51)]

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


figTD = plt.figure(figsize=(7,6))

ax = figTD.add_subplot(211)
ax.plot(T_d, Z_plot.real/max(abs(Z_plot))*1.1, '-', color='b',alpha=0.7)
ax.set_ylabel('intensity [a.u.]')
#ax.set_xlabel(r'$\tau$ [fs]')
ax.legend([r'$Re(S)$'])
ax.set_ylim(-1.3,1.3)
ax.set_xticklabels([])
ax.set_xlim(-data_window[1],-data_window[0])


ax = figTD.add_subplot(212)
ax.plot(T_d, abs(Z_plot)/max(abs(Z_plot))*1.1, '-', color='g',alpha=0.7)
ax.plot(0,0,color='r')
#if dev=='dev3265' and d=='1':
#    ax.plot(T_d,decay(T_d,popt[0],-popt[1])/max(abs(decay(T_d,popt[0],-popt[1]))), '--',color='black', linewidth=2)
#    ax.plot(t_fano_theo,decay_theo(t_fano_theo,popt_theo[0],-popt_theo[1])/max(abs(decay_theo(t_fano_theo,popt_theo[0],-popt_theo[1]))), '-',color='black', linewidth=2)
ax.set_ylabel('intensity [a.u.]')
ax.set_xlabel(r'$\tau$ [fs]')
ax.legend([r'$A(\tau)$',r'$\phi(\tau)$'],loc='center right')
ax.set_ylim(-.3,1.3)
ax.set_yticks([-0.2,0,0.2,0.4,0.6,0.8,1,1.2])
ax.set_yticklabels([-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2])
ax2 = ax.twinx()
ax2.plot(T_d, (phase-b+phase0)/(2*np.pi),color='r', alpha=0.7)
#ax2.plot(T_d, (a/(2*np.pi)*T_d+b/(2*np.pi)))
ax2.set_ylabel(r'phase [$2\pi$ rad]')
ax2.set_ylim(15,65)
ax.set_xlim(-data_window[1],-data_window[0])
plt.tight_layout()
plt.subplots_adjust(hspace=0.0)





