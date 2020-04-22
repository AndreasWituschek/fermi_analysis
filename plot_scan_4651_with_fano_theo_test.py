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
data_window = [-600,-150.]
gauss = 0 # gauss window on TD (true) or not (False)
modfreq = 3.028


'''Phase correction'''
phas_corr = 5.4 #phase correction for everything as obtained by 5H measurements
print('Total phase correction: {} degree'.format(phas_corr))
phas_corr += (180.*harmonic) % 360. #accounts for phase shift due to the fact that ref signal is other output of interferometer.
phas_corr -= 180. #to take into account that we measure negative signals by ion detection


""" Parameters for theoretical curve """
phi = 180. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
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
    low = int(len(t_fano_theo)/2)
else:
    low = np.max(np.argwhere(t_fano_theo<data_window[0]))
high = np.min(np.argwhere(t_fano_theo>data_window[1]))
t_fano_theo = t_fano_theo[low:high] 
t_fano_theo = t_fano_theo[::-1]
I_fano_theo = I_fano_theo[low:high]
I_fano_theo = I_fano_theo[::-1]
I_fano_theo = fk.Phasing_TD(I_fano_theo,data_window,E_r) #time domain phasing of theoretical dipole response so it starts with Im(d(t))=-1 at beginning of data_window in case of a lorenzian line shape

   

#==============================================================================
# ''' Resolution from Data window '''
#==============================================================================
print('The resolution from the data window length is {} eV'.format(spc.h/(abs(data_window[-1]-data_window[0])*1E-15*spc.e)))

#==============================================================================
# ''' i0 correction for different harmonic contents of i0 '''
#==============================================================================
i0_fft = data['LDM']['i0_fft']
i0_0H = np.sum(i0_fft[::,395:405],axis=1)
i0_1H = np.sum(i0_fft[::,443:453],axis=1)
i0_2H = np.sum(i0_fft[::,492:502],axis=1)
i0_3H = np.sum(i0_fft[::,540:550],axis=1)
i0_4H = np.sum(i0_fft[::,589:599],axis=1)
i0_5H = np.sum(i0_fft[::,637:647],axis=1)
i0_6H = np.sum(i0_fft[::,685:695],axis=1)
i0_6H /= np.max(i0_6H)



#==============================================================================
# Analysis
#==============================================================================
mfli_harmonic = data[device]['harmonic'][int(demod)]
i = str(mfli_harmonic) + 'H'
Z = (data[device]['x' + demod] + 1j * data[device]['y' + demod])# / data['LDM']['I0']
#Z *= np.exp(-1j*np.pi/180.*200)
#plt.plot(delay,Z.real/max(abs(Z)))
##plt.plot(delay, Z.imag/max(abs(Z)))
#plt.plot(delay, Z_uhlia.real,color='r', alpha=0.5)
##plt.plot(delay, Z_uhlia.imag,color='g', alpha=0.5)
#
#plt.plot(delay,np.unwrap(np.angle(Z)))
#plt.plot(delay,np.unwrap(np.angle(Z_uhlia)))

#Z = Z_uhlia
if i0_correction:
    Z = Z*(np.nanmean(data['LDM']['I0'])/data['LDM']['I0'])
if i0_6H_correction:
    Z = Z*i0_6H


Z *= np.exp(1j*np.pi/180.*phas_corr)  #phasing of dataset according to phase difference between reference and WPI signal
R = np.abs(Z)
Phi = np.angle(Z, deg=False)
#        plt.plot(np.unwrap(Phi))


# Create theoretical curve
Ttheo = np.linspace(delay[0],delay[-1], 30000)
Xtd,Ytd,Xt,Yt = fk.Curve(spc.h*spc.c/(spc.e*l_trans)*1e9, l_ref, harmonic, phi, a*max(abs(Z)), offset, delay[0], delay[-1], 30000) # theo curve for transition
Xtd,Ytd,Xt,Yt = fk.Curve(l_fel/harmonic, l_ref, harmonic, phi, a*max(abs(Z)), offset, delay[0], delay[-1], 30000) # theo curve for laser CWL

T_d,X_d,Y_d = fk.CutDataSet(delay, Z.real, Z.imag, data_window)
Z_cut = X_d + 1j*Y_d
Z_plot = Z_cut
Z_cut = fk.Phasing_TD(Z_cut,T_d,E_r) #time domain phasing of data set to take into account position of data window

scaling = 1 #max(Z.real)/max(Ttheo)

# Fourier traffo

Td = np.mean(np.unique(np.diff(delay)))
if gauss:
    Zg = fk.GaussWindow(T_d, Z_cut, suscept)
    wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
else:
    wn, dft = fk.DFT(T_d, Z_cut, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
#        dft = dft[::-1]
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

#            plt.plot(T_d,R_fit)
#            plt.plot(T_d, decay(T_d,popt[0],popt[1]))
#            plt.plot(t_fano_theo, decay(t_fano_theo,popt_theo[0],popt_theo[1]))
   
   
#==============================================================================
# ''' Plotting '''
#==============================================================================
figTD = plt.figure('comparison Fano from data and theo',figsize=(14, 10))

ax = figTD.add_subplot(221)
ax.set_title('comparison Fano from data and theo')
ax.plot(delay, Z.real, '-', color=color)
ax.plot(delay, Z.imag, '-', color='g')
#ax.plot(Ttheo, scaling*Theo, 'k-')
if plotTheo: ax.plot(Ttheo, Xt, 'b', alpha=0.3)
ax.set_ylabel('a.u.')
ax.set_xlabel('delay [fs]')
ax.grid()
ax.set_xlim(data_window[0],data_window[-1])
ax.legend(['X data','Y data'])

ax = figTD.add_subplot(222)
#ax.plot(Ttheo, scaling*Theo, 'k-')
ax.plot(delay, Z.real/max(abs(Z.real)), '-', color=color)
ax.plot(t_fano_theo, I_fano_theo.real/max(abs(I_fano_theo.real)), '-', color='g')
#        ax.plot(T_d, Z_cut.real/max(abs(Z_cut.real)), '-', color='g')
#ax.plot(t_fano_theo, I_fano_theo.imag/max(abs(I_fano_theo.real)), '-', color='g')
##        ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.imag)),color='black',alpha=0.3)
if plotTheo: ax.plot(Ttheo, Yt, 'r', alpha=0.3)
ax.set_ylabel('a.u.')
ax.set_xlabel('delay [fs]')
ax.grid()
ax.set_xlim(data_window[0],data_window[-1])
ax.legend(['X theo','Y theo'])

ax = figTD.add_subplot(223)
ax.plot(delay, R/max(R), '-', color=color)
ax.plot(t_fano_theo, abs(I_fano_theo)/max(abs(I_fano_theo)), '-', color='g')
#        ax.fill_between(t_theo_seed,sl_AC*np.max(R),color='black',alpha=0.3)#,'-', color='black')
#        if dev=='dev3265' and d=='1':
#            ax.plot(delay_fit,decay(delay_fit,popt[0],popt[1])/max(R), '--',color='black', linewidth=2)
#            ax.plot(t_fano_theo,decay(t_fano_theo,popt_theo[0],popt_theo[1])/max(abs(I_fano_theo)), '--',color='black', linewidth=2)
ax.set_ylabel('R')
ax.set_xlabel('delay [fs]')
ax.grid()
ax.set_xlim(data_window[0],data_window[-1])
ax.legend(['R data','R theo','fits'])


'''plot frequency domain'''

if plot_eV:
    factor = spc.h*spc.c/spc.e*100.
    wn *= factor 
    wn_theo *= factor
    wn_lim *= factor
    l_ref /= factor
    l_trans /= factor
    FWHM *= factor
    print('The resolution of the DFT (after zeropadding) is {} eV'.format(FWHM))
undersampled_freq = wn-harmonic*1E7/l_ref        
#        dft = dft**2 #squaring to get intensity insted of amplitude
#        dft_theo = dft_theo**2

  
ax = figTD.add_subplot(224)
ax.plot(wn, dft.real/max(abs(dft)), '-', color='b', linewidth = 2)
ax.plot(wn, dft.imag/max(abs(dft)), '-', color='g', linewidth = 2)
#        ax.plot(wn, dft.real/max(abs(dft)), '-', color='magenta', linewidth = 2)
#ax.plot(wn, abs(dft)/max(abs(dft)), '-', color='black', linewidth = 2)
ax2 = ax.twinx()
#ax2.plot(wn[fk.find_index(wn,28.45):fk.find_index(wn,28.57)],spec_phase[fk.find_index(wn,28.45):fk.find_index(wn,28.57)],color = 'k',linestyle='--',alpha=0.7)
ax2.set_ylabel(r'phase [rad]')
ax2.plot(wn,np.angle(dft), '--', color='k', linewidth = 2)
ax2.plot(wn_theo,np.angle(dft_theo), '--', color='magenta', linewidth = 2)
#        ax.plot(wn_theo, abs(dft_theo)/max(abs(dft_theo)), 'o-', color='g', linewidth = 2)
ax.plot(wn_theo, dft_theo.real/max(abs(dft_theo)), '-', color='r', linewidth = 2)
#        ax.plot(wn_theo, dft_theo.real/max(abs(dft_theo)), '-', color='b', linewidth = 2)
#        ax.plot(t_fano,(I_fano+0.03)/np.max(abs(I_fano+0.03)), color='black',linewidth = 2)
ax.plot(E_fano+E_undersamp,fano_resonance,color='grey',linewidth = 2)
#add second x axis to plot:
#        ax2 = ax.twiny()
#        ax2.set_xlabel('rotating frame energy [eV]')         
#        ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)
#        axs.plot(padres_roi,padres_spec,color='grey',linewidth = 2) #FEL spectrum
if plot_eV:
    ax.set_xlabel(r'energy [eV]')
else:
    ax.set_xlabel(r'wavenumber [cm$^{-1}$]')
ax.set_ylim(-1.1,1.1)
ax.set_ylabel('intensity [a.u.]')
ax.set_xlim(wn_lim[0],wn_lim[1])
ax.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
plot_x_range = np.diff(ax.get_xlim())[0]
plot_y_range = np.diff(ax.get_ylim())[0]
ax.axvspan(l_trans-FWHM/2, l_trans + FWHM/2, ymin=0.1, ymax=0.2,alpha=0.3)
ax.legend(['data.real','data.imag','theo.real'])
ax2.legend(['data.phase', 'theo.phase'],loc='lower right')
plt.tight_layout()

if not interactive_plots:
    plt.close('all')
else:
    plt.show()



if spectrogram:
    wn_lim_s = [226000,233000]
    l_ref *= factor
    delay = -delay[::-1]
    Z = np.conjugate(Z)[::-1]
    data_window = np.negative(data_window)
#    l_trans *= factor
    # plot spec
    figFT = plt.figure('scan_{0:03}_{1}_d{2}_SG_{3}'.format(run, 2, 1, i))
    slide_step = Td  #in fs
    FWHM_slide = 20.0  #in fs
    t_lim = data_window
    slide_positions = np.arange(delay[0],delay[-1], slide_step)
    S = np.zeros((np.size(delay), np.size(wn)))
    i = 0
    for slide_pos in slide_positions:
        Z_slide = fk.slide_window(delay, Z, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
        wn,  DFT_slide = fk.DFT(delay, Z_slide, Td, l_ref, harmonic, zeroPaddingFactor = zeroPaddingFactor)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian           
        # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
#                for i in range(np.size(T_d)):
#                    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))*fk.weighting_coeff(T_d[i], slide_pos, FWHM_slide)
        S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
        i += 1
    S[:,:] /= np.size(slide_positions)
    fk.plot_spectrogram(figFT, wn, delay, S, l_fel/6., 1E7/l_trans, wn_lim_s, t_lim)
#            plt.show()
    


#print delay

#figTD = plt.figure("zoom",figsize=(7,6))
#
#ax = figTD.add_subplot(111)
##ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
#ax.plot(delay, Z.real, 'o-', color='b')
#ax.plot(delay, Z.imag, 'o-', color='g')
#ax.plot(Ttheo, Xt*0.34, 'b', alpha=0.3)
#ax.legend(['real','imag','real_theo'])
#ax.set_xlim(-20,20)
#ax.grid()



#fig = plt.figure()
#
#ax = fig.add_subplot(111)
##ax.plot(abs(np.fft.fftshift(np.fft.fft(Z))))
#dft2 = abs(np.fft.fftshift(np.fft.fft(Z,n=256)))
#freq2 = 1e15*spc.h/spc.e*np.fft.fftshift(np.fft.fftfreq(len(dft2),delay[1]-delay[0]))
#ax.plot(freq2,dft2)


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




'''plot frequency domain'''
undersampled_freq = wn-harmonic*1E7/l_ref
spec_phase = np.angle(dft)
#spec_phase =  [x + 2.*np.pi if x <= 0 else x for x in spec_phase]
  
figFD = plt.figure(figsize=(9,4))

ax = figFD.add_subplot(111)
ax.plot(wn, dft.real/max(abs(dft)), '-', color='b', alpha=0.7)
ax.plot(wn, dft.imag/max(abs(dft)), '-', color='r', alpha=0.7)
ax.plot(E_fano+E_undersamp,fano_resonance,color='black', alpha=0.9)

ax2 = ax.twinx()
ax2.plot(wn,spec_phase,color = 'k',linestyle='--',alpha=0.7)
ax2.set_ylabel(r'phase [rad]')
ax2 = ax.twiny() #add second x axis to plot:
ax2.set_xlabel('rotating frame energy [eV]')         
ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)

ax.set_xlabel(r'energy [eV]')
ax.set_ylim(-1.3,1)
ax.set_ylabel('intensity [a.u.]')
ax.set_xlim(wn_lim[0],wn_lim[1])
plot_x_range = np.diff(ax.get_xlim())[0]
plot_y_range = np.diff(ax.get_ylim())[0]
ax.axvspan(l_trans-FWHM/2, l_trans + FWHM/2, ymin=0.1, ymax=0.2,alpha=0.3)
ax.plot(0,0,'o',color='k',alpha=0.7)
ax.legend([r'$Re(FT(S))$',r'$Im(FT(S))$',r'lit.',r'$\phi(FT(S))$'],ncol=2,loc='lower right')

plt.tight_layout()

