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
import pandas as pd


#plt.close('all')
""" Data parameters """
run = 893 #first run of delay scan
demod = ['1']
device = ['dev3265']
#demod = ['0','1', '2']
#device = ['dev3265', 'dev3269']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime1/combined/'
interactive_plots = 1
save_fig = 0
plotTheo = 1
spectrogram = 1
plot_eV = 1
plot_l_ref = 0

"""Experimental Parameters"""
# parameters for theoretical curve
l_trans = 52.2213067  # helium 1s-4p transition in nm 
#l_trans = 260.92/6.0 # fanoresonance in argon
#l_ref = 266.0234244226323  #265.98  # reference laser wavelenght in nm
fwhm_FEL = 0.075 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Scan_{0:03}/'.format(int(run))
color = 'b'
data_window = [-400.,200.0] 
gauss = 0 # gauss window on TD (true) or not (False)

""" Parameters for theoretical curve """
phi = 50. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 0.9   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
i0_correction = False
unwrap_phase = True
slide_step = 32
FWHM_slide = slide_step / np.sqrt(48)
zeroPaddingFactor = 2
suscept = False
#wn_lim = [217100,233800]
wn_lim = np.asarray([185000.,196000.])
#wn_lim = [189800,193100]

# Load preanalysed data
file_path = root_file_path + 'scan_{0:03}/'.format(int(run))
h5f = h5py.File(file_path + 'scan_{0:03}.h5'.format(int(run)), 'r')

# sort data by delay in ascending order
sort_inds = np.array(h5f.get('LDM/delay')).argsort()

data = {
    'run_numbers': np.array(h5f.get('run_numbers'))[sort_inds],
    'LDM': {
        'I0': np.array(h5f.get('LDM/I0'))[sort_inds],
        's_I0': np.array(h5f.get('LDM/s_I0'))[sort_inds],
        'delay': np.array(h5f.get('LDM/delay'))[sort_inds],
        's_delay': np.array(h5f.get('LDM/s_delay'))[sort_inds],
        'l_seed': np.array(h5f.get('LDM/l_seed'))[sort_inds],
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

l_ref = 266.003
#l_fel = 43.46*harmonic
l_fel = 261.0726087 #261.7 #261.0726087

T = data['LDM']['delay']
#data_window = [T[0]-5, T[-1]+5]



''' laser spectrum '''
points = np.linspace(wn_lim[0],wn_lim[1],1000)
fwhm_FEL_wn = 1E7*(fwhm_FEL/(l_fel/harmonic)**2) # wavenumber FWHM of FEL
l_FEL = 1E7/l_fel # wavenumber CWL of FE FalseL

''' Fano resonance simulation '''
Ttheo = np.linspace(T[0],T[-1],1000)
AC_FWHM = 100.0  # in fs
Fano_FWHM = 140.0  # in fs
AC = np.exp(-4.0*np.log(2)*((Ttheo)/AC_FWHM)**2)
Fano = np.exp(-4.0*np.log(2)*((Ttheo)/Fano_FWHM)**2)
Theo = 3*(3*Fano-7*AC)

''' Seed laser AC simulation '''
sl_fwhm = 99. # seed intensity fwhm in fs
t_theo_seed = np.linspace(-200,200,5000)
sl_AC = np.exp(-4.0*np.log(2)*(t_theo_seed/(sl_fwhm*np.sqrt(2)))**2) # seed laser AC
sl_AC /= np.max(sl_AC)


# loop over the mfli devices
for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = (data[dev]['x' + d] + 1j * data[dev]['y' + d])# / data['LDM']['I0']
#        Z *= np.exp(1j*40./360.*2.*np.pi)
        if i0_correction:
            Z = Z*(np.nanmean(data['LDM']['I0'])/data['LDM']['I0'])
        Z_s = data[dev]['s_x' + d] + 1j * data[dev]['s_y' + d]

        R = np.sqrt(Z.real**2 + Z.imag**2)
        R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)
        Phi = np.angle(Z, deg=False)


# Create theoretical curve
#        if draw_theory:
        #Ttheo = np.linspace(T[0],T[-1], 1000)
        Xtd,Ytd,Xt,Yt = fk.Curve(l_trans, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], 1000)


        delay = T
        X = Z.real
        X_s = Z_s.real
        Y = Z.imag
        Y_s = Z_s.imag

        T_d,X_d,Y_d = fk.CutDataSet(T, X, Y, data_window)
        Z = X_d + 1j*Y_d
        scaling = 1E-4 #max(Z.real)/max(Ttheo)

        # Fourier traffo
        Zg = fk.GaussWindow(T_d, Z, suscept)
        Td = T[1]-T[0]
        if gauss:
            wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        else:
            wn, dft = fk.DFT(T_d, Z, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        FWHM , peakFullWidth = fk.PeakWidth(T, 0.1, False, zeroPaddingFactor)
        FWHM = FWHM

        # theoretical spectrum
        center = 1E7*(1/l_fel- 1/l_ref)*mfli_harmonic + harmonic*1E7/l_ref
        spectrum = np.exp(-4.0*0.693147*(points-center)**2/fwhm_FEL_wn**2) # spectrum in wavenumber space
        spectrum*=(max(abs(dft))-min(abs(dft)))+min(abs(dft))
        
        ''' Plot time domain '''
        figTD = plt.figure('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i), figsize=(14, 10))

        ax = figTD.add_subplot(321)
        ax.set_title('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i))
#        ax.errorbar(delay, X, yerr=X_s, color=color, linestyle='')
        ax.plot(delay, X, '-', color=color)
        #ax.plot(Ttheo, scaling*Theo, 'k-')
        if plotTheo: ax.plot(Ttheo, Xt, 'r', alpha=0.3)
        ax.set_ylabel('X in a.u.')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(322)
        #ax.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
        #ax.plot(Ttheo, scaling*Theo, 'k-')
        ax.plot(delay, Y, '-', color=color)
        if plotTheo: ax.plot(Ttheo, Yt, 'r', alpha=0.3)
        ax.set_ylabel('Y in a.u.')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(323)
        #ax.errorbar(delay, R, yerr=R_s, color=color, linestyle='')
        ax.plot(delay, R, '-', color=color)
        ax.set_ylabel('R')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(324)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        if unwrap_phase:
            ax.plot(delay, np.unwrap(Phi), '-', color=color)
        else: 
            ax.plot(delay, Phi, '-', color=color)
        ax.set_ylabel('Phi (rad)')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        
        ##Plot I0
        ax = figTD.add_subplot(325)
#        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(delay, data['LDM']['I0'], '-', color= 'g')
        ax.set_ylabel('I0 in $\mu$J')
        ax.set_xlabel('Delay in fs')
#        ax.set_ylim(0,3)
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])
        
        ###Plot Delay
#        ax = figTD.add_subplot(515)
##        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
#        ax.plot(delay, '-', color=color)
#        ax.set_ylabel('Delay in fs')
#        ax.set_xlabel('Run')
##        ax.set_ylim(0,3)
#        ax.grid()
##        ax.set_xlim(data_window[0],data_window[-1])
#
        ##Plot Reflaser wavelenght
        if plot_l_ref:
            ax = figTD.add_subplot(326)
    #        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle=''
            ax.set_ylabel('wavelenght in nm')
            ax.set_xlabel('delay')
    #        ax.set_ylim(0,3)
            ax.grid()
            ax.set_xlim(data_window[0],data_window[-1])

        plt.tight_layout()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_TD_{3}.png'.format(run, dev, d, i), dpi=400)


        '''plot frequency domain'''
        if plot_eV:
            wn *= 1.239842E-4 
            wn_lim *= 1.239842E-4
            l_ref /= 1.239842E-4
            l_trans /= 1.239842E-4
                
        figFT = plt.figure('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs = figFT.add_subplot(111)
        axs.plot(wn, abs(dft), '-', color='r', linewidth = 2)
        axs.plot(points,spectrum,'-',color='grey',linewidth = 2)
        if plot_eV:
            axs.set_xlabel(r'energy [eV]', fontsize=14)
        else:
            axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
        axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
        axs.set_title('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
#        axs.set_xlim(180000.,195000.)
        axs.set_xlim(wn_lim[0],wn_lim[1])
        axs.axvline(1E7/l_trans,color = 'k' ,linestyle='-', linewidth= 2)
        axs.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(axs.get_xlim())[0]
        plot_y_range = np.diff(axs.get_ylim())[0]
#        axs.axhline(y= plot_y_range * 0.9 + axs.get_ylim()[0], xmin=0.1, xmax=0.1 + FWHM / plot_x_range, linewidth=3)
#        axs.text(0.1, plot_y_range * 0.9, str(FWHM))
        axs.grid()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_FD_{3}.png'.format(run, dev, d, i), dpi=400)

        if spectrogram:
            wn_lim /= 1.239842E-4
            l_ref *= 1.239842E-4
            l_trans *= 1.239842E-4
            # plot spec
            figFT = plt.figure('scan_{0:03}_{1}_d{2}_SG_{3}'.format(run, dev, d, i))
            print(T_d[0], T_d[-1])
            slide_step = Td  #in fs
            FWHM_slide = 20.0  #in fs
            t_lim = data_window
            slide_positions = np.arange(T_d[0],T_d[-1], slide_step)
            S = np.zeros((np.size(T_d), np.size(wn)))
            i = 0
            for slide_pos in slide_positions:
                Z_slide = fk.slide_window(T_d, Z, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
                wn,  DFT_slide = fk.DFT(T_d, Z_slide, Td, l_ref, harmonic, zeroPaddingFactor = 2)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian           
                # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
#                for i in range(np.size(T_d)):
#                    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))*fk.weighting_coeff(T_d[i], slide_pos, FWHM_slide)
                S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
                i += 1
            S[:,:] /= np.size(slide_positions)
            fk.plot_spectrogram(figFT, wn, T_d, S, l_fel/5, l_trans, wn_lim, t_lim)
#            plt.show()

if not interactive_plots:
    plt.close('all')
else:
    plt.show()
    
    
#''' Plots for Paper '''
#ticksize= 2.
#ticklength = 5.
#fontsize=16.
#plt.rcParams['xtick.labelsize'] = fontsize
#plt.rcParams['ytick.labelsize'] = fontsize
#plt.rcParams['axes.labelsize'] = fontsize
#plt.rcParams['xtick.major.width'] = ticksize
#plt.rcParams['xtick.major.size'] = ticklength
#plt.rcParams['ytick.major.width'] = ticksize
#plt.rcParams['ytick.major.size'] = ticklength
#plt.rcParams['axes.linewidth'] = ticksize
#plt.rcParams['lines.linewidth'] = 2.
#
##importing 5H Signal from AC scan with 5th harmonic and Argon (scan2828)
#df = pd.read_csv('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_2828/Z_5H_vs_delay.csv')
#T_2828 = np.array(df['delay'])
#X_5H_2828 = np.array(df['X'])
#
##straight line that guides the eye when looking at the phase
#line =  39.5/400.*Ttheo
#
#figTD = plt.figure(figsize=(15, 8))
#
#ax = figTD.add_subplot(111)
#ax.plot(delay, X/max(abs(X)), color=color,alpha=0.7)
#ax.plot(Ttheo, Xt/max(abs(Xt))*0.7, 'r', alpha=0.3)
##ax.plot(T_2828,Z_5H_2828.real/max(abs(Z_5H_2828)),'g',alpha=0.6)
#ax.fill_between(t_theo_seed,sl_AC*0.4,color='grey')
#ax.fill_between(t_theo_seed,-sl_AC*0.4,color='grey')
#ax.plot(500,5,'grey')
#ax.legend([r'X','X_theo','seed AC'],'upper center')
##ax.plot(Ttheo, scaling*Theo, 'k-')
#ax.set_xlabel(r'$\tau$ [fs]')
#ax.set_ylabel(r'intensity [a.u.]')
#ax2 = ax.twinx()
#ax2.plot(delay, (np.unwrap(np.angle(X+1j*Y))-245.5)/(2*np.pi), color = 'g', alpha=0.7)
#ax2.plot(Ttheo,line,color = 'black', alpha=0.5)
#ax2.set_ylabel(r'phase [$2\pi$ rad]')
##ax2.set_ylim(-20,20)
#ax2.legend([r'$\phi(S)$',r'$0.1 \tau$'],'lower center')
#ax.set_xlim(-410,220)
#ax.set_ylim(-1,1)
##ax.set_yticks([-0.2,-0.1,0.0,0.1,0.2])
##ax.legend(['data','theory'],ncol=2,fontsize=12)
#plt.tight_layout()