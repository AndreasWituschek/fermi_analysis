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


""" Data parameters """
run = 893
demod = ['1'] #  ['0', '1', '2']
device = ['DEV3265'] #  ['DEV3265', 'DEV3269']
root_file_path = 'C:/Users/FemtoMeasure/Desktop/setup/'
interactive_plots = True
save_fig = False

"""Experimental Parameters"""
# parameters for theoretical curve
l_He = 52.2186  # helium 1s-4p transition in nm
l_ref = 265.98  # reference laser wavelenght in nm
fwhm_FEL = 0.09 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Comparison'
color = 'r'
data_window = [-1000,1000]

""" Parameters for theoretical curve """
phi =-0#-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
slide_step = 100
FWHM = 100
zeroPaddingFactor = 2
suscept = False

# Load preanalysed data
file_path = root_file_path + '/scan_{0}/'.format(int(run))
h5f = h5py.File(file_path + 'scan_{0}.h5'.format(int(run)), 'r')

# sort data by delay in ascending order
sort_inds = np.array(h5f.get('LDM/delay')).argsort()

data = {
    'run_numbers': np.array(h5f.get('run_numbers'))[sort_inds],
    'LDM': {
        'I0': np.array(h5f.get('LDM/I0'))[sort_inds],
        's_I0': np.array(h5f.get('LDM/s_I0'))[sort_inds],
        'delay': np.array(h5f.get('LDM/delay'))[sort_inds],
        's_delay': np.array(h5f.get('LDM/s_delay'))[sort_inds],
        'l_fel': np.array(h5f.get('LDM/l_fel'))[sort_inds], },
    'DEV3265': {
        'harmonic': np.array(h5f.get('DEV3265/harmonic')),
        'x0': np.array(h5f.get('DEV3265/x0'))[sort_inds],
        'y0': np.array(h5f.get('DEV3265/y0'))[sort_inds],
        'x1': np.array(h5f.get('DEV3265/x1'))[sort_inds],
        'y1': np.array(h5f.get('DEV3265/y1'))[sort_inds],
        'x2': np.array(h5f.get('DEV3265/x2'))[sort_inds],
        'y2': np.array(h5f.get('DEV3265/y2'))[sort_inds],
        's_x0': np.array(h5f.get('DEV3265/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('DEV3265/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('DEV3265/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('DEV3265/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('DEV3265/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('DEV3265/s_y2'))[sort_inds], },
    'DEV3269': {
        'harmonic': np.array(h5f.get('DEV3269/harmonic')),
        'x0': np.array(h5f.get('DEV3269/x0'))[sort_inds],
        'y0': np.array(h5f.get('DEV3269/y0'))[sort_inds],
        'x1': np.array(h5f.get('DEV3269/x1'))[sort_inds],
        'y1': np.array(h5f.get('DEV3269/y1'))[sort_inds],
        'x2': np.array(h5f.get('DEV3269/x2'))[sort_inds],
        'y2': np.array(h5f.get('DEV3269/y2'))[sort_inds],
        's_x0': np.array(h5f.get('DEV3269/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('DEV3269/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('DEV3269/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('DEV3269/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('DEV3269/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('DEV3269/s_y2'))[sort_inds], }, }
h5f.close()

T = data['LDM']['delay']

''' laser spectrum '''
points = np.linspace(180000.0,195000.0,1000)
fwhm_FEL_wn = 1E7*(fwhm_FEL/data['LDM']['l_fel'][0]**2) # wavenumber FWHM of FEL
l_FEL = 1E7/data['LDM']['l_fel'][0] # wavenumber CWL of FEL
spectrum = np.exp(-4.0*0.693147*(points-l_FEL)**2/fwhm_FEL_wn**2) # spectrum in wavenumber space

# loop over the mfli devices
for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = data[dev]['x' + d] + 1j * data[dev]['y' + d]
        Z_s = data[dev]['s_x' + d] + 1j * data[dev]['s_y' + d]

        R = np.sqrt(Z.real**2 + Z.imag**2)
        R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)
        Phi = np.angle(Z, deg=True)

# Create theoretical curve
#        if draw_theory:
        Ttheo = np.linspace(T[0],T[-1], 1000)
        Xtd,Ytd,Xt,Yt = fk.Curve(l_He, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], 1000)

        delay = T
        X = Z.real
        X_s = Z_s.real
        Y = Z.imag
        Y_s = Z_s.imag

        T_d,X_d,Y_d = fk.CutDataSet(T, X, Y, data_window)
        Z = X_d + 1j*Y_d

        # Fourier traffo
        Zg = fk.GaussWindow(T_d, Z, suscept)
        Td = T[1]-T[0]
        wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        FWHM , peakFullWidth = fk.PeakWidth(T, 0.1, False, zeroPaddingFactor)
        FWHM = FWHM

        # Plot time domain
        figTD = plt.figure('scan_{0}_{1}_d{2}_TD_{3}'.format(run, dev, d, i), figsize=(9, 12))

        ax = figTD.add_subplot(511)
        ax.set_title('scan_{0}_{1}_d{2}_TD_{3}'.format(run, dev, d, i))
        ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
        ax.plot(delay, X, 'b-')
        ax.plot(Ttheo, Xt, 'b', alpha=0.3)
        ax.set_ylabel('Amplitude in a.u.')
        ax.grid()

        ax = figTD.add_subplot(512)
        ax.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
        ax.plot(delay, Y, '-', color=color)
        ax.plot(Ttheo, Yt, 'r', alpha=0.3)
        ax.set_ylabel('Amplitude in a.u.')
        ax.grid()

        ax = figTD.add_subplot(513)
        ax.errorbar(delay, R, yerr=R_s, color=color, linestyle='')
        ax.plot(delay, R, '-', color=color)
        ax.set_ylabel('R')
        ax.grid()

        ax = figTD.add_subplot(514)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(delay, Phi, '-', color=color)
        ax.set_ylabel('Phi')
        ax.grid()

        ax = figTD.add_subplot(515)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(delay, np.abs(np.append(np.diff(delay)[0], np.diff(delay))), '-', color=color)
        ax.set_ylabel('step sizes')
        ax.set_xlabel('Delay in fs')
        ax.set_ylim(0,3)
        ax.grid()

        plt.tight_layout()
        if save_fig:
            plt.savefig(file_path + 'scan_{0}_{1}_d{2}_TD_{3}.png'.format(run, dev, d, i), dpi=400)

        # plot frequency domain
        figFT = plt.figure('scan_{0}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs = figFT.add_subplot(111)
        axs.plot(wn, abs(dft), '-', color=color, linewidth = 2)
        axs.plot(points-(1-mfli_harmonic/5.)*1E7*(1/l_He-1/(l_ref/5.0)),spectrum*(max(abs(dft))-min(abs(dft)))+min(abs(dft)),'-',color='grey',linewidth = 2)
        axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
        axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
        axs.set_title('scan_{0}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs.set_xlim(180000.,195000.)
        axs.axvline(1E7/l_He,color = 'k' ,linestyle='-', linewidth= 2)
        axs.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(axs.get_xlim())[0]
        plot_y_range = np.diff(axs.get_ylim())[0]
        axs.axhline(y= plot_y_range * 0.9 + axs.get_ylim()[0], xmin=0.1, xmax=0.1 + FWHM / plot_x_range, linewidth=3)
#        axs.text(0.1, plot_y_range * 0.9, str(FWHM))
        axs.grid()
        if save_fig:
            plt.savefig(file_path + 'scan_{0}_{1}_d{2}_FD_{3}.png'.format(run, dev, d, i), dpi=400)
        
        # spectrogram
        figSG = plt.figure('scan_{0}_{1}_d{2}_SGn_{3}'.format(run, dev, d, i))
#        plt.specgram(Z, Fs=1./np.mean(np.diff(T_d))*1e15, NFFT=32, noverlap=16)
        f, t, specr = signal.spectrogram(np.real(Z), fs=1./np.mean(np.diff(T_d*1e-15)),
                                        nperseg=32, noverlap=31,
                                        return_onesided=True, # mode='complex',
                                        scaling='density', nfft=256,
                                        window=signal.get_window(('gaussian', int(32 / np.sqrt(48))), 32))
        f, t, speci = signal.spectrogram(np.real(Z), fs=1./np.mean(np.diff(T_d*1e-15)),
                                        nperseg=32, noverlap=31,
                                        return_onesided=True, # mode='complex',
                                        scaling='density', nfft=256,
                                        window=signal.get_window(('gaussian', int(32 / np.sqrt(48))), 32))
        f = f/100.0/spc.c + harmonic*10**7/l_ref

        spec = specr + 1j * speci
        print np.min(np.abs(spec))
        print np.max(np.abs(spec))
        print type(spec[0][0])
        print spec.shape
#        plt.imshow(np.real(spec))
        plt.pcolormesh(t, f, np.real(spec), vmin= np.min(np.abs(spec)), vmax=np.max(np.abs(spec))*1.0)
        plt.hlines(1E7/l_He, np.min(t), np.max(t))
        plt.show()

##""" slide fourier trafo """
#        S = np.zeros((np.size(T_d), np.size(wn)))
#        slide_positions = np.arange(T_d[0],T_d[-1], slide_step)
#        for slide_pos in slide_positions:
#            Z_slide = fk.slide_window(T_d, Z, slide_pos, FWHM)    # multiply gaussian onto data set which is centered at current sliding position
#            DFT_slide, wn = fk.DFT(T_d, Z_slide, Td, l_ref , harmonic, zeroPaddingFactor = 2)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
#            # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
#            for i in range(np.size(T_d)):
#                S[i,:] += DFT_slide*fk.weighting_coeff(T_d[i], slide_pos, FWHM)
#        # normalize to get proper average
#        S[:,:] /= np.size(slide_positions)
#        figSlide = plt.figure('slide')
#        plt.imshow(S.transpose(), origin='lower', aspect='auto')
#        plt.show()

if not interactive_plots:
    plt.close('all')