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
run = 31
demod = ['1']
device = ['dev3265']
#demod = ['0','1', '2']
#device = ['dev3265', 'dev3269']
root_file_path = '/home/ldm/ExperimentalData/Online4LDM/RBT-UOF_4/Data/combined/'
interactive_plots = True
save_fig = False
plotTheo = False
spectrogram = False

"""Experimental Parameters"""
# parameters for theoretical curve
#l_trans = 52.2186  # helium 1s-4p transition in nm
l_trans = 260.92/6.0 # Fano resonance in argon
l_ref = 265.98  # reference laser wavelenght in nm
fwhm_FEL = 0.09 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 6.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Comparison'
color = ['r','b']
data_window = [-25, 300] #[-200,370]
gauss = False # gauss window on TD (true) or not (False)

""" Parameters for theoretical curve """
phi =0.#-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
slide_step = 31
FWHM_slide = slide_step / np.sqrt(48)
zeroPaddingFactor = 2
suscept = True
t_corr =1  # correction on time zero
phi_sus = -150.0  # phase correction of data
wn_lim = [217100,233800]

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
l_fel = data['LDM']['l_seed'][0]
#l_fel = 262.2
l_fel = 43.48*6.
print l_fel
T = data['LDM']['delay']

''' laser spectrum '''
points = np.linspace(180000.0,195000.0,1000)
fwhm_FEL_wn = 1E7*(fwhm_FEL/(l_fel/5)**2) # wavenumber FWHM of FEL
l_FEL = 1E7/l_fel # wavenumber CWL of FEL
spectrum = np.exp(-4.0*0.693147*(points-l_FEL*5)**2/fwhm_FEL_wn**2) # spectrum in wavenumber space

Ttheo = np.linspace(T[0],T[-1], np.size(T))
AC_FWHM = 100.0  # in fs
Fano_FWHM = 140.0  # in fs
AC = np.exp(-4.0*np.log(2)*((Ttheo)/AC_FWHM)**2)
Fano = np.exp((Ttheo)/Fano_FWHM)


# loop over the mfli devices
for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = data[dev]['x' + d] + 1j * data[dev]['y' + d]
        Z_s = data[dev]['s_x' + d] + 1j * data[dev]['s_y' + d]

# Create theoretical curve
#        if draw_theory:
        Xtd,Ytd,Xt,Yt = fk.Curve(l_trans, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], np.size(T))
        Xtd_gaus,Ytd_gaus,Xt_gaus,Yt_gaus = fk.Curve(l_fel/harmonic, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0]-t_corr, T[-1]-t_corr, np.size(T))
        Xt *= Fano
        Yt *= Fano

#        Xt -= 1*Xt_gaus*abs(data[dev]['x0'] + 1j * data[dev]['y0'])/max(abs(data[dev]['x0'] + 1j * data[dev]['y0']))
#        Yt -= 1*Yt_gaus*abs(data[dev]['x0'] + 1j * data[dev]['y0'])/max(abs(data[dev]['x0'] + 1j * data[dev]['y0']))

        # correct time and phase for susceptibility
        T += t_corr
        Z *= np.exp(1j*np.pi*phi_sus/180.0) # correction of phase
        if plotTheo:
            T = Ttheo + t_corr
            Z = Xt + 1j*Yt

        # do Fourier transform
        T_d,X_d,Y_d = fk.CutDataSet(T, Z.real, Z.imag, data_window)
        Z_d = X_d + 1j*Y_d
        # Fourier traffo
        Zg = fk.GaussWindow(T_d, Z_d, suscept)
        Td = T[1]-T[0]
        if min(T_d) <0:
            T_d = -np.flipud(T_d)
            Z_d = np.flipud(Z_d)
            Td = -Td
        print T_d
        if gauss:
            T = T_d
            Z = Zg
            wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        else:
            wn, dft = fk.DFT(T_d, Z_d, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        FWHM , peakFullWidth = fk.PeakWidth(T, 0.1, False, zeroPaddingFactor)

        # calculate R and phi
        R = np.sqrt(Z.real**2 + Z.imag**2)
        R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)
        Phi = np.angle(Z, deg=False)

        # Plot time domain
        figTD = plt.figure('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i), figsize=(9, 12))

        ax = figTD.add_subplot(511)
        ax.set_title('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i))
#        ax.errorbar(T, Z.real, yerr=Z_s.real, color='b', linestyle='')
        ax.plot(T, Z.real, 'b-')
        if plotTheo: ax.plot(Ttheo, Xt, 'b', alpha=0.3)
        ax.set_ylabel('Amplitude in a.u.')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(512)
#        ax.errorbar(T, Z.imag, yerr=Z_s.imag, color=color, linestyle='')
        ax.plot(T, Z.imag, '-', color=color[0])
        if plotTheo: ax.plot(Ttheo, Yt, 'r', alpha=0.3)
        ax.set_ylabel('Amplitude in a.u.')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(513)
#        ax.errorbar(T, R, yerr=R_s, color=color, linestyle='')
        ax.plot(T, R, '-', color=color[0])
        ax.set_ylabel('R')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(514)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(T, np.unwrap(Phi), '-', color=color[0])
        ax.set_ylabel('Phi (degree)')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(515)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(T, np.abs(np.append(np.diff(T)[0], np.diff(T))), '-', color=color[0])
        ax.set_ylabel('step sizes (fs)')
        ax.set_xlabel('Delay in fs')
        ax.set_ylim(0,3)
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        plt.tight_layout()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_TD_{3}.png'.format(run, dev, d, i), dpi=400)

        # plot frequency domain
        figFT = plt.figure('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs = figFT.add_subplot(111)
        axs.plot(wn, abs(dft), '-', color='k', linewidth = 1)
        axs.plot(wn, dft.real, '-', color=color[0], linewidth = 2)
        axs.plot(wn, dft.imag, '-', color=color[1], linewidth = 2)
        axs.plot(points-(1-mfli_harmonic/5.)*1E7*(1/l_trans-1/(l_ref/5.0)),spectrum*(max(abs(dft))-min(abs(dft)))+min(abs(dft)),'-',color='grey',linewidth = 2)
        axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
        axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
        axs.set_title('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs.set_xlim(wn_lim[0],wn_lim[1])
        axs.axvline(1E7/l_trans,color = 'k' ,linestyle='-', linewidth= 2)
        axs.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(axs.get_xlim())[0]
        plot_y_range = np.diff(axs.get_ylim())[0]
        axs.axhline(y= plot_y_range * 0.9 + axs.get_ylim()[0], xmin=0.1, xmax=0.1 + FWHM / plot_x_range, linewidth=3)
#        axs.text(0.1, plot_y_range * 0.9, str(FWHM))
        axs.grid()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_FD_{3}.png'.format(run, dev, d, i), dpi=400)

        if spectrogram:
            # plot spec
            figFT = plt.figure('scan_{0:03}_{1}_d{2}_SG_{3}'.format(run, dev, d, i))
            print(T_d[0], T_d[-1])
            slide_step = 2.0  #in fs
            FWHM_slide = 20.0  #in fs
            t_lim = data_window
            slide_positions = np.arange(T_d[0],T_d[-1], slide_step)
            S = np.zeros((np.size(T_d), np.size(wn)))
            for slide_pos in slide_positions:
                Z_slide = fk.slide_window(T_d, Z, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
                wn,  DFT_slide = fk.DFT(T_d, Z_slide, Td, l_ref , harmonic, zeroPaddingFactor = 2)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
                # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
                for i in range(np.size(T_d)):
                    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))*fk.weighting_coeff(T_d[i], slide_pos, FWHM_slide)
            S[:,:] /= np.size(slide_positions)
            fk.plot_spectrogram(figFT, wn, T_d, S, l_fel, wn_lim, t_lim)
            plt.show()

if not interactive_plots:
    plt.close('all')
else:
    plt.show()