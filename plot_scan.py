# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:31:42 2019

@author: ldm
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 02:47:44 2019

@author: ldm
"""
import numpy as np
import fermi_analysis.functions as fk
import h5py
import matplotlib.pyplot as plt
import scipy.constants as spc


"""Experimental Parameters"""
# parameters for theoretical curve
l_He = 52.2186  # helium 1s-4p transition in nm
l_ref = 265.98  # reference laser wavelenght in nm
l_FEL = 52.22 # FEL wavelength in nm, used for calculation of FEL spectrum
fwhm_FEL = 0.09 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Comparison'
color = 'r'
data_window = [-1000,1000]


''' laser spectrum '''
points = np.linspace(180000.0,195000.0,1000)
fwhm_FEL = 1E7*(fwhm_FEL/l_FEL**2)
l_FEL = 1E7/l_FEL # wavenumber CWL of FEL
spectrum = np.exp(-4.0*0.693147*(points-l_FEL)**2/fwhm_FEL**2) # spectrum in wavenumber space



""" Parameters for theoretical curve """
phi =-0#-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
slide_step = 100
FWHM = 100
zeroPaddingFactor = 2

""" Experimental parameters """
run = 1453
#delays = [50, -450]
#demod_harmonic = ['5_0', '5_1', '5_2', '9_0', '9_1', '9_2']
#demod_harmonic = ['5_1']
demod = ['0', '1', '2']
device = ['DEV3265', 'DEV3269']
assignment = {
    'DEV3265': {
        '0': '4H',
        '1': '5H',
        '2': '10H',},
    'DEV3269': {
        '0': '1H',
        '1': '2H',
        '2': '3H',}}

file_path = '//online4ldm.esce.elettra.trieste.it/store/20149020/results/HeDroplet/scan_{0}/'.format(int(run))
h5f = h5py.File(file_path + 'scan_{0}.h5'.format(int(run)), 'r')


data = {
    'run_numbers': np.array(h5f.get('run_numbers')),
    'I0': np.array(h5f.get('I0')),
    's_I0': np.array(h5f.get('s_I0')),
    'delay': np.array(h5f.get('delay')),
    's_delay': np.array(h5f.get('s_delay')),
    'DEV3265': {
        'x0': np.array(h5f.get('DEV3265/x0')),
        'y0': np.array(h5f.get('DEV3265/y0')),
        'x1': np.array(h5f.get('DEV3265/x1')),
        'y1': np.array(h5f.get('DEV3265/y1')),
        'x2': np.array(h5f.get('DEV3265/x2')),
        'y2': np.array(h5f.get('DEV3265/y2')),
        's_x0': np.array(h5f.get('DEV3265/s_x0')),
        's_y0': np.array(h5f.get('DEV3265/s_y0')),
        's_x1': np.array(h5f.get('DEV3265/s_x1')),
        's_y1': np.array(h5f.get('DEV3265/s_y1')),
        's_x2': np.array(h5f.get('DEV3265/s_x2')),
        's_y2': np.array(h5f.get('DEV3265/s_y2')), },
    'DEV3269': {
        'x0': np.array(h5f.get('DEV3265/x0')),
        'y0': np.array(h5f.get('DEV3265/y0')),
        'x1': np.array(h5f.get('DEV3265/x1')),
        'y1': np.array(h5f.get('DEV3265/y1')),
        'x2': np.array(h5f.get('DEV3265/x2')),
        'y2': np.array(h5f.get('DEV3265/y2')),
        's_x0': np.array(h5f.get('DEV3265/s_x0')),
        's_y0': np.array(h5f.get('DEV3265/s_y0')),
        's_x1': np.array(h5f.get('DEV3265/s_x1')),
        's_y1': np.array(h5f.get('DEV3265/s_y1')),
        's_x2': np.array(h5f.get('DEV3265/s_x2')),
        's_y2': np.array(h5f.get('DEV3265/s_y2')), }}

T = data['delay']

for dev in device:
    for d in demod:
        Z = data[dev]['x' + d] + 1j * data[dev]['y' + d]
        Z_s = data[dev]['s_x' + d] + 1j * data[dev]['s_y' + d]

        R = np.sqrt(Z.real**2 + Z.imag**2)
        R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)
        Phi = np.angle(Z, deg=True)


        Ttheo = np.linspace(T[0],T[-1], 1000)
        Xtd,Ytd,Xt,Yt = fk.Curve(l_He, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], 1000)

        # Optimal phase
        Phi_theo = T * 1e-6 * (spc.c/l_ref * harmonic - spc.c/l_He) * 360.0
        Phi_theo -= (Phi_theo[0]- Phi[0])
        Phi -= Phi_theo

        delay = T
        X = Z.real
        X_s = Z_s.real
        Y = Z.imag
        Y_s = Z_s.imag

        #print T
        T_d = T[:]
        X_d = X[:]
        Y_d = Y[:]

        T_d,X_d,Y_d = fk.CutDataSet(T, X, Y, data_window)
        Z = X_d + 1j*Y_d


        # Fourier traffo
        Zg = fk.GaussWindow(T_d, Z, False)
        Td = T[1]-T[0]
        wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
#        FWHM , peakFullWidth = fk.PeakWidth(T, 0.1, False, zeroPaddingFactor)
#        FWHM = FWHM
#        print(FWHM)


        # Plot time domain
        figTD = plt.figure('scan_{0}_TD_{1}.png'.format(run, assignment[dev][d]))
        ax = figTD.add_subplot(411)
        ax.set_title('scan_{0}_TD_{1}'.format(run, assignment[dev][d]))

        ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
        ax.plot(delay, X, 'b-')
        ax.grid()
    #    ax.plot(Ttheo, Xt, 'k-')
        ax = figTD.add_subplot(412)
        ax.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
        ax.plot(delay, Y, '-', color=color)
        #ax.plot(Ttheo, -Yt, '-', color=color)
        ax.set_ylabel('Amplitude in a.u.')
        ax.grid()
        ##ax.set_xlim(T[0]+1,T[-1]-1)

        ax = figTD.add_subplot(413)
        ax.errorbar(delay, R, yerr=R_s, color=color, linestyle='')
        ax.plot(delay, R, '-', color=color)
        ax.set_ylabel('R')
        ax.grid()
        #ax.set_xlim(T[0],T[-1])

        ax = figTD.add_subplot(414)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(delay, Phi, '-', color=color)
        ax.set_ylabel('Phi')
        ax.set_xlabel('Delay in fs')
        ax.grid()
#        ax.set_xlim(T[0], T[-1])
        plt.tight_layout()
        plt.savefig(file_path + 'scan_{0}_TD_{1}.png'.format(run, assignment[dev][d]), dpi=400)
#        plt.close()


        # plot frequency domain
        figFT = plt.figure('scan_{0}_FD_{1}.png'.format(run, assignment[dev][d]))
        axs = figFT.add_subplot(111)
        axs.plot(wn, abs(dft), '-', color=color, linewidth = 2)
        axs.plot(points,spectrum*(max(abs(dft))-min(abs(dft)))+min(abs(dft)),'-',color='grey',linewidth = 2)
        axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
        axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
        axs.set_title('scan_{0}_FD_{1}'.format(run, assignment[dev][d]))
        axs.set_xlim(180000.,195000.)
        axs.axvline(1E7/l_He,color = 'k' ,linestyle='-', linewidth= 2)
        axs.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(axs.get_xlim())[0]
        plot_y_range = np.diff(axs.get_ylim())[0]
        axs.axhline(y= plot_y_range * 0.9 + axs.get_ylim()[0], xmin=0.1, xmax=0.1 + FWHM / plot_x_range, linewidth=3)
#        axs.text(0.1, plot_y_range * 0.9, str(FWHM))
        axs.grid()
        plt.savefig(file_path + 'scan_{0}_FD_{1}.png'.format(run, assignment[dev][d]), dpi=400)
#        plt.close()

        #np.savefig(figFT, dpi=600)

        """ slide fourier trafo """
        #S = np.zeros((np.size(T_d), np.size(wn)))
        #slide_positions = np.arange(T_d[0],T_d[-1], slide_step)
        #for slide_pos in slide_positions:
        #    Z_slide = fk.slide_window(T_d, Z, slide_pos, FWHM)    # multiply gaussian onto data set which is centered at current sliding position
        #    DFT_slide, wn = fk.DFT(T_d, Z_slide, Td, l_ref , harmonic, zeroPaddingFactor = 2)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
        #    # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
        #    for i in range(np.size(T_d)):
        #        S[i,:] += DFT_slide*fk.weighting_coeff(T_d[i], slide_pos, FWHM)
        ## normalize to get proper average
        #S[:,:] /= np.size(slide_positions)
        #figSlide = plt.figure('slide')
        #plt.imshow(S.transpose(), origin='lower', aspect='auto')
#        plt.show()