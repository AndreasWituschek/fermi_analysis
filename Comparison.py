# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 02:47:44 2019

@author: ldm
"""
import numpy as np
import os
import fermi_analysis.functions as fk
import h5py
import matplotlib.pyplot as plt
import scipy.constants as spc
import pickle
import fermi_analysis.functions as fk


"""Experimental Parameters"""
# parameters for theoretical curve
l_He = 52.22  # helium 1s-4p transition in nm
l_ref = 266.  # reference laser wavelenght in nm
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Comparison'
color = 'r'

""" Parameters for theoretical curve """
phi = 0.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
offset = 0.  # offset


imp_data = []
fd_list = ['/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod0_01',
           '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod0_09',
           '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod1_015',
           '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod1_09',
           '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod2_025',
           '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419_demod2_09']
#            '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Late/delayscan_run409-419']
#fd_list = ['/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/combined/Late/delayscan_run409-419_deomd2_025',
#            '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/combined/Late/delayscan_Run241-252',
#            '/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Morning/delayscan_run395-405']

plot_off = 0
color_idx = np.linspace(0,1, len(fd_list))


figTDY, axY = plt.subplots(1, 1,'TD Y')
axY.set_title('Comparison Y')
axY.set_ylabel('Amplitude in a.u.')

figTDX, axX = plt.subplots(1, 1,'TD X')
axX.set_title('Comparison X')
axX.set_ylabel('Amplitude in a.u.')

figFD, axFD = plt.subplots(1, 1, 'FD')
axFD.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
axFD.set_ylabel('spectral amp. [arb. u.]', fontsize=14)

#fX = plt.figure()

for fn in fd_list:
    imp_data = []
    fp = open(fn , 'rb')
    label = os.path.basename(os.path.normpath(fn))[10:]
    while 1:
        try:
            imp_data.append(pickle.load(fp))
        except EOFError:
            break
    fp.close()
    T = imp_data[0]
    Z = imp_data[1] + 1j*imp_data [3]
    Z_s = imp_data[2] + 1j*imp_data [4]
    R = np.sqrt(Z.real**2 + Z.imag**2)
    R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)


    Phi = np.angle(Z, deg=True)
    Ttheo = np.linspace(T[0],T[-1], 1000)
    Xtd,Ytd,Xt,Yt = fk.Curve(l_He, l_ref, harmonic, -90.0, max(abs(Z)), offset, T[0], T[-1], 1000)

    # Optimal phase
    Phi_theo = T * 1e-6 * (spc.c/l_ref * harmonic - spc.c/l_He) * 360.0
    Phi_theo -= (Phi_theo[0]- Phi[0])
    Phi -= Phi_theo

    delay = T
    X = Z.real
    X_s = Z_s.real
    Y = Z.imag
    Y_s = Z_s.imag

    # Fourier traffo
    Zg = fk.GaussWindow(T, Z, False)
    Td = T[1]-T[0]
    wn, dft = fk.DFT(T, Zg, Td, l_ref , harmonic, zeroPaddingFactor = 2)

    color = plt.cm.jet(color_idx[plot_off])  # get color for each iter

    # Plot X time domain
    axX.errorbar(delay, X, yerr=X_s, color=color, linestyle='o')
    axX.plot(delay, X, 'o', color=color, label=label)
    axX.plot(Ttheo, Xt, '-', color=color)

    # Plot Y time domain
    axY.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
    axY.plot(delay, Y, 'o', color=color, label=label)
    axY.plot(Ttheo, -Yt, '-', color=color)

    # plot frequency domain
    axFD.plot(wn, abs(dft), '-', color = color, linewidth = 2, label=label)
    axFD.axvline(191492.711,color = 'k' ,linestyle='-', linewidth= 2)
    if plot_off <= 0:
        axFD.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)

    plt.show()
    plot_off += 1
    print(plot_off, fn)

axX.legend()
axY.legend()
axFD.legend()
plt.show()