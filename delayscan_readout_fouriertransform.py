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
data_window = [-1000,1000]

""" Parameters for theoretical curve """
phi =-0#-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analzsis parameters """
slide_step = 100
FWHM = 100

""" Experimental parameters """
Runs = [1748, 1799]
delays = [240, 260]

imp_data = []
#fp = open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_1/combined/Night/delayscan_run38-144' , 'rb')
#fp = open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/combined/Late/delayscan_Run231-239' , 'rb')
#fp = open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/combined/Late/delayscan_Run220-230_demod2' , 'rb')
#fp = open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_2/combined/Late/delayscan_Run241-252' , 'rb')
#fp = open('/home/ldm/ExperimentalData/Online4LDM/20149020/Day_3/combined/Morning/delayscan_run395-405' , 'rb')
#fp = open('/home/ldm/Desktop/PM_Data/PM_Data_Day3/Late/delayscan_467-568_demod1' , 'rb')

file_path = '/home/ldm/ExperimentalData/Online4LDM/20149020/results/PM_Data/PM_Data_Day5/Late/'
file_path += 'run_{0}-{1}_ion_{2}-{3}fs/delayscan_run_{0}-{1}_ion_{2}-{3}fs_5_1'.format(Runs[0], Runs[1] - 1, delays[0], delays[1])

fp = open(file_path, 'rb')

while 1:
    try:
        imp_data.append(pickle.load(fp))
    except EOFError:
        break
T = imp_data[0]
Z = imp_data[1] + 1j*imp_data [3]
Z_s = imp_data[2] + 1j*imp_data [4]
T = T[1:]
Z = Z[1:]
Z_s = Z_s[1:]

R = np.sqrt(Z.real**2 + Z.imag**2)
R_s = np.sqrt(Z_s.real**2 + Z_s.imag**2)
Phi = np.angle(Z, deg=True)


#Phi = np.arctan(Z.imag/Z.real)/np.pi*180.0
#k = 0
#for n in range(len(Phi)-1):
#    if Phi[n] < Phi[n+1]:
#        k +=1
#    Phi[n] += k*180.0
#Phi[-1] += k*180.0

#Ttheo = np.linspace(T[0],T[-1], 1000)
#Xtd,Ytd,Xt,Yt = fk.Curve(l_He, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], 1000)

# Optimal phase
Phi_theo = T * 1e-6 * (spc.c/l_ref * harmonic - spc.c/l_He) * 360.0
Phi_theo -= (Phi_theo[0]- Phi[0])
Phi -= Phi_theo

delay = T
X = Z.real
X_s = Z_s.real
Y = Z.imag
Y_s = Z_s.imag

T_d,X_d,Y_d = fk.CutDataSet(T, X, Y, data_window)
Z = X_d + 1j*Y_d

# Fourier traffo
Zg = fk.GaussWindow(T_d, Z, False)
Td = T[1]-T[0]
wn, dft = fk.DFT(T_d, Zg, -Td, l_ref , harmonic, zeroPaddingFactor = 2)


# Plot time domain
figTD = plt.figure('TD')
ax = figTD.add_subplot(411)
ax.set_title(title)

ax.errorbar(delay, X, yerr=X_s, color='k', linestyle='')
ax.plot(delay, X, 'ko')
#ax.plot(Ttheo, Xt, 'k-')
ax = figTD.add_subplot(412)
ax.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
ax.plot(delay, Y, 'o', color=color)
#ax.plot(Ttheo, -Yt, '-', color=color)
ax.set_ylabel('Amplitude in a.u.')
##ax.set_xlim(T[0]+1,T[-1]-1)

ax = figTD.add_subplot(413)
ax.errorbar(delay, R, yerr=R_s, color=color, linestyle='')
ax.plot(delay, R, 'o-', color=color)
ax.set_ylabel('R')
#ax.set_xlim(T[0],T[-1])

ax = figTD.add_subplot(414)
#ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
ax.plot(delay, Phi, 'o-', color=color)
ax.set_ylabel('Phi')
ax.set_xlabel('Delay in fs')
ax.set_xlim(T[0], T[-1])
plt.savefig(file_path + '_TD.png', dpi=400)


# plot frequency domain
figFT = plt.figure('FT')
axs = figFT.add_subplot(111)
axs.plot(wn, abs(dft), '-', color=color, linewidth = 2)
axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
axs.axvline(191492.711,color = 'k' ,linestyle='-', linewidth= 2)
axs.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
plt.savefig(file_path + '_FT.png', dpi=400)

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

plt.show()


print(np.where(abs(dft) == np.max(abs(dft)), abs(dft)))