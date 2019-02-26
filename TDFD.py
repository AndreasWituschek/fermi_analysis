# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:07:20 2019

@author: Andreas
"""

import fermi_analysis.functions as fk
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')


"""Experimental Parameters"""
#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266. #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)
offset=0.3 #offset


"""File I/O"""
path = "//nanoserver/user/Projekte/BMBF/FERMI/fermi_analysis/masterData.csv" #file path of masterData.csv:
data=fk.MasterFileReader(path) #importin masterData.csv as dictionary


"""Plotting parameters"""
demod = 0 #in [0,1,2,3]. Choose which demodulator channel to analyse and plot
downshifted= False #plot downshifted spectrum
absolute=True #plot absorption spectrum
absorptive=False #plot absorptive and dispersive spectrum


"""TD Analysis & Plotting"""
#delay parameters
start=min(data['delay'])
stop=max(data['delay'])
length=len(data['delay'])

#calculation of theoretical curve
Xtd,Ytd,X,Y=fk.Curve(l_He,l_ref,h,phi,A,offset, start, stop, length) #Xtd/Ytd is time domain data, where 10fs are added/substracted to start/stop of dataset,X and Y start/stop at edges of data set and are used for DFT
#plot in phase and in quadrature components of theoretical curve and of data
plt.figure()
fk.PlotCurve(Xtd,Ytd,start, stop)
fk.PlotTdData(data,demod)


"""FD Analysis & Plotting"""
""""""""""""""""""
#discrete complex fourier transform of theoretical curve:
Z= X + 1j*Y
cdft=np.fft.fft(Z)
#discrete complex fourier transform of experimental data:
Z_d = np.asarray(data['mX%d' % demod])[np.argsort(data['delay'])] + 1j*np.asarray(data['mY%d' % demod])[np.argsort(data['delay'])] #arrays sorted by delay for DFT
cdft_d=np.fft.fft(Z_d)

#### plotting FD results. Maximum value of theoretical curve is scaled to maximum of experimental data.
#plot downshifted frequency spectrum
if downshifted: 
    fk.PlotFdCurve(stop,start,cdft,cdft_d,l_ref,l_He,h)

#plot absorption spectrum in eV
if absolute:
    fk.PlotFdCurveAbsEV(stop,start,cdft,cdft_d,l_ref,l_He,h)

#plot dispersive and absorptive part of spectrum
if absorptive:
    fk.PlotFdCurveEV(stop,start,cdft,cdft_d,l_ref,l_He,h)


