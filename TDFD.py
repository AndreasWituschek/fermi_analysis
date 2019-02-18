# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:07:20 2019

@author: Andreas
"""

import fermi_analysis.functions as fk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')

"""File I/O"""
""""""""""""""""""
#file path:
path = "//nanoserver/user/Projekte/BMBF/FERMI/fermi_analysis/exampleData.csv"

#Import file with data from MFLI:
dataDF = pd.read_csv(path, index_col=0)
#Convert the DataFrame to the dictionary format:
data = dataDF.to_dict("list")

"""Experimental Parameters"""
#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266. #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)

"""Plotting parameters"""
points=20000 # number of TD points for plotting of theroetical curve 
downshifted= False #plot downshifted spectrum
absolute=True #plot absorption spectrum
absorptive=False #plot absorptive and dispersive spectrum


"""TD Analysis & Plotting"""
""""""""""""""""""
#delay parameters
start=min(data['Delay'])
stop=max(data['Delay'])
length=len(data['Delay'])

#calculation of theoretical curve
Xtd,Ytd,X,Y=fk.curve(l_He,l_ref,h,phi,A,start, stop, points, length) #Xtd/Ytd is time domain data, where 10fs are added/substracted to start/stop of dataset, so one can see where next datapoint should lie
#plot in phase and in quadrature components of theoretical curve and of data
plt.figure()
fk.plotCurve(Xtd,Ytd,start, stop, points)
fk.plotTdData(data)


"""FD Analysis & Plotting"""
""""""""""""""""""
#discrete complex fourier transform of theoretical curve:
Z= X + 1j*Y
cdft=np.fft.fft(Z)
#discrete complex fourier transform of experimental data:
Z_d = np.asarray(data['mx'])+ 1j*np.asarray(data['my'])
cdft_d=np.fft.fft(Z_d)

#### plotting FD results. Maximum value of theoretical curve is scaled to maximum of experimental data.
#plot downshifted frequency spectrum
if downshifted: 
    fk.plotFdCurve(stop,start,cdft,cdft_d,l_ref,l_He,h)

#plot absorption spectrum in eV
if absolute:
    fk.plotFdCurveAbsEV(stop,start,cdft,cdft_d,l_ref,l_He,h)

#plot dispersive and absorptive part of spectrum
if absorptive:
    fk.plotFdCurveEV(stop,start,cdft,cdft_d,l_ref,l_He,h)






