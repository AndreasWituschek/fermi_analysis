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
#File I/O:
#file path:
path = "//nanoserver/user/Projekte/BMBF/FERMI/fermi_analysis/exampleData.csv"
##import of data
#data=fk.importFile(path)

#Read the CSV file as a DataFrame:
dataDF = pd.read_csv(path, index_col=0)
#Convert the DataFrame to the dictionary format:
data = dataDF.to_dict("list")

###### TD Analysis
#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266. #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)

start=min(data['Delay'])-10
stop=max(data['Delay'])+10
points=20000
plotRange = np.linspace(start,stop,points)

#calculation of theoretical curve
X,Y=fk.curve(l_He,l_ref,h,phi,A,plotRange)
#plot in phase and in quadrature components of theoretical curve and of data
plt.figure()
fk.plotCurve(X,Y,plotRange)
fk.plotTdData(data)


######### FD Analysis
#discrete fourier transform of theoretical curve:
Z= X + 1j*Y
cdft=np.fft.fft(Z)
#downshifted oscillation frequency in PHz:
fDown = abs(299792.4*(h*l_He-l_ref)/(l_ref*l_He))
#plot downshifted frequency spectrum
fk.plotFdCurve(stop,start,cdft,fDown)
#plot spectrum in eV
fk.plotFdCurveEV(stop,start,cdft,fDown,l_ref,l_He,h)


#discrete fourier transform of data points:
#Zd = data['mx']+ 1j*data['my']





