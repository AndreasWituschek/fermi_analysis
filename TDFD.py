# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:07:20 2019

@author: Andreas
"""

import fermi_analysis.functions as fk
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
#File I/O:
#file path:
path = "//nanoserver/user/Projekte/BMBF/FERMI/fermi_analysis/masterData.txt"
#import of data
data=fk.importFile(path)

###### TD Analysis
#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266.1 #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)

start=min(data['Delay'])-100
stop=max(data['Delay'])+100
points=20000
plotRange = np.linspace(start,stop,points)

#calculation of theoretical curve
X,Y=fk.curve(l_He,l_ref,h,phi,A,plotRange)

plt.figure()
plt.plot(X)

#plot in phase and in quadrature components of theoretical curve and of data
plt.figure()
fk.plotCurve(X,Y,plotRange)
fk.plotData(data)


######### FD Analysis
#discrete fourier transform of theoretical curve:
Z= X + 1j*Y
step=1/(2*points) #stepsize of theoretical curve in TD
cdft=np.fft.fft(Z)
#downshifted oscillation frequency in PHz:
fDown = abs(299792.4*(h*l_He-l_ref)/(l_ref*l_He))

plt.figure()
plt.plot([i* 10**3./((stop-start)) for i in range(len(cdft))],abs(cdft),linestyle='--')
plt.xlim([0,200])
plt.axvline(x=fDown, color='black', linestyle='--')


#plot



#fk.hello()