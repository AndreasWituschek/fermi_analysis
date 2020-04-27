# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 21:59:36 2019

@author: Andi

Analyzes data that have been preanalyzed with the program saturation_preanalysis.py 
This way saturation effects of HGHG process and the ion/electron detectors can be analyzed.
"""
import numpy as np
#import os
#import errno
import h5py
import matplotlib.pyplot as plt


scan = 5395

file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_{}/scan_{}_saturation_preanalysis.h5'.format(scan,scan)
data = h5py.File(file_path, 'r')

i0 = np.array(data['/LDM/I0'])
ion = np.array(data['/LDM/ion_avg'])
elc = np.array(data['/LDM/elc_avg'])
delay = np.array(data['/LDM/delay'])

data.close()

#n= 40
#test = np.max(np.reshape(i0,(len(i0)/n,n)),axis=1)
#ion_avg = np.max(np.reshape(-ion+442,(len(ion)/n,n)),axis=1)*4.
#test = np.min(np.reshape(i0,(len(i0)/n,n)),axis=1)
#plt.plot(test)
#plt.plot(ion_avg)

delay_plot = np.linspace(min(delay),max(delay),len(i0))[::-1]
cut = np.max(np.argwhere(delay_plot>-200))
cut=0
tt = np.linspace(0,70,1000)


plt.figure('i0 vs ion yield')
plt.title('i0 vs ion yield')
plt.plot(i0[cut::],ion[cut::],'o',alpha=0.04)
plt.plot(tt,-0.31*tt+400.3,color='k',alpha = 0.5)
plt.xlabel('i0 [uJ]')
plt.ylabel('ion yield [a.u.]')


plt.figure('i0 vs electron yield')
plt.title('i0 vs electron yield')
plt.plot(i0[cut::],elc[cut::],'o',alpha=0.01)
plt.xlabel('i0 [uJ]')
plt.ylabel('electron yield [a.u.]')


plt.figure(2)
plt.hist2d(i0,ion,bins=100,range=[[6,max(i0)],[370,max(ion)]],cmap ='Blues')
plt.plot(tt,-0.31*tt+400.3,color='k',alpha = 0.5)
plt.xlabel('i0 [uJ]')
plt.ylabel('ion yield [a.u.]')

plt.figure('ion yield vs delay')
plt.title('ion yield vs delay')
plt.plot(delay_plot,ion,alpha=0.7)
plt.xlabel('delay [fs]')
plt.ylabel('ion yield [a.u.]')


plt.figure('electron yield vs delay')
plt.title('electron yield vs delay')
plt.plot(delay_plot,ion,alpha=0.7)
plt.xlabel('delay [fs]')
plt.ylabel('electron yield [a.u.]')

plt.figure('i0 vs delay')
plt.title('i0 vs delay')
plt.plot(delay_plot,i0,alpha=0.7)
plt.xlabel('delay [fs]')
plt.ylabel('i0 [uJ]')

