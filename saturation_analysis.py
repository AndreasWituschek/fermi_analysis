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
import Tkinter, tkFileDialog


root = Tkinter.Tk()
root.withdraw()

file_path = tkFileDialog.askopenfilename()
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

plt.figure('i0 vs ion yield')
plt.title('i0 vs ion yield')
plt.plot(i0[cut::],ion[cut::],'o',alpha=0.01)
plt.xlabel('i0 [uJ]')
plt.ylabel('ion yield [a.u.]')


plt.figure('i0 vs electron yield')
plt.title('i0 vs electron yield')
plt.plot(i0[cut::],elc[cut::],'o',alpha=0.01)
plt.xlabel('i0 [uJ]')
plt.ylabel('electron yield [a.u.]')



plt.figure(2)
plt.hist2d(i0,elc,bins=100,range=[[5,max(i0)],[380,max(ion)]])

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

