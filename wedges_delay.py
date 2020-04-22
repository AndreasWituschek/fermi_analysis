# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 19:11:54 2019

@author: Andreas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.close('all')

#file_name='//10.5.71.1/User/Projekte/BMBF/FERMI/Fermi Beamtime/19-03-16_18-05_SEED_Elite1_THG2 pulses wedge calibration.txt'
#
file_name='//10.5.71.1/User/Projekte/BMBF/FERMI/Fermi Beamtime/19-03-22_09-55_SEED_Elite1_THGwedge calibration with both pulses.txt'

#file_name='//10.5.71.1/User/Projekte/BMBF/FERMI/Fermi Beamtime/19-03-09_09-28_SEED_Elite1_THG2-pulse seed 262nm 1.25nm wedge pulse 1st- it is the left.txt'


df = pd.read_csv(file_name,sep='\t',header=None)


do_fit =  True

    
def gauss(x,mean,fwhm,amp):
    return amp*np.exp(-4*np.log(2)*(x-mean)**2/(fwhm**2))
    
#plt.plot(df[0].values.astype(np.float),df[1].values.astype(np.float))
runs = np.asarray([0,2,6,8,10,12,14])
runs = np.asarray([0])

#for i in runs:
#    t=df[i].values.astype(np.float)
#
#    intensity=df[i+1].values.astype(np.float)
#
#    plt.plot(t,intensity)
#    plt.xlabel('time [ps]')
#    plt.ylabel('intensity [a.u.]')
#
#plt.legend(runs/2)

    
means = np.array([])
fwhms = np.array([])
for i in runs:
    t=df[i].values.astype(np.float)[:-1]
    b = [not j for j in np.isnan(t)]
    idx = [j for j, x in enumerate(b) if x]
    t = t[idx]

    intensity=df[i+1].values.astype(np.float)[:-1]
#    b = [not j for j in np.isnan(intensity)]
#    idx = [j for j, x in enumerate(b) if x]
    intensity = intensity[idx]
    if do_fit:
        popt, pcov = curve_fit(gauss, t, intensity, p0=[t[np.argmax(intensity)], 0.1,50000.])
        print np.sqrt((popt[1]*1000)**2-30.0**2)
        plt.plot(t,intensity,t,gauss(t,popt[0],popt[1],popt[2]))
        plt.xlabel('time [ps]')
        plt.ylabel('intensity [a.u.]')
        means = np.append(means,popt[0])
        fwhms =  np.append(fwhms,popt[1]*1000)
    else: 
        plt.plot(t,intensity)
        plt.xlabel('time [ps]')
        plt.ylabel('intensity [a.u.]')

#if do_fit:
#    zero_error = (means[0]-means[3])*1000
#    print('The delay zero error is {:.4} fs'.format(zero_error))
#    mean_stepsize = np.mean(np.diff(means[[1,2,3,4,5]]))*1000
#    mean_stepsize_error = np.std(np.diff(means[[1,2,3,4,5]]))*1000
#    mean_fwhm =  np.mean(fwhms)
#    print('The mean stepsize was {:.5}+-{:.5} fs instead of 200 fs. This is a factor of {:.5}'.format(mean_stepsize,mean_stepsize_error,mean_stepsize/200))
#    print('The mean pulse duration was {} fs'.format(mean_fwhm))


