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


file_name='//10.5.71.1/user/Projekte/BMBF/Data/Fermi/190122CrossCorrelations/19-01-23_11-45_SEED_Elite1_THG7 deg wedge calibration.txt'
df = pd.read_csv(file_name,sep='\t',header=None)

    
def gauss(x,mean,std,amp):
    return amp*np.exp(-(x-mean)**2/(2*std**2))
    
    
means = np.array([])
for i in [0,2,4,6,8,10]:
    t=df[i].values.astype(np.float)[:-1]
    b = [not j for j in np.isnan(t)]
    idx = [j for j, x in enumerate(b) if x]
    t = t[idx]

    intensity=df[i+1].values.astype(np.float)[:-1]
    b = [not j for j in np.isnan(intensity)]
    idx = [j for j, x in enumerate(b) if x]
    intensity = intensity[idx]
    print(sum(np.isnan(intensity)))


    plt.plot(t,intensity)
    plt.xlabel('time [ps]')
    plt.ylabel('intensity [a.u.]')
    
means=[5815,5556,5306]
    
#plt.plot(df[0].values.astype(np.float),df[1].values.astype(np.float))
#    
#means = np.array([])
#for i in [2,4,6,8,10,12]:
#    t=df[i].values.astype(np.float)[:-1]
#    b = [not j for j in np.isnan(t)]
#    idx = [j for j, x in enumerate(b) if x]
#    t = t[idx]
#
#    intensity=df[i+1].values.astype(np.float)[:-1]
#    b = [not j for j in np.isnan(intensity)]
#    idx = [j for j, x in enumerate(b) if x]
#    intensity = intensity[idx]
#
#    popt, pcov = curve_fit(gauss, t, intensity, p0=[t[np.argmax(intensity)], 0.1,50000.])
#
#    plt.plot(t,intensity,t,gauss(t,popt[0],popt[1],popt[2]))
#    plt.xlabel('time [ps]')
#    plt.ylabel('intensity [a.u.]')
#    means = np.append(means,popt[0])
#
#zero_error = (means[0]-means[3])*1000
#print('The delay zero error is {:.4} fs'.format(zero_error))
#mean_stepsize = np.mean(np.diff(means[[1,2,3,4,5]]))*1000
#mean_stepsize_error = np.std(np.diff(means[[1,2,3,4,5]]))*1000
#print('The mean stepsize was {:.5}+-{:.5} fs instead of 200 fs. This is a factor of {:.5}'.format(mean_stepsize,mean_stepsize_error,mean_stepsize/200))