# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:25:19 2019

@author: FemtoMeasure
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

##argon scan:
#runs = [4651, 5020] #first and last run of a delay scan
#run_remove = [4964]

#helium scan:
runs = [602, 1209] #first and last run of a delay scan
run_remove = [226,227,292,293,294,308]

run_list = range(runs[0], runs[1] + 1)
for rr in run_remove:
    if rr in run_list:
        run_list.remove(rr)

delay_zero_pos = 11025.66

sl_spec = [] #seed spectrum of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files

delay = np.array([])
dft = np.array([])
peak_pos =  np.array([])


#which portion of the seed laser spectrum we want to use:
lim = [160, 2020] #upper and lower bound
mid = np.diff(lim)/2

#==============================================================================
# Some functions
#==============================================================================
def Lin(x,a,b):
    return a*x+b
    
def Gauss(x,amp,mean,fwhm,offs):
    return amp*np.exp(-4*np.log(2)*(x-mean)**2/fwhm**2)+offs


#==============================================================================
# Run for center substraction in dft
#==============================================================================
run = 892

ldm_file_path = '//10.5.71.28/FermiServer/Beamtime1/Day_3/Run_' + str(run) +'/rawdata/'.format(int(run))
ldm_files = os.listdir(ldm_file_path)
center_sl = np.zeros((0, np.diff(lim))) #seed laser spectum

for ldm_file in ldm_files:
    ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
  
    #seed laser spectrum
    ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)[::,lim[0]:lim[1]]
    ldm_sl = ldm_sl-np.mean(ldm_sl,axis=1).reshape(300,1) #removing dc offset in seed laser spectrum
    center_sl = np.append(center_sl, ldm_sl,axis=0) 
    ldm_data.close
dft_center = np.mean(abs(np.fft.fftshift(np.fft.fft(center_sl))),axis=0)[mid:]
dft_center /= max(abs(dft_center))
dft_center = np.append(dft_center[:370],np.zeros(len(dft_center)-len(dft_center[:370]))+0.0015)


#==============================================================================
# ''' Loop over all runs '''
#==============================================================================
for run in run_list:
    print run

    if run<100:
        ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_0' + str(run) +'/rawdata/'.format(int(run))
    else:
        ldm_file_path = '//10.5.71.28/FermiServer/Beamtime1/Day_3/Run_' + str(run) +'/rawdata/'.format(int(run))
    ldm_files = os.listdir(ldm_file_path)
    sl = np.zeros((0, np.diff(lim))) #seed laser spectum
    
    try:
        for ldm_file in ldm_files:
            ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
            #delay stage position
            ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
           
            #seed laser spectrum
            ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)[::,lim[0]:lim[1]]
#            plt.plot(ldm_sl[1])
            ldm_sl = ldm_sl-np.mean(ldm_sl,axis=1).reshape(300,1) #removing dc offset in seed laser spectrum

            sl = np.append(sl, ldm_sl,axis=0) 
            ldm_data.close
        dft_sl = np.mean(abs(np.fft.fftshift(np.fft.fft(sl))),axis=0)
#        freq = np.fft.fftshift(np.fft.fftfreq(sl.shape[1],1))
        dft_sl = dft_sl[mid:]        
        dft_sl /= max(abs(dft_sl))
#        plt.plot(dft_sl)
        dft_sl = dft_sl-dft_center # removing center peak of the spectrum
        peak_pos_sl = np.nonzero(dft_sl==np.max(dft_sl))

        #fitting peak:
        try:
            popt_gauss, pcov_gauss = curve_fit(Gauss, range(len(dft_sl)), dft_sl, p0=[1.,peak_pos_sl[0],90.,0.0015], absolute_sigma=False)
#            perr = np.sqrt(np.diag(pcov))
            peak_pos_sl = popt_gauss[1]
        except:
            print('Fit did not converge')
        peak_pos = np.append(peak_pos,peak_pos_sl)
        delay =  np.append(delay,ldm_delay)
#        plt.plot(dft_sl)
    except:
        print('No seed laser spectrum')

delay = delay - delay_zero_pos

cut = 500
popt, pcov = curve_fit(Lin, delay[cut:], peak_pos[cut:], absolute_sigma=False)
perr = np.sqrt(np.diag(pcov))

a = popt[0]
b = popt[1]
sa = perr[0]
sb = perr[1]

print('Zero delay is {} +- {} fs'.format(-b/a, -b/a*(abs(sa/a)+abs(sb/b)) ))

figTD = plt.figure('Fringe Spacing Scan31',figsize=(8, 5))

ax = figTD.add_subplot(111)
ax.set_title('Fringe Spacing Scan31')
ax.plot(delay,peak_pos,'o-')
ax.plot(delay[cut:], peak_pos[cut:],'o-')
ax.plot(delay,Lin(delay,popt[0],popt[1]))
ax.grid()
ax.legend(['all vlaues','fitted values','fit'],loc='lower right')
ax.set_xlabel('delay [fs]')
ax.set_ylabel('fringe spacing [a.u.]')
