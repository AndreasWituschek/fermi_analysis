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
import pandas as pd
import fermi_analysis.functions as fk

#all scans of whole beamtime:
#runs = [10, 6749] #first and last run of a delay scan

#argon scan:
runs = [4651, 5021] #first and last run of a delay scan
#runs = [4651, 4653] #first and last run of a delay scan

#helium scan:
#runs = [31, 690] #first and last run of a delay scan

run_list = range(runs[0], runs[1] + 1)

delay_zero_pos = 11025.66

sl_spec = [] #seed spectrum of all runs
ldm_sl= np.array([]) #raw seed laser spectrum from ldm files
delay = np.array([]) #pump probe delay [fs]
#dft = np.array([])
peak_pos =  np.array([]) #peak position of the fringe spacing peak in the dft [a.u.]
runs_analyzed = np.array([]) #list of runs that actually have been analysed (and didnt have problems like missing spectrum etc.)

#which portion of the seed laser spectrum we want to use:
lim = [160, 2032] #upper and lower bound
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
run = 5000

ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_' + str(run) +'/rawdata/'.format(int(run))
ldm_files = os.listdir(ldm_file_path)
center_sl = np.zeros((0, np.diff(lim))) #seed laser spectum

for ldm_file in ldm_files:
    ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
  
    #seed laser spectrum
    ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)[::,lim[0]:lim[1]]
    ldm_sl = ldm_sl-np.mean(ldm_sl,axis=1).reshape(400,1) #removing dc offset in seed laser spectrum
    center_sl = np.append(center_sl, ldm_sl,axis=0) 
    ldm_data.close
dft_center = np.mean(abs(np.fft.fftshift(np.fft.fft(center_sl))),axis=0)[mid:]
dft_center /= max(abs(dft_center))
dft_center = np.append(dft_center[:200],np.zeros(len(dft_center)-len(dft_center[:200]))+0.0015)


#==============================================================================
# ''' Loop over all runs '''
#==============================================================================
for run in run_list:
    print run
    try:
        if run<100:
            ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_0' + str(run) +'/rawdata/'.format(int(run))
        else:
            ldm_file_path = '//10.5.71.28/FermiServer/Beamtime2/Run_' + str(run) +'/rawdata/'.format(int(run))
        ldm_files = os.listdir(ldm_file_path)
        sl = np.zeros((0, np.diff(lim))) #seed laser spectum
        
        try:
            for ldm_file in ldm_files:
                ldm_data = h5py.File(ldm_file_path + ldm_file, 'r')
                #delay stage position
                ldm_delay = np.array(ldm_data['photon_source']['SeedLaser']['trls_sl_03_pos'])
               
                #seed laser spectrum
                ldm_sl = np.array(ldm_data['photon_source']['SeedLaserSpectrum_FEL01']['Spectrum'],dtype=np.float32)[::,lim[0]:lim[1]]
                ldm_sl = ldm_sl-np.mean(ldm_sl,axis=1).reshape(400,1) #removing dc offset in seed laser spectrum
        
                sl = np.append(sl, ldm_sl,axis=0) 
                ldm_data.close
                
            peak_positions = np.array([])
            for i in range(6):
                dft_sl = abs(np.fft.fftshift(np.fft.fft(sl)))[i,::]
            #        freq = np.fft.fftshift(np.fft.fftfreq(sl.shape[1],1))
                dft_sl = dft_sl[mid:]        
                dft_sl /= max(abs(dft_sl))
            #        plt.plot(dft_sl)
                dft_sl = dft_sl-dft_center # removing center peak of the spectrum
                peak_pos_sl = np.nonzero(dft_sl==np.max(dft_sl))
                peak_positions = np.append(peak_positions,peak_pos_sl)
            #fitting peak:
#            try:
#                popt_gauss, pcov_gauss = curve_fit(Gauss, range(len(dft_sl)), dft_sl, p0=[1.,peak_pos_sl[0],90.,0.0015], absolute_sigma=False)
#        #            perr = np.sqrt(np.diag(pcov))
#                peak_pos_sl = popt_gauss[1]
#            except:
#                print('Fit did not converge')
            peak_pos = np.append(peak_pos,peak_positions)
            delay =  np.append(delay,ldm_delay)
        #        plt.plot(dft_sl)
            runs_analyzed = np.append(runs_analyzed, run)
        except:
            print('No seed laser spectrum')
    except: 
        print ('File not found')

delay = delay - delay_zero_pos
peak_pos = peak_pos.reshape((len(delay),len(peak_positions)))

#==============================================================================
# Fitting the linar part
#==============================================================================
zero_delay = np.array([])
#iterator =  np.arange(160,len(delay)-50,1) # fpr Scan31
iterator =  np.arange(100,len(delay)-50,1) #for scan4651
delay_range =  np.array([])
for cut_delay in iterator:
#    cut_delay = -270 #value where peak_pos still depends linearly on delay
#    idx = np.nonzero(delay>cut_delay)
    
    popt, pcov = curve_fit(Lin, delay[cut_delay:cut_delay+50], np.min(peak_pos,axis=1)[cut_delay:cut_delay+50], absolute_sigma=False)
    perr = np.sqrt(np.diag(pcov))
    
    a = popt[0]
    b = popt[1]
    sa = perr[0]
    sb = perr[1]
    
#    print('Zero delay is {} +- {} fs'.format(-b/a, -b/a*(abs(sa/a)+abs(sb/b)) ))
    delay_range = np.append(delay_range,delay[cut_delay])
    zero_delay = np.append(zero_delay,-b/a)
    
figTD = plt.figure('zero_delay_vs_fit_start_point',figsize=(8, 5))
ax = figTD.add_subplot(111)
ax.set_title('zero_delay_vs_fit_start_point')   
ax.plot(delay_range,zero_delay,'o')
ax.set_xlabel('starting point of fit [fs]')
ax.set_ylabel('zero_delay [fs]')


fit_range = [-200,-370] #for argon
mean_zero =  np.mean(zero_delay[fk.find_index(delay_range,fit_range[0]):fk.find_index(delay_range,fit_range[1])])
std_zero =  np.std(zero_delay[fk.find_index(delay_range,fit_range[0]):fk.find_index(delay_range,fit_range[1])])
print('Zero delay is {} +- {} fs'.format(mean_zero,std_zero))


#==============================================================================
# Plotting the results
#==============================================================================
figTD = plt.figure('Fringe Spacing Scan4651',figsize=(8, 5))
ax = figTD.add_subplot(111)
ax.set_title('Fringe Spacing Scan4651')
ax.plot(delay,peak_pos,'o')
#ax.plot(delay[idx], peak_pos[idx],'o')
ax.plot(delay,Lin(delay,popt[0],popt[1]))
ax.grid()
ax.legend(['all values','fitted values','fit'],loc='lower left')
ax.set_xlabel('delay [fs]')
ax.set_ylabel('fringe spacing [a.u.]')


##==============================================================================
## Exporting results to csv
##==============================================================================
#
#d = {'run': runs_analyzed,\
#    'delay': delay,\
#    'peak_pos': peak_pos
#    }
#df = pd.DataFrame.from_dict(d)
#df.to_csv("//nanoserver/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/zero_delay_calibration/test.csv")
#
##==============================================================================
## Importin results from csv
##==============================================================================
##df = pd.DataFrame.from_csv("//nanoserver/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/zero_delay_calibration/Scan4651/delay_vs_peak_pos_scan4651.csv")
#df = pd.DataFrame.from_csv("//nanoserver/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/zero_delay_calibration/test.csv")
#delay = np.array(df['delay'])
#peak_pos = np.array(df['peak_pos'])

