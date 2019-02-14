# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:22:02 2019

@author: Andreas
"""
import numpy as np
import matplotlib.pyplot as plt

def importFile(path):
    with open(path, 'r') as document:
        data = {}
        header = document.readline().split()
        for i in header:
            data[i]=[]
  
        for line in document:
            line = line.split()
            for idx,i in enumerate(line):
                data[header[idx]].append(float(i))
    return data
    
def curve(l_He,l_ref,h,phi,A,tau):
    X=np.cos(2.*np.pi *299.8*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    Y=np.cos(2.*np.pi *299.8*tau*(h*l_He-l_ref)/(l_ref*l_He) + (phi+90)/180*np.pi)
    return X,Y

def plotCurve(X,Y,plotRange):
    plt.plot(plotRange,X,'b--')
    plt.plot(plotRange,Y, 'r--')
    plt.xlabel('Delay [fs]')
    plt.ylabel('Amplitude [V]')

def plotTdData(data):
    #plot data points
    plt.errorbar(data['Delay'],data['mx'], yerr=data['sx'],color='b',linestyle='')
    plt.plot(data['Delay'],data['mx'],'bo')
    plt.errorbar(data['Delay'],data['my'], yerr=data['sy'],color='r',linestyle='')
    plt.plot(data['Delay'],data['my'],'ro')
    
def plotFdCurve(stop,start,cdft,fDown):
    plt.figure()
    plt.plot([i* 10**3./((stop-start)) for i in range(len(cdft))],abs(cdft),linestyle='--')
    plt.xlim([0,200])
    plt.axvline(x=fDown, color='black', linestyle='--')
    plt.xlabel('Downshifted quantum interference frequency [THz]')
    plt.ylabel('Intensity [a.U.]')

def plotFdCurveEV(stop,start,cdft,fDown,l_ref,l_He,h):
    transition=299.8/l_He*4.135 #helium transition energy
    plt.figure()
    plt.plot([(i* 1/((stop-start))+h*299.8/l_ref)*4.135 for i in range(len(cdft))],abs(cdft),linestyle='--')
    plt.xlim([20,28])
    plt.axvline(x=transition, color='black', linestyle='--')
    plt.text(transition+0.1,0.8*max(abs(cdft)),'He 4P',rotation=0)
    plt.xlabel('Energy [eV]')
    plt.ylabel('Intensity [a.U.]')
 