# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:22:02 2019

@author: Andreas
"""
import numpy as np
import matplotlib.pyplot as plt

def hello():
    print('hello')

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

def plotData(data):
    #plot data points
    plt.errorbar(data['Delay'],data['mx'], yerr=data['sx'],color='r',linestyle='')
    plt.plot(data['Delay'],data['mx'],'bo')
    plt.errorbar(data['Delay'],data['my'], yerr=data['sy'],color='b',linestyle='')
    plt.plot(data['Delay'],data['my'],'ro')
    
