# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 11:09:30 2019

@author: Andreas
"""


import pandas as pd
import numpy as np
import fermi_analysis.functions as fk


#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266. #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)


#scan parameters
step=4 #delay stepsize in fs
start=70 #delay scan start in fs
stop=1  #delay scans stop in fs
delay=np.arange(start,stop,step) #delay points in fs
n=len(delay) #number of time domain data points


#creating theoretical curve and adding phase and amplitude noise:
aNoise=0.2*A #amplitude noise in units of the "correct" amplitude
pNoise= 20. #phase noise in degree, both for X and Y
pNoise= pNoise/180*np.pi
mx,my= fk.CurveCreator(l_He,l_ref,h,phi,A,delay, pNoise, aNoise)

# incrementing run number at each delay step
run=np.array(range(1,n+1)) 
#error bars
s=0.05 #5% error on all data points
sx, sy = np.full(n,s), np.full(n,s)


#write all data in dictionary:
d = {'run': run, 'delay': delay,'mX': mx,'mY': my, 'sX': sx, 'sY': sy}
df = pd.DataFrame.from_dict(d)
#Create a CSV file:
df.to_csv("exampleData.csv")