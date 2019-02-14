# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 11:09:30 2019

@author: Andreas
"""


import pandas as pd
import numpy as np
import fermi_analysis.functions as fk
import matplotlib.pyplot as plt


#parameters for theoretical curve
l_He = 52.22 #helium 1s-4p transition in nm
l_ref =266. #reference laser wavelenght in nm
h = 5. #harmonic of FEL
phi = 0. #phase offset between reference and signal in degrees
A = 1. #amplitude (=R from MFLI)


#scan parameters
step=3 #delay stepsize in fs
start=70 #delay scan start in fs
stop= 570 #delay scans stop in fs
delay=np.arange(start,stop,step) #delay points in fs
n=len(delay) #number of time domain data points
run=np.array(range(1,n+1)) # incrementing run number at each delay step

#ideal data:
mx,my= fk.curve(l_He,l_ref,h,phi,A,delay)
#adding 50% noise on data
for m in [mx,my]:
    m = (np.random.rand(n)*1.5+1)*m

#error bars
s=0.05 #5% error on all data points
sx, sy = np.full(n,s), np.full(n,s)

#background runs:
b=np.round(np.random.rand(n),0)

#write all data in dictionary:
d = {'Run': run, 'Delay': delay,'mx': mx,'my': my, 'sx': sx, 'sy': sy, 'Background': b}
df = pd.DataFrame.from_dict(d)
#Create a CSV file:
df.to_csv("exampleData.csv")

##Read the CSV file back as a DataFrame:
#
#df = pd.read_csv("exampleData.csv", index_col=0)
#
##Convert the DataFrame back to the original dictionary format:
#
#data = df.to_dict("list")
#
#
#for key in data:
#    print key
#    
plt.plot(data['mx'])