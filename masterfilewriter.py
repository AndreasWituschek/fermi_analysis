# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:14:24 2019

@author: Andreas

Reads hd5 files from LDM data acquisition system and hd5 files from mfli.py program containing the demodulated data from the MFLI. 

"""

import h5py
import numpy as np
import fermi_analysis.functions as fk
import pandas as pd


''' FILE I/O '''
#importing LDM data:
file_path='//nanoserver/user/Projekte/BMBF/FERMI Beispieldaten/Nh3_044_72274751.h5'
ldm_data = h5py.File(file_path, 'r')
I0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a']) # I0 in uJ
delay = np.array(ldm_data['user_laser']['delay_line']['position']) #position of delay line
bunches= len(I0) # number of bunches per file

#importing MFLI data:
mfli_data = fk.ReadMfliData(bunches)

''' Data filtering '''
apply_filter = False #choose wether to apply following filters:
I0_threshold= np.mean(I0)-2*np.std(I0) #lower threshold for I0 power filtering
modfreq_set = 1000. # modulation frequency in Hz set by the user at the AD9959 board
modfreq_threshold= 0.9 #disregards data when not within +- modfreq_threshold from modfreq_set because then MFLI was not properly locked


''' processing data '''
if all(x == delay[0] for x in delay):
    delay=delay[0]
else: print "Delaystage moved during run!"

#average values and standard deviation of filtered demodulated values
mX,mY,sX,sY=fk.AveragingMfliData(mfli_data,I0,apply_filter,I0_threshold,modfreq_set,modfreq_threshold)

''' Writing to masterfile'''
#write all data in dictionary:
d = {'run': [mfli_data['run']], 
     'delay': [delay],
     'mX0': mX[0], 
     'mY0': mY[0],
     'mX1': mX[1],
     'mY1': mY[1],
     'mX2': mX[2],
     'mY2': mY[2],
     'mX3': mX[3], 
     'mY3': mY[3],     
     'sX0': sX[0], 
     'sY0': sY[0],
     'sX1': sX[1],
     'sY1': sY[1],
     'sX2': sX[2],
     'sY2': sY[2],
     'sX3': sX[3], 
     'sY3': sY[3],
     'bunchnumber': [mfli_data['bunchnumber']],
     'timestamp' : [mfli_data['timestamp']]}


df = pd.DataFrame.from_dict(d) #converting dict to data frame

##Create a CSV file:
#df.to_csv("masterData.csv")

#appending data to masterData.csv File
with open('masterData.csv', 'a') as f:
    df.to_csv(f, header=False)


