# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:14:24 2019

@author: Andreas

Reads hd5 files from LDM data acquisition system and hd5 files from mfli.py program containing the demodulated data from the MFLI. 

"""

import h5py
import numpy as np
import fermi_analysis.functions as fk

'''Experimental parameters '''
ri=1.6116 #refractive index of wedges at seed wavelenght
delay_zero_pos= -19.8 #position of delay stage in mm where we 0 delay between seed pulses
modfreq_set = 1000. # modulation frequency in Hz set by the user at the AD9959 board

''' FILE I/O '''
#importing LDM data:
ldm_file_path='//nanoserver/user/Projekte/BMBF/FERMI Beispieldaten/Nh3_044_72274751.h5'
mfli_file_path = "//nanoserver/user/Projekte/BMBF/FERMI Beispieldaten/test.h5"

ldm_data = h5py.File(ldm_file_path, 'r')

I0 = np.array(ldm_data['photon_diagnostics']['FEL01']['I0_monitor']['iom_sh_a']) # I0 in uJ
delay_pos = np.array(ldm_data['user_laser']['delay_line']['position']) #position of delay stage
bunches= len(I0) # number of bunches per file

#importing MFLI data:
mfli_data = fk.ReadMfliData(mfli_file_path,bunches)

''' Data filtering '''
apply_filter = True #choose wether to apply following filters:
ignore=5 #number of bunches to ignore after filter incident (depends on MFLI time constant)
I0_threshold= np.mean(I0)-2*np.std(I0) #lower threshold for I0 power filtering
modfreq_threshold= 2.9 #[Hz] Disregards data when not within modfreq_set +- modfreq_threshold because then MFLI was not properly locked

fk.DelayMove(delay_pos) #check if delaystage moved during run.

''' processing data '''
#average values and standard deviation of filtered demodulated values
mX,mY,sX,sY,I0_filter,modfreq_filter=fk.AveragingMfliData(mfli_data,I0,apply_filter,I0_threshold,modfreq_set,modfreq_threshold,ignore)
delay= (delay_pos[0] - delay_zero_pos)*(ri-1)*np.sin(7./180*np.pi)/0.000299792458 #convert position [mm] of delay line to fs
delay=-45

#delay=4

''' Writing to masterfile'''
#write all data to one dataframe:
df=fk.MasterFileWriter(mfli_data,delay,mX,mY,sX,sY)

###Create a CSV file. Do this when running this program for the first time. Otherwise this will overwrite existing masterData.csv file. 
#df.to_csv("masterData.csv")

#appending data to masterData.csv File
with open('masterData.csv', 'a') as f:
    df.to_csv(f, header=False)