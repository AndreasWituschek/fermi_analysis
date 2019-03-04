# -*- coding: utf-8 -*-
"""
Created on Sun Mar 03 13:17:52 2019

@author: FemtoMeasure
"""

import fermi_analysis.functions as fk
import numpy as np
#import matplotlib.pyplot as plt



timestamp = np.array([])
for i in np.arange(2, 29703, 300):
    mfli_data = fk.ReadMfliData('C:\\Users\\FemtoMeasure\\Documents\\FermiTest\\Setup\\Run_095\\Run_095_' + str(i) + '.h5')
#    print(type(mfli_data))
    timestamp = np.append(timestamp, mfli_data['timestamp'])

print(np.unique(np.diff(timestamp)))

print((timestamp[-1] - timestamp[0])/6e7)
