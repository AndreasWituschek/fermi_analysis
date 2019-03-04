# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 18:17:27 2019

@author: FemtoMeasure
"""
import fermi_analysis.functions as fk
import numpy as np
import matplotlib.pyplot as plt



mfli_data = fk.ReadMfliData('C:\\Users\\FemtoMeasure\\Documents\\FermiTest\\Setup\\Run_097\\Run_097_2.h5')
d = mfli_data['trigger']
d = np.append(0., d)
d = np.flatnonzero((d[:-1] == 0) & (d[1:] == 2))
print(np.unique(np.diff(d)))
print(mfli_data.keys())
print(mfli_data['harmonic'])
print(mfli_data['oscselect'])
print(mfli_data['timeconstant'])
print(mfli_data['order'])

print((mfli_data['timestamp'][-1] - mfli_data['timestamp'][0])/6e7)
fig, ax = plt.subplots(3,1)
ax[0].plot(mfli_data['timestamp'], mfli_data['x0'])
ax[0].plot(mfli_data['timestamp'], mfli_data['y0'])
ax[1].plot(mfli_data['timestamp'], mfli_data['trigger'])
ax[1].set_ylim((-0.1, 2.1))
ax[2].plot(d)
for a in ax:
    a.grid()
plt.show()