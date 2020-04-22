# -*- coding: utf-8 -*-
"""
Created on Fri May 31 10:07:34 2019

@author: LukasB
"""
import numpy as np
import fermi_analysis.functions as fk 
import scipy.constants as spc

factor = spc.h*spc.c/spc.e*100.
T = np.linspace(0.0, 550.0, 200)
peakMargin = 0.05
suscept = False
zeroPaddingFactor = 1
FWHM, Full = fk.PeakWidth(T, peakMargin, suscept, zeroPaddingFactor)
print FWHM, Full
FWHM *= factor
print FWHM