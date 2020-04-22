# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 19:11:54 2019

@author: Andreas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

n = 1.6226

steps = np.array([1,43.5634])
delay = steps/const.c*np.sin(7./180*np.pi)*(n-1)*1e12
print(delay)