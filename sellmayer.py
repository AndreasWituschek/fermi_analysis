# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 12:26:02 2019

@author: Andreas
"""

import numpy as np

import numpy as np
import matplotlib.pyplot as plt

lam=0.261

B1=6.694225753e-1
B2=4.34583937e-1
B3=8.71694723e-1
C1=4.48011239e-3
C2=1.32847049e-2
C3=9.53414824e1

nsq = np.sqrt(1+B1*lam**2/(lam**2-C1) + B2*lam**2/(lam**2-C2) + B3*lam**2/(lam**2-C3))
print(nsq)
