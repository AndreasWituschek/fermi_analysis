# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:10:41 2019

@author: ldm
"""

import os, shutil
for run in range(1038, 1134):
    print 'Run_{:03}'.format(run)
    folder = '/home/ldm/ExperimentalData/Online4LDM/RBT-UOF_4/Data/Run_{:03}/work'.format(run)
#    for df in ['/dev3265', '/dev3269', '/niDAQ6321']:
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)