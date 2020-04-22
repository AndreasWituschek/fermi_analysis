# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 10:51:35 2019

@author: FemtoMeasure
"""
import glob
import os
import shutil
import winsound

import pyttsx

def recursive_copy_files(source_path, destination_path, override=False):
    """ Function to copy the data file to the fermi server.
    """
    files_count = 0
    if not os.path.exists(destination_path):
        os.mkdir(destination_path)
    items = glob.glob(source_path + '\\*')
    for item in items:
        if os.path.isdir(item):
            path = os.path.join(destination_path, item.split('\\')[-1])
            files_count += recursive_copy_files(source_path=item,
                                                destination_path=path,
                                                override=override)
        else:
            f = os.path.join(destination_path, item.split('\\')[-1])
            if not os.path.exists(f) or override:
                shutil.copyfile(item, f)
                files_count += 1
    return files_count


def recursive_copy_folder(source_path, destination_path, dev, override=False):
    """ Function to copy the data file to the fermi server.
    """
    files_count = 0
    if not os.path.exists(destination_path):
        os.mkdir(destination_path)
    items = glob.glob(source_path + '\\*')
    for item in items:
        if os.path.isdir(item):
            path = os.path.join(destination_path, item.split('\\')[-1])
            files_count += recursive_copy_files(source_path=item,
                                                destination_path=path,
                                                override=override)
        else:
            f = os.path.join(destination_path, item.split('\\')[-1])
            if not os.path.exists(f) or override:
                shutil.copyfile(item, f)
                files_count += 1
    return files_count


# \\online4ldm.esce.elettra.trieste.it\store\20149020\Setup

if __name__ == '__main__':
    start_run = 694
    end_run = 707
#    engine = pyttsx.init()
#    engine.say('I am starting to shift the files to the other computer.')
#    engine.runAndWait()   
    experiment = 'Data'
    for run in range(start_run, end_run + 1):
        print('Run_{:03}\\'.format(run))
#        root_local = 'D:\\beamtime2\\'
        root_local = '\\10.5.71.28\mfli\Beamtime2\\'
#        root_server = '\\\\online4ldm.esce.elettra.trieste.it\\store\\RBT-UOF_4\\'
        root_server = '\\10.5.71.28\FermiServer\Beamtime2\\'
        for dev in ['dev3265', 'dev3269', 'niDAQ6321']: # RBT-UOF_4
            local_path = root_local + experiment + '\\Run_{0:03}\\{1}\\'.format(run, dev)
            server_path = root_server + experiment + '\\Run_{0:03}\\work\\{1}\\'.format(run, dev)
#            print(local_path)
#            print(server_path)
            recursive_copy_folder(local_path, server_path, dev, override=False)
    execfile('C:\\Users\\FemtoMeasure\\Documents\\Python\\fermi_analysis\\c_file.py')
    winsound.Beep(440, 1000)
#    engine.say('Finished! Juhu. Analyze the data!')
#    engine.runAndWait()

