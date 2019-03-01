# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:22:02 2019

@author: Andreas
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd



#physical constants:
c = 299.792458 #speed of light in nm/fs
hPlanck = 4.135667662 #plancks constant in 10e-15 eV/s


def import_file(path):
    with open(path, 'r') as document:
        data = {}
        header = document.readline().split()
        for i in header:
            data[i] = []
  
        for line in document:
            line = line.split()
            for idx, i in enumerate(line):
                data[header[idx]].append(float(i))
    return data
    

def curve(l_He, l_ref, h, phi, A, offset, start, stop, length):
    """ creates datapoints on cosine/sine curve with given parameters for
        frequency, ampliotude and phase.
    """
    points=30000 # number of TD points for plotting of theroetical curve 
    plotRange = np.linspace(start-10,stop+10,points)
    tau = np.linspace(start,stop,length)
    Xtd = offset+np.cos(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    Ytd = offset+np.sin(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    X = offset+np.cos(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    Y = offset+np.sin(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    return Xtd,Ytd,X,Y
    

def curve_creator(l_He, l_ref, h, phi, A, delay, pNoise, aNoise):
    """ creates datapoints on cosine/sine curve with given parameters for
        frequency, ampliotude and phase. Phase noise and amplitude noise can
        be imparted.
    """
    n=len(delay)
    delPhiX= np.random.rand(n)*2*np.pi*pNoise
    delPhiY= np.random.rand(n)*2*np.pi*pNoise
    X=np.cos(2.*np.pi *c*delay*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi + delPhiX)
    Y=np.sin(2.*np.pi *c*delay*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi + delPhiY)
    X=((np.random.rand(n)-.5)*aNoise+1)*X
    Y=((np.random.rand(n)-.5)*aNoise+1)*Y
    return X,Y


def plot_curve(X, Y, start, stop):
    """ Plots the theoretical curve of the time domain scan
    """
    points=30000 # number of TD points for plotting of theroetical curve 
    plotRange = np.linspace(start-10,stop+10,points)
    fig, ax = plt.subplots(1,1)
    ax.plot(plotRange,X,'b--', label="theoretical curve demodX")
    ax.plot(plotRange,Y, 'r--', label="theoretical curve demodY")
    ax.set_xlabel('Delay [fs]')
    ax.set_ylabel('Amplitude [V]')
    ax.set_title("Downshifted quantum interferences (TD)")
    ax.legend()
    return fig, ax


def plot_td_data(data, demod):
    """ Plots the acquiered time domain data
    """
    #plot data points
    fig, ax = plt.subplots(1,1)
    ax.errorbar(data['delay'],data['mX%d' % demod], yerr=data['sX%d' % demod],color='b',linestyle='')
    ax.plot(data['delay'],data['mX%d' % demod],'bo')
    ax.errorbar(data['delay'],data['mY%d' % demod], yerr=data['sY%d' % demod],color='r',linestyle='')
    ax.plot(data['delay'],data['mY%d' % demod],'ro')
    return fig, ax
    
    
def plot_fd_curve(stop, start, cdft, cdft_d, l_ref, l_He, h):
    """ plots the absoprtion spectrum of theoretical curve and data, x- axis
        is downshifted frequency.
    """
#    plt.figure()
    fig, ax = plt.subplots(1,1)
    fDown = abs(c*1000*(h*l_He-l_ref)/(l_ref*l_He)) # downshifted frequency of demodulated signal
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data   
    xaxis=[i* 10**3./(stop-start) for i in range(len(cdft))] #x-axis conversion in THz
    ax.plot(xaxis,abs(cdft),linestyle='-', label="theoretical curve")
    ax.plot(xaxis,abs(cdft_d),linestyle='-', label="experimental data")
    ax.set_xlim([0,200])
    ax.axvline(x=fDown, color='black', linestyle='--')
    ax.set_xlabel('Downshifted quantum interference frequency [THz]')
    ax.set_ylabel('Intensity [a.U.]')
    ax.legend()
    ax.set_title("Downshifted quantum interference")
    return fig, ax


def plot_fd_curve_abs_ev(stop, start, cdft, cdft_d, l_ref, l_He, h):
    """ Plots the absoprtion spectrum of the theoretical curve and data.
    """
#    plt.figure()
    fig, ax = plt.subplots(1,1)
    transition=c/l_He*hPlanck #helium transition energy
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data    
    xaxis=[(i* (1./(stop-start)*hPlanck)+h*c/l_ref*hPlanck) for i in range(len(cdft))] # x-axis, convert from downshifted frequency axis to eV
    ax.plot(xaxis,abs(cdft),linestyle='-', label="theoretical curve")
    ax.plot(xaxis,abs(cdft_d),linestyle='-', label="experimental data")
    ax.set_xlim([23,25])
    ax.axvline(x=transition, color='black', linestyle='--')
    ax.text(transition+0.02,0.8*max(abs(cdft)),'He 4P') #position of He4P transition
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Intensity [a.U.]')
    ax.legend()
    ax.set_title("Absorption spectrum")
    return fig, ax


def plot_fd_curve_ev(stop, start, cdft, cdft_d, l_ref, l_He, h):
    """ Plots the absorptive and dispersive part of the spectrum of theoretical
        curve and data
    """
#    plt.figure()
    fig, ax = plt.subplots(1,1)
    transition=c/l_He*hPlanck #helium transition energy
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data    
    xaxis=[(i* (1./(stop-start)*hPlanck)+h*c/l_ref*hPlanck) for i in range(len(cdft))] # x-axis, convert from downshifted frequency axis to eV
    ax.plot(xaxis,cdft.real,linestyle='-', label="theoretical curve absorptive")
    ax.plot(xaxis,cdft_d.real,linestyle='-', label="experimental data absorptive")
    ax.plot(xaxis,cdft.imag,linestyle='-', label="theoretical curve dispersive")
    ax.plot(xaxis,cdft_d.imag,linestyle='-', label="experimental data dispersive")
    ax.set_xlim([23,25])
    ax.axvline(x=transition, color='black', linestyle='--')
    ax.text(transition+0.02,0.8*max(abs(cdft)),'He 4P') #position of He4P transition
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Intensity [a.U.]')
    ax.legend()
    ax.set_title("Absorptive and dispersive spectra")
    return fig, ax


def read_mfli_mata(mfli_file_path, bunches): #bunches= number of shots in one file
    """ Reads data from files created with mfli.py containing the raw
        demodulated data. ONLY DUMMY TO CREATE SOME DATA RIGHT NOW.
    """
    X = [np.random.rand(bunches) for i in np.arange(0,4)] #demodulated X data for each bunch
    Y = [np.random.rand(bunches) for i in np.arange(0,4)] #demodulated Y data for each bunch
    modfreq= 1000+np.random.rand(bunches) # modulation frequency for each bunch in Hz
    run= 1 # run number
    timestamp = 10 # timestamp of first bunch in this run
    bunchnumber = 3 #bunchnumber of first bunch in this run
    mfli_data = {'run': run, 
                 'bunchnumber': bunchnumber,
                 'timestamp': timestamp, 
                 'modfreq': modfreq,
                 'X': X, 
                 'Y': Y,} #write all into one dictionary
    return mfli_data
    

def averaging_mfli_data(mfli_data, I0, apply_filter, I0_threshold, modfreq_set, modfreq_threshold, ignore):
    """ Applys filters (IO and modfreq) on mfli data and averages filtered
        values.
    """
    if apply_filter:
        I0_filter = np.asarray([i>I0_threshold for i in I0])*1 #I0 filter
        I0_filter = np.concatenate([I0_filter,np.ones(ignore)])
        a=np.zeros((ignore,len(I0_filter)))
        for i in range(ignore):
            a[i] = np.roll(I0_filter,i)
        I0_filter = (np.sum(a,0)>(ignore-1))[ignore-1:len(I0)+(ignore-1)] #I0 filter extended by number of ignored shots
        
        modfreq_filter = np.asarray([-modfreq_threshold < i-modfreq_set < modfreq_threshold for i in mfli_data['modfreq']]) #filter for modfreq
        modfreq_filter = np.concatenate([modfreq_filter,np.ones(ignore)])
        a=np.zeros((ignore,len(modfreq_filter)))
        for i in range(ignore):
            a[i] = np.roll(modfreq_filter,i)
        modfreq_filter=(np.sum(a,0)>(ignore-1))[ignore-1:len(I0)+ignore-1] #modfreq filter extended by number of ignored shots

        b = np.logical_and(I0_filter,modfreq_filter) #combining both filters above
                
        if sum(b)>0:
            print("# total bunches = %d" % len(I0))
            print("# filtered Bunches = %d" % sum(b))
            mX = [np.mean(mfli_data['X'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            mY = [np.mean(mfli_data['Y'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
            sX = [np.std(mfli_data['X'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            sY = [np.std(mfli_data['Y'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
        else:
            print('Filters to tight, no data to average over. Script stopped, no data written in masterData.csv!')
            sys.exit(1)
    else: 
        mX = np.mean((mfli_data['X']),1)
        mY = np.mean((mfli_data['Y']),1)
        sX = np.std((mfli_data['X']),1)
        sY = np.std((mfli_data['Y']),1)
        print("# averaged bunches = %d" % len(I0))
    return mX, mY, sX, sY, I0_filter, modfreq_filter
    

def master_file_writer(mfli_data, delay, mX, mY, sX, sY):
    """ Writes the preanalyzed time domain data to the so called master file.
    """
    d = {'run': [mfli_data['run']],
         'delay': [delay],
         'mY0': mY[0],'mX0': mX[0]
         'mY1': mY[1],'mX1': mX[1],
         'mY2': mY[2],'mX2': mX[2],
         'mY3': mY[3],'mX3': sX[3],
         'sY0': sY[0],'sX0': sX[0],
         'sY1': sY[1],'sX1': sX[1],
         'sY2': sY[2],'sX2': sX[2],
         'sY3': sY[3],'sX3': sX[3],
         'bunchnumber': [mfli_data['bunchnumber']],
         'timestamp' : [mfli_data['timestamp']]}
    df = pd.DataFrame.from_dict(d)    
    return df
    

def master_file_reader(path): 
    """ Reads the so called master file, which contains the preanalyzed data.
    """
    dataDF = pd.read_csv(path, index_col=0) #Import file with data from MFLI:   
    data = dataDF.to_dict("list") #Convert the DataFrame to the dictionary format:
    return data


def delay_move(delay_pos):
    """ What ever.
    """
    if all(x == delay_pos[0] for x in delay_pos):
        delay_pos = delay_pos[0]
    else: print "Delaystage moved during run!"
