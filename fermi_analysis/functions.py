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
c=299.792458 #speed of light in nm/fs
hPlanck=4.135667662 #plancks constant in 10e-15 eV/s

def importFile(path):
    with open(path, 'r') as document:
        data = {}
        header = document.readline().split()
        for i in header:
            data[i]=[]
  
        for line in document:
            line = line.split()
            for idx,i in enumerate(line):
                data[header[idx]].append(float(i))
    return data
    
def Curve(l_He,l_ref,h,phi,A,offset, start, stop, length):
    "creates datapoints on cosine/sine curve with given parameters for frequency, ampliotude and phase."
    points=30000 # number of TD points for plotting of theroetical curve 
    plotRange = np.linspace(start-10,stop+10,points)
    tau = np.linspace(start,stop,length)
    Xtd=offset+np.cos(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    Ytd=offset+np.cos(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + (phi+90)/180*np.pi)
    X=offset+np.cos(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi)
    Y=offset+np.cos(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + (phi+90)/180*np.pi)
    return Xtd,Ytd,X,Y
    
def CurveCreator(l_He,l_ref,h,phi,A, delay,pNoise,aNoise):
    "creates datapoints on cosine/sine curve with given parameters for frequency, ampliotude and phase. Phase noise and amplitude noise can be imparted."
    n=len(delay)
    delPhiX= np.random.rand(n)*2*np.pi*pNoise
    delPhiY= np.random.rand(n)*2*np.pi*pNoise
    X=np.cos(2.*np.pi *c*delay*(h*l_He-l_ref)/(l_ref*l_He) + phi/180*np.pi + delPhiX)
    Y=np.cos(2.*np.pi *c*delay*(h*l_He-l_ref)/(l_ref*l_He) + (phi+90)/180*np.pi + delPhiY)
    X=((np.random.rand(n)-.5)*aNoise+1)*X
    Y=((np.random.rand(n)-.5)*aNoise+1)*Y
    return X,Y

def PlotCurve(X,Y,start, stop):
    points=30000 # number of TD points for plotting of theroetical curve 
    plotRange = np.linspace(start-10,stop+10,points)
    plt.plot(plotRange,X,'b--', label="theoretical curve demodX")
    plt.plot(plotRange,Y, 'r--', label="theoretical curve demodY")
    plt.xlabel('Delay [fs]')
    plt.ylabel('Amplitude [V]')
    plt.title("Downshifted quantum interferences (TD)")
    plt.legend()

def PlotTdData(data,demod):
    #plot data points
    plt.errorbar(data['delay'],data['mX%d' % demod], yerr=data['sX%d' % demod],color='b',linestyle='')
    plt.plot(data['delay'],data['mX%d' % demod],'bo')
    plt.errorbar(data['delay'],data['mY%d' % demod], yerr=data['sY%d' % demod],color='r',linestyle='')
    plt.plot(data['delay'],data['mY%d' % demod],'ro')
    
    
def PlotFdCurve(stop,start,cdft,cdft_d,l_ref,l_He,h):
    "plots the absoprtion spectrum of theoretical curve and data, x- axis is downshifted frequency"
    plt.figure()
    fDown = abs(c*1000*(h*l_He-l_ref)/(l_ref*l_He)) # downshifted frequency of demodulated signal
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data   
    xaxis=[i* 10**3./(stop-start) for i in range(len(cdft))] #x-axis conversion in THz
    plt.plot(xaxis,abs(cdft),linestyle='-', label="theoretical curve")
    plt.plot(xaxis,abs(cdft_d),linestyle='-', label="experimental data")
    plt.xlim([0,200])
    plt.axvline(x=fDown, color='black', linestyle='--')
    plt.xlabel('Downshifted quantum interference frequency [THz]')
    plt.ylabel('Intensity [a.U.]')
    plt.legend()
    plt.title("Downshifted quantum interference")


def PlotFdCurveAbsEV(stop,start,cdft,cdft_d,l_ref,l_He,h):
    "plots the absoprtion spectrum of theoretical curve and data"
    plt.figure()
    transition=c/l_He*hPlanck #helium transition energy
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data    
    xaxis=[(i* (1./(stop-start)*hPlanck)+h*c/l_ref*hPlanck) for i in range(len(cdft))] # x-axis, convert from downshifted frequency axis to eV
    plt.plot(xaxis,abs(cdft),linestyle='-', label="theoretical curve")
    plt.plot(xaxis,abs(cdft_d),linestyle='-', label="experimental data")
    plt.xlim([23,25])
    plt.axvline(x=transition, color='black', linestyle='--')
    plt.text(transition+0.02,0.8*max(abs(cdft)),'He 4P') #position of He4P transition
    plt.xlabel('Energy [eV]')
    plt.ylabel('Intensity [a.U.]')
    plt.legend()
    plt.title("Absorption spectrum")


def PlotFdCurveEV(stop,start,cdft,cdft_d,l_ref,l_He,h):
    "plots the absorptive and dispersive part of the spectrum of theoretical curve and data"
    plt.figure()
    transition=c/l_He*hPlanck #helium transition energy
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data    
    xaxis=[(i* (1./(stop-start)*hPlanck)+h*c/l_ref*hPlanck) for i in range(len(cdft))] # x-axis, convert from downshifted frequency axis to eV
    plt.plot(xaxis,cdft.real,linestyle='-', label="theoretical curve absorptive")
    plt.plot(xaxis,cdft_d.real,linestyle='-', label="experimental data absorptive")
    plt.plot(xaxis,cdft.imag,linestyle='-', label="theoretical curve dispersive")
    plt.plot(xaxis,cdft_d.imag,linestyle='-', label="experimental data dispersive")
    plt.xlim([23,25])
    plt.axvline(x=transition, color='black', linestyle='--')
    plt.text(transition+0.02,0.8*max(abs(cdft)),'He 4P') #position of He4P transition
    plt.xlabel('Energy [eV]')
    plt.ylabel('Intensity [a.U.]')
    plt.legend()
    plt.title("Absorptive and dispersive spectra")

def ReadMfliData(mfli_file_path,bunches): #bunches= number of shots in one file
    "reads data from files created with mfli.py containing the raw demodulated data. ONLY DUMMY TO CREATE SOME DATA RIGHT NOW"
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
    
def AveragingMfliData(mfli_data,I0,apply_filter,I0_threshold,modfreq_set,modfreq_threshold,ignore):
    "Applys filters (IO and modfreq) on mfli data and averages filtered values"
    if apply_filter:
        I0_filter=np.asarray([i>I0_threshold for i in I0])*1 #I0 filter
        I0_filter=np.concatenate([I0_filter,np.ones(ignore)])
        a=np.zeros((ignore,len(I0_filter)))
        for i in range(ignore):
            a[i]=np.roll(I0_filter,i)
        I0_filter=(np.sum(a,0)>(ignore-1))[ignore-1:len(I0)+(ignore-1)] #I0 filter extended by number of ignored shots
        
        modfreq_filter= np.asarray([-modfreq_threshold < i-modfreq_set < modfreq_threshold for i in mfli_data['modfreq']]) #filter for modfreq
        modfreq_filter=np.concatenate([modfreq_filter,np.ones(ignore)])
        a=np.zeros((ignore,len(modfreq_filter)))
        for i in range(ignore):
            a[i]=np.roll(modfreq_filter,i)
        modfreq_filter=(np.sum(a,0)>(ignore-1))[ignore-1:len(I0)+ignore-1] #modfreq filter extended by number of ignored shots

        b = np.logical_and(I0_filter,modfreq_filter) #combining both filters above
                
        if sum(b)>0:
            print("# total bunches = %d" % len(I0))
            print("# filtered Bunches = %d" % sum(b))
            mX= [np.mean(mfli_data['X'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            mY= [np.mean(mfli_data['Y'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
            sX= [np.std(mfli_data['X'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            sY= [np.std(mfli_data['Y'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
        else:
            print('Filters to tight, no data to average over. Script stopped, no data written in masterData.csv!')
            sys.exit(1)
    else: 
        mX= np.mean((mfli_data['X']),1)
        mY= np.mean((mfli_data['Y']),1)
        sX= np.std((mfli_data['X']),1)
        sY= np.std((mfli_data['Y']),1)
        print("# averaged bunches = %d" % len(I0))
    return mX,mY,sX,sY,I0_filter,modfreq_filter
    
def MasterFileWriter(mfli_data,delay,mX,mY,sX,sY):
    d = {'run': [mfli_data['run']], 'delay': [delay],'mX0': mX[0], 'mY0': mY[0],'mX1': mX[1],'mY1': mY[1],'mX2': mX[2],'mY2': mY[2],'mX3': mX[3],'mY3': mY[3],'sX0': sX[0],'sY0': sY[0],'sX1': sX[1],'sY1': sY[1],'sX2': sX[2],'sY2': sY[2],'sX3': sX[3],'sY3': sY[3],'bunchnumber': [mfli_data['bunchnumber']],'timestamp' : [mfli_data['timestamp']]}
    df = pd.DataFrame.from_dict(d)    
    return df
    
def MasterFileReader(path):   
    dataDF = pd.read_csv(path, index_col=0) #Import file with data from MFLI:   
    data = dataDF.to_dict("list") #Convert the DataFrame to the dictionary format:
    return data

def DelayMove(delay_pos):
    if all(x == delay_pos[0] for x in delay_pos):
        delay_pos=delay_pos[0]
    else: print "Delaystage moved during run!"