# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:22:02 2019

@author: Andreas
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pandas as pd
import scipy as sp
import scipy.constants as spc
from scipy import fftpack
import h5py

#physical constants:
c = spc.speed_of_light / 1000000.0  # speed of light in nm/fs
hPlanck = 4.135667662  # plancks constant in 10e-15 eV*s


def importFile(path):
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


def Curve(l_He, l_ref, h, phi, A, offset, start, stop, length):
    "creates datapoints on cosine/sine curve with given parameters for frequency, ampliotude and phase."
    points=30000 # number of TD points for plotting of theroetical curve
    plotRange = np.linspace(start-10,stop+10,points)
    tau = np.linspace(start,stop,length)
    Xtd=offset+A*np.cos(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + phi/180.0*np.pi)
    Ytd=offset+A*np.sin(2.*np.pi *c*plotRange*(h*l_He-l_ref)/(l_ref*l_He) + phi/180.0*np.pi)
    X=offset+A*np.cos(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180.0*np.pi)
    Y=offset+A*np.sin(2.*np.pi *c*tau*(h*l_He-l_ref)/(l_ref*l_He) + phi/180.0*np.pi)
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
    #plt.xlim([23,25])
    plt.axvline(x=transition, color='black', linestyle='--')
    plt.text(transition+0.02,0.8*max(abs(cdft)),'He 4P') #position of He4P transition
    plt.xlabel('Energy [eV]')
    plt.ylabel('Intensity [a.U.]')
    plt.legend()
    plt.title("Absorption spectrum")


def PlotFdCurveEV(xaxis,cdft,cdft_d,l_ref,l_He,h):
    "plots the absorptive and dispersive part of the spectrum of theoretical curve and data"
    plt.figure()
    transition=c/l_He*hPlanck #helium transition energy
    cdft= cdft/max(abs(cdft))*max(abs(cdft_d)) #normalizing theoretical curve on experimental data
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

def ReadMfliData_fake(mfli_file_path,bunches): #bunches= number of shots in one file
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

def ReadMfliData(mfli_file_path): #bunches= number of shots in one file
    """ reads data from files created with mfli.py containing the raw
        demodulated data.
    """
    h5f = h5py.File(mfli_file_path, 'r')
    # Create dict with:
    mfli_data = {}
    for key in h5f.keys():
        mfli_data.update({key: np.array(h5f[key])})
    h5f.close()
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

        modfreq_filter= np.asarray([-modfreq_threshold < i-modfreq_set < modfreq_threshold for i in mfli_data['frequency']]) #filter for modfreq
        modfreq_filter=np.concatenate([modfreq_filter,np.ones(ignore)])
        a=np.zeros((ignore,len(modfreq_filter)))
        for i in range(ignore):
            a[i]=np.roll(modfreq_filter,i)
        modfreq_filter=(np.sum(a,0)>(ignore-1))[ignore-1:len(I0)+ignore-1] #modfreq filter extended by number of ignored shots

        b = np.logical_and(I0_filter,modfreq_filter) #combining both filters above

        if sum(b)>0:
            print("# total bunches = %d" % len(I0))
            print("# filtered Bunches = %d" % sum(b))
            mX= [np.mean(mfli_data['x2'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            mY= [np.mean(mfli_data['y2'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
            sX= [np.std(mfli_data['x2'][i][np.where(b)]) for i in range(len(mfli_data['X']))]
            sY= [np.std(mfli_data['y2'][i][np.where(b)]) for i in range(len(mfli_data['Y']))]
        else:
            print('Filters to tight, no data to average over. Script stopped, no data written in masterData.csv!')
            sys.exit(1)
    else:
        mX= np.mean((mfli_data['x2']))
        mY= np.mean((mfli_data['y2']))
        sX= np.std((mfli_data['x2']))
        sY= np.std((mfli_data['y2']))
        print("# averaged bunches = %d" % len(I0))
        I0_filter=[]
        modfreq_filter=[]
    return mX,mY,sX,sY,I0_filter,modfreq_filter


#def MasterFileWriter(mfli_data, delay, mX, mY, sX, sY):
#    d = {'run': [mfli_data['run']],\
#        'delay': [delay],\
#        'mX0': [mX[0]],'mY0': [mY[0]],'mX1': [mX[1]],'mY1': [mY[1]],'mX2': [mX[2]],'mY2': [mY[2]],'mX3': [mX[3]],'mY3': [mY[3]],\
#        'sX0': [sX[0]],'sY0': [sY[0]],'sX1': [sX[1]],'sY1': [sY[1]],'sX2': [sX[2]],'sY2': [sY[2]],'sX3': [sX[3]],'sY3': [sY[3]],\
#        'bunchnumber': [mfli_data['bunchnumber']],\
#        'timestamp' : [mfli_data['timestamp']]\
#        }
#    df = pd.DataFrame.from_dict(d)
#    return df


#def MasterFileReader(path):
#    data = np.genfromtxt(path, delimiter='\t', dtype='|S')
#    a = np.array(data[1:]).astype(np.float)
#    sorted_data = np.array(sorted(a, key=lambda a_entry: a_entry[2]))
#    data_dict = dict(zip(data[0][1:], sorted_data.astype(np.float).T[1:]))
##    print(sorted_data)
#    print('dict anfang')
#    print(data_dict['delay'])
#    print('dict ende')
#    return data_dict


def DelayMove(delay_pos):
    if all(x == delay_pos[0] for x in delay_pos):
        delay_pos=delay_pos[0]
    else: print "Delaystage moved during run!"


def GaussWindow(T, X, suscept):
    """
    Applies Gaussian window to data set which decays to 5% at edges of data set
    . If suscept==False, the window function is centered with respect to the
    data set. Otherwise, it is centered at T[0] which is 0fs if susceptibility
    is wanted.
    """
    start = T[0]
    stop = T[-1]
    if suscept:
        mu = 0.0
        std = 2*(abs(stop-start))/sp.sqrt(48)  # so Gaussian will be dropped to 5% (e.g. 1/e^3) at edges of dataset
    else:
        mu = (stop+start)/2
        std = (abs(stop-start))/sp.sqrt(48)  # so Gaussian will be dropped to 5% (e.g. 1/e^3) at edges of dataset
    return X*sp.exp(-((T-mu)/(2*std))**2)


def zero_padding(T, factor):
    return 2**(int(round(np.log2(T.size)+0.5))+factor)  # Anzahl Datenpunkte auf die zur Not mit Zero-Padding aufgefuellt wird


def find_index(array, value):
    return np.argmin(abs(array-value))  # loops through array from index zero and returns first index that is closest to value


def redistribute(array, index):
    return np.concatenate((array[index:],array[0:index]))

def CutDataSet(T, X, Y, data_window):
    # T and data_window should be of same unit, e.g. femtoseconds
    # sort start,stop according to increasing/decreasing numbers in T
    if T[-1]>T[0]:
        start_user = min(data_window)
        stop_user = max(data_window)
    else:
        start_user = max(data_window)
        stop_user = min(data_window)
    # cut data
    lower_index = np.argmin(np.abs(T-start_user))
    upper_index = np.argmin(np.abs(T-stop_user))
    T = T[lower_index : upper_index]
    X = X[lower_index : upper_index]
    Y = Y[lower_index : upper_index]
    return (T,X,Y)

def DFT(T, Z, dT, l_ref, harmonic, zeroPaddingFactor):
    """
    Calculates DFT of complex valued input array Z. Applays zero padding as spe
    cified by input factor. Establishes correct monotonoic order of frequency
    axis. Returns wavenumber array and complex valued DFT.
    """
    print(np.diff(T))
    step = dT
    wn_Mono = 10**7/l_ref
    # determine zero padding
    N = zero_padding(T, zeroPaddingFactor)
    # calculate wavenumber axis
    freq = fftpack.fftfreq(N, d = step*10**(-15))  # frequency axis in Hz
    if step<0:
        freq = freq[::-1]   # reverse array
    cut_index = np.argmin(freq)
    freq = redistribute(freq, cut_index)    # sortiere wn Array so, dass Wellenzahlen monoton steigen
    wn = freq/100.0/spc.c + harmonic*wn_Mono # wn wird um Mono geshiftet wie in rpt File gegeben
    # calculate DFT
    DFT = redistribute(fftpack.fft(Z, n=N), cut_index)    # sortiere wn Array so, dass Wellenzahlen monoton steigen
    if step < 0:
        DFT = DFT[::-1]   # reverse array
    # correct for delay offset error
#    offset = T[find_index(T,0.0)]*10**(-15)   # in s
#    phaseCorrection = np.array([np.exp(-1j*2*np.pi*f*offset) for f in freq])
#    DFT *= phaseCorrection
    return wn, DFT

def slide_window(T, Z, t_center, FWHM):
    # multiplies a gaussion onto the dataset and returns the data set of the whole range of T
    # T, t_center, FWHM in ps
    # Z may be complex
    return Z*sp.exp(-4*sp.log(2)*((T-t_center)/FWHM)**2)

def weighting_coeff(t, t_center, FWHM): # calculates amplitude of gaussian at time t
    return sp.exp(-4*sp.log(2)*((t-t_center)/FWHM)**2)

def plot_spectrogram(fig, wn_lim, t_lim):
    ax = fig.add_axes([Offx, Offy, Rsubx, 2*Rsuby])  # define current axes
    xmin_index = np.argmin(abs(wn-wn_lim[0]))
    xmax_index = np.argmin(abs(wn-wn_lim[1]))
    tmin_index = np.argmin(abs(T-t_lim[0]))
    tmax_index = np.argmin(abs(T-t_lim[1]))
    clim= [0.0, 0.1*S.max()] # 0.65
    plt.imshow(S.transpose()[xmin_index:xmax_index, tmin_index:tmax_index], clim=clim, origin='lower', aspect='auto', extent=[T[tmin_index],T[tmax_index], wn[xmin_index],wn[xmax_index]])
    ax.set_ylabel('frequency (cm$^{-1}$)', fontsize = fs)
    ax.set_xlabel('pump-probe delay (ps)', fontsize = fs)
    ax.tick_params(labelsize=fs)
    locs, labels = plt.yticks()
    plt.setp(labels, rotation=90)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().set_tick_params(which='both', direction='out', pad=-5, length=5, width=1.5)
    ax.get_yaxis().set_tick_params(which='both', direction='out', pad=-5, length=5, width=1.5)
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax.get_yaxis().tick_left()
    ax.get_yaxis().set_tick_params(which='both', length=5, width=1.5)
    ax.get_xaxis().set_tick_params(which='both',  length=5, width=1.5)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    #ax.axvspan(-overlap,overlap, color='grey', alpha=0.9)
    ax.axhline(12816, linestyle='--', color='w', linewidth=lw)
    ax.axhline(12579, linestyle='--', color='w', linewidth=lw)
#    ax.axhline(11685, linestyle='--', color='w', linewidth=lw)
    ax.axhline(mono_calib, linestyle='--', color='k', linewidth=lw)
    ax.set_xlim(t_lim[0], t_lim[1])
#    ax.text(0.9, 0.9, '(b)', transform=ax.transAxes, ha='left', fontweight='bold',fontsize = 22, color='w')

    #annotate()
    plotbar()