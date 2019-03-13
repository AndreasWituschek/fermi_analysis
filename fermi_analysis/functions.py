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
from scipy import optimize
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

def PeakWidth(T, peakMargin, suscept, zeroPaddingFactor):
    """
    Calculates the FWHM of spectral peaks as well as the full width at which
    spectral peaks will be dropped to peakMargin (fraction of peak maximum)
    of their amplitude. For this purpose a simulated perfect cosine signal is
    Fourier transformed with Gaußian window and zero padding and its FWHM and
    peak margin (peakWidth) is measured and returned.
    """
    # simulate perferct cosine signal
    length = abs(T[-1]-T[0]) # in fs
    step = length/(len(T)-1)
    w = 2*np.pi*1/(20*step)   # choose frequency factor 10 below nyquist
#    paras = {'Schrittweite' : step, 'Monochromator Wavenumber' : 0 }
    sim = np.exp(1j*w*T)    # define simulated signal as in-phase + 1j*in-quad
    sim = GaussWindow(T, sim, suscept)
    wn, dft = DFT(T, sim, step, 1, 1, zeroPaddingFactor)
    if suscept:
        dft = dft.real
    else:
        dft = abs(dft)
    # normalize dft to 1
    dft /= dft.max()
    wnFitted, fFitted, fitParas = PeakFit("simulation", 1.0, wn, dft, np.abs(wn[-1]-wn[0]), 0.)

    # determine peak width of simulated signal
    peakFullWidth = abs(wn[findPeakIndexes(dft, peakMargin)[1]]-wn[findPeakIndexes(dft, peakMargin)[2]])
#    FWHM = abs(wn[findPeakIndexes(dft, 0.5)[1]]-wn[findPeakIndexes(dft, 0.5)[2]])
    # calculate widths from theoretical derivation of gaussian window width
#    delta = length*np.sqrt(np.log(2)/3)
#    FWHM = 8.0*np.log(2.0)/delta# FWHM of Gaussian Window
#    peakFullWidth = 4*np.sqrt(-np.log(2)*np.log(peakMargin))/delta
    return fitParas['FWHM'], peakFullWidth

def findPeakIndexes(Y, peakMargin):
    """
    Returns index of peak maximum, left and right index where dropped to
    fraction defined by peakMargin.
    """
    indexMax = np.argmax(abs(Y))    # find peak maximum
    indexLow = np.argmin(abs(Y[:indexMax+1]-peakMargin*Y.max())) # search index in first half of peak
    indexHigh = np.argmin(abs(Y[indexMax:]-peakMargin*Y.max())) + indexMax # search in second half of peak
    return  indexMax, indexLow, indexHigh

def PeakFit(peakName, peakWaveNo, wn, spectrum, windowLength, guessNoise):
    """
    Performes a gaussian fit of a spectral peak with no constant offset.
    windowLength : full width of fit window, usually windowLength=x*PeakWidth(T)
    spectrum : should be ABSOLUTE value of spectrum
    guessNoise : typical a definition like this is choosen: noiseFloor*(1+(Ham-1)*5)
    """
    totalError = False
    # construct fit window
#    wnFit, spectrumFit = FitData(wn, spectrum, peakWaveNo, windowLength)
    wnFit = wn
    spectrumFit = spectrum
    # guess fit paras
    pguess = GuessPeakParameter(wnFit, spectrumFit)
    sigmas = np.zeros(len(wnFit)) + guessNoise   # assume every data point has same errorbars
    # proceed fitting
    popt, pcov = optimize.curve_fit(f_fixoff, wnFit, spectrumFit, pguess[0:-1]) #, sigma=sigmas) #, absolute_sigma=True)
    #popt, pcov = optimize.curve_fit(f, wnFit, spectrumFit, pguess[0:-1], sigma=sigmas)
    errors = sp.sqrt(pcov.diagonal())
    if totalError:
        # calculate total error of fit function evaluated at maximum of fit function
        err_fmax = GaussFitError(wn,popt,pcov)
        # use total fit error as estimate for Amp_err
        errors[0] = err_fmax
        AmpErrText = 'Error Amp tot'
        AreaErrText = 'Error Area tot'
    else:
        AmpErrText = 'Error Amp'
        AreaErrText = 'Error Area'
    # calculate peak area:
    area = popt[0]*popt[2]*np.sqrt(np.pi/np.log(16))
    err_area = abs(area)*np.sqrt((errors[0]/popt[0])**2 + (errors[2]/popt[2])**2) # Gaussian error propagation

    wnFitted = np.linspace(wnFit[0], wnFit[-1], num=3000)
    #f_fitted = f(wn_fitted, popt[0], popt[1], popt[2], popt[3])
    fFitted = f_fixoff(wnFitted, popt[0], popt[1], popt[2])
    print('peak_fit: ' + peakName,popt[0],popt[1],popt[2],area)
    fitParas = {'Peak': peakName , \
                'Amp': round_sig(popt[0],3) , \
                AmpErrText: round_sig(errors[0],1) , \
                'Center': round(popt[1],3) , \
                'Error Center':round_sig(errors[1],1) , \
                'FWHM':round(popt[2],2) , \
                'Error FWHM':round_sig(errors[2],1) ,\
                'Area':round_sig(area,3) , \
                AreaErrText:round_sig(err_area,1)}
    return wnFitted, fFitted, fitParas

def FitData(wn, spectrum, peakWaveNo, windowLength):
    # construct fit window
    lowIndex = find_index(wn, peakWaveNo-0.5*windowLength)
    highIndex = find_index(wn, peakWaveNo+0.5*windowLength)
    wnFit = wn[lowIndex:highIndex+1].copy() # fit window of the wavenumber array
    spectrumFit = spectrum[lowIndex:highIndex+1].copy()  # fit window of the DFT amplitudes array
    return wnFit, spectrumFit

def GuessPeakParameter(Xfit, Yfit):
    """
    Estimates peak parameters for peak fit. XFit and Yfit should be restricted
    to a small area around the peak of interest. If more than one peak is within
    Xfit and Yfit, wrong parameters will be estimated. Provide ABSOLUTE value
    of spectrum.
    """
    pguess = np.zeros(4)
    pguess[0] = Yfit.max()   # Amp
    pguess[1] = Xfit[np.argmax(Yfit)]   # center
    pguess[2] = abs(Xfit[findPeakIndexes(Yfit, 0.5)[1]]-Xfit[findPeakIndexes(Yfit, 0.5)[2]]) # FWHM
    pguess[3] = 0.0 # offset
    return pguess

def GaussFitError(X,popt,pcov):
    """
    Calculates the total error of fit function for a Gaussian fit with no offset evaluated at maximum of fit function
    """
    a = popt[0]
    b = popt[1]
    c = popt[2]
    Y = f_fixoff(X, a, b, c)
    i = find_index(X,popt[1])   # find index of peak maximum
    # berechne derivatives
    da = Y[i]/a
    db = Y[i]*(8*np.log(2)/(c**2)*(X[i]-b))
    dc = Y[i]*(8*np.log(2)/(c**3)*(X[i]-b)**2)
    dadb = db/a
    dadc = dc/a
    dbdc = db*dc/Y[i]-2*db*Y[i]/c
    var_f = da**2*pcov[0][0] + db**2*pcov[1][1] + dc**2*pcov[2][2] +\
            2*(dadb*pcov[0][1] + dadc*pcov[0][2] + dbdc*pcov[1][2])
    sigma_f = np.sqrt(abs(var_f))
    return sigma_f

def f_fixoff(wn, Amp, wn_center, FWHM):
    """Gaussian function without vertical offset"""
    return Amp*sp.exp(-4*sp.log(2)*((wn-wn_center)/FWHM)**2)

def round_sig(x, sig):
   return round(x, sig-int(np.floor(np.log10(x)))-1)

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))