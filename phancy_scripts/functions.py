# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 16:48:46 2016

Known bugs:

TODO:
    * Example in class Peak does not work anymore and remove is missing.

@author: LukasB
"""
import scipy as sp
import scipy.constants as const
from scipy import fftpack
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
from math import log10, floor
from Tkinter import Tk
from tkFileDialog import askopenfilename
import csv
from scipy.integrate import simps   # method to integrate fixed sample which is unequally spaced
import ast
from scipy.odr import *

# File IO =====================================================================
def OpenDialog (initDir, title):
    """
    Creates open file dialog.
    Returns user specified file path.
    """
    # define Masterwindow/Canva and close hide it
    root = Tk()
    root.withdraw()
    # show an "Open" dialog box and return the path to the selected file
    filePath = askopenfilename(initialdir= initDir, defaultextension='.dat', title=title) # returns complete file path
    return filePath
    
def OpenMultipleFiles(initDir, command):
    """
    Dispalys open file dialog several times until user press cancel on open 
    dialog. Then all user-chosen file locations are stored in a list and returned.
    """
    files = [OpenDialog(initDir, command)]  # first file
    while files[-1] != '':
        files.append(OpenDialog(initDir, command))
    return files[0:-1]

def GetFileDictName(file_path):
    """
    Returns file directory and file name without ending '.dat'
    """
    index1 = string.rfind(file_path, '/') #python 2.7
    index2 = string.rfind(file_path, '\\') #python 2.7
    index = max([index1, index2])
#    index = file_path.rfind('/') #python 3.*
    return file_path[:index+1], file_path[index+1:-4]    #truncate '.dat'

def ImportReport(reportFileName):
    """
    Opens report file specified by ReportFileName located in current directory.
    Returns dictionary with all keywords and values listed in report file.
    """
    report = np.genfromtxt(reportFileName + '.rpt', delimiter = ':', skip_header = 12,\
                    usecols=(0,1), dtype=[('keys','S40'), ('values','<f8')])    # wichtig dtype
    return dict(zip(report['keys'],report['values']))

def ReadParameter(filePath, X, key='' , listIndex=None):
    """
    Opens rpt-File and reads parameter specified by key. If key contains a list
    one can optionally provide for listIndex to choose an individual vlue from 
    list. Otherwise the hole list is taken. This/these values are appended to X.
    """
    with open(filePath[:-4] + ".rpt", "rb") as myfile:
        csv_read = csv.reader(myfile, delimiter=':')    # csv_read is a list of rows in the report file
        for row in csv_read:
            if len(row)>1:
                if row[0]==key: 
                    row[1] = ast.literal_eval(row[1])
                    if listIndex != None:
                        X.append(row[1][listIndex])
                    else: 
                        X.append(row[1])    
    
def ImportLaserSpectrum(filePath, boxcarWindow, boxcarLength, skipHeader=None, delimiter=None):
    """
    Imports spectrum data acquired with OceanView software. 
    Applies boxcar smoothing:
    boxcarWindow: window-type ('flat', 'hanning', 'hamming', 'bartlett', 'blackman')
    boxcarLength: No of sample points for each smoothing step
    No boxcar averaging is performed for boxcarLength < 3
    if keyword agr 'skipHeader' is supplied, code will skip header lines as defined. 
    Otherwise header lines will be skipped until '>>>>>' string appears.  
    Returns wavelength, wavenumber and power spectral density
    For Avantes spectrometer specifiy delimiter=';' an skipHeader=8
    """
    if delimiter == None:
        delimiter = '\t'
    # open data
    with open(filePath, "r") as data:
    # import header
        header = [] 
        if skipHeader == None: skipHeader = 100    # if user did not specify skip header, then skip header until '>>>>' string appears
        while len(header)<skipHeader:    # to make sure, that no endless loop
            line = data.readline()
            header.append(line)
            if line.startswith('>>>>>'): break            
    # import data 
        spec = np.genfromtxt(data, delimiter = delimiter, comments = '#', skip_header = skipHeader, skip_footer = 1, usecols=(0,1), dtype='<U40')
    # convert comma to decimal point
    wavelength = np.array(list(map(lambda valstr: float(valstr.replace(',','.')),spec.T[0])))
    PSD = np.array(list(map(lambda valstr: float(valstr.replace(',','.')),spec.T[1])))
    # apply boxcar smoothing 
    wavelength = BoxcarSmoothing(wavelength, boxcarLength, boxcarWindow)   
    PSD = BoxcarSmoothing(PSD, boxcarLength, boxcarWindow)
    # calculate wavenumber array
    wavenumber = 10**7/wavelength
    return wavelength, wavenumber, PSD

def BoxcarSmoothing(x,window_len=11,window='flat'):
    """
    Performs boxcar smoothing of data array x.
    window: defines window-type for smoothing ('flat', 'hanning', 'hamming', 'bartlett', 'blackman')
    No boxcar averaging is performed for window_len < 3
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def DataScaling(report,X,Y,No):
    amp = report[b'Lock-In_' + No + b' Amplification']
    sensitivity = report[b'Lock-In_' + No + b' Sensitivity /V,A ']
    X *= sensitivity/amp/10.0  # orignial signal magnitude in Volts
    Y *= sensitivity/amp/10.0  # orignial signal magnitude in Volts
    return (X,Y)

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

def ImportWPIData(filePath, noChannels, monoCorrWn, refLD, data_window):
    """
    Imports WPI data and report file.
    Scales data according to LockIn sensitivity and amplification.
    Reads monochromator wavelength and applies monochromator calibration. 
    Takes user specified monochromator correction and correction specified in 
    report file, if listed.
    Returns data and dictionary of all experimental parameters.
    """
    fileName = GetFileDictName(filePath)[1]    # extracts filename from file path, i.e. removes directory and data type
    # import parameters
    paras = ImportReport(filePath[:-4])
    paras['File Path'] = filePath
    title = fileName
    paras['Title'] = title  # add title to paras dictionary
    lambda_L = paras['Chameleon Spektrum']
    paras['Laser Wavenumber'] = 10**7/lambda_L
    lambda_M_scale = paras['Monochromator Wavelength /nm ']
    if 'Monochromator Correction /cm-1' in paras.keys():
        monoCorrWn += paras['Monochromator Correction /cm-1']
    if refLD:
        lambda_M = lambda_M_scale   # no calibartion correction has to be applied
    else:
        # if Mono is used for referencening and not laser diode, than apply calibration
        lambda_M = lambda_M_scale*1.004-49.872  # Eichung vom 22.01.14
    paras['Monochromator Wavenumber'] = 10**7/lambda_M + monoCorrWn

    # import WPI data
    if noChannels==0:
        data = np.genfromtxt(filePath, delimiter = '\t', skip_header = 1 , \
        usecols = (2))
        T = data.T
        X = np.zeros(len(T))
        Y = np.zeros(len(T))
    elif noChannels==1:
        data = np.genfromtxt(filePath, delimiter = '\t', skip_header = 1 , \
        usecols = (2,3,4))
        T,X1,Y1 = data.T
        X = [X1]
        Y = [Y1]
    elif noChannels==2:
        data = np.genfromtxt(filePath, delimiter = '\t', skip_header = 1 , \
        usecols = (2,3,4,5,6))
        T,X1,Y1,X2,Y2, = data.T
        X = [X1,X2]
        Y = [Y1,Y2]
    elif noChannels==3:
        data = np.genfromtxt(filePath, delimiter = '\t', skip_header = 1 , \
        usecols = (2,3,4,5,6,7,8))
        T,X1,Y1,X2,Y2,X3,Y3 = data.T
        X = [X1,X2,X3]
        Y = [Y1,Y2,Y3]
    elif noChannels > 3:
        data = np.genfromtxt(filePath, delimiter = '\t', skip_header = 1 , \
        usecols = (2,3,4,5,6,7,8,9,10))
        T,X1,Y1,X2,Y2,X3,Y3,X4,Y4 = data.T
        X = [X1,X2,X3,X4]
        Y = [Y1,Y2,Y3,Y4]

    LockInNo=['1','2','3','4']
    T_dummy = T
    for i in range(noChannels):
        X[i],Y[i] = DataScaling(paras,X[i],Y[i],LockInNo[i])
        T, X[i], Y[i] = CutDataSet(T_dummy, X[i], Y[i], data_window)
    # return scaled LockIn data and acquisition parameters
    return T,X,Y,paras
    
def AddLeakPeaks(leakingPeaks, Peak, paras, Harmonics, annotate, leakColor, analysis, monoName):
    """
    Adds leak peaks to the peak register. As origin for leak peaks in higher harmonics
    will be considered the Peaks defined in leakingPeaks. If certain leak Peaks 
    have been added already to the Peak register, their wavenumber will be adjusted
    according to the current calibrated monochromator wavenumber.
    """
    calMono = paras['Monochromator Wavenumber']    # get current calibrated mono wavenumber
    for h in Harmonics:       
        for peak in leakingPeaks:  # loop through all peaks that cause leaking
            peakName = peak.get_name()
            peakWn = peak.get_wavenumber()
            peakHarmonic = peak.get_harmonic()
            if peakName != (str(peakHarmonic)+monoName) and h != peakHarmonic:  # exclude mono peaks in leakingPeaks and don't add leak peaks to harmonic from which leak is originating
                leakPeak = Peak.select(name=peakName+'_leak', harmonic=h)  # ist leak peak in harmonic h bereits eingetragen?
                if len(leakPeak)>0: # falls leak Peak bereits vorhanden, überschreibe Wellenlänge mit neuer Monochromatorkalibration, da jeder Datensatz eigene Mono-Kalib hat
                    leakPeak[0].set_wavenumber(peakWn+(h-1)*calMono)
                else:   # falls leak Peak noch nicht vorhanden lege neu an
                    Peak.add_line(name=peakName+'_'+str(h)+'Hleak',wavenumber=peakWn+(h-1)*calMono,color=leakColor,annotate=annotate,analysis=analysis,window=0,harmonic=h)


def InitReport(filePath,analysis,scriptAcronym, reportDict):
    """
    Initializes report if global variable report=TRUE. Report is appended to current .rpt file
    """    
    with open(filePath[:-4] + ".rpt", "a") as myfile:
        csv_write = csv.writer(myfile, delimiter=':')
        csv_write.writerow([])
        csv_write.writerow([])
        csv_write.writerow(['#---------------------------------------------------------------------'])
        csv_write.writerow(['# '+analysis])
        csv_write.writerow(['# Skript ' + scriptAcronym + '.py'])
        csv_write.writerow(['#---------------------------------------------------------------------'])
        for key, value in reportDict.iteritems():
            csv_write.writerow([key] + [value])

def ReportLine(filePath, key, value):
    """
    Appends new line in .rpt file if global variable report=TRUE
    """
    with open(filePath[:-4] + ".rpt", "a") as myfile:
        csv_write = csv.writer(myfile, delimiter=':')
        csv_write.writerow([key] + [value])   # write new line
            
def ReportComment(filePath, comment):
    """
    Appends new comment line in .rpt file if global variable report=TRUE
    """
    with open(filePath[:-4] + ".rpt", "a") as myfile:
        csv_write = csv.writer(myfile, delimiter=':')
        csv_write.writerow(['# ' + comment])   # write new line

def ReportEmptyLine(filePath):
    """
    Appends empty line in .rpt file if global variable report=TRUE
    """
    with open(filePath[:-4] + ".rpt", "a") as myfile:
        csv_write = csv.writer(myfile, delimiter=':')
        csv_write.writerow([])   # write new line
            
def ReportNoise(filePath, Harmonics, reportNoise, reportNoiseExcluded):
    """
    Appends calculated noise for all harmonics to .rpt file if global variable report=TRUE
    """
    for i in range(len(Harmonics)):
        ReportLine(filePath, 'noise'+str(Harmonics[i])+'H', reportNoise[i][:])
    for i in range(len(Harmonics)):
        ReportLine(filePath, 'noise'+str(Harmonics[i])+'H'+'excluded', reportNoiseExcluded[i][:])

def ReportPeaks(filePath, Harmonics, reportPeakList):
    """
    Appends fit results for all peaks to .rpt file if global variable report=TRUE
    """
    for h in Harmonics:
        ReportEmptyLine(filePath)
        ReportComment(filePath, str(h)+'H peaks') # jede Harmonische wird mit Kommentar gekennzeichnet
        ReportComment(filePath, 'Peak fit list [fit type, Center, Center_err, Amp, Amp_err, Area, Area_err, FWHM, FWHM_err]')    
        for fitParas in reportPeakList:
            # check if total fit error taken for amp Error estimate            
            if 'Error Amp tot' in fitParas:        
                keyAmpErr = 'Error Amp tot' 
                keyAreaErr = 'Error Area tot'
            else: 
                keyAmpErr = 'Error Amp'
                keyAreaErr = 'Error Area'
            for peak in Peak.select(harmonic=h):
                if fitParas['Peak'] == peak.get_name():
                    peakReport = [
                    fitParas['type'],\
                    fitParas['Center'], fitParas['Error Center'] ,\
                    fitParas['Amp'], fitParas[keyAmpErr] ,\
                    fitParas['Area'], fitParas[keyAreaErr] ,\
                    fitParas['FWHM'], fitParas['Error FWHM'] ,\
                    ]
                    ReportLine(filePath,fitParas['Peak'], peakReport)
                        
def init_csv(filePath, scriptAcronym, logParas):
    with open(filePath+".csv", "wb") as myfile:
        csv_write = csv.writer(myfile, delimiter='\t', dialect='excel')
        csv_write.writerow(['#---------------------------------------------------------------------'])
        csv_write.writerow(['# Skript ' + scriptAcronym + '.py'])
        csv_write.writerow(['#---------------------------------------------------------------------'])
        for key, value in logParas.iteritems():
            csv_write.writerow([key] + [value])
        csv_write.writerow([])

def write_2_csv(filePath, line):
    with open(filePath+".csv", "a") as myfile:
        csv_write = csv.writer(myfile, delimiter='\t', dialect='excel')
        csv_write.writerow(line)   # write new line

def SimWPISignal(T, peakDict, harmonic, paras):
    """
    Simulate cos+isin signal für a list of peaks. List should be delivered as 
    dict: peakDict = {'peak1' : [wn, amp, phase], ...}
    T in fs, phase in degree
    """ 
    simSig = np.zeros(len(T))+1j*np.zeros(len(T))
    for name, p in peakDict.iteritems():
        # wichtig dass man p[i] an lokale variablen übergibt und nicht überschreibt,
        # ansosnten sind diese Werte auch im Hauptskript in peakDict überschrieben
        freq = (p[0]-harmonic*paras['Monochromator Wavenumber'])*(100.0*const.c)
        phi = p[2]/360.0*2*np.pi # conversion to radians 
        amp = p[1]
        simSig += amp*(np.cos(2*np.pi*freq*T*10**-15+phi)+ \
    1j*np.sin(2*np.pi*freq*T*10**-15+phi))
    return simSig

# Data Analysis ===============================================================
def round_sig(x, sig):
   return round(x, sig-int(floor(log10(x)))-1)

def redistribute(array, index):
    return np.concatenate((array[index:],array[0:index]))

def find_index(array, value):
    return np.argmin(abs(array-value))  # loops through array from index zero and returns first index that is closest to value

def Dict2List(d):
    l = []
    for key, value in d.iteritems():
        l.append(key)
        l.append(value)
    return l

def findPeakIndexes(Y, peakMargin):
    """
    Returns index of peak maximum, left and right index where dropped to
    fraction defined by peakMargin.
    """
    indexMax = np.argmax(abs(Y))    # find peak maximum
    indexLow = np.argmin(abs(Y[:indexMax+1]-peakMargin*Y.max())) # search index in first half of peak
    indexHigh = np.argmin(abs(Y[indexMax:]-peakMargin*Y.max())) + indexMax # search in second half of peak
    return  indexMax, indexLow, indexHigh

def zero_padding(T, factor):
    return 2**(int(round(np.log2(T.size)+0.5))+factor)  # Anzahl Datenpunkte auf die zur Not mit Zero-Padding aufgefuellt wird

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
        mu = (stop+start)/2.
        std = (abs(stop-start))/sp.sqrt(48)  # so Gaussian will be dropped to 5% (e.g. 1/e^3) at edges of dataset
    return X*sp.exp(-((T-mu)/(2*std))**2)

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
    paras = {'Schrittweite' : step, 'Monochromator Wavenumber' : 0 }
    sim = np.exp(1j*w*T)    # define simulated signal as in-phase + 1j*in-quad
    sim = GaussWindow(T, sim, suscept)
    wn, dft = DFT(T, sim, paras, 1, zeroPaddingFactor)
    if suscept:
        dft = dft.real
    else:
        dft = abs(dft)
    # normalize dft to 1
    dft /= dft.max()
    # determine peak width of simulated signal
    peakFullWidth = abs(wn[findPeakIndexes(dft, peakMargin)[1]]-wn[findPeakIndexes(dft, peakMargin)[2]])
    FWHM = abs(wn[findPeakIndexes(dft, 0.5)[1]]-wn[findPeakIndexes(dft, 0.5)[2]])
    # calculate widths from theoretical derivation of gaussian window width
#    delta = length*np.sqrt(np.log(2)/3)
#    FWHM = 8.0*np.log(2.0)/delta# FWHM of Gaussian Window
#    peakFullWidth = 4*np.sqrt(-np.log(2)*np.log(peakMargin))/delta
    return FWHM , peakFullWidth

def DFT(T, Z, paras, harmonic, zeroPaddingFactor):
    """
    Calculates DFT of complex valued input array Z. Applays zero padding as spe
    cified by input factor. Establishes correct monotonoic order of frequency
    axis. Returns wavenumber array and complex valued DFT.
    """
    step = paras[b'Schrittweite']
    wn_Mono = paras['Monochromator Wavenumber']
    # determine zero padding
    N = zero_padding(T, zeroPaddingFactor)
    # calculate wavenumber axis
    freq = fftpack.fftfreq(N, d = step*10**(-15))  # frequency axis in Hz
    if step<0:
        freq = freq[::-1]   # reverse array
    cut_index = np.argmin(freq)
    freq = redistribute(freq, cut_index)    # sortiere wn Array so, dass Wellenzahlen monoton steigen
    wn = freq/100.0/const.c + harmonic*wn_Mono # wn wird um Mono geshiftet wie in rpt File gegeben
    # calculate DFT
    DFT = redistribute(fftpack.fft(Z, n=N), cut_index)    # sortiere wn Array so, dass Wellenzahlen monoton steigen
    if step<0:
        DFT = DFT[::-1]   # reverse array
    # correct for delay offset error
    offset = T[find_index(T,0.0)]*10**(-15)   # in s
    phaseCorrection = np.array([np.exp(-1j*2*np.pi*f*offset) for f in freq])
    DFT *= phaseCorrection    
    return wn, DFT
   
def Normalize(normalize, T, X, Y, paras, zeroPaddingFactor, Harmonics):
    """
    Returnes normalized data. normalize: 'none', 'delay', 'nH'
    'delay' : normalize to pump probe delay, 1ps delay correspnds to factor of 1
    'nH' : defines to which harmonic should be normalized, e.g. '1H'
    T should be in fs, and should be already cut to used data window
    """    
    if normalize =='none':
        normFactor = 1.0        
    elif normalize == 'delay':  #Fourier amplitudes scale linear with delay length
        length = abs(T[-1]-T[0])    #T in fs
        normFactor = length*1E-3    # 1ps corresponds to normFactor=1
    else: 
        normalizeH = int(filter(str.isdigit, normalize))    # extract Harmonic             
        index = Harmonics.index(normalizeH)  # index of harmonic to which should be normalized
        Znorm = X[index] + 1j*Y[index]    
        Zg = GaussWindow(T, Znorm, suscept)
        wn, dft = DFT(T, Zg, paras, 0, zeroPaddingFactor)   # harmonic is irrelevant, since wavenumber array is not needed, only dft amplitudes are needed
        normFactor = max(abs(dft))
    return normFactor

def NoiseFloor(wn, spectrum, Peaks, windowLength, noiseDef):
    """
    Extracts the noise floor of a given spectrum. It excludes spectral peaks.
    Peaks is a list of peak objects. This should include all spectral peaks that
    should not be included in noise evaluation (i.e. including mono peak, leak
    peaks). WindowLength is the full length around a peak to be excluded.
    A wavenumber vector and a noise vector is returned, containig the data that
    should be used for noise evaluation.
    """
    # construct noise vector by excluding all peaks
    for peak in Peaks:
        indexPeakLow = find_index(wn, peak.get_wavenumber()-0.5*windowLength)   # will return the last index of wn if peak not in wn
        indexPeakHigh = find_index(wn, peak.get_wavenumber()+0.5*windowLength)
        if indexPeakLow < indexPeakHigh: 
            wn = np.delete(wn,range(indexPeakLow,indexPeakHigh+1))   # remove peak values from wn vector        
            spectrum = np.delete(spectrum,range(indexPeakLow,indexPeakHigh+1))   # remove peak values from spectrum
        else:
            print peak.get_name() + ' is not in data range'
    # calculate noise
    avgNoise = np.mean(spectrum)
    stdNoise = np.std(spectrum)
    noiseLevel = avgNoise + noiseDef*stdNoise
    return avgNoise, stdNoise, noiseLevel, wn, spectrum

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

def f_fixoff(wn, Amp, wn_center, FWHM):
    """Gaussian function without vertical offset"""
    return Amp*sp.exp(-4*sp.log(2)*((wn-wn_center)/FWHM)**2)
    
def f_fixoffSum(wn, *p):
    """
    Sum of Gaussians. *p is Pointer on a list with following structure: [Amp1, 
    wn_center1, FWHM1, Amp2, wn_center2, FWHM2, Amp3,...]. Always pass pointer 
    *p instead of p!    
    """    
    # deconvolve list p    
    Amp = p[0::3]
    wn_center = p[1::3]
    FWHM = p[2::3]
    # construct sum of gaussians. It is important to explicitely use namespaces for numpy here
    # otherwise curve_fit does not work with this function
    f = sum([ Amp[i]*np.exp(-4*np.log(2)*((wn-wn_center[i])/FWHM[i])**2) for i in range(len(Amp)) ])
    return f

def f_SumFixed(wn, *amps):
    """
    Sum of Gaussians.   
    """    
    # construct sum of gaussians. It is important to explicitely use namespaces for numpy here
    # otherwise curve_fit does not work with this function
    f = sum([amps[i]*np.exp(-4*np.log(2)*((wn-fixedCenters[i])/fixedFWHMs[i])**2) for i in range(len(fixedCenters)) ])
    return f
    
def f_SumFixedPara(wn, fixedCenters, fixedFWHMs, *amps):
    f = sum([amps[i]*np.exp(-4*np.log(2)*((wn-fixedCenters[i])/fixedFWHMs[i])**2) for i in range(len(amps))])
    return f
    
def LinearFitLogLog(X,Y,Xerr,Yerr,key, omittedPoints):
    """
    Performs linear fit in a log-log scale. Considers X- and Yerrors
    Attention: Xerr and Yerr should not contain elements equal to zero!!
    """
    # remove omitted points
    # linear fit data in log-log
    Xfit = []
    Yfit = []
    Xerr_fit = []
    Yerr_fit = []
    for i in range(len(X)):
        if i not in omittedPoints : Xfit.append(np.log(X[i]))
        if i not in omittedPoints : Yfit.append(np.log(Y[i]))
        # to avoid devision by zero, replace in X and Y all zero elements with smallest element
        X_min = min([x for x in X if x>0])
        Y_min = min([y for y in Y if y>0])
        if i not in omittedPoints : 
            if X[i]==0.0:
                Xerr_fit.append(Xerr[i]/X_min)   # avoid devision by zero
            else:
                Xerr_fit.append(Xerr[i]/X[i])   # proper error propagation of log(x)                
        if i not in omittedPoints : 
            if Y[i]==0.0:
                Yerr_fit.append(Yerr[i]/Y_min)   # avoid devision by zero
            else:
                Yerr_fit.append(Yerr[i]/Y[i])   # proper error propagation of log(x)
    Xfit = np.array(Xfit)
    Yfit = np.array(Yfit)
    Xerr_fit = np.array(Xerr_fit)
    Yerr_fit = np.array(Yerr_fit)
    mGuess = (Yfit[0]-Yfit[-1])/(Xfit[0]-Xfit[-1])
    cGuess = Yfit[0]
    def f_lin(p, x):
        m, c = p
        return m*x+c    
    linear_model = Model(f_lin) # Create a model for fitting.
    # Create a RealData object using our initiated data from above.
    data = RealData(Xfit, Yfit, sx=Xerr_fit, sy=Yerr_fit)
    # Set up ODR with the model and data.
    odr = ODR(data, linear_model, beta0=[mGuess, cGuess])
    # Run the regression.
    out = odr.run()
    # Use the in-built pprint method to give us results.
    print key + ' ----------------------------------------------------------------' 
    out.pprint()
    print '-----------------------------------------------------------------------'
    
    legend_fit = 'lin. fit: m='+str(round(out.beta[0],2))+ '$\pm$' + str(round(out.sd_beta[0],2))
    fFitted = f_lin(out.beta, np.log(X))
    return X, np.exp(fFitted), legend_fit

def PeakFit(peakName, peakWaveNo, wn, spectrum, windowLength, guessNoise):
    """
    Performes a gaussian fit of a spectral peak with no constant offset.
    windowLength : full width of fit window, usually windowLength=x*PeakWidth(T)
    spectrum : should be ABSOLUTE value of spectrum
    guessNoise : typical a definition like this is choosen: noiseFloor*(1+(Ham-1)*5)
    """
    totalError = False
    # construct fit window
    wnFit, spectrumFit = FitData(wn, spectrum, peakWaveNo, windowLength)
    # guess fit paras
    pguess = GuessPeakParameter(wnFit, spectrumFit)
    sigmas = np.zeros(len(wnFit)) + guessNoise   # assume every data point has same errorbars
    # proceed fitting
    popt, pcov = optimize.curve_fit(f_fixoff, wnFit, spectrumFit, pguess[0:-1], sigma=sigmas, absolute_sigma=True)
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

def getFitPeaks(peakCenterList, peaks):
    peakList = []
    for wn in peakCenterList:
        for peak in peaks:
            if wn==peak.get_wavenumber():
                peakList.append(peak)
    return peakList

def BundelList(peakList, peakWidth):
    """
    Splitts up peakList in bundles of overlapping peaks. Returns a list of lists
    conatining bundles of Peak objects that have overlapping peak centers. The
    overlap criterium is peakCenter+/-peakWidth
    Possible issue: The peakList is transfered to a wavenumber list and 
    transfered back. It could be that due to float precission the back transform
    is in some case not unique.
    """
    tempList = [] # list of peak wavenumbers
    bundeledList = []
    # create list ordered with increasing wavenumbers
    centers=[]    
    for peak in peakList:
        centers.append(peak.get_wavenumber())
    centers.sort()   # sort produces a list of monotonic increasing wavenumbers
    # loop through centers and perform fitting
    for i in range(0, len(centers)):
        tempList.append(centers[i])
        if i==len(centers)-1:   # aboard criterium of loop
            bundeledList.append(tempList)
            break
        # next peak overlapping with current peak?                   
        if abs(tempList[-1]-centers[i+1])>peakWidth:
            bundeledList.append(tempList) # last peak is not overlapping with adjecent peak
            tempList = [] # start new bundle
    return bundeledList
    
def MultiPeakFit(Peaks, wn, spectrum, peakFWHM, singlePeakMargin, guessNoise):
    """
    Fits a sum of peaks with a sum of gaussians, having no vertical offset.
    peakFWHM : expected FWHM of individual deconvolved peaks
    singlePeakMargin : expected full width of individual deconvolved peak
    spectrum : should be ABSOLUTE value of spectrum
    guessNoise : typical a definition like this is choosen: noiseFloor*(1+(Ham-1)*5)
    """
    totalError = False
    # guess fit parameters
    pguess = []
    for i in range(len(Peaks)):
        peakWaveNo = Peaks[i].get_wavenumber()        
        pguess.append(FitData(wn, spectrum, peakWaveNo, peakFWHM)[1].max())   # Amp
        pguess.append(peakWaveNo)   # center
        pguess.append(peakFWHM) # FWHM
    # construct fit window
    window_center = np.mean(pguess[1::3])
    wnFit, spectrumFit = FitData(wn, spectrum, window_center, singlePeakMargin*len(Peaks)) # singlePeakMargin ist standard Breite für single Peak
    # estimate error
    sigmas = np.zeros(len(wnFit)) + guessNoise   # assume every data point has same errorbars
    # perform fit
    popt, pcov = optimize.curve_fit(f_fixoffSum, wnFit, spectrumFit, p0=pguess, sigma=sigmas, absolute_sigma=True)  # if absolute_sigma=True, dann werden Fehelr als abs Fehler gewertet und nicht rel Gewichtung
    errors = sp.sqrt(pcov.diagonal())
    wnFitted = np.linspace(wnFit[0], wnFit[-1], num=1000)
    fFitted = f_fixoffSum(wnFitted, *popt)
    # construct list of fitParas dictionaries    
    i=0    
    fitParasList = []
    for i in range(len(Peaks)):
        if totalError:
            # calculate total error of fit function evaluated at maximum of fit function
            err_fmax = GaussFitError(wn,popt[i*3:i*3+3],pcov[i*3:i*3+3])
            # use total fit error as estimate for Amp_err
            errors[i*3+0] = err_fmax
            AmpErrText = 'Error Amp tot'
            AreaErrText = 'Error Area tot'
        else:
            AmpErrText = 'Error Amp'
            AreaErrText = 'Error Area'
        # calculate peak area:
        area = popt[i*3+0]*popt[i*3+2]*np.sqrt(np.pi/np.log(16))
        err_area = abs(area)*np.sqrt((errors[i*3+0]/popt[i*3+0])**2 + (errors[i*3+2]/popt[i*3+2])**2) # Gaussian error propagation
        fitParas = {'Peak': Peaks[i].get_name() , \
                'Amp': round_sig(popt[i*3+0],3) , \
                AmpErrText: round_sig(errors[i*3+0],1) , \
                'Center': round(popt[i*3+1],3) , \
                'Error Center':round_sig(errors[i*3+1],1) , \
                'FWHM':round(popt[i*3+2],2) , \
                'Error FWHM':round_sig(errors[i*3+2],1) ,\
                'Area':round_sig(area,3) , \
                AreaErrText:round_sig(err_area,1)}
        print('multipeak_fit: ' + fitParas['Peak'],fitParas['Amp'],fitParas['Center'],fitParas['FWHM'],fitParas['Area'])
        fitParasList.append(fitParas)
    return wnFitted, fFitted, fitParasList
 
def MultiPeakFitFixed(Peaks, wn, spectrum, peakFWHM, singlePeakMargin, guessNoise):
    """
    Fits a sum of peaks with a sum of gaussians, having no vertical offset.
    Here, FWHM and peakCenter are kept fixed and only peak amplitudes are fitted. 
    peakFWHM : expected FWHM of individual deconvolved peaks
    singlePeakMargin : expected full width of individual deconvolved peak
    spectrum : should be ABSOLUTE value of spectrum
    guessNoise : typical a definition like this is choosen: noiseFloor*(1+(Ham-1)*5)
    """
    # guess fit parameters
    pguess = []
    fixedCenters = []
    fixedFWHMs = []
    for i in range(len(Peaks)):
        peakWaveNo = Peaks[i].get_wavenumber()        
        pguess.append(FitData(wn, spectrum, peakWaveNo, peakFWHM)[1].max())   # Amp
        fixedCenters.append(peakWaveNo)   # center
        fixedFWHMs.append(peakFWHM) # FWHM
    # construct fit window
    window_center = np.mean(fixedCenters[1::3])
    wnFit, spectrumFit = FitData(wn, spectrum, window_center, singlePeakMargin*len(Peaks)) # singlePeakMargin ist standard Breite für single Peak
    # estimate error
    sigmas = np.zeros(len(wnFit)) + guessNoise   # assume every data point has same errorbars
    # perform fit
    popt, pcov = optimize.curve_fit(lambda wn, *amps: f_SumFixedPara(wn, fixedCenters, fixedFWHMs, *amps), wnFit, spectrumFit, p0=pguess, sigma=sigmas, absolute_sigma=True)  # if absolute_sigma=True, dann werden Fehelr als abs Fehler gewertet und nicht rel Gewichtung
    errors = sp.sqrt(pcov.diagonal())
    # add fixed parameters to fit result
    amps = popt
    ampsErr = errors
    popt = np.zeros(3*len(amps))
    errors = np.zeros(3*len(amps))
    for i in range(len(amps)):
        popt[i*3+0]=amps[i]
        popt[i*3+1]=fixedCenters[i]
        popt[i*3+2]=fixedFWHMs[i]
        errors[i*3+0]=ampsErr[i]
        errors[i*3+1]=0.0
        errors[i*3+2]=0.0
    wnFitted = np.linspace(wnFit[0], wnFit[-1], num=1000)
    fFitted = f_fixoffSum(wnFitted, *popt)
    # construct list of fitParas dictionaries    
    i=0    
    fitParasList = []
    for i in range(len(Peaks)):
        AmpErrText = 'Error Amp'
        AreaErrText = 'Error Area'
        # calculate peak area:
        area = popt[i*3+0]*popt[i*3+2]*np.sqrt(np.pi/np.log(16))
        err_area = abs(area)*np.sqrt((errors[i*3+0]/popt[i*3+0])**2 + (errors[i*3+2]/popt[i*3+2])**2) # Gaussian error propagation
        fitParas = {'Peak': Peaks[i].get_name() , \
                'Amp': round_sig(popt[i*3+0],3) , \
                AmpErrText: round_sig(errors[i*3+0],1) , \
                'Center': round(popt[i*3+1],3) , \
                'Error Center':errors[i*3+1], \
                'FWHM':round(popt[i*3+2],2) , \
                'Error FWHM':errors[i*3+2] ,\
                'Area':round_sig(area,3) , \
                AreaErrText:round_sig(err_area,1)}
        print('multipeak_fit: ' + fitParas['Peak'],fitParas['Amp'],fitParas['Center'],fitParas['FWHM'],fitParas['Area'])
        fitParasList.append(fitParas)
    return wnFitted, fFitted, fitParasList
    

def PeakRawAnalysis(peakName, peakWaveNo, wn, spectrum, windowLength, noiseLevel):
    """
    Performes a crude analysis of a peak by evaluating the data at peakWaveNo
    +/- windowLength/2
    - Estimates peak centre as the balanced point
    - Amplitude as max value in data window
    - FWHM as the two values at half the Amplitude
    - Area as integral of data window
    - most errors are very crude approximations!
    windowLength : full width of fit window, usually windowLength=x*PeakWidth(T)
    spectrum : should be ABSOLUTE value of spectrum
    noiseLevel : typical a definition like this is choosen: noiseFloor*(1+(Ham-1)*5)
    """
    # construct fit window
    wnFit, spectrumFit = FitData(wn, spectrum, peakWaveNo, windowLength)
    # guess fit paras
    popt = GuessPeakParameter(wnFit, spectrumFit)
    popt[1] = np.average(wnFit,weights=spectrumFit/spectrumFit.max())   # center = gew. Mittel
    errors = np.zeros(len(popt))
    # Fehler auf Amplitude:
    # Fehler auf Amplitude ist 2*Amp falls nur ein Datenpunkt im Fenster
    # Fehler nimmt assymptotisch ab je mehr Datenpunkte im Fenster
    # Fehler konvergiert asyptotisch gegen Noisefloor
    errors[0] = 2*popt[0]/len(wnFit)+noiseLevel  # Grobe Abschätzung des Fehlers
    # Fehler auf Center = StdDev des Gewichtetes Mittel
    errors[1] = sp.sqrt(np.average((wnFit-popt[1])**2, weights=spectrumFit/spectrumFit.max()))
    # Fehler auf FWHM:
    errors[2] = errors[1]
    # calculate peak area:
    if wnFit[0] > wnFit[1]:
        X = wnFit[::-1]
        Y = spectrumFit[::-1]
    else:
        X = wnFit
        Y = spectrumFit
    # falls wnFit nicht exakt länge des PeakWindows hat, dann wird link und rechts noch ein ungerader wn-Punkt
    # ergänzt mit dem dann PeakWindow exakte Länge hat
    if peakWaveNo-0.5*windowLength < X[0]:
        X = np.insert(X,0,peakWaveNo-0.5*windowLength)  
        Y = np.insert(Y,0,0)
    if peakWaveNo+0.5*windowLength > X[-1]:
        X = np.append(X,peakWaveNo+0.5*windowLength)
        Y = np.append(Y,0)
    area = simps(Y,X)   # numeric integration using Simpson's rule
    # Fehler auf Fläche
    err_area = area*errors[0]/popt[0] # prozentualer Fehler der Amplitude
    print('peak_raw: ' + peakName,popt[0],popt[1],popt[2],area)
    fitParas = {'Peak': peakName , \
                'Amp': round_sig(popt[0],3) , \
                'Error Amp': round_sig(errors[0],1) , \
                'Center': round(popt[1],3) , \
                'Error Center':round_sig(errors[1],1) , \
                'FWHM':round(popt[2],2) , \
                'Error FWHM':round_sig(errors[2],1) ,\
                'Area':round_sig(area,3) , \
                'Error Area':round_sig(err_area,1)}
    return wnFit, spectrumFit, fitParas

def MonoCalibration(T, X, Y, Harmonics, paras, peakName, peakMargin, guessMonoCorr, noiseDef, zeroPaddingFactor, plotCal, plotColor, data_window):
    targetPeak = Peak.select_name(peakName)[0]
    peakWaveNo = targetPeak.get_wavenumber()
    harmonic = targetPeak.get_harmonic()
    append = False
    lw=2

    # generate calibration data
    index = Harmonics.index(harmonic)
    Tcalib, Xcalib, Ycalib = CutDataSet(T.copy(), X[index], Y[index], data_window) # for high resolution scans it is helpful to reduce resolution for mono peak fitting
    Zcalib = Xcalib + 1j*Ycalib
    # calculate DFt spectrum
    Zg = GaussWindow(Tcalib, Zcalib, False)
    wn, dft = DFT(Tcalib, Zg, paras, harmonic, zeroPaddingFactor)
    # guess fit parameters
    guess_peak_pos = peakWaveNo-harmonic*guessMonoCorr  # man muss also den Fehler des Monos berücksichtigen wenn Peak auch wirklich im Fit WIndow sein soll
    guess_noise = NoiseFloor(wn, abs(dft), Peak.select(harmonic = harmonic), PeakWidth(Tcalib,peakMargin,False, zeroPaddingFactor)[1], noiseDef)[1]  # take StdDev instead of Avg+StdDev, since Avg is just background signal but not Noise fluctuations
    # do Peak Fit with factor 1.5 enlarged fit window
    wnFitted, fFitted, fitParas = PeakFit(peakName, guess_peak_pos, wn, abs(dft), 2.5*PeakWidth(Tcalib, peakMargin, False, zeroPaddingFactor)[1], guess_noise)
    # calculate calibrated monochromator wavenumber
    mono_corr = (peakWaveNo - fitParas['Center'])/harmonic # in cm^-1, the corrected Mono is then : para_dict['wn_M']+mono_corr
    error_mono_corr = fitParas['Error Center']
    calibrated_mono = paras['Monochromator Wavenumber'] + mono_corr
    print '\n' + 'mono correction = ' + str(round(mono_corr, 3)) + ' cm^-1'
    print 'calibrated mono = ' + str(round(calibrated_mono, 3)) + ' +/- ' + str(round_sig(error_mono_corr,1)) + ' cm^-1'
    # create fit label
    fit_label = 'LS-fit:\n' + \
                'Amp= '+ str(round_sig(fitParas['Amp'], 2)) + r'$\pm$' + str(round_sig(fitParas['Error Amp'],1))+ '\n'+ \
                'FWHM= ' + str(round(fitParas['FWHM'], 3)) + r'$\pm$' + str(round_sig(fitParas['Error FWHM'],1)) + 'cm$^{-1}$\n' + \
                r'Mono$_{corr}$= ' + str(round(mono_corr, 3)) + r'$\pm$' + str(round_sig(error_mono_corr,1)) + r'cm$^{-1}$' +'\n'+\
                r'calibr. Mono= ' + str(round(calibrated_mono, 3)) + r'$\pm$' + str(round_sig(error_mono_corr,1)) + r'cm$^{-1}$'
    if plotCal:
        # plot calibration
        figName = 'MonoCal-' + paras['Title']
        plt.close(figName)  # close previous calibration figure
        fig = plt.figure(figName, figsize=(10,5))
        axCalib = fig.add_subplot(111)
        c = plotColor[harmonic-1]
        label = peakName + ' data'
        xlim = [min([wnFitted[0]-5, peakWaveNo-5]) , max([wnFitted[-1]+10.0, peakWaveNo+10])]
        axCalib.set_xlim(xlim)
        axCalib.plot(wn, abs(dft), '-', label=label, color = c, linewidth = lw)
        axCalib.plot(wnFitted, fFitted, '-', label= fit_label, color = 'k', linewidth = lw)
        axCalib.set_xlabel('$\\nu$ [cm$^{-1}$]')
        axCalib.set_ylabel('abs. amplitude [a.u.]')
        axCalib.set_title(figName, loc = 'center')
        axCalib.legend(loc='best')
    # write calibrated mono in rpt file
    if append:
        # write in rpt-file
        with open(para_dict['file path'][:-4] + ".rpt", "a") as myfile:
            myfile.write('\n\n' + 'Calibrated Mono /cm^-1 : ' + str(round(calibrated_mono, 3)) + \
            '\nCalibration Error /cm^-1 : ' + str(round_sig(error_mono_corr,1)) + \
            '\n' + 'Monochromator Corr /nm : ' + str(round(1e7*(1/calibrated_mono - 1/para_dict['wn_M']),3)) + ' nm')
    # return values
    return calibrated_mono, error_mono_corr, mono_corr


def SpectralPhaseOffset(T, X, Y, phaseStep, phaseAcc, Harmonics, paras, peakName, peakMargin, zeroPaddingFactor, plotCal, plotColor):
    """
    Find global phase correction (phase shift between 1H WPI signal and 1H Ref 
    signal) after delay-offset error has been corrected in DFT(). 
    Phase is found by symmetrize dispersion curve at the peak position of peakName.
    """
    # It would be more precise to fit theoretical dispersion curve according to
    # expected dft resolution and a global phase as free fit parameter. The negative
    # value of this fit parameter would then be the correction phase.
    
    # simulate signal for debugging
#    sstart = 0#-3.35
#    sstep = 53.0
#    sstop = 100000
#    paras[b'Schrittweite'] = sstep
#    T = np.linspace(sstart,int((sstop-sstart)/sstep)*sstep, num=int((sstop-sstart)/sstep), endpoint=True)
#    T += sstart
#    peakDict = {'peak1': [12985.19, 1.0, 0.0], 'peak2': [13042.90, 1.0, 0.0]} # [wn amp phase]
#    Z = SimWPISignal(T, peakDict, 1, paras)    

    # start algorithm
    targetPeak = Peak.select_name(peakName)[0]
    peakWaveNo = targetPeak.get_wavenumber()
    harmonic = targetPeak.get_harmonic()
    lw=2
    # generate calibration data
    index = Harmonics.index(harmonic)
#    Tcalib, Xcalib, Ycalib = CutDataSet(T.copy(), Z.real, Z.imag, [0.0,100000.0]) # for high resolution scans it is helpful to reduce resolution for mono peak fitting    
    Tcalib, Xcalib, Ycalib = CutDataSet(T.copy(), X[index], Y[index], [0.0,100000.0]) # for high resolution scans it is helpful to reduce resolution for mono peak fitting
    Zcalib = Xcalib + 1j*Ycalib
    peakWidth = 2.0*PeakWidth(Tcalib, peakMargin, False, zeroPaddingFactor)[1]  #calculate for abs SPectrum but use factor of 2.0    
    # calculate DFt spectrum
    Zg = GaussWindow(Tcalib, Zcalib, True)  # important that suscept=True
    wn, dft = DFT(Tcalib, Zg, paras, harmonic, zeroPaddingFactor)
    wnFit, dftFit = FitData(wn, dft, peakWaveNo, peakWidth)
    # find phase
    criterion = max(abs(dftFit))/phaseAcc
    phi = -180-phaseStep # Startphase für iteration
    dftFit *= np.exp(-1j*2*np.pi*phi/360.0)
    while (abs(max(dftFit.imag)+min(dftFit.imag))>criterion):
        dftFit *= np.exp(-1j*2*np.pi*phaseStep/360)
        phi += phaseStep
        if phi>180:
            print 'Phasenanpassung fehlgeschlagen'
            phi = 0
            break
    if dftFit.real[find_index(wnFit,peakWaveNo)]<0:
        phi += 180
    dft *= np.exp(-1j*2*np.pi*phi/360) 
    # create fit label
    label = 'Corr. phase /deg: ' + str(-phi)
    # plot calibration    
    if plotCal:
        c = plotColor[harmonic-1]
        xlim = [min([wnFit[0]-5, peakWaveNo-5]) , max([wnFit[-1]+5, peakWaveNo+5])]
        ylim = [-1.2, 1.2]
        figName = 'PhaseCal-' + paras['Title']
        plt.close(figName)  # close previous calibration figure
        fig = plt.figure(figName, figsize=(10,5))
        # plot imag
        axCalib = fig.add_subplot(121)
        axCalib.plot(wn, dft.imag/max(abs(dftFit.imag)), '-', label=label, color = c, linewidth = lw)
        axCalib.set_xlim(xlim)
        axCalib.set_ylim(ylim)
        axCalib.grid(which='both')
        axCalib.set_xlabel('$\\nu$ [cm$^{-1}$]')
        axCalib.set_title(figName, loc = 'left')
        axCalib.get_xaxis().get_major_formatter().set_useOffset(False)
        # plot real
        axCalib = fig.add_subplot(122)
        axCalib.plot(wn, dft.real/max(abs(dftFit.real)), '-', label=label, color = c, linewidth = lw)
        axCalib.set_xlim(xlim)
        axCalib.set_ylim(ylim)
        axCalib.grid(which='both')
        axCalib.set_xlabel('$\\nu$ [cm$^{-1}$]')
        axCalib.legend(loc='best')
        axCalib.get_xaxis().get_major_formatter().set_useOffset(False)
        axCalib.yaxis.set_ticklabels([])
    # return values
    return -phi #correction phase, i.e. data was shifted by +phi, in degree



def Interpolate(X, Y, pos):
    """
    The Y-value at the specified X-position (pos) is reurned. This amplitude is 
    determined by linear interpolation of the Y-data.
    """    
    index = find_index(X,pos)
    if X[index] < pos:
        index1=index
        index2=index+1
    else:
        index1=index-1
        index2=index
    m = (Y[index2]-Y[index1])/(X[index2]-X[index1])   #Steigung der Geraden
    return m*(pos-X[index1])+Y[index1]    # uebergirbt Interpolation der PSD an Pos des Peaks
    
def SpectralAmplitude(Peaks, wnPSD, PSD):
    """
    Calculates spectral amplitudes at all peak positions. Spectral amplitudes 
    are normalized to the PSD peak amplitude. Returns dictionary of spectral 
    amplitudes. 
    """    
    spectralAmplitudesDict = {}    
    for peak in Peaks:
        peakName = peak.get_name()        
        peakWaveNo = peak.get_wavenumber()     
        spectralAmplitudesDict[peakName] = Interpolate(wnPSD, PSD, peakWaveNo)/PSD.max()
    return spectralAmplitudesDict
    
def AnalyseLaserSpectrum(X,Y):
    """
    Calculates moments of laser spectrum. Background has to be subtracted.
    Problem: calculated values are extremely sensitive to residual background.
    Large errors are already obtained for residual background >0.1%
    """    
    # Normierungsfaktor des Spektrums    
    Y /= sum(Y)
    # Erwartungswert
    E = sum(X[i]*Y[i] for i in range(len(Y)))
    # Varianz
    Var = sum(X[i]**2*Y[i] for i in range(len(Y))) - E**2
    # Schiefe
    v = ( sum(X[i]**3*Y[i] for i in range(len(Y))) - 3*Var*E - E**3 )/Var**(1.5)
    # calculating parameters under assumption of a Gaussian distribution
    center = E
    FWHM = np.sqrt(8*np.log(2)*Var)
    skewness = v    # v<0 : Übergewicht rechts (>0 links), v=0 : symmetrisch
    return center, FWHM, skewness
    
def FitLaserSpectrum(Wn, PSD, SpectrumPath, plotFit):
    fileName = GetFileDictName(SpectrumPath)[1]    # extracts filename from file path, i.e. removes directory and data type
    peakName = 'laser'
    lw=2
   # guess fit parameters
    guess_peak_pos = Wn[np.argmin(np.abs(PSD-PSD.max()))]
    peakMargin = 0.02
    peakWidth = abs(Wn[findPeakIndexes(PSD, peakMargin)[1]]-Wn[findPeakIndexes(PSD, peakMargin)[2]])
    # guess noise
    indexPeakLow = find_index(Wn, guess_peak_pos-4*peakWidth)   # will return the last index of wn if peak not in wn
    indexPeakHigh = find_index(Wn, guess_peak_pos+4*peakWidth)
    noiseVector = np.delete(PSD,range(indexPeakLow,indexPeakHigh+1))   # remove peak values from spectrum    
    guess_noise = np.mean(noiseVector) + 3*np.std(noiseVector)
    # do Peak Fit with factor 1.5 enlarged fit window
    wnFitted, fFitted, fitParas = PeakFit(peakName, guess_peak_pos, Wn[::-1], PSD[::-1], 1.5*peakWidth, guess_noise)
    # add peak position to fitParas
    fitParas['Peak'] = guess_peak_pos
    # create fit label
    fit_label = 'LS-fit:\n' + \
                'Peak= ' + str(round(fitParas['Peak'], 3)) + 'cm$^{-1}$'+ '\n'+ \
                'Center= '+ str(round(fitParas['Center'], 3)) + r'$\pm$' + str(round_sig(fitParas['Error Center'],1)) + 'cm$^{-1}$'+ '\n'+ \
                'FWHM= ' + str(round(fitParas['FWHM'], 3)) + r'$\pm$' + str(round_sig(fitParas['Error FWHM'],1)) + 'cm$^{-1}$'
    if plotFit:
        # plot laser fit
        figName = 'laserFit-' + fileName
        plt.close(figName)  # close previous calibration figure
        figLaser = plt.figure(figName, figsize=(10,5))
        axCalib = figLaser.add_subplot(111)
        c = 'r'
        label = peakName + ' data'
        xlim = [fitParas['Center']-3*fitParas['FWHM'] , fitParas['Center']+10*fitParas['FWHM']]
        axCalib.set_xlim(xlim)
        axCalib.plot(Wn, PSD, '-', label=label, color = c, linewidth = lw)
        axCalib.plot(wnFitted, fFitted, '-', label= fit_label, color = 'k', linewidth = lw)
        axCalib.set_xlabel('$\\nu$ [cm$^{-1}$]')
        axCalib.set_ylabel('abs. amplitude [a.u.]')
        axCalib.set_title(figName, loc = 'center')
        axCalib.legend(loc='best')
    return wnFitted, fFitted, fitParas
    


# plot data====================================================================

def PlotData(ax, X, Y, color, xlim, ylim, legend='', legendLoc='best', legendFrame=True, plainLegend=False, annotatePeaks=[], annotateColor='k'):  
    lw = 2.5 # linewidth in DFTs
    fs = 22 # font size for annotations
    texty = 0.99 # in %
    textx = 0.7 # in cm^-1
    ax.plot(X, Y, '-', label= legend, color= color, linewidth = lw)
    if len(xlim) >0:    
        ax.set_xlim(xlim)
    if len(ylim) >0:
        ax.set_ylim(ylim)
    # annotate peaks
    for peak in annotatePeaks:
        if peak.get_annotate():          
            ax.axvline(peak.get_wavenumber(), linestyle='--', color = peak.get_color(), linewidth=1.8, zorder=0)
            ax.annotate(peak.get_name() , (peak.get_wavenumber()+textx, Y.max()*texty), color = annotateColor, backgroundcolor='w' , va='top', fontsize = fs, zorder=0)    
    # legend
    # define loc for plainLegend
    xoff = 0.02
    yoff= 0.04
    loc = {'upper left' : [xoff, 1-yoff, 'left','top'],'upper right' : [1-xoff, 1-yoff,'right','top'],\
            'lower left' : [xoff, yoff,'left','bottom'], 'lower right' : [1-xoff, yoff,'right','bottom'],\
            'upper center' : [0.5, 1-yoff,'center','top'],'lower center' : [0.5, yoff,'center','bottom'],\
            'best' : [0.5, 0.5,'center','center']}
    if len(legend)>0:
        if plainLegend:
            ax.text(loc[legendLoc][0], loc[legendLoc][1], legend, transform=ax.transAxes, fontweight='bold', ha=loc[legendLoc][2], va=loc[legendLoc][3], fontsize = fs, color=color)#,bbox=dict(fc='w',ec='w'))
        else:
            ax.legend(loc=legendLoc, fontsize=18, frameon=legendFrame,numpoints=1)

def FigureLayout(size,numSubPlots,title,xlabel,ylabel,majorLocX=None,majorLocY=None,layout='vertical',tightLayout=False,commonAxis=False):
    """
    size = (length,height)
    majorLoc/minorLoc : defines increment of tick labels. If minorLoc=None the 
    increment for minor Ticks will be half of major ticks
    layout= 'vertial' or 'horizontal' defines arrangement of subplots, vertial is default 
    (horizontal ist not implemented yet)
    tightLayout : if True will merge subplots to common axis
    """  
    lw = 2 # linewidth
    fs = 20 # font size 
    majorLen = 7
    minorLen = 5
    paddingX = 5    # distence between axes and tick labels
    paddingY = 6    
    fig = plt.figure(title, figsize=(size)) # create figure canvas
    ax = [None]*numSubPlots     # create empty array of axes
    # make subplots for vertical/horizontal alignment
    if layout == 'vertical' or layout == None:    # vertical alignment
        for i in range(numSubPlots):
            ax[i] = fig.add_subplot(numSubPlots,1,numSubPlots-i) # define axes vertically aligned
            # change axes layout
            ax[i].get_yaxis().get_major_formatter().set_powerlimits((0,0))    # use always scientific format
            ax[i].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False)) # change style of scientifc offset
            ax[i].get_xaxis().get_major_formatter().set_useOffset(False)   # never use offset notation
            # format xaxis
            ax[i].get_xaxis().set_tick_params(which='major', pad=paddingX, length=majorLen, width=lw)
            ax[i].get_xaxis().set_tick_params(which='minor', pad=paddingX, length=minorLen, width=lw)
            if majorLocX != None:    
                majxLoc = plt.MultipleLocator(majorLocX[i])
                minxLoc = plt.MultipleLocator(majorLocX[i]/2)
                ax[i].xaxis.set_major_locator(majxLoc)
                ax[i].xaxis.set_minor_locator(minxLoc)
            # format yaxis
            ax[i].get_yaxis().set_tick_params(which='major', pad=paddingY, length=majorLen, width=lw)
            ax[i].get_yaxis().set_tick_params(which='minor', pad=paddingY, length=minorLen, width=lw)   
            if majorLocY != None:  
                majyLoc = plt.MultipleLocator(majorLocY[i])
                minyLoc = plt.MultipleLocator(majorLocY[i]/2)
                ax[i].yaxis.set_major_locator(majyLoc)
                ax[i].yaxis.set_minor_locator(minyLoc)
            # draw thick frame
            for axis in ['top','bottom','left','right']:
              ax[i].spines[axis].set_linewidth(lw) 
            # label fontsize
            ax[i].tick_params(labelsize=fs)
            ax[i].yaxis.get_offset_text().set_fontsize(fs) # font size of scientific format offset 
            # hide first x tick because of common 0 at corner
            #xticks = ax[i].xaxis.get_major_ticks()
            #xticks[1].label1.set_visible(False)
            # remove common x labels if wished
            if commonAxis == True:
                if i>0: ax[i].xaxis.set_ticklabels([])            
            # set title            
            if i==numSubPlots-1: ax[i].set_title(title)              
            # set x label            
            if i==0:
                ax[i].set_xlabel(xlabel, fontsize = fs)  
            # set y label in center of subplots
            if len(ylabel)>0:
                if numSubPlots % 2 == 0:
                    if i == numSubPlots/2-1 : 
                        ax[i].get_yaxis().set_label_coords(-0.1,1.0)
                        ax[i].set_ylabel(ylabel, fontsize = fs) 
                elif i == int(numSubPlots/2):
                    ax[i].get_yaxis().set_label_coords(-0.1,0.5)
                    ax[i].set_ylabel(ylabel, fontsize = fs) 
    return fig, ax
 

class Peak:
    """ Class that represents peaks as objects and stores them in a list.

    Attributes:
    name       -- name of the peak in latex notation (string)
    wavenumber -- wavenumber of the peak in cm^{-1} (float)
    annotate   -- can be 'True' or 'False' depending on whether the peak show be annotated in the plot
    analysis   -- determines the analysis type, can be
                    * 'fit': the peak will be fitted
                    * 'raw': the values of the array in of the peak will be summed up
                    * 'None': no analysis of the peak will happend
                  (string)
    color      -- contains the color of the shown fit or annotation in the plot (color notation that matplotlib can work with)
    window     -- contains the length of the window that is used for the analysis in cm^{-1} (distance to the upper and lower border)

    Methods:
    |Class:
    ||add_line       -- adds one peak to the peak-list
    ||add_lines      -- adds multiple peaks to the peak-list.
    ||select         -- returns a list of peak objects with specific attributes
    ||remove_line    -- removes and deletes an peak-object from the list
                        Be carfull this method removes every line that fits the given parameter
    |Object:
    ||get_name       -- returns the name of the peak (string)
    ||get_wavenumber -- returns the wavenumber of the peak in cm^{-1} (float)
    ||get_annotate   -- returns the whether the peak gets an annotation in the plot ('True' or 'False')
    ||get_analysis   -- returns type of analysis ('fit','raw' or 'None')
    ||get_color      -- returns the color used in the plot of the peak (matplotlib notation)
    ||get_window     -- returns analysis 'radius' in cm^{-1} (float)
    ||get_harmonic   -- returns the harmonic of the peak (int)
    ||set_name       -- sets the name of the peak (string)
    ||set_wavenumber -- sets the wavenumber of the peak in cm^{-1} (float)
    ||set_annotate   -- sets the whether the peak gets an annotation in the plot ('True' or 'False')
    ||set_analysis   -- sets type of analysis ('fit','raw' or 'None')
    ||set_color      -- sets the color used in the plot of the peak (matplotlib notation)
    ||set_window     -- sets analysis 'radius' in cm^{-1} (float)
    ||set_harmonic   -- sets the harmonic of the peak (int)


    Examples:
    Add peaks to the list:
    >>> Peak.add_line(name='Userpeak',wavenumber=12345,color='k',annotate=False,analysis='raw',window=2.4,harmonic=1)
    >>> Peak.add_lines(wavenumbers().one_photon(),color='k',annotate=False,analysis='fit',window=2.4,harmonic=1)

    Select peaks with the attribute harmonic of 1 and gives back all the names
    >>> for peak in Peak.select(harmonic = 1):
    ...    print peak.get_name()
    Userpeak
    Dd
    D$_1$

    Change the color a peaks of first harmonic that get fitted to green ('g')
    >>> peak = Peak.select(harmonic = 1, analysis = 'fit')
    >>> peak[1].set_color('g')
    >>> for peak in Peak.select(harmonic = 1, color = 'k'):
    ...    print peak.get_name() + ' : ' + peak.get_color()
    Userpeak : k
    D$_2$ : k

    """
    __author__  = "Daniel Uhl"
    __email__   = "daniel.uhl@mars.uni-freiburg.de"
    __version__ = "0.8"
    __status__  = "Prototype"

    peaks = []
    def __init__(self,name,wavenumber,color,annotate,analysis,window,harmonic):
        if (analysis != 'raw' and analysis != 'none' and analysis != 'fit'):
            print "Analysis has to be \'none\', \'raw\' or \'fit\'. Your entry was ", analysis ,"."
        self.__name = name
        self.__wavenumber = wavenumber
        self.__color = color
        self.__annotate = annotate
        self.__analysis = analysis
        self.__window = window
        self.__harmonic = harmonic

    def __repr__(self):
        return "<Peak %s>" % repr(self.properties)

    @staticmethod
    def add_lines(d, **kwds):
        for name in d:
            Peak.add_line(name = name,wavenumber = d[name], **kwds)

    @staticmethod
    def add_line(**kwds):
        Peak.peaks.append(Peak(**kwds))

    @staticmethod
    def remove_line(name=None,wavenumber=None,color=None,annotate=None,
                                analysis=None,window=None,harmonic=None):
        savety = False
        for attr in [name,wavenumber,color,annotate,analysis,window,harmonic]:
            savety = attr != None or savety
        if savety:
            result =Peak.select(name,wavenumber,color,annotate,analysis,window,harmonic)
            for unneeded in result:
                Peak.peaks.remove(unneeded)
            print str(len(result)) + 'line(s) has been deleted.'
        else:
            print 'No line has been deleted.'

    @staticmethod
    def select(name=None,wavenumber=None,color=None,annotate=None,
                    analysis=None,window=None,harmonic=None):
        result = list(Peak.peaks)
        if name != None:
            result = filter(lambda x: x.get_name()==name, result)
        if wavenumber != None:
            result = filter(lambda x: x.get_wavenumber()==wavenumber, result)
        if color != None:
            result = filter(lambda x: x.get_color()==color, result)
        if window != None:
            result = filter(lambda x: x.get_window()==window, result)
        if analysis != None:
            result = filter(lambda x: x.get_analysis()==analysis, result)
        if annotate != None:
            result = filter(lambda x: x.get_annotate()==annotate, result)
        if harmonic != None:
            result = filter(lambda x: x.get_harmonic()==harmonic, result)
        return result

    @staticmethod
    def system_select(key):
        result = []
        for p in Peak.peaks:
            if key(p): result.append(p)
        return result

    @staticmethod
    def select_harmonic(value):
        return Peak.system_select(lambda x: x.get_harmonic() == value)

    @staticmethod
    def select_name(value):
        return Peak.system_select(lambda x: x.get_name() == value)

    @staticmethod
    def select_wavenumber(value):
        return Peak.system_select(lambda x: x.get_wavenumber() == value)

    @staticmethod
    def select_color(value):
        return Peak.system_select(lambda x: x.get_color() == value)

    @staticmethod
    def select_window(value):
        return Peak.system_select(lambda x: x.get_window() == value)

    @staticmethod
    def select_analysis(value):
        return Peak.system_select(lambda x: x.get_analysis() == value)

    @staticmethod
    def select_annotate(value):
        return Peak.system_select(lambda x: x.get_annotate() == value)

    def get_name(self):
        return self.__name

    def get_wavenumber(self):
        return self.__wavenumber

    def get_annotate(self):
        return self.__annotate

    def get_analysis(self):
        return self.__analysis

    def get_color(self):
        return self.__color

    def get_window(self):
        return self.__window

    def get_harmonic(self):
        return self.__harmonic

    def get_info(self):
        return (self.__name,self.__wavenumber,self.__annotate,self.__analysis,
                self.__color,self.__window,self.__harmonic)

    def set_name(self,name):
        self.__name = name

    def set_wavenumber(self,wavenumber):
        self.__wavenumber = wavenumber

    def set_annotate(self,annotate):
        self.__annotate = annotate

    def set_analysis(self,analysis):
        self.__analysis = analysis

    def set_color(self,color):
        self.__color = color

    def set_window(self,window):
        self.__window = window

    def set_harmonic(self,harmonic):
        self.__harmonic = harmonic