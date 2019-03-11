 # -*- coding: utf-8 -*-
"""
Program to evaluate WPI data:
imports WPI data
calculates DFTs
plots time domain and frequency domain data
performs analysis of peaks and noise floor
cannot handle overlapping peaks
creats report file with peak analysis
optional: plots laser spectrum of 1st harmonic
for susceptibility: use WPIsus.py

WPImultifit : performs multipeak fits of overlapping peaks

@author: Luki
"""
from inputs import *
import sys
sys.path.append('//10.5.71.1/user/Projekte/2D-Spectro/Data Analysis Stuff/python/Programs')
from functions import *

plt.close('all')
scriptAcronym = 'WPIfit'
suscept = False # specifies wether susceptibility or absolute value of spectra

#Script Start==================================================================

# import WPI data
FilePath = OpenDialog(initDir, 'Open WPI data set')
#FilePath = u'//10.5.71.1/user/Projekte/2D-Spectro/Data/WPI-SpectroCell/2017-12-07/8kHz_780nm_782mono.dat'
#FilePath =u'//10.5.71.1/user/Projekte/2D-Spectro/Data/Multi-QCs/2016-06-06/2016-06-06_vapor_778nm_LD_2.dat'
#FilePath = '\\\\10.5.71.1\\user\\Projekte\\2D-Spectro\\Data\\Multi-QCs\\2016-05-19\\2016-05-19_vapor_769nm_81525A_11.dat'
T, X, Y, paras = ImportWPIData(FilePath, len(Harmonics), guessMonoCorr, refLD, data_window)

# add spectral peaks to peak register
# here only peaks are will be fitted that are not overlaping with other peaks
Peak.add_lines(Rb1HLines,color='k',annotate=True,analysis='fit',window=0,harmonic=1)
Peak.add_lines(Rb2HLines,color='k',annotate=True,analysis='fit',window=0,harmonic=2)
Peak.add_lines(Rb3HLines,color='k',annotate=True,analysis='fit',window=0,harmonic=3)
Peak.add_lines(Rb4HLines,color='k',annotate=True,analysis='fit',window=0,harmonic=4)

# initialize report
reportDict = {
    'dataWindow' : [x/1000.0 for x in data_window] ,\
    'zeroPaddingFactor' : zeroPaddingFactor ,\
    'windowFunction' : 'Gaussian' ,\
    'normalize' : normalize ,\
    'AF' : AF ,\
    }
if report: InitReport(paras['File Path'],'Peak Analysis',scriptAcronym, reportDict)

# calibrate monochromator and add mono peaks to peak register
#calibrated_mono, error_calibration, mono_correction = MonoCalibration(T, X, Y, Harmonics, paras, calibPeak, peakMargin, guessMonoCorr, noiseDef, zeroPaddingFactor, plotCal, plotColor, data_window)        
calibrated_mono = 12719#37594 #12721 #12719#
paras['Monochromator Wavenumber'] = calibrated_mono
if report: ReportLine(paras['File Path'], 'calibration', [calibPeak, calibrated_mono, error_calibration, mono_correction])   #log
# add mono peaks to peak register
for h in Harmonics:
    Peak.add_line(name=str(h)+monoName,wavenumber=h*calibrated_mono,color='#777777',annotate=True,analysis='none',window=2.4,harmonic=h)

# add leak peaks to peak register
leakingHarmonic = 1 # consider only 1H peaks causing leaking into higher harmonics
leakingPeaks = Peak.select(harmonic=leakingHarmonic)    # list of all peaks that generate leaking
AddLeakPeaks(leakingPeaks, Peak, paras, Harmonics, True, leakColor, 'none', monoName)


# main loop starts ============================================================  
# calculate peak margin
peakFWHM, peakWidth = PeakWidth(T, peakMargin, suscept, zeroPaddingFactor) # keep same also for suscept=True since then peaks have broad tails
if report: ReportLine(paras['File Path'],'resolution',[peakMargin, peakFWHM, peakWidth]) #log

# calculate normalization factor for spectral amplitudes
normFactor = Normalize(normalize, T, X, Y, paras, zeroPaddingFactor, Harmonics)

# initialize log variables (necessary for ordering in log file)
reportNoise=[]
reportNoiseExcluded=[]
reportPeakList=[]

# make figures and axis
if plotTD: figTD, axTD = FigureLayout((10,10),len(Harmonics),scriptAcronym+'_'+paras['Title']+'_TD','pump-probe delay [ps]','amplitude [arb. u.]',majorLocX=majorLocXTD,majorLocY=majorLocYTD,commonAxis=True)
if plotFD: figFD, axFD = FigureLayout((10,10),len(Harmonics),scriptAcronym+'_'+paras['Title']+'_FDRe','wavenumber [cm$^{-1}$]','amplitude [arb. u.]',majorLocX=majorLocXFD,majorLocY=majorLocYFD)

# looping through harmonics====================================================    
for i in range(len(Harmonics)):  
    # select current data set 
    Z = (X[i] + 1j*Y[i])
    # correct frequency function of electronics
    if Harmonics[i] > 1:    # not first Harmonic, relevant for frequency function
        Z *= AF[i-1]
   
# calculate DFT and properties of spectrum=====================================
    Zg = GaussWindow(T, Z, suscept)
    wn, dft = DFT(T, Zg, paras, Harmonics[i], zeroPaddingFactor)
    dft = abs(dft)  # for suscept use WPIsus.py
    # normalize dft according to delay, harmonic...
    dft /= normFactor
#    if i==0:
#        for n in range(len(wn)):
#            print str(wn[n]) + '\t' + str(dft[n])

# plot data ===================================================================  
    if plotFD:    
        PlotData(axFD[i], wn, dft, plotColor[i], [xlim[2*i], xlim[2*i+1]], [], legend=legend[i], legendLoc=legendLoc, plainLegend=plainLegend, annotatePeaks=Peak.select(harmonic = Harmonics[i]), annotateColor=annotateColor)
    if plotTD:
        PlotData(axTD[i], T, Z.real, plotColor[i], [], [], legend=legend[i], legendLoc=legendLoc, plainLegend=plainLegend)
 
# analyse data set ============================================================  
    # analyse laser spectrum
    if loadSpectrum and Harmonics[i]==1:    # if spectrum für SHG then use here ==2
        # import laser spectrum and apply boxcar smoothing       
        SpectrumPath = OpenDialog(GetFileDictName(FilePath)[0], 'Open laser spectrum')
        wlLaserSpectrum, wnLaserSpectrum, PSDLaserSpectrum = ImportLaserSpectrum(SpectrumPath, boxcarWindow, boxcarLength)
        # normalize spectrum to WPI data
        PSDLaserSpectrum = PSDLaserSpectrum/PSDLaserSpectrum.max()*abs(dft).max()  
        # calculate parameters of laser spectrum
        wnLaserFitted, PSDLaserFitted, laserFitParas = FitLaserSpectrum(wnLaserSpectrum, PSDLaserSpectrum, SpectrumPath, plotLaserFit)
        # calculate spectral amplitudes at peak wavenumbers as a fraction of peak amplitude of laser spectrum
        spectralAmplitudesDict = SpectralAmplitude(Peak.select(harmonic = Harmonics[i]), wnLaserSpectrum, PSDLaserSpectrum) # determine spectral amlitude with raw data
        #spectralAmplitudesDict = SpectralAmplitude(Peak.select(harmonic = Harmonics[i]), wnLaserFitted, PSDLaserFitted) # determine spectral amlitude with fitted data
        # plot spectum
        if plotFD: axFD[i].fill_between(wnLaserSpectrum, PSDLaserSpectrum, facecolor='#BBBBBB', label = "PSD laser", alpha=0.5)
        # log values
        if report: ReportLine(paras['File Path'],'laserSpectrum',[GetFileDictName(SpectrumPath)[1], boxcarWindow, boxcarLength, laserFitParas['Peak'], laserFitParas['Center'],laserFitParas['FWHM']]) #log              
        if report: ReportLine(paras['File Path'],'PSDamps', Dict2List(spectralAmplitudesDict))
        
    # determine noise floor of n'th harmonic, all peaks including Mono and Leaks are excluded for noise determination
    excludedNoisePeaks = Peak.select(harmonic = Harmonics[i])
    avgNoise, stdNoise, noiseLevel, wnNoise, noiseVector = NoiseFloor(wn, dft, excludedNoisePeaks, peakWidth, noiseDef)
    # plot noise level
    if plotFD: axFD[i].axhline(y=noiseLevel, linestyle='--', color= 'k', linewidth = lw)
    # add to log message
    reportNoise.append([avgNoise, stdNoise, noiseLevel])
    print "Noise level of " +str(Harmonics[i])+'H = ' + str(noiseLevel)
    names = []   
    for peak in excludedNoisePeaks:
        names.append(peak.get_name())
    reportNoiseExcluded.append(names)

    # analyse peaks
    for peak in Peak.select(harmonic = Harmonics[i]):
        peakAnalysis = peak.get_analysis()
        peakName = peak.get_name()
        peakWaveNo = peak.get_wavenumber()
#        if peakAnalysis == 'fit':
#            # do Peak Fit
#            wnFitted, fFitted, fitParas = PeakFit(peakName, peakWaveNo, wn, dft, peakWidth, noiseLevel)
#            if plotFD: axFD[i].plot(wnFitted, fFitted, '-', color=annotateColor, linewidth=lw)
#            #if plotFD: axFD[i].axvspan(wnFitted[0],wnFitted[-1], color='#999999', zorder=0)
#            fitParas['type']='fit'            
#            reportPeakList.append(fitParas)    # for log file
#        if peakAnalysis == 'raw':
#            wnFitted, fFitted, fitParas = PeakRawAnalysis(peakName, peakWaveNo, wn, dft, peakWidth, noiseLevel)
#            if plotFD: axFD[i].plot(wnFitted, fFitted, 'o', color=annotateColor, markersize=2*lw)    
#            if plotFD: axFD[i].axvspan(wnFitted[0],wnFitted[-1], color='#999999', zorder=0)
#            fitParas['type']='raw' 
#            reportPeakList.append(fitParas) # for log file
## finish log
#if report:
#    ReportEmptyLine(paras['File Path'])
#    ReportComment(paras['File Path'], 'Noise')
#    ReportComment(paras['File Path'], 'Noise list [avgNoise, stdNoise, noiseLevel]')
#    ReportLine(paras['File Path'], 'noiseDef', noiseDef)
#    ReportNoise(paras['File Path'], Harmonics,reportNoise, reportNoiseExcluded)
#    ReportPeaks(paras['File Path'], Harmonics, reportPeakList)
