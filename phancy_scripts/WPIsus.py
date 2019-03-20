 # -*- coding: utf-8 -*-
"""
Program to evaluate susceptibility of WPI data:
imports WPI data
performs phase correction
calculates DFTs
determines susceptibility
plots time domain and frequency domain data
optional: plots laser spectrum of 1st harmonic
for peak analysis: use WPIfit.py

@author: Luki
"""

from inputs import *
import sys
sys.path.append('//10.5.71.1/user/Projekte/2D-Spectro/Data Analysis Stuff/python/Programs')
from functions import *

plt.close('all')
scriptAcronym = 'WPIsus'
suscept = True # specifies wether susceptibility or absolute value of spectra
data_window = [t-tCorr for t in data_window]     # To obtain T as exactly specified in data_window, data_window has to be shifted by tCorr, since later in the code T is shifted by tCorr

#Script Start==================================================================

# import WPI data
#FilePath = OpenDialog(initDir, 'Open WPI data set')
#FilePath = u'//10.5.71.1/user/Projekte/2D-Spectro/Data/WPI-SpectroCell_Legend/2016-07-12/2016-07-12_vapor_778_1-4H_8.dat'
FilePath =u'//10.5.71.1/user/Projekte/2D-Spectro/Data/Multi-QCs/2016-06-06/2016-06-06_vapor_778nm_LD_2.dat'
T, X, Y, paras = ImportWPIData(FilePath, len(Harmonics), guessMonoCorr, refLD, data_window)

# add spectral peaks to peak register
Peak.add_lines(Rb1HLines,color='k',annotate=True,analysis='none',window=0,harmonic=1)
Peak.add_lines(Rb2HLines,color='k',annotate=True,analysis='none',window=0,harmonic=2)
Peak.add_lines(Rb3HLines,color='k',annotate=True,analysis='none',window=0,harmonic=3)
Peak.add_lines(Rb4HLines,color='k',annotate=True,analysis='none',window=0,harmonic=4)

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
calibrated_mono, error_calibration, mono_correction = MonoCalibration(T, X, Y, Harmonics, paras, calibPeak, peakMargin, guessMonoCorr, noiseDef, zeroPaddingFactor, plotCal, plotColor)        
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

# find correction phase for 1H data
T += tCorr
phi1H = SpectralPhaseOffset(T, X, Y, phaseStep, acc, Harmonics, paras, calibPeak, peakMargin, zeroPaddingFactor, plotCal, plotColor)# returns phase for 1H (in degree) 

# make figures and axis
if plotTD: figTD, axTD = FigureLayout((10,10),len(Harmonics),scriptAcronym+'_'+paras['Title']+'_TD','pump-probe delay [ps]','amplitude [arb. u.]',majorLocX=majorLocXTD,majorLocY=majorLocYTD,commonAxis=True)
if plotFD: figFD, axFDRe = FigureLayout((10,10),len(Harmonics),scriptAcronym+'_'+paras['Title']+'_FDRe','wavenumber [cm$^{-1}$]','amplitude [arb. u.]',majorLocX=majorLocXFD,majorLocY=majorLocYFD)
if plotFD: figFD, axFDIm = FigureLayout((10,10),len(Harmonics),scriptAcronym+'_'+paras['Title']+'_FDIm','wavenumber [cm$^{-1}$]','amplitude [arb. u.]',majorLocX=majorLocXFD,majorLocY=majorLocYFD)
# looping through harmonics==================================================== 
for i in range(len(Harmonics)):  
    # select current data set    
    Z = (X[i] + 1j*Y[i])
    # correct frequency function of electronics
    if Harmonics[i] > 1:    # not first Harmonic, relevant for frequency function
        Z *= AF[i-1]
    # adjust phase if susceptibility is evaluated
    if Harmonics[i]==1:    
        phi = phi1H
    else: # use 1H phase to determine higher harmonic phases
        phi = phi1H+PF[i-1]   # correcting for phase funciton of electronics (in degree)
    phi += phiCorr[i]   # user correction (in degree) for all harmonics (including 1H)
    print str(i+1) + 'H correction phase (degree): ' + str(phi)
    if report: ReportLine(paras['File Path'], 'corrPhase'+str(Harmonics[i])+'H', phi)
    phi = 2*np.pi*phi/360   # conversion to radians
    Z *= np.exp(1j*phi) # correction of phase
   
# calculate DFT and properties of spectrum=====================================
    Zg = GaussWindow(T, Z, suscept)
    wn, dft = DFT(T, Zg, paras, Harmonics[i], zeroPaddingFactor)
    # normalize dft according to delay, harmonic...
    dft /= normFactor

# plot data ===================================================================
    # plot FD data
    if plotFD:    
        # plot FD data
        PlotData(axFDRe[i], wn, dft.real, plotColor[i], [xlim[2*i], xlim[2*i+1]], [], legend=legend[i], legendLoc=legendLoc, plainLegend=plainLegend, annotatePeaks=Peak.select(harmonic = Harmonics[i]), annotateColor=annotateColor)
#        axFDRe[i].axvline(paras['Laser Wavenumber']*Harmonics[i], linestyle='-', color = '#888888', linewidth=lw, zorder=0)
        PlotData(axFDIm[i], wn, dft.real, plotColor[i], [xlim[2*i], xlim[2*i+1]], [], legend=legend[i], legendLoc=legendLoc, plainLegend=plainLegend, annotatePeaks=Peak.select(harmonic = Harmonics[i]), annotateColor=annotateColor)
#        axFDIm[i].axvline(paras['Laser Wavenumber']*Harmonics[i], linestyle='-', color = '#888888', linewidth=lw, zorder=0)
#        axFDIm.grid(which='both')
    if plotTD:
        # plot TD data
        PlotData(axTD[i], T, abs(Z), plotColor[i], [], [], legend=legend[i], legendLoc=legendLoc, plainLegend=plainLegend)