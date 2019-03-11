# -*- coding: utf-8 -*-

# user inputs =================================================================


# Data Analysis ===============================================================
data_window = [-500 , 500]   # in fs
peakMargin = 0.01   # factor to which peak amplitude has to be dropped for peak width calculation
normalize =  'none' #'delay'   # none : keine Normalisierung, delay : 1ps entspricht Faktor 1, 1H : 1'te Harmonische wird auf 1 normiert und alle anderen Harmonischen relativ dazu 
calibPeak = 'Rb_D2'   # peak used for mono calibration
guessMonoCorr = 0  # nur nötig wenn Calibration sonst nicht Peak findet
refLD = True #True  # wurde laserdiode als ref benutzt?
Harmonics = [1,2,3]   # define which Harmonic is applied to which input channel, list can be shorter or longer, but no longer than 4 since only 4 input channels exist at ADC board
#AF = [1.02,1.05,1.09] # Amplituden Fkt PD M1599 @PM=1kHz und 1M Wiederstand, Werte für 2H und 4H wurden interpoliert da nicht gemessen
#AF = [1.13,1.27,1.44] # Amplituden Fkt PD M1599 @PM=3kHz und 1M Wiederstand, Array sollte len(Harmonics)-1, Faktor sollte relativ zu 1H data sein!
AF = [1.21,1.47,1.80] # Amplituden Fkt PD M1599 @PM=5kHz und 1M Wiederstand, Array sollte len(Harmonics)-1, Faktor sollte relativ zu 1H data sein!
noiseDef = 3    #noise level is defined as noiseFloor= avgNoise + noiseDef*stdNoise
zeroPaddingFactor = 1   # how many powers of 2 should data set be filled with zeros

# Phasenanpassung H1 wird auf optimales T-signal angepasst durch iterativen Vorgang mit folgenden Anfangsbedigungen
tCorr = 0.0  # correction of time axis in fs, if WPI scan started at positive delay then use tCorr > 0. E.g. AC yielded zero delay at -30fs
phaseStep= 0.01 # Phasenschritte bei interativer automatischer Phasenanpassung, in degree
acc = 10000.0 #Mass für Genauigkeit der Phasenanpassung ergibt sich durch Y-max/acc.
#PF =[5.9, 11, 15.7] # Phasenfkt PD M1599 @PM=1kHz, in degree, Array sollte len(Harmonics)-1, Phasenshift sollte relativ zu 1H data sein!
#PF =[12.8,22.6,35.8] # Phasenfkt PD M1599 @PM=3kHz, in degree, Array sollte len(Harmonics)-1, Phasenshift sollte relativ zu 1H data sein!
#PF =[14.9,26,34.2] # Phasenfkt PD M1599 @PM=4kHz, in degree, Array sollte len(Harmonics)-1, Phasenshift sollte relativ zu 1H data sein!
PF =[16.6, 28, 36.3] # Phasenfkt PD M1599 @PM=5kHz, in degree, Array sollte len(Harmonics)-1, Phasenshift sollte relativ zu 1H data sein!
#PF = [2.5, 4.0, 5.0] # Phasenfkt PMT @ PM=1kHz, in degree
#PF = [3.2, 5.5, 7.3] # Phasenfkt PMT @ PM=5kHz, in degree
phiCorr = [-90*n for n in [0,1,2,3]]     #manuelle Nachkorrektur der Phase (z.B. Offset Phase des Ref Signals), in degree, Phasenshift sollte relativ zu 1H data sein! 
#phiCorr = [90,0,-90,-180]
#phiCorr = [p-15 for p in phiCorr]

# File IO =====================================================================
report = False  # should be written to rpt file?
initDir = '\\\\10.5.71.1\\user\Projekte\\2D-Spectro\\Dat' # directory for "open file"-dialog
loadSpectrum = False
boxcarWindow = 'flat'   #window type for boxcar averaging. Supported types: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman' 
boxcarLength = 5    # if < 3, no boxcar averaging is performed


# Plot parameters =============================================================
plotCal = False  # plot mono calibration
plotLaserFit = False # plot laser spectrum fit
plotTD = True
plotFD = True
saveFig = False
plotColor = ['r', 'b', 'g', 'm']    # colors for plots, should have same length as Harmonics
annotateColor = 'k' # color for annotating labels
leakColor = '#999999' # color for annotating leak
monoName = 'M'  # necessary to define globally
legend = ['1H', '2H', '3H', '4H']   # legend text, no legend displayed if empty list
legendLoc = 'upper right'   # location of legend, auto location if empty string
plainLegend = True  # box around legend?
majorLocXTD=None #[4,4,4,4]   # specifies for each harmonic  separation between xticks, if =None ticks are optimized by python
majorLocYTD=None#[5e-8,5e-10,5e-11,0.4e-10]   # specifies for each harmonic  separation between yticks, if =None ticks are optimized by python
majorLocXFD=None #[4,4,4,4]   # specifies for each harmonic  separation between xticks, if =None ticks are optimized by python
majorLocYFD=None#[5e-8,5e-10,5e-11,0.4e-10]   # specifies for each harmonic separation between yticks, if =None ticks are optimized by python
lw = 2
#xlim = [12980,13050, 2*12980, 2*13050, 38930, 39150, 51900, 52200]  # K D1&D2
#xlim = [12800,12900, 25600, 25800, 38400, 38700, 51200, 51600] # Rb D2
#xlim = [12550,12620, 25140, 25450, 38400, 38700, 51250, 51600]  # Rb D1
#xlim = [12530,12950, 25600, 26600, 39200, 39800, 52400, 52900] # Rb 7s+5D
#xlim = [12100,12600, 24400, 25000, 36300, 37800, 48600, 50200] # K SHG
#xlim = [11500,12500, 23500, 24200, 35000, 36500, 47100, 48100] # Rb D2 SHG
#xlim = [11650, 12150, 23400, 24100, 35000, 36500, 47300, 47700] # Rb 6P SHG
#xlim = [12775,13085, 25600, 26010, 38000, 39500, 50500, 52000] # Cs D3 SHG
#xlim = [12750,13100, 25600, 25900, 38500, 38750, 51900, 52200]  # K,Rb,Cs SHG+Fund
#xlim = [12750,13100, 25600, 26050, 38930, 39150, 51900, 52200]  # K+Rb mixedSpecies
#xlim = [12985.19-25,12985.19+25, 2*12980, 2*13050, 38930, 39150, 51900, 52200]  # K D1&D2
#xlim = [12350,13150, 25000, 25750, 12350,13150, 51200, 51600] # Rb NOPA
xlim=[12200,13300,24900, 26000, 37600,38600,0,50000]
#xlim= [10000,15000, 22500, 27500, 10000,15000, 22500, 27500]
center=38360
width=1000
xlim=[center-width/2.0,center+width/2.0 , 2*center-width/2.0,2*center+width/2.0 , 3*center-width/2.0,3*center+width/2.0]




# define spectral lines =======================================================
# potassium:
K1HLines = { \
    'K_D1': 12985.19, \
    'K_D2': 13042.90, \
#    '1H_K_5P_1/2': 0.5*24701.382, \
#    '1H_K_5P_3/2': 0.5*24720.139, \
    }
K2HLines = { \
    'K_2D1': 2*12985.19, \
    'K_2D2': 2*13042.90, \
    'K_D1D2': 12985.19 + 13042.90, \
#    'K_5P_1/2': 24701.382, \
#    'K_5P_3/2': 24720.139, \
    }
K3HLines = { \
    'K_3D1': 3*12985.19, \
    'K_3D2': 3*13042.90, \
    'K_2D1D2': 2*12985.19 + 13042.90, \
    'K_D12D2': 12985.19 + 2*13042.90, \
#    '3H_K_5P_1/2': 1.5*24701.382, \
#   '3H_K_5P_3/2': 1.5*24720.139, \
    }
K4HLines = { \
    'K_4D1': 4*12985.19, \
    'K_4D2': 4*13042.90, \
    'K_3D1D2': 3*12985.19 + 13042.90, \
    'K_D13D2': 12985.19 + 3*13042.90, \
    'K_2D12D2': 2*12985.19 + 2*13042.90, \
#    '4H_K_5P_1/2': 2*24701.382, \
#    '4H_K_5P_3/2': 2*24720.139, \
    }


#rubidium
Rb1HLines = { \
    'Rb_D1': 12578.95, \
    'Rb_D2': 12816.55, \
    'Rb_p3d5' : 12886.95, \
    'Rb_p3d3' : 12883.99, \
#    'Rb/2': 0.5*23792.591, \
#    'Rb/2': 0.5*23715.081,\
    'Rb_p1d3' : 25700.54-12578.95, \
    'Rb_7s-pd' : 26311.437-(25700.54-12578.95) ,\
    }
Rb2HLines = { \
#    'Rb_6p32': 23792.591, \
#    'Rb_6p12': 23715.081,\
    'Rb_2D1': 2*12578.95, \
    'Rb_2D2' : 2*12816.55, \
#    'Rb_2p3d5' : 2*12886.95, \
#    'Rb_2p3d3' : 2*12883.99, \
##    'Rb_d3+d5' : 12886.95 + 12883.99 ,\
    'Rb_s1d3' : 25700.54 ,\
    'Rb_s1d5' : 25703.5 ,\
    'Rb_5s7s' : 26311.437 ,\
    }
Rb3HLines = { \
    'Rb_3D1': 3*12578.95, \
    'Rb_3D2' : 3*12816.55, \
#    'Rb_D2+2p3d3' : 12816.55+2*12883.99 ,\
#    'Rb_D2+2p3d5' : 12816.55+2*12886.95 ,\
#    'Rb_D2+d3+d5' : 12816.55+12883.99+12886.95 ,\
#    'Rb_D2+s1d3' : 12816.55+25700.54 ,\
#    'Rb_D2+s1d5' : 12816.55+25703.5 ,\
    }
Rb4HLines = { \
#    'Rb_4D1': 4*12578.95, \
#    '2Rb_6P_1/2': 2*23715.081,\
#    '2Rb_6P_3/2': 2*23792.591,\
    'Rb_4D2' : 4*12816.55 ,\
#    'Rb_4p3d3' : 4*12883.99 ,\
#    'Rb_4p3d5' : 4*12886.95 ,\
#    'Rb_3D2+p3d3' : 3*12816.55+12883.99 ,\
#    'Rb_3D2+p3d5' : 3*12816.55+12886.95 ,\
#    'Rb_D2+3p3d3' : 12816.55+3*12883.99 ,\
#    'Rb_D2+3p3d5' : 12816.55+3*12886.95 ,\
###    'Rb_2D2+2p3d3' : 2*12816.55+2*12883.99 ,\  # same as 'Rb_2s1d3'
###    'Rb_2D2+2p3d5' : 2*12816.55+2*12886.95 ,\  # same as 'Rb_2s1d5'
###    'Rb_2D2+s1d3' : 2*12816.55+25700.54 ,\ # same as'Rb_3D2+p3d3'
###    'Rb_2D2+s1d5' : 2*12816.55+25703.5, \  # same as 'Rb_3D2+p3d5'
#    'Rb_2s1d3' : 2*25700.54 ,\
#    'Rb_2s1d5' : 2*25703.5, \
    '2Rb_5s7s' : 2*26311.437 ,\
    }
 
Cs1HLines = { \
    'Cs_D1': 11178.268, \
    }
   
Cs2HLines = { \
    'Cs_2D1': 2*11178.268, \
    'Cs_8p32': 25791.508, \
    'Cs_8p12': 25708.83548,\
    'Rb+Cs': 0.5*25791.508 + 12816.55,\
    }
Cs3HLines = { \
    'Rb+Cs': 25791.508 + 12816.55,\
    }
    
Na1HLines = { \
    'Na_3s6p' : 37297.61 \
    }
    
#rubidium
Rbn1HLines = { \
    'Rb_D1': 12578.95, \
    'Rb_D2': 12816.55, \
    'Rb_p3d5' : 12886.95, \
    'Rb_p1d3' : 25700.54-12578.95, \
    'Rb2_T' : 13500 ,\
    'Rb2_S' : 10800 ,\
    'Rb3_Q' : 11760 ,\
    }