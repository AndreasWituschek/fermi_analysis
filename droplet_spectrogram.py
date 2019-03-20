# -*- coding: utf-8 -*-
"""
Created on Thu Nov 07 10:47:33 2013

@author: Luki



"""

# to do:
# error handling für import_report()
# deutsche umlaute in report datei


import scipy as sp
import scipy.constants as const
from scipy import fftpack
import numpy as np
import matplotlib.pyplot as plt
import string
from math import log10, floor
from Tkinter import Tk
from tkFileDialog import askopenfilename
import os

#=========user inputs==========================================================
script_acronym = 'es'
init_dir = '\\\\10.5.71.1\\user\\Projekte\\2D-Spectro\\Data\\WPI-simple_PE_detector'
os.chdir('\\\\NANOSERVER\\Marcelb\\PhD-stuff\\DataAnalysisStuff\\Python\\Droplet\\RbHeN_dissociation_DPG18')
data_window = [-3,3]   # which part of the interferogram should be considered in ps
t_lim = [0,.6]
wn_lim = [12530,12950]
mono_calib = 13068.815  #calibrated mono value in cm^-1, if this is 0.0, sktript will use value from rpt-file + mono_corr
zero_delay = 0.0 # in ps
slide_step = 0.01 # in ps
FWHM = .150 # in ps
t_desorbed = 1.8
overlap = 4.0#0.8*4.65#4.65 # half width in ps
clim = [0,0]#[0.0, 0.6E-14]   # limit of color range
mono_corr = -0.19 # monochromator correction in nm: mono_actual = mono_scale + mono_corr
lw = 2 # linewidth in DFTs
fs = 14 # font size for annotations
Ham =1     # 1, 2, 3, 4te Harmonische, vorrausgesetzt LockIn 1 = Ham1,...

#init_dir = '\\\\10.5.71.1\\user\\Projekte\\2D-Spectro\\Data\\WPI-simple_PE_detector'
#os.chdir('\\\\10.5.71.1\\user\\Projekte\\2D-Spectro\\Data\\WPI-simple_PE_detector')
#file_path = '\\\\nanoserver\\user\\Projekte\\2D-Spectro\\Data\\QMS_2014\\2014-10-02\\2014-10-02_Rb87He_QMS_776nm_8215A_1+2H.dat'
file_path = '\\\\NANOSERVER\\Marcelb\\PhD-stuff\\DataAnalysisStuff\\Python\\Droplet\\RbHeN_dissociation_DPG18\\2017-08-07_1D_droplet_3.dat'
#file_name = '2014-10-02_Rb87He_QMS_776nm_8215A_1+2H'
file_name = '2017-08-07_1D_droplet_3'

num_sub_plots = 3
Rx = 0.8    # % der Gesamtgröße der fig
Ry = 0.8    # % der Gesamtgröße der fig
Rsubx = Rx  # breite des Subplots in %
Rsuby = Ry/num_sub_plots   # Höhe des Subplots in %
Offx = (1-Rx)*3/4.0 # Offset-Koordinate in %
Offy = (1-Ry)/2.0   # Offset-Koordinate in %



#Definig some functions========================================================

def import_report (keyword_substr, report_name):
    """
    Searches for the first occurence of keyword_substr in report file and
    returns value associated with the corresponding keyword. If keyword_substr
    occures several times in the report file, only the value of first occurence
    is returned, no error will be thrown.
    Function can't handle german vocals ä,ö,ü, in this case truncate
    keyword_substr just before any of these characters.
    """
    report = np.genfromtxt(report_name+'.rpt', delimiter = ':', skip_header = \
    12 , usecols = (0,1), dtype=str)    # wichtig hier den dtype auf str festzulegen sonst macht genfromtxt probleme mit colum 2
    keyword_list, values = report.T
    keyword_index = -1
    for i in np.arange(0, keyword_list.size):
        if string.find(keyword_list[i], keyword_substr) >= 0 :
            keyword_index = i
            break
    return float(values[keyword_index])

def cut_dataset (T, X, Y, window):
    """
    new method!!!, 02.07.14
    """
    # sort start,stop according to increasing/decreasing numbers in T
    if T[-1]>T[0]:
        start_user = min(window)
        stop_user = max(window)
    else :
        start_user = max(window)
        stop_user = min(window)
    # cut data
    lower_index = np.argmin(np.abs(T-start_user))
    upper_index = np.argmin(np.abs(T-stop_user))
    T = T[lower_index : upper_index]
    X = X[lower_index : upper_index]
    Y = Y[lower_index : upper_index]
    return (T,X,Y)

def round_sig(x, sig=2):
   return round(x, sig-int(floor(log10(x)))-1)

def get_file_name (file_path):
    index = string.rfind(file_path, '/')
    return file_path[index+1:-4]    #truncate '.dat'

def get_export_folder (file_path):
    index = string.rfind(file_path, '/')T_d
    return file_path[0:index]    #truncate '.dat'

def open_file (init_folder):
    # define Masterwindow/Canva and close hide it
    root = Tk()
    root.withdraw()
    # show an "Open" dialog box and return the path to the selected file
    file_path = askopenfilename(initialdir= init_dir, defaultextension='.dat') # returns complete file path
    #file_path = 'D://Dokumente//Diplomarbeit//2D-FT-Spectro//Data//WPI//07.03.13//WPI_07.03.13_6.dat'
    return file_path

def alkali_one_photon_wavenumbers():
    # Nist wavenumbers /cm^-1:
    wn_Rb_D1 = 12578.95
    wn_Rb_D2 =  12816.55
    wn_Rb_d3p = 12883.99    # 5pP_3/2 -> 5dD_3/2
    wn_Rb_d5p = 12886.95    # 5pP_3/2 -> 5dD_5/2
    wn_dict = {'Rb_D1': wn_Rb_D1, 'Rb_D2': wn_Rb_D2, \
    'Rb_d3p': wn_Rb_d3p, 'Rb_d5p': wn_Rb_d5p}
    return wn_dict

def mono_frequency(lambda_M):
    # diffraction index for lambda_M (formula Zula SabrinaL)
    lambda_M = lambda_M*10**-9  # conversion to meter
    n = 1 + \
    8.06051*10**(-5) + \
    1.87564*10**(-4)*lambda_M**2/(lambda_M**2-7.56006*10**(-15)) + \
    4.43831*10**(-6)*lambda_M**2/(lambda_M**2-2.54262*10**(-15))
    # calculate monochromator frequency in vacuum
    return const.c/lambda_M/n  # in Hz

def data_mirrow_extend(T,X,Y):
    T_neg = -T[::-1]
    T= -np.append(T_neg, T)
    X_neg = X[::-1]
    X= np.append(X_neg, X)
    Y_neg = -Y[::-1]
    Y= np.append(Y_neg, Y)
    Z = X +1j*Y
    return T, Z

def import_data2():    # import data of two lock-ins
    # ceate file strings
    #file_path = open_file('')
    file_name = get_file_name(file_path)
    title = script_acronym + '_' + file_name +'_spectrogram'
    # import data
    data = np.genfromtxt(file_path, delimiter = '\t', skip_header = 1 , usecols = (2,3,4,5,6))
    T,X1,Y1,X2,Y2 = data.T
    # import experiment parameters
    start = import_report('Startposition', file_path[:-4])   #in fs
    stop = import_report('Endposition', file_path[:-4])   #in fs
    step = import_report('Schrittweite', file_path[:-4])   #in fs
    amp_1 = import_report('Lock-In_1 Amplification', file_path[:-4])
    sensitivity_1 = import_report('Lock-In_1 T_dSensitivity', file_path[:-4])
    amp_2 = import_report('Lock-In_2 Amplification', file_path[:-4])
    sensitivity_2 = import_report('Lock-In_2 Sensitivity', file_path[:-4])
    lambda_L = import_report('Chameleon Spektrum', file_path[:-4])
    wn_L = 10**7/lambda_L   # wavenumbers in cm^-1
    lambda_M_scale = import_report('Monochromator Wavelength', file_path[:-4])
    #lambda_M_corr = import_report('Monochromator Correction', file_path[:-4])
    print 'laser=' + str(lambda_L) +' - ' + file_name + '\n'
    lambda_M_corr = mono_corr
    #Prepare Data==================================================================
    # converting wavelength to actual air wavelength
    lambda_M = lambda_M_scale*1.004-49.872 + lambda_M_corr  # Eichung vom 22.01.14
    freq_Mono = mono_frequency(lambda_M)    # in Hz
    wn_M = Ham*freq_Mono/100.0/const.c # wavenumbers in cm^-1
    # scaling signal, neglecting 2Ham signal
    X1 *= sensitivity_1/amp_1/10.0  # orignial signal magnitude in Volts
    Y1 *= sensitivity_1/amp_1/10.0  # orignial signal magnitude in Volts
    X2 *= sensitivity_2/amp_2/10.0  # orignial signal magnitude in Volts
    Y2 *= sensitivity_2/amp_2/10.0  # orignial signal magnitude in Volts
    # cut dataset
    T1, X1, Y1 = cut_dataset(T, X1, Y1, data_window)
    T2, X2, Y2 = cut_dataset(T, X2, Y2, data_window)
    #Z1 = X1 + 1j*Y1
    Z1 = X1
    Z2 = X2 + 1j*Y2
    para_dict = {'start' : start, 'stop' : stop, 'step' : step, 'lambda_L' : lambda_L, \
    'wn_L' : wn_L, 'lambda_M' : lambda_M, 'wn_M' : wn_M, 'title' : title}
    return T1, Z1, Z2, para_dict

def slide_window(T, Z, t_center, FWHM):
    # multiplies a gaussion onto the dataset and returns the data set of the whole range of T
    # T, t_center, FWHM in ps
    # Z may be complex
    return Z*sp.exp(-4*sp.log(2)*((T-t_center)/FWHM)**2)

def weighting_coeff(t, t_center, FWHM): # calculates amplitude of gaussian at time t
    return sp.exp(-4*sp.log(2)*((t-t_center)/FWHM)**2)

def zero_padding(T):
    n = int(round(np.log2(T.size)+0.5))  #potenz um datensatz auf nächste 2^N zu erweitern
    n += 1  # falls zero-padding erwünscht
    N = 2**n    # Anzahl Datenpunkte auf die zur Not mit Zero-Padding aufgefüllt wird
    return N

def redistribute(array, index):
    array_neg = array[index:]
    array_pos = array[0:index]
    return np.concatenate((array_neg,array_pos))

def plotbar():
    ##cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    #cax = fig.add_axes([0.85, 0.1, 0.2, 0.8])
    #cax.get_xaxis().set_visible(False)
    #cax.get_yaxis().set_visible(False)
    #cax.patch.set_alpha(0)
    #cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')

def annotate():
    ax = plt.gca()
    ax.axvline(para_dict['wn_M'] , linestyle='--', color = 'w', linewidth=lw)
    ax.axvline(para_dict['wn_L'] , linestyle='--', color = '#ff0000', linewidth=lw)
    ax.axhline(zero_delay , linestyle='--', color = 'w', linewidth=lw)
    ax.axvline(12883.99 , linestyle='--', color = '#888888', linewidth=lw)
    ax.axvline(12886.95 , linestyle='--', color = '#888888', linewidth=lw)

def plot_spectrogram(fig, wn_lim, t_lim):
    ax = fig.add_axes([Offx, Offy, Rsubx, 2*Rsuby])  # define current axes
    xmin_index = np.argmin(abs(wn-wn_lim[0]))
    xmax_index = np.argmin(abs(wn-wn_lim[1]))
    tmin_index = np.argmin(abs(T-t_lim[0]))
    tmax_index = np.argmin(abs(T-t_lim[1]))
    clim= [0.0, 0.9*S.max()] # 0.65
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
    ax.set_xlim(t_lim[0], t_lim[1])
#    ax.text(0.9, 0.9, '(b)', transform=ax.transAxes, ha='left', fontweight='bold',fontsize = 22, color='w')

    #annotate()
    #plotbar()

def find_index(array, value):
    return np.argmin(abs(array-value))

#Script Start==================================================================
## weight DFT with window size
## works only if start = 0ps!!!
#if slide_pos < HWFP:
#    # slide window is cut at start/end of data set
#    DFT_slide *= (slide_pos + HWFP)/(2*HWFP)

data_window[0] *= 1000.0   # conversion to fs
data_window[1] *= 1000.0   # conversion to fs
T, Z1, Z2, para_dict = import_data2() # returns T in fs
T /= 1000.0 # convert to ps
if Ham ==1:
    Z = Z1
else:
    Z = Z2
if mono_calib != 0.0:
    para_dict['wn_M'] = mono_calib# measured wavenumbers in cm^-1

# zero-padding and calulate frequency axis
N = zero_padding(T)
freq = fftpack.fftfreq(N, d = para_dict['step']*10**(-15))  # frequency axes in Hz
if para_dict['step']<0:
    freq = freq[::-1]   # reverse array
wn = freq/100.0/const.c + para_dict['wn_M'] # measured wavenumbers in cm^-1
cut_index = np.argmin(wn)
wn = redistribute(wn, cut_index)
# create spectrogram-matrix
S = np.zeros((np.size(T), np.size(wn)))
# calculate start position
data_window[0] /= 1000.0   # back to ps
data_window[1] /= 1000.0   # back to ps
HWFP = sp.sqrt(3.0/sp.log(2))*FWHM/2.0 # half width where gaussian has dropped to 5% i.e. 1/e^3
slide_pos = data_window[0] - HWFP
n=0

while slide_pos <= data_window[1] + HWFP :
    Z_slide = slide_window(T, Z, slide_pos, FWHM)    # multiply gaussian onto data set which is centered at current sliding position
    DFT_slide = redistribute(abs(fftpack.fft(Z_slide, n=N)),cut_index)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
    n+=1
    if para_dict['step']<0:
        DFT_slide = DFT_slide[::-1] # reverse array
    # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
    for i in range(np.size(T)):
        S[i,:] += DFT_slide*weighting_coeff(T[i], slide_pos, FWHM)
    slide_pos += slide_step  # move to next slide position
# normalize to get proper average
S[:,:] /=n


# plot spectrogram
#plt.close('all')
fig = plt.figure('spectrogram_paper', figsize=(7,7))
plot_spectrogram(fig, wn_lim, t_lim)
Z *= 10**9
ax2 = fig.add_axes([Offx, Offy+2*Rsuby, Rsubx, Rsuby])  # define current axes
ax2.plot(T,abs(Z), '-', color = 'k', linewidth = 0.7)
ax2.set_xlim(t_lim[0],t_lim[1])
ax2.set_ylim(-0.1,abs(Z).max()*1.1)
desorb_index = find_index(T, t_desorbed)
max_amp = np.max(abs(Z[desorb_index:]))
ax2.axhline(max_amp, linestyle='--', color='#666666', linewidth=lw)
ax2.set_xticklabels([])
ax2.get_yaxis().set_tick_params(which='both', length=5, width=1.5)
ax2.get_xaxis().set_tick_params(which='both',  length=5, width=1.5)
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(1.5)
ax2.set_ylabel('ion yield (arb.u.)', fontsize=fs)
ax2.tick_params(labelsize=fs)
yticks = ax2.yaxis.get_major_ticks()
ax2.axvspan(-0.4,0.4, color='grey', alpha=0.5)
ax2.text(0.9, 0.84, '(a)', transform=ax2.transAxes, ha='left', fontweight='bold',fontsize = 22, color='k')

#AC_FWHM = T(find_index(abs(Z), 1.0)[1])-T(find_index(abs(Z), 1.0)[0])
#print AC_FWHM

ax3 = fig.add_axes([Offx+Rsubx/3.5, Offy+2*Rsuby+Rsuby-Rsuby/1.8, Rsubx/2.5, Rsuby/1.8])
ax3.plot(T,abs(Z), '-', color = 'k', linewidth = 1.2)
ax3.axhline(max_amp, linestyle='--', color='#666666', linewidth=lw)
ax3.set_xlim(15,20)
ax3.set_ylim(-0,0.3)
ax3.set_yticklabels([])
ax3.tick_params(labelsize=fs)
ax3.get_yaxis().set_tick_params(which='both', length=5, width=1.5)
ax3.get_xaxis().set_tick_params(which='both', length=5, width=1.5)
ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
ax3.set_yticks([])
ticks = ax3.xaxis.get_major_ticks()
#xticks[0].label1.set_visible(False)
#xticks[-1].label1.set_visible(False)
#for axis in ['bottom','left', 'right']:
    #ax3.spines[axis].set_linewidth(1.5)
ax2.annotate("", xy=(17.5, 0.3), xycoords='data',xytext=(25.7, 0.95), textcoords='data', \
arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.3", linewidth=lw))

plt.show()