# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:46:56 2019

@author: FemtoMeasure
"""

import numpy as np
import fermi_analysis.functions as fk 
import h5py
from scipy import signal
import matplotlib.pyplot as plt
import scipy.constants as spc
from matplotlib.gridspec import GridSpec

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

#plt.close('all')
""" Data parameters """
run = 31 #first run of delay scan
demod = ['0']
device = ['dev3265']
#demod = ['0','1', '2']
#device = ['dev3265', 'dev3269']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
root_file_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
interactive_plots = 1
save_fig = 0
plotTheo = 0
spectrogram = 1
plot_eV = 1

"""Experimental Parameters"""
# parameters for theoretical curve
E_r = 23.74207019 # helium 1s-4p transition in eV 
l_trans =  spc.h*spc.c/(E_r*spc.e)*1e9  # helium 1s-4p transition in nm 
#l_trans = 260.92/6.0 # fanoresonance in argon
#l_ref = 266.0234244226323  #265.98  # reference laser wavelenght in nm
fwhm_FEL = 0.075 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 5.  # harmonic of FEL
delay_zero_pos = 11025.66
title = 'Scan_{0:03}/'.format(int(run))
color = 'r'
data_window = [150,700] #[-200,370]
gauss = 1 # gauss window on TD (true) or not (False)
modfreq = 3.96


'''Phase correction'''
phas_corr = harmonic*8. #phase of reference transmitter [degree] (gamma in Lukas phasing routine)
phas_corr -= 0.01*(harmonic*modfreq)*360 #phase of boxcar integrator [degree] (beta in Lukas phasing routine)
print('A phase correction of {} degrees due to electronics was applied'.format(phas_corr))
phas_res = 30. #more or less arbitrary phasing factor that makes real part absorptive and imaginary part dispersive
print('A phase correction of {} degrees was applied to yield a perfect absorptive line shape'.format(phas_res))
phas_corr += phas_res
print('Total phase correction: {} degree'.format(phas_corr))
phas_corr -= 180. #to take into account that we measure negative signals by ion detection
phas_corr += (180.*harmonic) % 360. #accounts for phase shift due to the fact that ref signal is other output of interferometer.

""" Parameters for theoretical curve """
phi = 0. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
i0_correction = False
i0_5H_correction = False
unwrap_phase = True
slide_step = 32
FWHM_slide = slide_step / np.sqrt(48)
zeroPaddingFactor = 1
suscept = False
peakMargin = 0.05
#wn_lim = [217100,233800]
wn_lim = np.asarray([185000.,196000.])
#wn_lim = [189800,193100]

# Load preanalysed data
data = fk.ImportPreanalysedData(root_file_path,run)

l_fel = data['LDM']['l_seed'][0]
ldm_l_ref = data['LDM']['l_ref']
l_ref = 266.003
l_fel = 261.0726087 #261.7 #261.0726087
T = data['LDM']['delay']-5.38


E_undersamp = harmonic*spc.h*spc.c/(l_ref*spc.e)*1e9 # undersampling energy
E_r = E_r - E_undersamp  #eV. undersampled resonance energy



''' theoretically calculated laser spectrum '''
points = np.linspace(wn_lim[0],wn_lim[1],1000)
fwhm_FEL_wn = 1E7*(fwhm_FEL/(l_fel/harmonic)**2) # wavenumber FWHM of FEL
l_FEL = 1E7/l_fel # wavenumber CWL of FEL


''' Laser spectrum from Run30 '''
#padres_roi_5H, padres_spec_5H =  fk.ImportPadresRun30()

''' Fano resonance simulation '''
Ttheo = np.linspace(data_window[0],data_window[1],1000)
AC_FWHM = 100.0  # in fs
Fano_FWHM = 140.0  # in fs
AC = np.exp(-4.0*np.log(2)*((Ttheo)/AC_FWHM)**2)
Fano = np.exp(-4.0*np.log(2)*((Ttheo)/Fano_FWHM)**2)
Theo = 3*(3*Fano-7*AC)

''' Seed laser AC simulation '''
sl_fwhm = 99. # seed intensity fwhm in fs
t_theo_seed = np.linspace(-250,250,5000)
sl_AC = np.exp(-4.0*np.log(2)*(t_theo_seed/(sl_fwhm*np.sqrt(2)))**2) # seed laser AC
sl_AC /= np.max(sl_AC)


''' Resolution from Data window '''
print('The resolution from the data window length is {} meV = bin size of the DFT'.format(1000.*spc.h/(abs(data_window[-1]-data_window[0])*1E-15*spc.e)))


''' i0 correction for different harmonic contents of i0 '''
i0_fft = data['LDM']['i0_fft']
i0_0H = np.sum(i0_fft[::,395:405],axis=1)
i0_1H = np.sum(i0_fft[::,458:468],axis=1)
i0_2H = np.sum(i0_fft[::,522:532],axis=1)
i0_3H = np.sum(i0_fft[::,585:595],axis=1)
i0_4H = np.sum(i0_fft[::,648:658],axis=1)
i0_5H = np.sum(i0_fft[::,711:721],axis=1)

#i0_5H = i0_5H-np.mean(i0_5H[200:])
i0_5H /= np.max(i0_5H)
#i0_5H = abs(i0_5H)


# loop over the mfli devices
for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = (data[dev]['x' + d] + 1j * data[dev]['y' + d])# / data['LDM']['I0']
        if i0_correction:
            Z = Z*(np.nanmean(data['LDM']['I0'][228:303])/data['LDM']['I0'])
        if i0_5H_correction:
            Z = Z*i0_5H

#        Z *= np.exp(1j*(5.*8./360.*2.*np.pi)) #phasing of reference offset
        Z *= np.exp(1j*np.pi/180.*phas_corr)
        R = np.sqrt(Z.real**2 + Z.imag**2)
        Phi = np.angle(Z, deg=False)


# Create theoretical curve
#        if draw_theory:
        #Ttheo = np.linspace(T[0],T[-1], 1000)
        Xtd,Ytd,Xt,Yt = fk.Curve(l_trans, l_ref, harmonic, phi, a*max(abs(Z)), offset, T[0], T[-1], 1000)

        delay = T

        T_d,X_d,Y_d = fk.CutDataSet(T, Z.real, Z.imag, data_window)
        Z_cut = X_d + 1j*Y_d
        Z_cut = fk.Phasing_TD(Z_cut,T_d,E_r)
        scaling = 1E-4 #max(Z.real)/max(Ttheo)

        # Fourier traffo
        Td = np.mean(np.unique(np.diff(T)))
        if gauss:
            Zg = fk.GaussWindow(T_d, Z_cut, suscept)
            wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        else:
            wn, dft = fk.DFT(T_d, Z_cut, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        FWHM , peakFullWidth = fk.PeakWidth(T_d, peakMargin, suscept, zeroPaddingFactor)



        # theoretical spectrum
        center = 1E7*(1/l_fel- 1/l_ref)*mfli_harmonic + harmonic*1E7/l_ref
        spectrum = np.exp(-4.0*0.693147*(points-center)**2/fwhm_FEL_wn**2) # spectrum in wavenumber space
        spectrum*=(max(abs(dft))-min(abs(dft)))+min(abs(dft))
        
        ''' Plot time domain '''
        figTD = plt.figure('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i), figsize=(14, 10))

        ax = figTD.add_subplot(321)
        ax.set_title('scan_{0:03}_{1}_d{2}_TD_{3}'.format(run, dev, d, i))
        #ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
        ax.plot(delay, Z.real, color=color)
        ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.real))*0.8,color='grey')
        ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z.real))*0.8,color='grey')
        if plotTheo: ax.plot(np.linspace(T[0]-10,T[-1]+10,len(Xtd)), Xtd, 'b', alpha=0.3)
        ax.set_ylabel('X [a.u.]')
        ax.set_xlabel('delay [fs]')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(322)
        #ax.errorbar(delay, Y, yerr=Y_s, color=color, linestyle='')
        #ax.plot(Ttheo, scaling*Theo, 'k-')
        ax.plot(delay, Z.imag, '-', color=color)
        ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.imag)),color='b',alpha=0.3)
        if plotTheo: ax.plot(np.linspace(T[0]-10,T[-1]+10,len(Ytd)), Ytd, 'r', alpha=0.3)
        ax.set_ylabel('Y in a.u.')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(323)
        #ax.errorbar(delay, R, yerr=R_s, color=color, linestyle='')
        ax.plot(delay, R, '-', color=color)
        ax.fill_between(t_theo_seed,sl_AC*np.max(R),color='black',alpha=0.3)
        ax.set_ylabel('R')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        ax = figTD.add_subplot(324)
        #ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        if unwrap_phase:
            ax.plot(delay, np.unwrap(Phi), '-', color=color)
        else: 
            ax.plot(delay, Phi, '-', color=color)
        ax.set_ylabel('Phi (rad)')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        
        ##Plot I0
        ax = figTD.add_subplot(325)
#        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        if i0_5H_correction:
            ax.plot(delay, i0_5H, '-', color= color)
            ax.fill_between(t_theo_seed,sl_AC,color='black',alpha=0.3)
            ax.set_ylabel('I0 5H [a.u.]')
        else:
            ax.plot(delay, data['LDM']['I0'], '-', color= color)
            ax.fill_between(t_theo_seed,sl_AC*np.max(data['LDM']['I0']),color='black',alpha=0.3)
            ax.set_ylabel('I0 in $\mu$J')        
        
        
        ax.set_xlabel('Delay in fs')
#        ax.set_ylim(0,3)
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])
        
        ###Plot Delay
#        ax = figTD.add_subplot(515)
##        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
#        ax.plot(delay, '-', color=color)
#        ax.set_ylabel('Delay in fs')
#        ax.set_xlabel('Run')
##        ax.set_ylim(0,3)
#        ax.grid()
##        ax.set_xlim(data_window[0],data_window[-1])
#
        ##Plot Reflaser wavelenght
        ax = figTD.add_subplot(326)
#        ax.errorbar(delay, Phi, yerr=R_s, color='k', linestyle='')
        ax.plot(delay, ldm_l_ref, '-', color = 'g')
        ax.set_ylabel('wavelenght in nm')
        ax.set_xlabel('delay')
        ax.axhline(l_ref, color='k')
#        ax.set_ylim(0,3)
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])

        plt.tight_layout()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_TD_{3}.png'.format(run, dev, d, i), dpi=400)


        '''plot frequency domain'''
        if plot_eV:
            factor = spc.h*spc.c/spc.e*100.
            wn *= factor 
            wn_lim *= factor
            l_ref /= factor
            l_trans /= factor
            FWHM *= factor
            print('The frequency resoulution is {} meV FWHM'.format(FWHM*1000))
#            print('The FWHM of the peak is {} meV'.format(1000*(wn[fk.findPeakIndexes(abs(dft)/max(abs(dft)),0.5)[2]] - wn[fk.findPeakIndexes(abs(dft)/max(abs(dft)),0.5)[1]])))
                
        figFT = plt.figure('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
        axs = figFT.add_subplot(111)
        axs.plot(wn, abs(dft)/max(abs(dft)), '-', color='black', linewidth = 2)
        axs.plot(wn, dft.real/max(abs(dft)), '-', color=color, linewidth = 2)
#        axs.plot(wn, dft.imag/max(abs(dft)), '-', color='r', linewidth = 2)
#        axs.plot(points,spectrum,'-',color='grey',linewidth = 2)
        if plot_eV:
            axs.set_xlabel(r'energy [eV]', fontsize=14)
        else:
            axs.set_xlabel(r'wavenumber [cm$^{-1}$]', fontsize=14)
        axs.set_ylabel('spectral amp. [arb. u.]', fontsize=14)
        axs.set_title('scan_{0:03}_{1}_d{2}_FD_{3}'.format(run, dev, d, i))
#        axs.set_xlim(180000.,195000.)
        axs.set_xlim(wn_lim[0],wn_lim[1])
        axs.set_ylim(-1.1,1.2)
        axs.axvline(23.74207019, color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(axs.get_xlim())[0]
        plot_y_range = np.diff(axs.get_ylim())[0]
        axs.legend(['abs','real'])
#        axs.axhline(y= plot_y_range * 0.9 + axs.get_ylim()[0], xmin=0.1, xmax=0.1 + FWHM / plot_x_range, linewidth=3)
#        axs.text(0.1, plot_y_range * 0.9, str(FWHM))
        axs.grid()
        if save_fig:
            plt.savefig(file_path + 'scan_{0:03}_{1}_d{2}_FD_{3}.png'.format(run, dev, d, i), dpi=400)

#        if spectrogram:
#            # plot spec
#            figFT = plt.figure('scan_{0:03}_{1}_d{2}_SG_{3}'.format(run, dev, d, i))
#            print(T_d[0], T_d[-1])
#            slide_step = Td  #in fs
#            FWHM_slide = 20.0  #in fs
#            t_lim = data_window
#            slide_positions = np.arange(T_d[0],T_d[-1], slide_step)
#            S = np.zeros((np.size(T_d), np.size(wn)))
#            i = 0
#            for slide_pos in slide_positions:
#                Z_slide = fk.slide_window(T_d, Z, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
#                wn,  DFT_slide = fk.DFT(T_d, Z_slide, Td, l_ref , harmonic, zeroPaddingFactor = 2)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
#                # add DFT to spectrum-matrix while weighting with temporal gaussian envelope
##                for i in range(np.size(T_d)):
##                    S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))*fk.weighting_coeff(T_d[i], slide_pos, FWHM_slide)
#                S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
#                i += 1
#            S[:,:] /= np.size(slide_positions)
#            fk.plot_spectrogram(figFT, wn, T_d, S, l_fel, wn_lim, t_lim)

        if spectrogram:
            wn_lim_s = np.asarray([188000.,194500.])
            l_ref *= factor
            l_trans*=factor
            figFT = plt.figure('scan_{0:03}_{1}_d{2}_SG_{3}'.format(run, 2, 1, i),figsize=(7, 4))
            slide_step = Td  #in fs
            FWHM_slide = 20.0  #in fs
            slide_positions = np.arange(delay[0],delay[-1], slide_step)
            S = np.zeros((np.size(delay), np.size(wn)))
            i = 0
            for slide_pos in slide_positions:
                Z_slide = fk.slide_window(delay, Z, slide_pos, FWHM_slide)    # multiply gaussian onto data set which is centered at current sliding position
                wn,  DFT_slide = fk.DFT(T_d, Z_slide, Td, l_ref, harmonic, zeroPaddingFactor = zeroPaddingFactor)   # FFT of interferogram which has been truncated by mulitplying wiht gaussian
                S[i,:] += abs(DFT_slide)/max(abs(DFT_slide))
                i += 1
            S[:,:] /= np.size(slide_positions)
            fk.plot_spectrogram(figFT, wn, delay, S, l_fel/5., l_trans, wn_lim_s, data_window)
            plt.tight_layout()


if not interactive_plots:
    plt.close('all')
else:
    plt.show()

#i0_list = [i0_0H,i0_1H,i0_2H,i0_3H,i0_4H,i0_5H,i0_0H+i0_1H+i0_2H+i0_3H+i0_4H+i0_5H] 
#for ii in i0_list:
#    plt.plot(delay,ii)
#    plt.xlabel('delay [fs]')
#    plt.ylabel('i0 [a.u.]')
#plt.legend(['0H','1H','2H','3H','4H','5H','sum'])

#    
''' Plots for PM_XUV paper '''

ticksize= 2.
ticklength = 5.
fontsize=16.
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.major.width'] = ticksize
plt.rcParams['xtick.major.size'] = ticklength
plt.rcParams['ytick.major.width'] = ticksize
plt.rcParams['ytick.major.size'] = ticklength
plt.rcParams['axes.linewidth'] = ticksize
plt.rcParams['lines.linewidth'] = 2.


wn_lim = [harmonic*1E7/l_ref,24.2]

## pld version
fig_paper = plt.figure(figsize=(7,8.3))
gs = GridSpec(12, 10)
ax = fig_paper.add_subplot(gs[0:5,:])

ax.plot(delay, Z.real, 'b-',alpha=0.7)
ax.plot(5,5,'grey')
ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.real))*0.8,color='grey')
ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z.real))*0.8,color='grey')
ax.legend([r'Re(S)','seed AC'],fontsize=12)
ax.set_ylabel(r'intensity [a.u.]')
ax.set_xlabel(r'$\tau$ [fs]')
#ax.grid()
ax.set_xlim(0,data_window[-1])
ax.set_ylim(-0.35,0.35)
#plt.tight_layout()

ax1 = fig_paper.add_subplot(gs[7:,2:])
ax1.plot(wn, abs(dft)/max(abs(dft)), '-', color='b', alpha=0.7)
#ax1.plot(padres_roi_5H, padres_spec_5H,'black', alpha=0.38)
ax1.set_xlabel(r'energy [eV]')
ax1.set_ylabel(r'intensity [a.u.]')
ax1.legend([r'$\vert$FT(S)$\vert$'],fontsize=12)
#ax1.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
ax1.axvline(E_r + E_undersamp, color = 'k',linestyle='--', linewidth= 2)
ax1.set_xlim(wn_lim[0],wn_lim[1])
ax1.set_xlim(23.31,24.2)
ax1.set_ylim(-0.05,1.15)
ax1.set_xticks([23.4,23.6,23.8,24.0,24.2])
ax1.set_xticklabels([23.4,23.6,23.8,24.0,24.2])
#add second x axis to plot:
ax2 = ax1.twiny()
ax2.set_xlabel('rotating frame energy [eV]')         
ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)


#new version
plt.figure(figsize=(7, 3.45))
plt.plot(delay, Z.real, 'b-',alpha=0.7)
plt.plot(5,5,'grey')
plt.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.real))*0.8,color='grey')
plt.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z.real))*0.8,color='grey')
plt.legend([r'Re(S)','seed AC'],fontsize=12)
plt.ylabel(r'intensity [a.u.]')
plt.xlabel(r'$\tau$ [fs]')
#ax.grid()
plt.xlim(0,data_window[-1])
plt.ylim(-0.35,0.35)
plt.tight_layout()


fig_paper1 = plt.figure(figsize=(5,3.9))
ax1 = fig_paper1.add_subplot(111)
ax1.plot(wn, abs(dft)/max(abs(dft)), '-', color='b', alpha=0.7)
#ax1.plot(padres_roi_5H, padres_spec_5H,'black', alpha=0.38)
ax1.set_xlabel(r'energy [eV]')
ax1.set_ylabel(r'intensity [a.u.]')
ax1.legend([r'$\vert$FT(S)$\vert$'],fontsize=12)
#ax1.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
ax1.axvline(E_r + E_undersamp, color = 'k',linestyle='--', linewidth= 2)
ax1.set_xlim(wn_lim[0],wn_lim[1])
ax1.set_xlim(23.31,24.2)
ax1.set_ylim(-0.05,1.19)
ax1.set_xticks([23.4,23.6,23.8,24.0,24.2])
ax1.set_xticklabels([23.4,23.6,23.8,24.0,24.2])
ax1.set_yticks([0,0.5,1])
ax1.set_yticklabels([0.0,0.5,1.0])
#add second x axis to plot:
ax2 = ax1.twiny()
ax2.set_xlabel('rotating frame energy [eV]')         
ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)
plt.tight_layout()


#==============================================================================
# Signal to Noise:
#==============================================================================
#mean = np.mean([np.mean(abs(dft[np.argwhere(wn<23.7)])/max(abs(dft))),np.mean(abs(dft[np.argwhere(wn>23.8)])/max(abs(dft)))])
#std = np.mean([np.std(abs(dft[np.argwhere(wn<23.7)])/max(abs(dft))),np.std(abs(dft[np.argwhere(wn>23.8)])/max(abs(dft)))])
#noise = mean+3*std
#print 1/noise

#==============================================================================
# Short part of Scan
#==============================================================================
#Xtd,Ytd,Xt,Yt = fk.Curve(52.2186, 266.003, 5, 0, max(abs(Z)), offset, delay[0], delay[-1], 1000)
#
#figTD = plt.figure(figsize=(4.5, 3.3))
#
#ax = figTD.add_subplot(111)
#ax.plot(delay, 1.4*Z.real/max(Z.real), 'o', color='b',alpha=0.7)
##ax.plot(Ttheo, scaling*Theo, 'k-')
#ax.plot(np.linspace(delay[0]-10,delay[-1]+10,30000), 1*Xtd/max(Xtd), '-', color='black', alpha=0.7)
#ax.set_xlabel(r'$\tau$ [fs]')
#ax.set_ylabel(r'$Re(S)$ [a.u.]')
#ax.set_xlim(220,280)
#ax.set_ylim(-1.2,2)
##ax.set_yticks([-0.2,-0.1,0.0,0.1,0.2])
#ax.legend(['data','theory'],ncol=2)
#plt.tight_layout()

''' Plots for Fermi_pulse_overlap paper '''

#fig_paper = plt.figure(figsize=(7,6))
#
#ax = fig_paper.add_subplot(211)
#ax.plot(delay, Z.real, 'b-',alpha=0.7)
#ax.plot(5,5,'silver')
#ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.real))*0.8,color='silver')
#ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(Z.real))*0.8,color='silver')
#ax.legend([r'Re(S)','seed AC'],'upper left')
#ax.set_ylabel(r'intensity [a.u.]')
#ax.set_xlabel(r'$\tau$ [fs]')
#ax.set_xlim(-310,310)
#ax.set_ylim(-0.3,0.4)
#
#ax = fig_paper.add_subplot(212)
#ax.plot(delay, data['LDM']['I0'], 'b-',alpha=0.7)
#ax.set_ylabel(r'FEL pulse energy [$\mu$J]')
#ax.set_xlabel(r'$\tau$ [fs]')
#ax.set_xlim(-310,310)
#ax.set_ylim(20,80)
#plt.tight_layout()