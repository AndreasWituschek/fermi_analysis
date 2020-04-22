# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:46:56 2019

@author: FemtoMeasure
"""

import numpy as np
import fermi_analysis.functions as fk 
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as spc


#plt.close('all')
""" Data parameters """
run = 4651 #first run of delay scan
demod = ['1']
device = ['dev3265']
#demod = ['0','1', '2']
#device = ['dev3265', 'dev3269']
#root_file_path = '//10.5.71.28/FermiServer/Beamtime2/combined/'
root_file_path = '//10.5.71.1/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/'
interactive_plots = 1
save_fig = 0
plotTheo = False
spectrogram = 0
plot_eV = 1

"""Experimental Parameters"""
# parameters for theoretical curve
#l_trans = 52.2186  # helium 1s-4p transition in nm 
l_trans =  28.510223  # [eV] fanoresonance in argon
#l_ref = 266.0234244226323  #265.98  # reference laser wavelenght in nm
fwhm_FEL = 0.075 # FEL fwhm, used for calculation of FEL spectrum
harmonic = 6.  # harmonic of FEL
title = 'Scan_{0:03}/'.format(int(run))
color = 'b'
data_window = [152,600.] #values have to be multiples of delay step size (2fs in this case)
gauss = 0 # gauss window on TD (true) or not (False)

""" Parameters for theoretical curve """
phi = 0. #-190.0 #-100 #-20.0  # phase offset between reference and signal in degrees
A = 1.  # amplitude (=R from MFLI)
a = 1.   # amplitude scaling factor
offset = 0.  # offset

""" analysis parameters """
i0_correction = False
i0_6H_correction = False
unwrap_phase = True
slide_step = 32
FWHM_slide = slide_step / np.sqrt(48)
zeroPaddingFactor = 2
suscept = True
wn_lim = np.asarray([228900.,231080.])

# Load preanalysed data
file_path = root_file_path + 'scan_{0:03}/'.format(int(run))
h5f = h5py.File(file_path + 'scan_{0:03}.h5'.format(int(run)), 'r')

# sort data by delay in ascending order
sort_inds = np.array(np.negative(h5f.get('LDM/delay'))).argsort()


data = {
    'run_numbers': np.array(h5f.get('run_numbers'))[sort_inds],
    'LDM': {
        'I0': np.array(h5f.get('LDM/I0'))[sort_inds],
        's_I0': np.array(h5f.get('LDM/s_I0'))[sort_inds],
        'delay': np.negative(np.array(h5f.get('LDM/delay'))[sort_inds]),
        's_delay': np.array(h5f.get('LDM/s_delay'))[sort_inds],
        'l_seed': np.array(h5f.get('LDM/l_seed'))[sort_inds],
        'l_ref': np.array(h5f.get('LDM/l_ref'))[sort_inds],
        'i0_fft': np.array(h5f.get('LDM/i0_fft'))[sort_inds],
            },
    'dev3265': {
        'harmonic': np.array(h5f.get('dev3265/harmonic')),
        'x0': np.array(h5f.get('dev3265/x0'))[sort_inds],
        'y0': np.array(h5f.get('dev3265/y0'))[sort_inds],
        'x1': np.array(h5f.get('dev3265/x1'))[sort_inds],
        'y1': np.array(h5f.get('dev3265/y1'))[sort_inds],
        'x2': np.array(h5f.get('dev3265/x2'))[sort_inds],
        'y2': np.array(h5f.get('dev3265/y2'))[sort_inds],
        's_x0': np.array(h5f.get('dev3265/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('dev3265/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('dev3265/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('dev3265/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('dev3265/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('dev3265/s_y2'))[sort_inds], },
    'dev3269': {
        'harmonic': np.array(h5f.get('dev3269/harmonic')),
        'x0': np.array(h5f.get('dev3269/x0'))[sort_inds],
        'y0': np.array(h5f.get('dev3269/y0'))[sort_inds],
        'x1': np.array(h5f.get('dev3269/x1'))[sort_inds],
        'y1': np.array(h5f.get('dev3269/y1'))[sort_inds],
        'x2': np.array(h5f.get('dev3269/x2'))[sort_inds],
        'y2': np.array(h5f.get('dev3269/y2'))[sort_inds],
        's_x0': np.array(h5f.get('dev3269/s_x0'))[sort_inds],
        's_y0': np.array(h5f.get('dev3269/s_y0'))[sort_inds],
        's_x1': np.array(h5f.get('dev3269/s_x1'))[sort_inds],
        's_y1': np.array(h5f.get('dev3269/s_y1'))[sort_inds],
        's_x2': np.array(h5f.get('dev3269/s_x2'))[sort_inds],
        's_y2': np.array(h5f.get('dev3269/s_y2'))[sort_inds], }, }
h5f.close()
l_fel = data['LDM']['l_seed'][0]
ldm_l_ref = data['LDM']['l_ref']
l_ref = np.mean(ldm_l_ref)
l_ref = 266.032
print l_ref

l_fel = 261.0726087 #261.7 #261.0726087


delay = data['LDM']['delay']
Ttheo = np.linspace(data_window[0],data_window[1],1000)

''' Seed laser AC simulation '''
sl_fwhm = 99. # seed intensity fwhm in fs
t_theo_seed = np.linspace(0,200,5000)
sl_AC = np.exp(-4.0*np.log(2)*(t_theo_seed/(sl_fwhm*np.sqrt(2)))**2) # seed laser AC
sl_AC /= np.max(sl_AC)


'''' Fano Resonance profile from Paper '''
#fano parameters
q = -0.17
E_r = 28.509 #eV
gamma = 0.0126 #eV

E_undersamp = harmonic*spc.h*spc.c/(l_ref*spc.e)*1e9 # undersampling energy
E_r = E_r - E_undersamp  #eV. undersampled fano resonance
E_fano = np.linspace(0,E_r+30,500000)

fano_resonance = fk.fano(E_fano,E_r,q,gamma,1.,-1.)
I_fano_theo = -np.conjugate(np.fft.fftshift(np.fft.fft(fano_resonance))*np.exp(1j*np.pi/2))
t_fano_theo = 1e15*spc.h/spc.e* np.fft.fftshift(np.fft.fftfreq(len(fano_resonance),E_fano[1]-E_fano[0]))
#adjusting data window
if data_window[0]==0.:
    low = len(t_fano_theo)/2
else:
    low = max(np.argwhere(t_fano_theo<data_window[0]))
high = min(np.argwhere(t_fano_theo>data_window[1]))
t_fano_theo = t_fano_theo[low:high] 
I_fano_theo = I_fano_theo[low:high]
I_fano_theo = fk.Phasing_TD(I_fano_theo,data_window,E_r) #time domain phasing of theoretical dipole response so it starts with Im(d(t))=-1 at beginning of data_window in case of a lorenzian line shape

   

''' Resolution from Data window '''
print('The resolution from the data window length is {} eV'.format(spc.h/(abs(data_window[-1]-data_window[0])*1E-15*spc.e)))




''' i0 correction for different harmonic contents of i0 '''
i0_fft = data['LDM']['i0_fft']
i0_0H = np.sum(i0_fft[::,395:405],axis=1)
i0_1H = np.sum(i0_fft[::,443:453],axis=1)
i0_2H = np.sum(i0_fft[::,492:502],axis=1)
i0_3H = np.sum(i0_fft[::,540:550],axis=1)
i0_4H = np.sum(i0_fft[::,589:599],axis=1)
i0_5H = np.sum(i0_fft[::,637:647],axis=1)
i0_6H = np.sum(i0_fft[::,685:695],axis=1)
i0_6H /= np.max(i0_6H)



# loop over the mfli devices
for dev in device:
    # loop over the demodulators
    for d in demod:
        mfli_harmonic = data[dev]['harmonic'][int(d)]
        i = str(mfli_harmonic) + 'H'
        Z = (data[dev]['x' + d] + 1j * data[dev]['y' + d])# / data['LDM']['I0']
        if i0_correction:
            Z = Z*(np.nanmean(data['LDM']['I0'])/data['LDM']['I0'])
        if i0_6H_correction:
            Z = Z*i0_6H
        

        Z *= np.exp(1j*63./360.*2.*np.pi) #phasing of dataset according to reference transmitter phase function
        R = np.sqrt(Z.real**2 + Z.imag**2)
        Phi = np.angle(Z, deg=False)
#        plt.plot(np.unwrap(Phi))


# Create theoretical curve
#        if draw_theory:
        Ttheo = np.linspace(delay[0],delay[-1], 30000)
        Xtd,Ytd,Xt,Yt = fk.Curve(spc.h*spc.c/(spc.e*l_trans)*1e9, l_ref, harmonic, phi, a*max(abs(Z)), offset, delay[0], delay[-1], 30000) # theo curve for transition
#        Xtd,Ytd,Xt,Yt = fk.Curve(l_fel/harmonic, l_ref, harmonic, phi, a*max(abs(Z)), offset, delay[0], delay[-1], 30000) # theo curve for laser CWL
        ''' Fitting the dephasing for 6th harmonic R'''
        if dev=='dev3265' and d=='1':
            def decay(t,amp,tau):
                return 0.002+amp*np.exp(-t/tau)
            
            R_fit = R[np.argmin(np.abs(delay-data_window[0])):np.argmin(np.abs(delay-data_window[1]))] #region of R that will be fitted
            delay_fit = delay[np.argmin(np.abs(delay-data_window[0])):np.argmin(np.abs(delay-data_window[1]))]
            popt, pcov = curve_fit(decay, delay_fit, R_fit, p0=[0.8,140.])
            
            #fitting theoretical curve:
            popt_theo, pcov_theo = curve_fit(decay, t_fano_theo, abs(I_fano_theo), p0=[0.8,140.])
            
            print('The time constant of the data fit is {} fs'.format(popt[1]))
            print('The time constant of the theo fit is {} fs'.format(popt_theo[1]))
        
        
        
#        X = Z.real
#        Y = Z.imag

        T_d,X_d,Y_d = fk.CutDataSet(delay, Z.real, Z.imag, data_window)
        Z_cut = X_d + 1j*Y_d
        Z_cut = fk.Phasing_TD(Z_cut,data_window,E_r) #time domain phasing of data set to take into account position of data window

        scaling = 1 #max(Z.real)/max(Ttheo)

        # Fourier traffo
        Td = np.mean(np.unique(np.diff(delay)))
#        print Td
        if gauss:
            Zg = fk.GaussWindow(T_d, Z_cut, suscept)
            wn, dft = fk.DFT(T_d, Zg, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
        else:
            wn, dft = fk.DFT(T_d, Z_cut, Td, l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
#        dft = dft[::-1]
        FWHM , peakFullWidth = fk.PeakWidth(T_d, 0.1, False, zeroPaddingFactor)
        
        #Fouriert Trafo of the theoretical dataset
        wn_theo, dft_theo = fk.DFT(t_fano_theo, I_fano_theo, t_fano_theo[1]-t_fano_theo[0], l_ref , harmonic, zeroPaddingFactor = zeroPaddingFactor)
    
       
        ''' Plotting '''
        figTD = plt.figure('comparison Fano from data and theo',figsize=(14, 10))

        ax = figTD.add_subplot(221)
        ax.set_title('comparison Fano from data and theo')
        ax.plot(delay, Z.real, '-', color=color)
        ax.plot(delay, Z.imag, '-', color='g')
        ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.real)),color='black',alpha=0.3)
        #ax.plot(Ttheo, scaling*Theo, 'k-')
        if plotTheo: ax.plot(Ttheo, Xt, 'b', alpha=0.3)
        ax.set_ylabel('a.u.')
        ax.set_xlabel('delay [fs]')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])
        ax.legend(['X data','Y data'])

        ax = figTD.add_subplot(222)
        #ax.plot(Ttheo, scaling*Theo, 'k-')
        ax.plot(t_fano_theo, I_fano_theo.real/max(abs(I_fano_theo.real)), '-', color=color)
        ax.plot(T_d, Z_cut.real/max(abs(Z_cut.real)), '-', color='g')
#        ax.plot(t_fano_theo, I_fano_theo.imag, '-', color='g')
        ax.fill_between(t_theo_seed,sl_AC*np.max(abs(Z.imag)),color='black',alpha=0.3)
        if plotTheo: ax.plot(Ttheo, Yt, 'r', alpha=0.3)
        ax.set_ylabel('a.u.')
        ax.set_xlabel('delay [fs]')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])
        ax.legend(['X theo','Y theo'])

        ax = figTD.add_subplot(223)
        ax.plot(delay, R/max(R), '-', color=color)
        ax.plot(t_fano_theo, abs(I_fano_theo)/max(abs(I_fano_theo)), '-', color='g')
        ax.fill_between(t_theo_seed,sl_AC*np.max(R),color='black',alpha=0.3)#,'-', color='black')
        if dev=='dev3265' and d=='1':
            ax.plot(delay_fit,decay(delay_fit,popt[0],popt[1])/max(R), '--',color='black', linewidth=2)
            ax.plot(t_fano_theo,decay(t_fano_theo,popt_theo[0],popt_theo[1])/max(abs(I_fano_theo)), '--',color='black', linewidth=2)
        ax.set_ylabel('R')
        ax.set_xlabel('delay [fs]')
        ax.grid()
        ax.set_xlim(data_window[0],data_window[-1])
        ax.legend(['R data','R theo','fits'])

        
        '''plot frequency domain'''
        
        if plot_eV:
            factor = spc.h*spc.c/spc.e*100.
            wn *= factor 
            wn_theo *= factor
            wn_lim *= factor
            l_ref /= factor
            l_trans /= factor
            FWHM *= factor
            print('The resolution of the DFT (after zeropadding) is {} eV'.format(FWHM))
        undersampled_freq = wn-harmonic*1E7/l_ref        
#        dft = dft**2 #squaring to get intensity insted of amplitude
#        dft_theo = dft_theo**2
        
      
        ax = figTD.add_subplot(224)
        ax.plot(wn, dft.real/max(abs(dft)), '-', color='b', linewidth = 2)
        ax.plot(wn, dft.imag/max(abs(dft)), '-', color='g', linewidth = 2)
#        ax.plot(wn, dft.real/max(abs(dft)), '-', color='magenta', linewidth = 2)
        ax.plot(wn, abs(dft)/max(abs(dft)), '-', color='black', linewidth = 2)
#        ax.plot(wn_theo, abs(dft_theo)/max(abs(dft_theo)), 'o-', color='g', linewidth = 2)
        ax.plot(wn_theo, dft_theo.imag/max(abs(dft_theo)), '-', color='r', linewidth = 2)
#        ax.plot(wn_theo, dft_theo.real/max(abs(dft_theo)), '-', color='b', linewidth = 2)
#        ax.plot(t_fano,(I_fano+0.03)/np.max(abs(I_fano+0.03)), color='black',linewidth = 2)
        ax.plot(E_fano+E_undersamp,fano_resonance,color='grey',linewidth = 2)
        #add second x axis to plot:
#        ax2 = ax.twiny()
#        ax2.set_xlabel('rotating frame energy [eV]')         
#        ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)
#        axs.plot(padres_roi,padres_spec,color='grey',linewidth = 2) #FEL spectrum
        if plot_eV:
            ax.set_xlabel(r'energy [eV]')
        else:
            ax.set_xlabel(r'wavenumber [cm$^{-1}$]')
        ax.set_ylim(-1.1,1.1)
        ax.set_ylabel('intensity [a.u.]')
        ax.set_xlim(wn_lim[0],wn_lim[1])
        ax.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
        plot_x_range = np.diff(ax.get_xlim())[0]
        plot_y_range = np.diff(ax.get_ylim())[0]
        ax.axvspan(l_trans-FWHM/2, l_trans + FWHM/2, ymin=0.1, ymax=0.2,alpha=0.3)
        ax.legend(['data.real','data.imag','data.abs','theo.dft.imag','theo.fd'])

        
        
        plt.tight_layout()

        

        

if not interactive_plots:
    plt.close('all')
else:
    plt.show()
    

figTD = plt.figure("zoom",figsize=(7,6))

ax = figTD.add_subplot(111)
#ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
ax.plot(delay, Z.real, 'o-', color='b')
ax.plot(delay, Z.imag, 'o-', color='g')
ax.plot(Ttheo, Xt*0.34, 'b', alpha=0.3)
ax.legend(['real','imag','real_theo'])
ax.set_xlim(-20,20)
ax.grid()



#fig = plt.figure()
#
#ax = fig.add_subplot(111)
##ax.plot(abs(np.fft.fftshift(np.fft.fft(Z))))
#dft2 = abs(np.fft.fftshift(np.fft.fft(Z,n=256)))
#freq2 = 1e15*spc.h/spc.e*np.fft.fftshift(np.fft.fftfreq(len(dft2),delay[1]-delay[0]))
#ax.plot(freq2,dft2)


''' plots for paper '''
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


figTD = plt.figure(figsize=(7,6))

ax = figTD.add_subplot(211)
#ax.errorbar(delay, X, yerr=X_s, color='b', linestyle='')
ax.plot(delay, Z.imag, '-', color=color,alpha=0.7)
#ax.plot(t_fano_theo, I_fano_theo.imag/max(I_fano_theo.real)*max(abs(X)), '-', color='g',alpha=0.7)
#ax.plot(delay, Y, '-', color='g')
#ax.fill_between(t_theo_seed,sl_AC*np.max(abs(X))*0.8,color='grey')
#ax.fill_between(t_theo_seed,-sl_AC*np.max(abs(X))*0.8,color='grey')
#ax.plot(Ttheo, scaling*Theo, 'k-')
if plotTheo: ax.plot(Ttheo, Xt, 'b', alpha=0.3)
ax.set_ylabel('Im(d) [a.u.]')
ax.set_xlabel(r'$\tau$ [fs]')
ax.set_xlim(data_window[0],data_window[-1])


'''plot frequency domain'''
undersampled_freq = wn-harmonic*1E7/l_ref        
  
ax = figTD.add_subplot(212)
ax.plot(wn, dft.real/max(abs(dft)), '-', color=color, alpha=0.7)
#ax.plot(wn_theo, abs(dft_theo)/max(abs(dft_theo)), '-', color='g',alpha=0.7)
ax.plot(E_fano+E_undersamp,fano_resonance,color='black', alpha=0.7)
ax2 = ax.twiny() #add second x axis to plot:
ax2.set_xlabel('rotating frame energy [eV]')         
ax2.set_xlim(wn_lim[0]-harmonic*1E7/l_ref,wn_lim[1]-harmonic*1E7/l_ref)
ax.set_xlabel(r'energy [eV]')
ax.set_ylim(-1.1,0.2)
ax.set_ylabel('intensity [a.u.]')
ax.set_xlim(wn_lim[0],wn_lim[1])
ax.axvline(harmonic*1E7/l_ref,color = 'k',linestyle='--', linewidth= 2)
plot_x_range = np.diff(ax.get_xlim())[0]
plot_y_range = np.diff(ax.get_ylim())[0]
ax.axvspan(l_trans-FWHM/2, l_trans + FWHM/2, ymin=0.1, ymax=0.2,alpha=0.3)
#ax.legend(['exp','theo'],loc=4)



plt.tight_layout()