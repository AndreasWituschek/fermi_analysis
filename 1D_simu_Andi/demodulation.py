# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:29:08 2020

@author: andreas

This code takes simulation from Matlab script "num_sim_1D". It demodulates the
phase cycled population at each delay step and gives the demodulated amplitude 
and phase at the harmonics of the modulation frequency.
"""

import scipy.io as spio
import numpy as np
import matplotlib.pyplot as plt

ticksize= 2.
ticklength = 5.
fontsize=13.
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.major.width'] = ticksize
plt.rcParams['xtick.major.size'] = ticklength
plt.rcParams['ytick.major.width'] = ticksize
plt.rcParams['ytick.major.size'] = ticklength
plt.rcParams['axes.linewidth'] = ticksize
plt.rcParams['lines.linewidth'] = 1.8
plt.rcParams['legend.fontsize'] = 13


# =============================================================================
# parameters
# =============================================================================
harm = 5 #FEL harmonic
tile = 100 # how many cycles of 2pi are tiled 
plot_pop =0 #plot populations?
plot_paper = 1
plot_R0_5H = 0 #plot the R values from 0-5H vs tau
# =============================================================================
# #importing multiple population scans for saturated HGHG
# =============================================================================
#
base_path = '//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/'
base_name = '2020-04-03_amp105_'

file_pulse = base_name + 'pulse'

tau_steps_nr = len(spio.loadmat(base_path + base_name + '1.mat', squeeze_me = True)['Popi']) #number of tau steps used in matlab script (=max_tau12 + 1 )

scan_nr = np.arange(0,20) #number pf phase cycling steps (defined in matlab script)
scans = np.zeros((len(scan_nr),tau_steps_nr)) #contains all the delay scans for all phase cycling steps
for n in scan_nr:
    mat = spio.loadmat(base_path + base_name + str(n) + '.mat', squeeze_me = True) # data population
    if n == 0: tau = mat['ttau12']
    scans[n] = mat['Popi']

mat2 = spio.loadmat(base_path + file_pulse + '.mat', squeeze_me = True) # simulated laser pulse
lambda_0 = mat2['WDp']# carrier wavelength in nm


if len(base_name) == 17: 
    amp = int(base_name[-3:-1])/100
else: amp = int(base_name[-4:-1])/100
print(amp)

# =============================================================================
# Demodulation of phase cycling
# =============================================================================
RvsTau = np.zeros((harm+1,len(tau))) #this vector contains amp of dft at every delay step for each harmonic demodulation channel
XvsTau = np.zeros((harm+1,len(tau)))
for i in np.arange(len(tau)):
    dft_t = np.fft.rfft(np.tile(scans[::,i],tile))#dft of tiled phase modulated data
    dft_x = np.angle(dft_t)
    dft_t = np.abs(dft_t)
    for h in np.arange(0,harm+1):
        RvsTau[h,i] = dft_t[h*tile] #picking dft amp at the harmonics
        XvsTau[h,i] = dft_x[h*tile]

# =============================================================================
# #importing multiple population scans for unsaturated HGHG
# =============================================================================

base_name_u = '2020-02-17_amp199_'

#file_pulse = base_name_u + 'pulse'


tau_steps_nr = len(spio.loadmat(base_path + base_name_u + '1.mat', squeeze_me = True)['Popi'])

scan_nr = np.arange(0,20) #number pf phase cycling steps (defined in matlab script)
scans_unsat = np.zeros((len(scan_nr),tau_steps_nr)) #contains all the delay scans for all phase cycling steps
for n in scan_nr:
    mat = spio.loadmat(base_path + base_name_u + str(n) + '.mat', squeeze_me = True) # data population
    if n == 0: tau = mat['ttau12']
    scans_unsat[n] = mat['Popi']

#mat2 = spio.loadmat(base_path + file_pulse + '.mat', squeeze_me = True) # simulated laser pulse
#lambda_0 = mat2['WDp']# carrier wavelength in nm
#
#amp_u = int(base_name[-4:-1])/100
##amp = 1/100

# =============================================================================
# Demodulation of phase cycling
# =============================================================================
RvsTau_unsat = np.zeros((harm+1,len(tau))) #this vector contains amp of dft at every delay step for each harmonic demodulation channel

for i in np.arange(len(tau)):
    dft_t = np.fft.rfft(np.tile(scans_unsat[::,i],tile))#dft of tiled phase modulated data
    dft_t = np.abs(dft_t)
    for h in np.arange(0,harm+1):
        RvsTau_unsat[h,i] = dft_t[h*100] #picking dft amp at the harmonics
        


# =============================================================================
# importing data from Scan_31
# =============================================================================
X = np.zeros((harm,657))
Y = np.zeros((harm,657))
for i in [0,1,2,3,4]:
    XX = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/X_vs_delay_' +str(i+1) + '.txt').transpose()
    YY = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_31/Y_vs_delay_' +str(i+1) + '.txt').transpose()
    X[i] = XX
    Y[i] = YY
Z= X + 1j*Y
tau_data = np.linspace(-310.0393465165588,1001.9601165780342,657)


# =============================================================================
# plotting
# =============================================================================
plt.close('all')
mean_pop = np.mean(scans,axis=0)
if plot_pop:
    plt.figure('populations')
    for i in scan_nr:
        plt.plot(tau,scans[i],label=i)
    plt.plot(tau,mean_pop,label = 'mean',linewidth = 3,color = 'k')
    plt.legend()
    plt.xlabel('delay [fs]')
    plt.ylabel('population')
    plt.title('amp= {}, l0= {}nm'.format(amp,lambda_0))
    plt.tight_layout()
    
#    plt.figure('population at delay vs. laboratory time')
#    plt.plot(np.tile(scans[::,80],4),'o-') # one random delay time
#    plt.xlabel('laboratory time [arb. u.]')
#    plt.ylabel('amplitude [arb. u.]')
#    plt.title('amp = {}'.format(amp))
#    plt.tight_layout()


if plot_R0_5H:
    plt.figure('R_0H-5H vs delay' +'_amp= {}, l0= {}nm'.format(amp,lambda_0))
    for i in np.roll(np.arange(0,harm+1),-1):
        plt.plot(tau,RvsTau[i],label='{}H'.format(i))
    plt.legend(loc=1)
    plt.xlabel('delay [fs]')
    plt.ylabel('modulation amplitude [arb. u.]')
    plt.title('amp= {}, l0= {}nm'.format(amp,lambda_0))
    plt.xlim([-300,300])
    plt.tight_layout()

# =============================================================================
# plotting data and simulations alongside
# =============================================================================

if plot_paper:
    colorcycle = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
    norm = max(abs(Z[-1]))*0.55 #normation on 5H amplitude
    Z /= norm
    norm_sim = max(abs(RvsTau[-1])) #normation on 5H amplitude
    RvsTau /= norm_sim
    norm_sim_unsat = max(abs(RvsTau_unsat[-1])) #normation on 5H amplitude
    RvsTau_unsat /= norm_sim_unsat
    
    nrow = int(harm)
    ncol = 1
    fig, axs = plt.subplots(nrow, ncol,figsize=(6,10.5), sharex = 'col',num=base_name[:-1])
#    ylim = [[-1.1,1.1],[-0.25,.25],[-0.15,0.15],[-0.11,0.11],[-0.04,0.04],[-0.03,0.03]]/(0.7*norm_theo)
    xlim =[-300,300]
    mult = [160,120,60,15,1]
    print(RvsTau.shape)
    for i in range(harm):
        axs[0].set_title('file= {}'.format(base_name[:-1]))
        ax = axs[i]
        XvsTau[i+1] = np.exp(1j*tau*0.15*(i+1)+1j*XvsTau[i+1])*RvsTau[i+1]
        ax.plot(tau_data,np.abs(Z[i]), label='{}H'.format(i+1), color=colorcycle[i])
        ax.plot(tau, RvsTau[i+1],color='k',linestyle='dotted')#,label = 'sat. calc.')
#        ax.plot(tau,XvsTau[i+1],color='blue')
        ax.plot(tau, RvsTau_unsat[i+1]/mult[i],color='k',linestyle = 'dashed')#,label = 'unsat. calc./{}'.format(mult[i]))
        ax.legend(loc = 'upper left')
        ax.set_xlim(xlim)
#       ax.set_ylim(ylim[i])
        if i == range(harm)[-1]: ax.set_xlabel('delay [fs]') 
#        if i==0: ax.set_title('amp= {}, l0= {}nm'.format(amp,lambda_0))
    plt.tight_layout()



# =============================================================================
# tests            
# =============================================================================

interp = np.interp(np.linspace(0,len(mean_pop),14001),np.arange(len(mean_pop)),mean_pop)

plt.savefig('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/'+ base_name)
