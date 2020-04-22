# -*- coding: utf-8 -*-
"""
Created on Tue Oct 01 13:59:45 2019

@author: andreas

Simulation of the nth order autocorrelation with saturation in the harmonic process.

"""

#plt.close('all')
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz, savgol_filter
from scipy.special import jv    
import time


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
plt.rcParams['lines.linewidth'] = 1.5


#fundamental seed pulse E field with unity amplitude (1H)
# freq = frequency in Hz, fwhm = intensity fwhm
def pulse(t,freq,cen,fwhm,amp,phase):
    return amp*np.exp(-2*np.log(2)*(t-cen)**2/fwhm**2 + 1j*(2.*np.pi*freq*(t-cen) + phase))
    
def pulse_steady(t,freq,cen,fwhm,chirp,amp,phase):
    return amp*np.exp(-2*np.log(2)*(t-cen)**2/fwhm**2 + 1j*(2.*np.pi*freq*(t-cen)+phase) + 1j*chirp*(t-cen)**2)

def pulse_scanned(t,freq,cen,fwhm,chirp,amp,phase,freqshift):
    return amp*np.exp(-2*np.log(2)*(t-cen)**2/fwhm**2 + 1j*(2.*np.pi*(freq + cen*freqshift)*(t-cen)+phase) + 1j*chirp*cen*(t-cen)**2)

#OLD SIGMPOID MODE DONT USE ANYMORE
#def hg(pulse,harmonic):
#    return np.real(pulse**harmonic/(1+np.abs(pulse**harmonic))), pulse**harmonic


#BESSEL MODEL
def hg(pulse,harmonic):
    seedenvelope = np.abs(pulse)
    seedphase = np.unwrap(np.angle(pulse))
    return jv(harmonic,harmonic* seedenvelope)*np.exp(1j*seedphase*harmonic), pulse**harmonic

#NEW SIGMOID MODEL
#def hg(pulse,harmonic):
#    xmax = (1+np.sqrt(2/3)*harmonic**(-2/3))
#    jmax = jv(harmonic,harmonic*xmax)
#    f= np.sqrt(2)/xmax
#    return np.real((f*pulse)**harmonic/(1+np.abs((f*pulse)**harmonic)))*jmax, pulse**harmonic
#

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def idft_window(dft_ac,cent,sigma,harm,points,order):
    points  = int(points)
    b, a = butter_bandpass(harm*cent-sigma, harm*cent+sigma, points, order=order)
    w, h = freqz(b, a, worN= int(points/2))
    window = np.abs(np.append(np.roll(np.flip(h),1),h))
    return np.fft.ifft(np.fft.ifftshift(dft_ac*window)), window

def idft_window_undersamp(dft_ac,cent,sigma,harm,points,order):
    points  = int(len(tau))
    b, a = butter_bandpass(harm*cent-sigma, harm*cent+sigma, points, order=order)
    w, h = freqz(b, a, worN= int(points/2))
    window = np.abs(np.append(np.roll(np.flip(h),1),h)) #window over respective harmonic
    
    dft_ac_window = dft_ac*abs(window) # selects one harmonic
    daw_low = dft_ac_window[:int(len(dft_ac_window)/2)] #divides dft_ac_window in a positive, zero and negative frequency part
    daw_zero = dft_ac_window[int(len(dft_ac_window)/2)]
    daw_high = dft_ac_window[int(len(dft_ac_window)/2+1):]
    idx = np.argmin(np.abs(freqs_ac-harm))-len(dft_ac_window)/2 #finds index in frequency axis where freqquency equals 1*freq
#    print('idx high ={}'.format(idx*(1-1/undersamp)))
    shift = int(idx*(1-1/undersamp)) 

    phase = np.exp(-1j*2*np.pi*shift/1.99987) #compensating phase shift originating from shifting the components to lower frequency
#    print('angle = {}'.format(np.angle(phase,deg=1)))
    daw_high_roll = np.roll(daw_high,-shift) * phase #shifts the positive and negative components to lower frequencies (=undersampling)
    daw_low_roll = np.roll(daw_low,shift) *phase.conjugate()
    daw_roll = np.append(np.append(daw_low_roll,daw_zero),daw_high_roll)
    return np.fft.fft(np.fft.ifftshift(daw_roll)), window
#==============================================================================
# 
#==============================================================================

plot_ac = 0
plot_efield = 0
plot_components = 0
plot_data_sim = 0
plot_coherent_addup = 0
plot_hgfunction = 0
plot_paper = 1

#parameters
t = np.linspace(-600,600,4000) #time variable for pulses (ca. 40000 is good to reslove fast oscillations)
tau=np.linspace(-400,400,2**14) #delay range for autocorrelation

amp = 0.83 #amplitude of seed field BESSEL
#amp = 0.979 #amplitude of seed field SIGMOID

freq = 1.1492 #[1/fs] seed laser frequency
undersamp = 52.23 #undersampling factor. If 0, no undersampling will be applied
chirp = 0# -freq*0.000005 #chirp parameter
phase = 0.0 #rad
harm = 6 #order of the harmonic generation process
fwhm = 98 #pulse intensity fwhm in fs
#fwhm = 120
freqshift =0# freq*0.00001 #frequency shift of the scanned pulse while scanning in Hz/fs


#==============================================================================
# simulation
#==============================================================================

start = time.time()
# ====== integrating AC integral for different tau vals
ac, ac_sat = np.zeros_like(tau), np.zeros_like(tau)

steady = pulse_steady(t,freq,0,fwhm,0,amp,phase)

for n,step in enumerate(tau):
    pulse_5H_sat, pulse_5H  = hg(steady +pulse_scanned(t,freq,step,fwhm,chirp,amp,phase,freqshift),harm) #5H fields   
#    ac[n] = np.sum(np.real(pulse_5H)**2)
#    ac_sat[n] = np.sum(np.real(pulse_5H_sat)**2)
    ac[n] = np.sum(np.abs(pulse_5H)**2)
    ac_sat[n] = np.sum(np.abs(pulse_5H_sat)**2)
print('loop duration {}s'.format(time.time()-start))


# ======== DFT of AC traces:
dft_ac_sat = np.fft.fftshift(np.fft.fft(ac_sat))
dft_ac_sat /= np.max(abs(dft_ac_sat))
dft_ac = np.fft.fftshift(np.fft.fft(ac))
dft_ac /= np.max(abs(dft_ac))
freqs_ac = np.fft.fftshift(np.fft.fftfreq(len(ac),d=tau[1]-tau[0]))/freq #DFT frequency axis


#creating parameters for butterworth filter
polyfit = np.polyfit(np.arange(len(freqs_ac)),freqs_ac-1,deg=1)
cent = -polyfit[1] / polyfit[0]-len(freqs_ac)/2 
sigma = cent*.02
#sigma = cent*.32
order = 4

print('This value must be < 1: {}'.format((harm*cent+sigma)/(0.5*len(tau))))

# ======== IDFT of AC to get harmonic components
Z_theo_sat = np.array([]) #contains the individual saturated harmonic demodulated AC Z components
Z_theo = np.array([]) #contains the individual UNsaturated harmonic demodulated AC Z components
windows = np.array([]) #windows used for IDFT in each harmoinic
harmonics = [1,2,3,4,5,6]


for n in harmonics:
    if undersamp != 0: idft, window = idft_window_undersamp(dft_ac,cent,sigma,n,len(tau),order)
    else: idft, window = idft_window(dft_ac,cent,sigma,n,len(tau),order)
    Z_theo = np.append(Z_theo,idft)
    if undersamp != 0: idft_sat, window  = idft_window_undersamp(dft_ac_sat,cent,sigma,n,len(tau),order)
    else: idft_sat, window = idft_window(dft_ac_sat,cent,sigma,n,len(tau),order)
    Z_theo_sat = np.append(Z_theo_sat,idft_sat)
    windows = np.append(windows,window)
Z_theo = Z_theo.reshape((-1,int(len(Z_theo)/len(harmonics)))).transpose()
Z_theo_sat = Z_theo_sat.reshape((-1,int(len(Z_theo_sat)/len(harmonics)))).transpose()
windows = windows.reshape((-1,int(len(windows)/len(harmonics)))).transpose()  

#==============================================================================
# importing X values from all 6 harmonics from scan5234
#==============================================================================
X = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/X_vs_delay.txt').transpose()
Y = np.loadtxt('//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/ScansOverZero/scan_5234/Y_vs_delay.txt').transpose()
Z = X + 1j*Y
delay_min = -316.66025 #fs
delay_max = 299.3396
delay = np.linspace(delay_min,delay_max,len(X[1]))

#phasing of experimental data to fit the phase of the simulated AC traces
phas = np.array([1.1,-0.45,0.25,0.95,1.5,0.])*np.pi
phas = np.exp(1j*phas)
for i, z in enumerate(Z):
    Z[i] = phas[i]*z

#==============================================================================
# evaluating frequency of data to obtain frequency for simulation
#==============================================================================
#freq_guess = 0.022 #guessing frequency of demodulated AC 1H signal
#plt.close('DFT of X data')
#figDFTX = plt.figure('DFT of X data',figsize=(9,10))
#for i, xval in enumerate(X):
#    dft_xval = np.fft.fftshift(np.fft.fft(xval))
#    freqs_xval = np.fft.fftshift(np.fft.fftfreq(len(xval),d=delay[1]-delay[0]))
#    ax = figDFTX.add_subplot(111)
#    ax.plot(freqs_xval,abs(dft_xval)/max(abs(dft_xval)))
#    ax.set_xlim([0,0.25])
#    ax.axvline(freq_guess*(i+1),color='black')

#==============================================================================
# Plotting
#==============================================================================

#time and frequency domain reperesentations of the harmonic autocorrelation
if plot_ac:
    figAC = plt.figure('Autocorrelation',figsize=(8,9))
    ax = figAC.add_subplot(211)
    ax.set_title('Time Domain')
    ax.plot(tau,ac_sat,alpha=0.7)
    ax.plot(tau,ac,alpha=0.7)
    #ax.plot(tau,ac_sat_conv)
    ax.set_xlabel('delay')
    ax.set_ylabel('intensity [a.u.]')
    ax.legend(['saturated','non-saturated'])
    ax.set_xlim([min(tau),max(tau)])
    
    ax2 = figAC.add_subplot(212)
    ax2.set_title('Frequency Domain')
    ax2.plot(freqs_ac,abs(dft_ac_sat),alpha=0.7)
    ax2.plot(freqs_ac,abs(dft_ac),alpha=0.7)
    ax2.plot(freqs_ac,windows,color='silver')
    ax2.set_xlabel('freq [f_0]')
    ax2.set_ylabel('intensity [a.u.]')
    ax2.legend(['saturated','non-saturated'])
    ax2.set_xlim([-0.5,8])
    plt.tight_layout()

#plotting the electric fields
if plot_efield:
    pulse_1H_1 = pulse_steady(t,freq,0,fwhm,chirp,amp,phase)
    pulse_1H_2 = pulse_scanned(t,freq,max(tau),fwhm,chirp,amp,phase,freqshift)
    plt.figure('E-Fields')
    plt.plot(t,(pulse_1H_1+pulse_1H_2),alpha=0.7)
    plt.plot(t,pulse_5H_sat,alpha=0.7)
    plt.plot(t,pulse_5H,alpha=0.7)
    plt.xlabel('t')
    plt.ylabel('E-field [a.u.]')
    plt.legend(['seed','saturated','non-saturated'])
    plt.tight_layout()

#plot individual harmonic TD AC components for saturated and unsaturated case
if plot_components:
    figIDFT = plt.figure('inverse DFT of AC components',figsize=(9,10))
    ax = figIDFT.add_subplot(211)
    ax.plot(tau,Z_theo,alpha=0.7,label='{}H'.format(n))
    ax.legend()
    ax.set_title('non-saturated HGHG')
    ax.set_xlim([min(tau),max(tau)])
    ax.axvline(0,color='black')
    
    ax2 = figIDFT.add_subplot(212)
    ax2.plot(tau,Z_theo_sat,alpha=0.7,label='{}H'.format(n))
    ax2.legend()
    ax2.set_xlabel('delay [a.u.]')
    ax2.set_title('saturated HGHG')
    ax2.set_xlim([min(tau),max(tau)])
    ax2.axvline(0,color='black')
    plt.tight_layout()
    

# plotting simulation vs. data
if plot_data_sim:
    colorcycle = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
    #data from scan5234
    nrow = int(harm)
    ncol = 1
    fig, axs = plt.subplots(nrow, ncol,figsize=(6,9))
    norm_data = max(abs(Z[0])) #normation on 1H
    for i, ax in enumerate(fig.axes):
        ax.plot(delay,Z[i]/norm_data,label='{}H'.format(i+1),color=colorcycle[i])
        ax.get_xaxis().set_visible(False)
        if i==nrow-1:
            ax.get_xaxis().set_visible(True)
        if i == 0:
            ax.set_title('data')
        ax.legend()
        ax.set_xlim([-300,300])
    #    ax.axvline(0,color='black')
        ax.set_xlabel('delay [fs]')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
#        
    #plotting saturated AC components
    fig, axs = plt.subplots(nrow, ncol,figsize=(6,9))
    norm_theo = max(abs(Z_theo_sat[::,0])) #normation on 1H
    for i, ax in enumerate(fig.axes):
        ax.plot(tau,Z_theo_sat[::,i]/norm_theo,label='{}H'.format(i+1),color=colorcycle[i])
        ax.get_xaxis().set_visible(False)
        if i==nrow-1:
            ax.get_xaxis().set_visible(True)
        if i == 0:
            ax.set_title('saturated amp= {}'.format(amp))
        ax.legend()
        ax.set_xlim([-300,300])
    #    ax.axvline(0,color='black')
        ax.set_xlabel('delay [fs]')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    #plotting unsaturated AC components
    fig, axs = plt.subplots(nrow, ncol,figsize=(6,9))
    norm_theo = max(abs(Z_theo[::,0])) #normation on 1H
    for i, ax in enumerate(fig.axes):
        ax.plot(tau,Z_theo[::,i]/norm_theo,label='{}H'.format(i+1),color=colorcycle[i])
        ax.get_xaxis().set_visible(False)
        if i==nrow-1:
            ax.get_xaxis().set_visible(True)
        if i == 0:
            ax.set_title('non-saturated')
        ax.legend()
        ax.set_xlim([-300,300])
    #    ax.axvline(0,color='black')
        ax.set_xlabel('delay [fs]')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

#for i in range(harm):
#    plt.plot(tau,Z_theo[::,i]/max(abs(Z_theo[::,i])))


# =============================================================================
# plots for paper 
# =============================================================================
if plot_paper:
    colorcycle = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']

    Z_theo /= max(abs(Z_theo[::,0]))
    Z_theo_sat /= max(abs(Z_theo_sat[::,0]))
    Z /= max(abs(Z[0])) #normation of data on data_1H

     #normation on 1H non-saturated simulation
    ylim = [[-1.1,1.1],[-0.3,.3],[-0.4,0.4],[-0.11,0.11],[-0.12,0.12],[-0.025,0.025]] #for BESSEL
    ylim = [[-1.1,1.1],[-0.3,.3],[-0.2,0.2],[-0.11,0.11],[-0.04,0.04],[-0.03,0.03]] #for SIGMOID
    xlim =[-200,200]

    nrow = int(harm)
    ncol = 3
    fig, axs = plt.subplots(nrow, ncol,figsize=(11,10), sharex = 'col')
   
    for col in range(ncol):
        if col == 0:
            for i in range(harm):
                ax = axs[i,col]
                ax.plot(tau, Z_theo[::,i],label='{}H'.format(i+1),color=colorcycle[i])
                ax.set_xlim(xlim)
                if i == 0: ax.set_ylim([-1.1,1.1])
                ax.legend(loc=1)
                if i == harm -1: ax.set_xlabel('delay [fs]')
#                ax.axvline(0,color='black')

        if col == 1:
            for i in range(harm):
                ax = axs[i,col]
                ax.plot(tau,Z_theo_sat[::,i], label='{}H'.format(i+1), color=colorcycle[i])
#                ax.legend(loc=1)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim[i])
                if i == harm -1: ax.set_xlabel('delay [fs]') 
                
#                ax.axvline(0,color='black')
        if col == 2:
            for i in range(harm):
                ax = axs[i,col]
                ax.plot(delay,Z[i],label='{}H'.format(i+1),color=colorcycle[i])     #data from scan5234
#                ax.legend(loc=1)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim[i])
                if i == harm -1: ax.set_xlabel('delay [fs]')
                ax.set_yticklabels([])

         
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)
    plt.subplots_adjust(wspace=0.15)

# =============================================================================
# 
# =============================================================================
# coherently adding up all non-DC AC components of data and simulation   
if plot_coherent_addup:
    sum_exp = np.sum(Z,axis=0)
    sum_theo = np.sum(Z_theo_sat,axis=1)
    plt.figure('ac experimental')
    plt.plot(delay,np.abs(sum_exp)**2/np.max(np.abs(sum_exp)**2))
    plt.plot(tau,np.abs(sum_theo)**2/np.max(np.abs(sum_theo)**2))
    plt.xlim([-300,300])

#plotting the HG function
if plot_hgfunction:
    plt.figure('HG function')
    xaxis = np.linspace(-2,2,1000)
    plt.plot(xaxis,hg(xaxis,harm)[0],label='sat HG function')
#    plt.plot(xaxis,bessel(xaxis,harm)[0],label ='bessel')
    plt.axvline(x=amp,color='k')
    #plt.plot(xaxis, hg(xaxis,harm)[1],label='x^6')
    plt.legend()
        


#==============================================================================
# What happens to a sinusoidal signal when saturating it?
#==============================================================================

#def sinus(t,amp,freq,phase):
#    return amp*np.exp(1j*(freq*t+phase))
#
#t = np.linspace(-10,10,200000)
#freq = 1  
#sin1 = sinus(t,.4,freq,0)
#sin2 = sinus(t,.4,freq*1.11,0)
#sin = sin1+sin2
##sin = 2*pulse_chirp(t,freq,0,3,0)
#
#sin_sat5H, sin_5H = hg(sin,5.)
##plotting the electric fields
#plt.figure('Sinus Saturation')
#plt.plot(t,sin_5H,alpha=0.7)
#plt.plot(t,sin_sat5H,alpha=0.7)
#
#dft_sin_sat = np.fft.fftshift(np.fft.fft(sin_sat5H))
#dft_sin_sat /= np.max(abs(dft_sin_sat))
#dft_sin = np.fft.fftshift(np.fft.fft(sin_5H))
#dft_sin /= np.max(abs(dft_sin))
##dft_ac_sat_conv = np.fft.fftshift(np.fft.fft(ac_sat_conv))
##dft_ac_sat_conv /= np.max(abs(dft_ac_sat_conv))
#freqs_sin = np.fft.fftshift(np.fft.fftfreq(len(sin_sat5H),d=t[1]-t[0]))/(freq/(2.*np.pi))
#
#plt.figure('DFT sinus saturated')
#plt.plot(freqs_sin,abs(dft_sin))
#plt.plot(freqs_sin, abs(dft_sin_sat))
#plt.xlim([-0.5,40])


# =============================================================================
#testing heterodyned detection     
# =============================================================================

#def sinus(t,amp,freq,phase):
#    return amp*np.sin((freq*t+phase))
#
#t = np.linspace(-10,30,2000)
#freq = 1  
#sin1 = sinus(t,1,freq,0)
#sin2 = sinus(t,1,freq*1.11,0)
#mult = sin1*sin2
#
#plt.figure()
#plt.plot(t,sin1,label = freq)
#plt.plot(t,sin2, label = freq*1.11)
#plt.plot(t,mult, label = 'mult')
#plt.legend()
#
#dft_mult = np.fft.fftshift(np.fft.fft(mult))
#
#
#plt.figure()
#plt.plot(np.abs(dft_mult))
#
#
#
#dft_ac_sat /= np.max(abs(dft_ac_sat))
#dft_ac = np.fft.fftshift(np.fft.fft(ac))
#dft_ac /= np.max(abs(dft_ac))
#freqs_ac = np.fft.fftshift(np.fft.fftfreq(len(ac),d=tau[1]-tau[0]))

# =============================================================================
# bessel function tests
# =============================================================================
#    
    
#def bsl(pulse,harmonic):
##    return 0.354141 * jv(harmonic,harmonic*pulse), pulse**harmonic
#    max_x_jv= harmonic*(1+np.sqrt(2/3)*harmonic**(-2/3))
#    max_jv = jv(harmonic,max_x_jv)
#    ones =np.ones_like(pulse)
#    ones[np.argwhere(pulse<0)] = -1
#    return jv(harmonic,-harmonic* pulse), np.real(pulse**harmonic)
#
#    
#h=5 
#    
#time = np.linspace(0,2*np.pi,100)
#exp = np.exp(1j*time)
#bessel = jv(h,h*exp)
#bslown, _ = bsl(exp,h)
#plt.figure()
##plt.plot(time,exp.real)
##plt.plot(time,exp.imag)
#plt.plot(time,bessel.real)
#plt.plot(time,bessel.imag)
#plt.plot(time,bslown.real,label = 'own')
#plt.plot(time,bslown.imag,label = 'own')
#plt.legend()
