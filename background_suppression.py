import numpy as np
import matplotlib.pyplot as plt
import h5py

#run_remove = [226, 227, 292, 293, 294, 308, ]
#h5fd = h5py.File('/home/daniel/Documents/Datasets/Fermi2019/scan_031/scan_031_correct_delay.h5', 'r')
# scan 4651 -----------------
#run = 4881
# ---------------------------

#h5fd = h5py.File('//10.5.71.1/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/scan_4651_digitizer.h5', 'r')
#h5fd = h5py.File('//10.5.71.28/FermiServer/Beamtime2/Run_4781/rawdata/Run_4781_274161561.h5', 'r') #200fs delay
h5fd = h5py.File('//10.5.71.28/FermiServer/Beamtime2/Run_4881/rawdata/Run_4881_274315077.h5', 'r') #400fs delay
signal = np.array(h5fd['/digitizer/channel1'])
h5fd.close()

baseline = np.mean(signal[::,20000:])
signal = signal - baseline
signal = signal [::,15800:16300]
plt.plot(signal[1])
integral = np.sum(signal,axis=1)
#plt.plot(integral)

#signal = rt.butter_highpass_filter(signal, 6.0, 50., order=2)

print(integral.shape)

x = np.arange(400)

fig, ax = plt.subplots(2, 1, sharex=False, figsize=(9,12))

ax[0].plot(x, integral)
ax[0].set_xlim(-1, 401)
ax[0].set_xlabel('Shots Count')
ax[0].set_ylabel('Sum over the window')

fx = np.fft.fftfreq(integral.size, d=1./50.)
fy = np.fft.fft(integral)
fx = np.fft.fftshift(fx)
fy = np.fft.fftshift(fy)
ax[1].plot(fx, np.abs(fy))
#ax[1].set_ylim(-2500, 600000)

for a in ax:
    a.grid()
plt.show()
