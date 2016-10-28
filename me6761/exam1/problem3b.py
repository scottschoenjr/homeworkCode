# -*- coding: utf-8 -*-
"""
Problem 3b
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz

# Define bandpass filter (Warren Weckesser on Stack Exchange) based on Python
# Cookbook example
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
#
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
#

# Define signal paraneters
Fs = 80E3; # Sampling frequency [Hz]
f0 = 2E3; # Pulse cycle frequency
N = 1; # Number of cycles
delay = 50E-3;
cpl = 500; # [m/s]
r = 2; # [m]
filterOrder = 8;

# Constants and initializations
omega0 = 2*np.pi*f0;
dt = 1./Fs;
endTime = 1;
tVector = np.arange(0, endTime, dt);
signal = 0*tVector;

# Define pulse signals
T = 1/f0;
samplesPerPeriod = int(T/dt);
delaySamples = int(delay/dt);
tVector_temp = np.arange(0, N*T, dt);
pulse = np.sin( omega0*tVector_temp );
ind0 = delaySamples;
ind1 = delaySamples + N*samplesPerPeriod;
signal[ind0:ind1] = pulse;

# Filter signal
fLow = 1E3;
fHigh = 10E3;
filteredSignal = \
    butter_bandpass_filter( signal, fLow, fHigh, Fs, order=filterOrder);
receivedSignal = (1/np.sqrt(r))*filteredSignal;
shiftedTime = tVector + r/cpl;


# Plot Signal
plt.plot( 1000*tVector, signal, '--k', label='Transmitted' );
plt.plot( 1000*shiftedTime, receivedSignal, 'k', label='Received' );

plt.xlim([49,58]);
plt.xlabel('Time [ms]', fontsize=14, family='serif');
plt.xticks( fontsize=12, family='serif' );

plt.ylabel('Normalized Amplitude',fontsize=14, family='serif');
plt.yticks( fontsize=12, family='serif' );

plt.ylim([-1, 1]);
plt.legend(loc='best', frameon=False, fontsize=12, \
    prop={'family':'serif'});

# Plot Filter
plt.figure()
b, a = butter_bandpass(fLow, fHigh, Fs, order=filterOrder)
omegaNorm, h = freqz(b, a, worN=2000)
fVector = Fs*omegaNorm/(2*np.pi)
plt.plot(fVector/1E3, 20*np.log10(abs(h)), 'k', label="order = %d" % filterOrder)

plt.xlabel('Frequency [Hz]');
plt.ylabel('Gain [dB]');

plt.ylim([-60, 0.01]);
plt.xlim([0, 12]);


