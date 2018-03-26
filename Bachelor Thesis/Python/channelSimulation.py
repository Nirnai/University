# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 01:52:47 2015

@author: Nirnai
"""
from __future__ import division



import numpy as np


def delayChannel(delay, SNR, data, fs, fc):
    
    r"""
    Delays the Data by Phasemodulation of a complex exponential function, 
    wich coresponds to a Windowed time delayed sinc.
    
    Parameters
    ----------
    delay : delay in s e.g. 1e-6s
    SNR : Signal to Noise ratio Es/N0 in dB
    data : Data to be delayed (need Padding before calling this funktion)
    fs : Sample frequency in Hz
    fc : Carrier frequency in Hz
    """
    
    # Calculates the std from the given SNR
    sigma = np.sqrt(np.sum(data**2)/(10**(SNR/10)*data.size))    
    
    # Creates a Phase shift out of a delay t wich has double the length of the data
    H = np.exp(-2 * np.pi * 1j *delay*(fc + fs * np.r_[:1 + len(data), 1 - len(data):0] / 2 / len(data)))
    
    #Convolves the Signal with the FIR of the channel and cuts of half the length
    s_delay = np.fft.ifft(H * np.fft.fft(data,2*len(data)))[:len(data)]
    
    # Creates a random Noise Signal with given Sigma
    noise = (np.random.normal(0,sigma,s_delay.size) + 1j*np.random.normal(0,sigma,s_delay.size))/np.sqrt(2)
    
    # Recieved Signal with white noise
    r = s_delay + noise
    
    return r 




def twowayChannel(t1,t2,alpha,phi,SNR,data,fs):
    
    r"""
    Simulates a two way channel with the secand path scaled 
    with the factor alpha and randomly turned by phi. Noise is calculated from the direct path and addid to both pathes
    so that the indirect path has have the SNR of the direkt path
    
    Parameters
    ----------
    t1: Line of sight
    t2: indirect path
    alpha: Factor to scale indirect path between 0 and 1
    phi : Randam phaseturn by a reflector
    SNR : Signal to noise Ratio
    data: signal
    fs: Sampling frequency
    """
    sigma = np.sqrt(np.sum(data**2)/(10**(SNR/10)*data.size))
    noise = (np.random.normal(0,sigma,data.size) + 1j*np.random.normal(0,sigma,data.size))/np.sqrt(2)
    randomPhase = np.exp(-1j*phi)    
    
    r = delayChannel(t1,1000,data,fs,868e6) + randomPhase * delayChannel(t2,1000,alpha*data,fs,868e6)+noise
    
    return r    
    







def SNR(t):
    
    r"""
    calculates a fixed value for the SNR from the delay t.
    
    """
    r = t*3e8 + 1
    #r = np.linspace(0,150,151)
     
    # SNR over the distance r 
    # Range from 0m to 150m     
    # Noise Figure from the reciever in dB
    NF = 8
    # Transmitted Power of -20dB
    PtxdB = -20  
    # Recieved Power with pathloss and vegation and attenuation
    PrxdB = PtxdB - 20 * np.log10((4*np.pi*r*868e6)/3e8) - 0.25*r  
    # Noisepower in dB
    NdB = 10*np.log10(1.38e-23 * 290 * 4e6)
    
    SNR = 10**((PrxdB - NF - NdB)/10)  
        
    return SNR


def twowayChannelfor3D(t1,t2,alpha,phi,data,fs):
    
    r"""
    Simulates a two way channel with the secand path scaled 
    with the factor alpha and randomly turned by phi. The SNR is calculated for both Paths depending on the delay of each.
    
    Parameters
    ----------
    t1: Line of sight
    t2: indirect path
    alpha: Factor to scale indirect path between 0 and 1
    phi : Randam phaseturn by a reflector
    data: signal
    fs: Sampling frequency
    """
    N0 = 1 / data.size
    noise = (np.random.normal(0,N0,data.size) + 1j*np.random.normal(0,N0,data.size))/np.sqrt(2)
    randomPhase = np.exp(-1j*phi)    
    SNR1 = np.sqrt(SNR(t1))
    SNR2 = np.sqrt(SNR(t2))
    
    r = SNR1*delayChannel(t1,SNR,data,fs,868e6) + SNR2*randomPhase * delayChannel(t2,SNR,data,fs,868e6) + noise
    
    return r       







