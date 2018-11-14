# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:12:11 2016

@author: Nirnai
"""
from __future__ import division
import numpy as np


def LundR(s, r, fs, typ):
    
    r"""
    Luise and Reggianni Spektrum estimator for delay estimation
    
    s : generated Signal (needed for Phase correction)
    r : Recieved Signal
    typ: number of subcarriers
    fs : Samplefrequency
    """
    
    S = np.fft.fftshift(np.fft.fft(s)[::5])
    R = np.fft.fftshift(np.fft.fft(r)[::5])
    
    # Frequencyaxis is normalized to 1MHz
    freq = np.fft.fftshift(np.fft.fftfreq(S.size,1/(fs/1e6)))
        
    # Zurückdrehen der Phase beim erzeugen des Signals 
    lower_bandlimit = np.where(freq>-1.001)[0][0]
    upper_bandlimit = np.where(freq>0.999)[0][0]
    spec = (np.conj(S) * R)[lower_bandlimit:upper_bandlimit]
    
    # Filtern der Frequenzbins
    freqbins = spec[np.where(np.abs(spec) > 1e-5)]
    
    # Damit APpromixmation gilt muss der Betrag der Zeiger auf die gleiche höhe
    freqbins = freqbins/(np.abs(freqbins))
    df = fs/(2*typ) 
    
        
    # Schätzer nach Luise & Reggiannini (Siehe Skript Prof. Koch Empfängersynchronisation)    
    L0 = freqbins.size
    R = np.zeros(L0, dtype=complex)
    
    # L0 -1 in arange um andere abstände zu nehmen
    for k in np.arange(L0-1)+1:    
        for i in np.arange(L0-k):
            R[k-1] += (np.conj(freqbins[i+k])*(freqbins[i]))/(L0-k)
        
    # 2 ist L0 um andere abstände zu nehmen
    t = np.angle(np.sum(R))/(np.pi*(L0)*df)
    
    return t
    
    
    
    
    
    
    
def gemittelte_Phasendifferenz(s,r,fs,typ):
    
    r"""
    avaraging over Subcarriers to estimate Delay
    
    s : generated Signal (needed for Phase correction)
    r : Recieved Signal
    typ: number of subcarriers
    fs : Samplefrequency
    """
        
    
    # Spektren 
    S = np.fft.fftshift(np.fft.fft(s)[::5])
    R = np.fft.fftshift(np.fft.fft(r)[::5])
    freq = np.fft.fftshift(np.fft.fftfreq(S.size,0.25))
    
    lower_bandlimit = np.where(freq>-1.001)[0][0]
    upper_bandlimit = np.where(freq>0.999)[0][0]
    
     # Zurückdrehen der Phase beim erzeugen des Signals 
    spec = (np.conj(S) * R)[lower_bandlimit:upper_bandlimit]
    
    
    # Filtern der Frequenzbins
    freqbins = spec[np.where(np.abs(spec) > 1e-5)]
    # abstand der Subträger    
    df = fs/(2*typ)
        
    # gemittelte Phasendifferenzschätzer     
    N = freqbins.size
    dPhi = 0  
    
    for i in np.arange(N-1):
        dPhi += np.angle(np.conj(freqbins[i+1])*(freqbins[i]))
        
    Phasedifference = dPhi/(N-1)
    
    t = Phasedifference/(2*np.pi*df)
    
    return t
    