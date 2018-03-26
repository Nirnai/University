# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:24:56 2015

@author: Nirnai Rao
"""
from __future__ import division
import mseq
import numpy as np
from scipy.linalg import hadamard



def generate_Signal(typ,l,fs,f):  
    
    r"""
    Generates a Signal with multiple subcarriers, depending on the parameter typ. The Signalenergy= 1, 
    The Signal is Padded to five times its length l
    
    Parameters
    ----------
    
    typ: can be 4 for a Hadamard sequence with 4 subcarriers, 8 for
         a Hadamard sequence with 8 subcarriers or 0 for a m-sequence
    l : defines the length of the Burstsignal
        (m-sequneces can only have length of power 2 and due to the implementation
        have a max. length of 1024)
    fs : defines the Samplingfrequency of the Signal
    f : defines the frequency of the Signal
    """

    samplerate = int(fs/f)
        
           
    if(typ == 2):
        s = hadamard_sequence(1,8,l,samplerate)           
       
    if(typ == 4):
        
        s = hadamard_sequence(5,16,l,samplerate)   
        
    if(typ == 8):
        
        s = hadamard_sequence(9,32,l,samplerate)
        
    if(typ == 16):
        
        s = hadamard_sequence(25,32,l,samplerate)
        
    if(typ == 1023): 
        
        
        s = mseq.m_sequence(int(np.log2(l)))
        for i in range(1,samplerate):    
            s = np.insert(s, slice(i,None,i),0, axis=0)
            s = np.lib.pad(s, (0,1), 'constant', constant_values=(0, 0))
        rect = ([1]*samplerate)
        s =  np.lib.pad(np.convolve(s, rect, mode='same'),(2*samplerate*1023,2*samplerate*1023),'constant',constant_values = (0,0))
            
        s = s/np.sqrt(np.sum(s**2))
       #s = s/np.sqrt(np.sum(s**2))
        
    

    if(typ == 7):
        
        s = np.tile(mseq.m_sequence(3), np.ceil(l/7))
        for i in range(1,samplerate):    
            s = np.insert(s, slice(i,None,i),0, axis=0)
            s = np.lib.pad(s, (0,1), 'constant', constant_values=(0, 0))
        rect = ([1]*samplerate)
        
        s =  np.lib.pad(np.convolve(s, rect, mode='same'),(2*samplerate*1022,2*samplerate*1022),'constant',constant_values = (0,0))
        s = s/np.sqrt(np.sum(s**2))         
        
    return s





def hadamard_sequence(n, N, l, fs):
    
    r"""
    Takes a Sequence from the generate_sequence funktion and 
    Modulates it on a Rect-Signal
    
    Parameters
    ----------
    
    n : Sequence number
    N : Shape of Hadarmard Matrix
    l : length of the Sequence
    fs : Samplerate
    """
    
        # Generates an NxN Hadamard Matrix
    H = hadamard(N)
    
    # Number of Repetitions needed to Generate M sampled sequences
    r = np.ceil(l/N)
    
    # Repeats Hadarmard Sequence r times 
    s = np.tile(H, (1, r))
    
    # Zero Stuffing to adjust the Symbolrate to the Samplingrate
    for i in range(1,fs):    
        s = np.insert(s, slice(i,None,i),0, axis=1)
        s = np.lib.pad(s, (0,1), 'constant', constant_values=(0, 0))
    
    # Rect-signal sampled with f values
    rect = ([1]*fs)
    
    # Modulate n. sequence with the rect.signal and normalize the result to a signalenergy 1
    # Padding needed, so that signal can't leave window when delayed in channel
    s_mod =  np.lib.pad(np.convolve(s[n], rect, mode='full')[0 : -(fs-1)],(2*fs*l,2*fs*l),'constant',constant_values = (0,0))
    
    # Normalized Signal with energy 1
    return s_mod/np.sqrt(np.sum(s_mod**2))
        
   
   
   



def beta_rms(s,fs):
    
    S = np.abs(np.fft.fftshift(np.fft.fft(s)[::5]))
    f = np.fft.fftshift(np.fft.fftfreq(S.size,1/fs))
    lower_bandlimit = np.where(f>-1.001e6)[0][0]
    upper_bandlimit = np.where(f>0.999e6)[0][0]
    f = f[lower_bandlimit:upper_bandlimit]
    S = S[lower_bandlimit:upper_bandlimit]
    E = np.sum(S**2)
    
    b_rms = np.sqrt(np.sum((f**2*S**2))/E)
    
    return b_rms




def b_rms(s, fs):
    
    # RMS Bandwidth
    # Spectrum and Frequency axes
    spec = np.abs(np.fft.fft(s)[::5]) 
    f = np.fft.fftfreq(spec.size,1/fs)
    
    #plt.plot(f,spec)
    # Signal Energy
    E = np.sum(spec**2)
    # Calculates b_rms with given formula
    b_rms = np.sqrt(np.sum((f**2*spec**2))/E)
    
    
    return b_rms #* fs / spec.size    

  






  
    

    
