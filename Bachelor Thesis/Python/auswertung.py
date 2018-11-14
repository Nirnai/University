# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 15:19:18 2016

@author: Nirnai
"""

from __future__ import division
from mseq import m_sequence
from fractions import Fraction
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes, mark_inset
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import SignalGeneration as s
import channelSimulation as ch
import PhaseEstimation as est
import SubspaceEstimation as sub
from algos import totlstsq
from scipy.stats import norm
from scipy.linalg import hadamard

def AWGN_KanalPlot(SNR):
    #Signal mit Signalenergie 1
    s = np.ones(1000)
    s = s/np.sqrt(np.sum(s**2))
    
    # Erzeuge Weißes Rauschen
    sigma = np.sqrt(np.sum(s**2)/(10**(SNR/10)*s.size))
    noise = np.random.normal(0,sigma,s.size)
    N= np.fft.fftshift(np.fft.fft(noise))
    
    # Delay Kanal
         
    
    freq = np.fft.fftshift(np.fft.fftfreq(s.size, 1/4))
    
    
    plt.figure(1)
    
    plt.subplot(121)
    plt.plot(freq,np.abs(H+N))
    plt.ylim(-5,5) 
    plt.ylabel('|H(f)|')
    plt.xlabel('normierte Frequenz f')  
    
    plt.title('Betragsgang' )
    plt.subplot(122)
    plt.plot(freq,np.angle(H+N))
    plt.ylabel('arg(H(f))')
    plt.xlabel('normierte Frequenz f')
    plt.title('Phasenverlauf' )
    
    
    
    
def multipath_KanalPlot():
    
    #Signal mit Signalenergie 1
    s = np.ones(1000)
    s = s/np.sqrt(np.sum(s**2))
    
    H1 = np.exp(2 * np.pi * 1j * 0.1e-6 *(1e9 + 2e6 * np.r_[1 - len(s):0,:1 + len(s)] / 2 / len(s))) + 0.5*np.exp(2 * np.pi * 1j * 1e-6 *(1e9 + 2e6 * np.r_[1 - len(s):0,:1 + len(s)] / 2 / len(s)))
    H2 = np.exp(2 * np.pi * 1j * 0.1e-6 *(1e9 + 2e6 * np.r_[1 - len(s):0,:1 + len(s)] / 2 / len(s))) + np.exp(2 * np.pi * 1j * 1e-6 *(1e9 + 2e6 * np.r_[1 - len(s):0,:1 + len(s)] / 2 / len(s)))
    H3 = np.exp(2 * np.pi * 1j * 0.1e-6 *(1e9 + 2e6 * np.r_[1 - len(s):0,:1 + len(s)] / 2 / len(s)))
    freq = np.fft.fftshift(np.fft.fftfreq(H1.size,1/2))
    
    plt.figure(2)
    plt.subplot(121)
    plt.plot(freq, np.angle(H3),'k--',label = 'Line of Sight')
    plt.plot(freq,np.angle(H1), label = 'Multipath')
    plt.ylim(-2,2)
    plt.ylabel('arg(H(f))')
    plt.xlabel('normierte Frequenz f')
    plt.title(r'Zweipfad Phasenverlauf mit $\alpha$ = 0.5' ) 
    plt.legend()    
    
    plt.subplot(122)
    plt.plot(freq, np.angle(H3),'k--', label = 'Line of Sight')
    plt.plot(freq,np.angle(H2), label = 'Multipath')
    plt.ylabel('arg(H(f))')
    plt.xlabel('normierte Frequenz f')
    plt.title(r'Zweipfad Phasenverlauf mir $\alpha=$ 1' )     
    plt.legend(loc = 2)
    
    plt.figure(3)
    plt.plot(freq,np.abs(H1), label = r'$\alpha$ = 0.5')
    plt.plot(freq,np.abs(H2), label = r'$\alpha$ = 1')
    plt.ylabel('|H(f)|')
    plt.xlabel('normierte Frequenz f')
    plt.title('Betragsgang eines Zweipfadkanals' ) 
    plt.legend()
    
def FD_FilterPlot():
    
    x1 = np.linspace(-5,9,1000)
    x2 = np.linspace(-5,9,15)
    x3 = np.linspace(-5,9,141)
    s1 = np.sinc(x1-2)
    s2 = np.sinc(x2-2)
    s3 = np.sinc(x3-2)[::13]

    ax1 = plt.subplot2grid((1,1),(0,0))
    ax1.plot(x1,s1, color = 'k')
    ax1.spines['top'].set_color('none')
    ax1.spines['right'].set_color('none')
    ax1.spines['bottom'].set_position('zero')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.stem(x2,s2, markerfmt = 'go', linefmt = 'g', basefmt = 'k', label = 'D = 2')
    ax1.stem(x3[::13],s3, markerfmt = 'ro', linefmt = 'r', basefmt = 'k', label = 'D = 2,3' )
    ax1.legend()
    plt.ylim(-0.5,1.5)
    
    
    plt.show()
    
def ML_estimatorPlot():
    
    x = np.linspace(0,6,10000)
    x2 = np.linspace(-6,6,10000)
    p1 = norm.pdf(x,3,0.3)
    p2 = norm.pdf(x,3,0.6)
    l1 = np.log(p1)
    l2 = np.log(p2)
    
    plt.figure()
    plt.subplot(121)
    plt.plot(x,p1)
    plt.plot(x,p2)
    plt.xlabel(r'Parameter $\theta$')
    plt.ylabel(r'p(x;$\theta$)')
    plt.title('Likelihood Funktion')
    plt.subplot(122)
    plt.plot(x2,l1)
    plt.plot(x2,l2)
    plt.xlabel(r'Parameter $\theta$')
    plt.ylabel(r'$\ln$(p(x;$\theta$))')
    plt.title('Log Likelihood Funktion')    
    

def hadamard_spectrums(): 
    
    H = hadamard(8)
    
    f = np.fft.fftshift(np.fft.fftfreq(2001,1/2))
    
    s1 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[0],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s2 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[1],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s3 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[2],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s4 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[3],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s5 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[4],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s6 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[5],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s7 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[6],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    s8 = np.abs(np.fft.fftshift(np.fft.fft(np.lib.pad(np.tile(H[7],250),(0,1),'constant',constant_values=(0,0)))))/f.size
    
    plt.figure()

    ax1 = plt.subplot2grid((2,4),(0,0))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax2 = plt.subplot2grid((2,4),(0,1))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax3 = plt.subplot2grid((2,4),(0,2))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax4 = plt.subplot2grid((2,4),(0,3))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax5 = plt.subplot2grid((2,4),(1,0))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax6 = plt.subplot2grid((2,4),(1,1))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax7 = plt.subplot2grid((2,4),(1,2))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    ax8 = plt.subplot2grid((2,4),(1,3))
    plt.ylim(0,1.1)
    plt.xlim(-1.1,1.1)
    
    ax1.plot(f,s1)
    ax2.plot(f,s2)
    ax3.plot(f,s3)
    ax4.plot(f,s4)
    ax5.plot(f,s5)
    ax6.plot(f,s6)
    ax7.plot(f,s7)
    ax8.plot(f,s8)
    
    
    


def m_sequencePlot():
    
    sequence = m_sequence(3)
    #signal = np.tile(sequence,1000)
    
    f1 = np.fft.fftshift(np.fft.fftfreq(sequence.size,1/2))
    spec1 = np.abs(np.fft.fftshift(np.fft.fft(sequence)))    
    
    plt.figure()
    ax1 = plt.subplot2grid((1,1),(0,0))
    ax1.stem(f1,spec1, markerfmt = 'ko', linefmt = 'k', basefmt = 'k')

    
    #plt.stem(f1,spec1, markerfmt = 'ko', linefmt = 'k', basefmt = 'k')
   
   
def Eigenwertzerlegung():
    signal = s.generate_Signal(8,1024,2e6,1e6)
    r = ch.twowayChannel(1.3e-6,2.7e-6,1,0,20,signal,2e6)
    
    R = sub.ChannelCovMatrix(signal,r,8)
    sigma, U = np.linalg.eigh(R)
    
    plt.figure()
    plt.stem(np.sort(sigma)[::-1])
    #plt.ylim(0,np.argmax(sigma)+1)
    plt.xlim(-1,15)
    plt.xlabel('M-Eigenwerte')
    


def Simulationsspektren(): 
    
    f = np.fft.fftshift(np.fft.fftfreq(1024,0.5))
    
    s4 = np.abs(np.fft.fftshift(np.fft.fft(np.tile(hadamard(8)[5], 128))))/f.size
    s8 = np.abs(np.fft.fftshift(np.fft.fft(np.tile(hadamard(32)[9], 32))))/f.size
    s16 = np.abs(np.fft.fftshift(np.fft.fft(np.tile(hadamard(32)[25], 32))))/f.size
    
    
    
    plt.figure(figsize = (10,3))
    #plt.subtitle('normiertes Betragsspektrum')
    
    ax1  = plt.subplot2grid((1,3),(0,0))
    plt.ylim(0,1)
    ax1.set_title('4-Ton')
    plt.xlabel('normierte Frequenz')
    plt.ylabel('|S(f)|')
    #ax1.set_yticks([])
    ax2 = plt.subplot2grid((1,3),(0,1))
    plt.ylim(0,1)
    ax2.set_title('8-Ton')
    plt.xlabel('normierte Frequenz')
    plt.ylabel('|S(f)|')
    #plt.xlabel('normierte Frequenz')
    #ax2.set_yticks([])
    ax3 = plt.subplot2grid((1,3),(0,2))
    plt.ylim(0,1)
    ax3.set_title('16-Ton')
    plt.xlabel('normierte Frequenz')
    plt.ylabel('|S(f)|')
    #ax3.set_yticks([])
   
    #fig.tight_layout()
    #plt.subplots_adjust(top=0.85)
    ax1.plot(f,s4)
    ax2.plot(f,s8)
    ax3.plot(f,s16)
    
    
def mSeqPlot():
    
    f1 = np.fft.fftshift(np.fft.fftfreq(1022,0.5))
    f2 = np.fft.fftshift(np.fft.fftfreq(1023,0.5))
    s1 = np.abs(np.fft.fftshift(np.fft.fft(np.tile(m_sequence(3),146))))/f1.size
    s2 = np.abs(np.fft.fftshift(np.fft.fft(m_sequence(10))))/f2.size
    
    plt.figure(figsize = (10,3))
    ax1 = plt.subplot2grid((1,2),(0,0))
    plt.ylim(0,0.5)
    ax1.set_title('7-Ton')
    ax2 = plt.subplot2grid((1,2),(0,1))
    plt.ylim(0,0.1)
    ax2.set_title('weisses Spektrum')
    
    ax1.plot(f1,s1)
    ax2.plot(f2,s2)
  

def Modspektrum():

    
    f1= np.fft.fftshift(np.fft.fftfreq(2048,0.25))
    signal1 = np.abs(np.fft.fftshift(np.fft.fft(s.generate_Signal(4,1024,4e6,2e6))[::5]))/f1.size
    
    f2 = np.fft.fftshift(np.fft.fftfreq(2046,0.25))    
    signal2 = np.abs(np.fft.fftshift(np.fft.fft(s.generate_Signal(1023,1024,4e6,2e6))[::5]))/f2.size
    
    
    
    plt.figure(figsize=(8,3))
    ax1 = plt.subplot2grid((1,2),(0,0))
    ax1.set_title('Modulierter 4-Ton')
    ax1.set_ylim(0,1)
    ax2 = plt.subplot2grid((1,2),(0,1))
    ax2.set_title('Modulierte m-Sequenz')
    ax2.set_ylim(0,0.1)

    ax1.plot(f1,signal1)
    ax2.plot(f2,signal2)
    
    
def SNR():
    
    r = np.linspace(0,150,151)
     
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
    
    
    
    def CRB_range():
    
    r"""
    Calculates the CRB for the given signal
    
    Parameters:
    -----------
    
    s : Signal to be evaluated
    """
    
    s = generate_Signal(2,1024,4e6,2e6) 
    s1 = generate_Signal(4,1024,4e6,2e6) 
    s2 = generate_Signal(8,1024,4e6,2e6)
    s3 = generate_Signal(16,1024,4e6,2e6)
    s4 = generate_Signal(7,1022,4e6,2e6)
    s5 = generate_Signal(1023,1024,4e6,2e6)
    
    # Parameters needed for CRB
    # Time of signal in seconds
    Ts = 0.5e-3
    # Bandwidth in Hz
    B = 4e6    
    
    r = np.linspace(0,150,151)

    # CRB with all given parameters 
    #SNR = np.arange(101)
    sigma_r = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s, B)**2 * 10**(SNR()/10)))
    sigma_r1 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s1, B)**2 * 10**(SNR()/10)))  
    sigma_r2 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s2, B)**2 * 10**(SNR()/10)))  
    sigma_r3 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s3, B)**2 * 10**(SNR()/10))) 
    sigma_r4 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s4, B)**2 * 10**(SNR()/10))) 
    sigma_r5 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*B*b_rms(s5, B)**2 * 10**(SNR()/10))) 
    
    
    plt.figure(1)
    plt.plot(r,sigma_r, color='k', label = '2-Ton')
    plt.plot(r,sigma_r1, color='b', label = '4-Ton')
    plt.plot(r,sigma_r2, color='g', label = '8-Ton')
    plt.plot(r,sigma_r3, color='r', label = '16-Ton')
    plt.plot(r,sigma_r4, color='c', label = '7-Ton')
    plt.plot(r,sigma_r5, color='m', label = '1023-Ton')
    plt.legend(loc = 2)
    plt.ylabel(r'$\mathbf{\sigma_{range}}$ in m', fontsize = 15)
    plt.xlabel('Distanz r in m', fontsize = 15)
    plt.xlim(0,150)
    plt.show()
    
    
    
    def CRLB_SNR():
    
    s = generate_Signal(2,1024,4e6,2e6) 
    s1 = generate_Signal(4,1024,4e6,2e6) 
    s2 = generate_Signal(8,1024,4e6,2e6)
    s3 = generate_Signal(16,1024,4e6,2e6)
    s4 = generate_Signal(7,1022,4e6,2e6)
    s5 = generate_Signal(1023,1024,4e6,2e6)
    
    SNR = np.linspace(0,100,101)
    # Time of signal in seconds
    Ts = 0.5e-3
    # Bandwidth in Hz
    fs= 4e6 
    #B = 2e6
    
    sigma_r = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s, fs))**2 * 10**(SNR/10)))
    sigma_r1 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s1, fs))**2 * 10**(SNR/10)))  
    sigma_r2 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s2, fs))**2 * 10**(SNR/10)))  
    sigma_r3 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s3, fs))**2 * 10**(SNR/10))) 
    sigma_r4 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s4, fs))**2 * 10**(SNR/10))) 
    sigma_r5 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(b_rms(s5, fs))**2 * 10**(SNR/10))) 
    
    plt.figure(10)
    plt.plot(SNR+10*np.log10(4e6),sigma_r, color='k', label= '2-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_r1, color='b', label = '4-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_r2, color='g', label = '8-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_r3, color='r', label = '16-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_r4, color='c', label = '7-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_r5, color='m', label = '1023-Ton')
    plt.legend()
    plt.ylabel(r'$\mathbf{\sigma_{range}}$ in m', fontsize = 15)
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dBHz', fontsize = 15)
    plt.show()















def varianz_Phasendifferenz():
    
    # Hadamard-sequenzen 
    signal1 = s.generate_Signal(4,1024,4e6,2e6)
    signal2 = s.generate_Signal(8,1024,4e6,2e6)
    signal3 = s.generate_Signal(16,1024,4e6,2e6)
    
    
    
    signal4 = s.generate_Signal(7,1022,4e6,2e6)
    #signal5 = s.generate_Signal(1023,1024,4e6,2e6)
    
    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6
    
    
    t_fehler1 = np.zeros(100)
    t_fehler2 = np.zeros(100)
    t_fehler3 = np.zeros(100)
    t_fehler4 = np.zeros(100)
    #t_fehler5 = np.zeros(1)
    sigma_t1 = np.zeros(SNR.size)
    sigma_t2 = np.zeros(SNR.size)
    sigma_t3 = np.zeros(SNR.size)
    sigma_t4 = np.zeros(SNR.size)
    #sigma_t5 = np.zeros(SNR.size)        
    # Referenz Empfangssignale ohne Rauschen 
    r_ref1 = ch.delayChannel(tau,1000,signal1,4e6,868e6)
    r_ref2 = ch.delayChannel(tau,1000,signal2,4e6,868e6)
    r_ref3 = ch.delayChannel(tau,1000,signal3,4e6,868e6)
    r_ref4 = ch.delayChannel(tau,1000,signal4,4e6,868e6)
    #r_ref5 = ch.delayChannel(tau,1000,signal5,4e6,868e6)
        
    # Referenz Laufzeitschätzungen
    t_ref1 = est.gemittelte_Phasendifferenz(signal1,r_ref1,4e6,4)
    t_ref2 = est.gemittelte_Phasendifferenz(signal2,r_ref2,4e6,8)
    t_ref3 = est.gemittelte_Phasendifferenz(signal3,r_ref3,4e6,16)
    t_ref4 = est.gemittelte_Phasendifferenz(signal4,r_ref4,4e6,7)
    #t_ref5 = est.gemittelte_Phasendifferenz(signal5,r_ref5,4e6,1023)
        
    for j in SNR:
        print j
        np.random.seed(3)
        # Würfelanzahl für Rasuchen 
        for x in np.arange(100):
            # Empfangssignale mit Rauschen
            r1 = ch.delayChannel(tau,j,signal1,4e6,868e6)
            r2 = ch.delayChannel(tau,j,signal2,4e6,868e6)
            r3 = ch.delayChannel(tau,j,signal3,4e6,868e6)
            r4 = ch.delayChannel(tau,j,signal4,4e6,868e6)
            #r5 = ch.delayChannel(tau,j,signal5,4e6,868e6)
                
            # Fehler in der Schätzung
            t_fehler1[x] = (est.gemittelte_Phasendifferenz(signal1,r1,4e6,4) - t_ref1) * c0
            t_fehler2[x] = (est.gemittelte_Phasendifferenz(signal2,r2,4e6,8) - t_ref2) * c0
            t_fehler3[x] = (est.gemittelte_Phasendifferenz(signal3,r3,4e6,16) - t_ref3) * c0
            t_fehler4[x] = (est.gemittelte_Phasendifferenz(signal4,r4,4e6,7) - t_ref4) * c0
            #t_fehler5[x] = (est.gemittelte_Phasendifferenz(signal5,r5,4e6,1023) - t_ref5) * c0
                
        # Varianz des Fehlers
        sigma_t1[j] = np.sqrt(np.mean(t_fehler1**2))
        sigma_t2[j] = np.sqrt(np.mean(t_fehler2**2))
        sigma_t3[j] = np.sqrt(np.mean(t_fehler3**2))
        sigma_t4[j] = np.sqrt(np.mean(t_fehler4**2))
            #sigma_t5[j] = np.sqrt(np.mean(t_fehler5**2))
    
    fig1 = plt.figure(1)
    plt.plot(SNR+10*np.log10(4e6),sigma_t1, label = '4-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_t2, label = '8-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_t3, label = '16-Ton')
    plt.plot(SNR+10*np.log10(4e6),sigma_t4, label = '7-Ton')
    #plt.plot(SNR,sigma_t5, label = '1023-Ton')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.legend()   
    plt.title('Gemittelte Phasendifferenz')
    fig1.show()
    
    return sigma_t1, sigma_t2, sigma_t3, sigma_t4    
    




def varianz_mseq_Phasendifferenz():
    
    signal5 = s.generate_Signal(1023,1024,4e6,2e6)

    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6

    t_fehler = np.zeros(10)
    sigma_5 = np.zeros(SNR.size)
    
    r_ref5 = ch.delayChannel(tau,1000,signal5,4e6,868e6)
    t_ref5 = est.gemittelte_Phasendifferenz(signal5,r_ref5,4e6,1023)
    
    for i in SNR:
        print i
        np.random.seed(0)
        for x in np.arange(10):
            r5 = ch.delayChannel(tau,i,signal5,4e6,868e6)
            t_fehler[x] = (est.gemittelte_Phasendifferenz(signal5,r5,4e6,1023) - t_ref5) * c0
            
        sigma_5[i] = np.sqrt(np.mean(t_fehler**2))
    
    plt.figure(1)
    plt.plot(SNR + 10*np.log10(4e6), sigma_5, color='m', label = '1023-Ton')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.legend()
    plt.title('gemittelte Phasendifferenz')
    plt.show()
    
    return sigma_5
    
    











def varianz_LundR():
    
    # Hadamard-sequenzen 
    signal1 = s.generate_Signal(4,1024,4e6,2e6)
    signal2 = s.generate_Signal(8,1024,4e6,2e6)
    signal3 = s.generate_Signal(16,1024,4e6,2e6)
    
    
    # m-Sequenzen
    signal4 = s.generate_Signal(7,1022,4e6,2e6)

    
    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6
    
    
    t_fehler1 = np.zeros(100)
    t_fehler2 = np.zeros(100)
    t_fehler3 = np.zeros(100)
    t_fehler4 = np.zeros(100)

    sigma_t1 = np.zeros(SNR.size)
    sigma_t2 = np.zeros(SNR.size)
    sigma_t3 = np.zeros(SNR.size)
    sigma_t4 = np.zeros(SNR.size)

    
        
    # Referenz Empfangssignale ohne Rauschen 
    r_ref1 = ch.delayChannel(tau,1000,signal1,4e6,868e6)
    r_ref2 = ch.delayChannel(tau,1000,signal2,4e6,868e6)
    r_ref3 = ch.delayChannel(tau,1000,signal3,4e6,868e6)
    r_ref4 = ch.delayChannel(tau,1000,signal4,4e6,868e6)

        
    # Referenz Laufzeitschätzungen
    t_ref1 = est.LundR(signal1,r_ref1,4e6,4)
    t_ref2 = est.LundR(signal2,r_ref2,4e6,8)
    t_ref3 = est.LundR(signal3,r_ref3,4e6,16)
    t_ref4 = est.LundR(signal4,r_ref4,4e6,7)

        
    for j in SNR:
        print j
        np.random.seed(0)
        # Würfelanzahl für Rasuchen 
        for x in np.arange(100):
            # Empfangssignale mit Rauschen
            r1 = ch.delayChannel(tau,j,signal1,4e6,868e6)
            r2 = ch.delayChannel(tau,j,signal2,4e6,868e6)
            r3 = ch.delayChannel(tau,j,signal3,4e6,868e6)
            r4 = ch.delayChannel(tau,j,signal4,4e6,868e6)

                
            # Fehler in der Schätzung
            t_fehler1[x] = (t_ref1 - est.LundR(signal1,r1,4e6,4) ) * c0
            t_fehler2[x] = (t_ref2 - est.LundR(signal2,r2,4e6,8) ) * c0
            t_fehler3[x] = (t_ref3 - est.LundR(signal3,r3,4e6,16)) * c0
            t_fehler4[x] = (t_ref4 - est.LundR(signal4,r4,4e6,7) ) * c0

            
        # rmse
        sigma_t1[j] = np.sqrt(np.mean(t_fehler1**2))
        sigma_t2[j] = np.sqrt(np.mean(t_fehler2**2))
        sigma_t3[j] = np.sqrt(np.mean(t_fehler3**2))
        sigma_t4[j] = np.sqrt(np.mean(t_fehler4**2))

        
    fig2 = plt.figure(2)
    plt.plot(SNR + 10*np.log10(4e6),sigma_t1, label = '4-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t2, label = '8-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t3, label = '16-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t4, label = '7-Ton')

    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.legend()   
    plt.title('L$&$R')
    fig2.show()
    
    return sigma_t1, sigma_t2, sigma_t3, sigma_t4  
    
    
    
def varianz_mseq_LundR():
    
    signal5 = s.generate_Signal(1023,1024,4e6,2e6)

    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6

    t_fehler = np.zeros(10)
    sigma_5 = np.zeros(SNR.size)
    
    r_ref5 = ch.delayChannel(tau,1000,signal5,4e6,868e6)
    t_ref5 = est.LundR(signal5,r_ref5,4e6,1023)
    
    for i in SNR:
        print i
        np.random.seed(0)
        for x in np.arange(10):
            r5 = ch.delayChannel(tau,i,signal5,4e6,868e6)
            t_fehler[x] = (est.LundR(signal5,r5,4e6,1023) - t_ref5) * c0
            
        sigma_5[i] = np.sqrt(np.mean(t_fehler**2))
    
    plt.figure(2)
    plt.plot(SNR + 10*np.log10(4e6), sigma_5, color='m', label = '1023-Ton')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.title(r'L$&$R')
    plt.show()
    
    return sigma_5
        








    
    

def varianz_ESPRIT_MUSIC():
    # Hadamard-sequenzen 
    signal1 = s.generate_Signal(4,1024,4e6,2e6)
    signal2 = s.generate_Signal(8,1024,4e6,2e6)
    signal3 = s.generate_Signal(16,1024,4e6,2e6)
    
    
    # m-Sequenzen
    signal4 = s.generate_Signal(7,1022,4e6,2e6)

    
    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6
    
    t_fehler1 = np.zeros(100)
    t_fehler2 = np.zeros(100)
    t_fehler3 = np.zeros(100)
    t_fehler4 = np.zeros(100)
    
    t_fehler5 = np.zeros(100)
    t_fehler6 = np.zeros(100)
    t_fehler7 = np.zeros(100)
    t_fehler8 = np.zeros(100)

    sigma_t1 = np.zeros(SNR.size)
    sigma_t2 = np.zeros(SNR.size)
    sigma_t3 = np.zeros(SNR.size)
    sigma_t4 = np.zeros(SNR.size)
    
    sigma_t5 = np.zeros(SNR.size)
    sigma_t6 = np.zeros(SNR.size)
    sigma_t7 = np.zeros(SNR.size)
    sigma_t8 = np.zeros(SNR.size)

    # Referenz Empfangssignale ohne Rauschen 
    r_ref1 = ch.delayChannel(tau,1000,signal1,4e6,868e6)
    r_ref2 = ch.delayChannel(tau,1000,signal2,4e6,868e6)
    r_ref3 = ch.delayChannel(tau,1000,signal3,4e6,868e6)
    r_ref4 = ch.delayChannel(tau,1000,signal4,4e6,868e6)
    #r_ref5 = ch.delayChannel(tau,1000,signal5,4e6,868e6)
    
    # Referenz Laufzeitschätzungen
    t_ref1 = sub.esprit(signal1,r_ref1,4e6,4)
    t_ref2 = sub.esprit(signal2,r_ref2,4e6,8)
    t_ref3 = sub.esprit(signal3,r_ref3,4e6,16)
    t_ref4 = sub.esprit(signal4,r_ref4,4e6,7)
    
    t_ref5 = sub.rtmusic(signal1,r_ref1,4e6,4)
    t_ref6 = sub.rtmusic(signal2,r_ref2,4e6,8)
    t_ref7 = sub.rtmusic(signal3,r_ref3,4e6,16)
    t_ref8 = sub.rtmusic(signal4,r_ref4,4e6,7)

    
    for j in SNR:
        print j
        np.random.seed(5)
        # Würfelanzahl für Rasuchen 
        for x in np.arange(100):
            # Empfangssignale mit Rauschen
            r1 = ch.delayChannel(tau,j,signal1,4e6,868e6)
            r2 = ch.delayChannel(tau,j,signal2,4e6,868e6)
            r3 = ch.delayChannel(tau,j,signal3,4e6,868e6)
            r4 = ch.delayChannel(tau,j,signal4,4e6,868e6)
                
            # Fehler in der Schätzung
            t_fehler1[x] = (sub.esprit(signal1,r1,4e6,4) - t_ref1) * c0
            t_fehler2[x] = (sub.esprit(signal2,r2,4e6,8) - t_ref2) * c0
            t_fehler3[x] = (sub.esprit(signal3,r3,4e6,16) - t_ref3) * c0
            t_fehler4[x] = (sub.esprit(signal4,r4,4e6,7) - t_ref4) * c0
            
            t_fehler5[x] = (sub.rtmusic(signal1,r1,4e6,4) - t_ref1) * c0
            t_fehler6[x] = (sub.rtmusic(signal2,r2,4e6,8) - t_ref2) * c0
            t_fehler7[x] = (sub.rtmusic(signal3,r3,4e6,16) - t_ref3) * c0
            t_fehler8[x] = (sub.rtmusic(signal4,r4,4e6,7) - t_ref4) * c0
            
            
            
        # Varianz des Fehlers
        sigma_t1[j] = np.sqrt(np.mean(t_fehler1**2))
        sigma_t2[j] = np.sqrt(np.mean(t_fehler2**2))
        sigma_t3[j] = np.sqrt(np.mean(t_fehler3**2))
        sigma_t4[j] = np.sqrt(np.mean(t_fehler4**2))
        
        sigma_t5[j] = np.sqrt(np.mean(t_fehler5**2))
        sigma_t6[j] = np.sqrt(np.mean(t_fehler6**2))
        sigma_t7[j] = np.sqrt(np.mean(t_fehler7**2))
        sigma_t8[j] = np.sqrt(np.mean(t_fehler8**2))

    fig1 = plt.figure(3)
    plt.plot(SNR + 10*np.log10(4e6),sigma_t1, label = '4-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t2, label = '8-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t3, label = '16-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t4, label = '7-Ton')

    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.legend()   
    plt.title('ESPRIT')
    fig1.show()
    
    fig2 = plt.figure(4)
    plt.plot(SNR + 10*np.log10(4e6),sigma_t5, label = '4-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t6, label = '8-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t7, label = '16-Ton')
    plt.plot(SNR + 10*np.log10(4e6),sigma_t8, label = '7-Ton')

    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.legend()   
    plt.title('MUSIC')
    fig2.show()
    
    
    return  sigma_t5, sigma_t6, sigma_t7, sigma_t8, sigma_t1, sigma_t2, sigma_t3, sigma_t4,
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def varianz_mseq_ESPRIT():

    signal5 = s.generate_Signal(1023,1024,4e6,2e6)

    c0 = 3e8
    SNR = np.linspace(0,100,101)    
    tau = 0.1e-6

    t_fehler1 = np.zeros(10)
    #t_fehler2 = np.zeros(100)
    sigma_5 = np.zeros(SNR.size)
    #sigma_6 = np.zeros(SNR.size)
    
    r_ref5 = ch.delayChannel(tau,1000,signal5,4e6,868e6)
    t_ref5 = sub.esprit(signal5,r_ref5,4e6,1023)
    #t_ref6 = sub.rtmusic(signal5,r_ref5,4e6,1023)
    
    for i in SNR:
        print i
        np.random.seed(5)
        for x in np.arange(10):
            r5 = ch.delayChannel(tau,i,signal5,4e6,868e6)
            t_fehler1[x] = (sub.esprit(signal5,r5,4e6,1023) - t_ref5) * c0
            #t_fehler2[x] = (sub.rtmusic(signal5,r5,4e6,1023)- t_ref6) * c0
        sigma_5[i] = np.sqrt(np.mean(t_fehler1**2))
        #sigma_6[i] = np.sqrt(np.mean(t_fehler2**2))
    
    plt.figure(3)
    plt.plot(SNR + 10*np.log10(4e6), sigma_5, color='m', label = '1023-Ton')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.title('ESPRIT')
    plt.show()
    
#    plt.figure(4)
#    plt.plot(SNR + 10*np.log10(4e6), sigma_6, color='m', label = '1023-Ton')
#    plt.ylabel('rmse in m')
#    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
#    plt.title('MUSIC')
#    plt.show()
    
    return sigma_5
    
    
    
    
    
    
    
    
    
    
    
    

def Plotvergleich(sigma1,sigma2,sigma3,sigma4, sigma5, label='estimatortype'):
    
    s1 = s.generate_Signal(4,1024,4e6,2e6) 
    s2 = s.generate_Signal(8,1024,4e6,2e6)
    s3 = s.generate_Signal(16,1024,4e6,2e6)
    s4 = s.generate_Signal(7,1022,4e6,2e6)
    s5 = s.generate_Signal(1023,1024,4e6,2e6)
    
    SNR = np.linspace(0,100,101) 
    
    Ts = 0.5e-3
    # Bandwidth in Hz
    fs= 4e6 
    #B = 2e6
    
    
    CRLB1 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s1, fs))**2 * 10**(SNR/10)))  
    CRLB2 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s2, fs))**2 * 10**(SNR/10)))  
    CRLB3 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s3, fs))**2 * 10**(SNR/10))) 
    CRLB4 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s4, fs))**2 * 10**(SNR/10)))  
    CRLB5 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s5, fs))**2 * 10**(SNR/10)))
    
    
    fig1 = plt.figure(5)
    plt.plot(SNR + 10*np.log10(4e6),sigma1, color = 'b', label= label)
    plt.plot(SNR + 10*np.log10(4e6),CRLB1, color='k', label = 'CRLB')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()       
    
    
    fig2 = plt.figure(6)
    plt.plot(SNR + 10*np.log10(4e6),sigma2, color = 'g' ,label=label)
    plt.plot(SNR + 10*np.log10(4e6),CRLB2, color='k', label = 'CRLB')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()   
    
    fig3 = plt.figure(7)
    plt.plot(SNR + 10*np.log10(4e6),sigma3, color='r' ,label=label)
    plt.plot(SNR + 10*np.log10(4e6),CRLB3, color='k', label = 'CRLB')
    plt.ylabel('rmse im m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()
    
    fig4 = plt.figure(8)
    plt.plot(SNR + 10*np.log10(4e6),sigma4*1.2, color= 'c', label=label)
    plt.plot(SNR + 10*np.log10(4e6),CRLB4, color='k', label = 'CRLB')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()   
    
    fig5 = plt.figure(9)
    plt.plot(SNR + 10*np.log10(4e6),sigma5, color= 'm', label=label)
    plt.plot(SNR + 10*np.log10(4e6),CRLB5, color='k', label = 'CRLB')
    plt.ylabel('rmse in m')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()  
    
    
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    fig5.show()

    
def estim_efficiency(sigma1,sigma2,sigma3,sigma4,sigma5):
    
    s1 = s.generate_Signal(4,1024,4e6,2e6) 
    s2 = s.generate_Signal(8,1024,4e6,2e6)
    s3 = s.generate_Signal(16,1024,4e6,2e6)
    s4 = s.generate_Signal(7,1022,4e6,2e6)
    s5 = s.generate_Signal(1023,1024,4e6,2e6)
    
    SNR = np.linspace(0,100,101) 
    
    Ts = 0.5e-3
    # Bandwidth in Hz
    fs= 4e6 
    #B = 2e6
    
    
    CRLB1 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s1, fs))**2 * 10**(SNR/10)))  
    CRLB2 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s2, fs))**2 * 10**(SNR/10)))  
    CRLB3 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s3, fs))**2 * 10**(SNR/10))) 
    CRLB4 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s4, fs))**2 * 10**(SNR/10)))  
    CRLB5 = np.sqrt((3e8**2)/((2*np.pi)**2*Ts*fs*(s.b_rms(s5, fs))**2 * 10**(SNR/10)))
    
    fig1 = plt.figure(5)
    plt.plot(SNR + 10*np.log10(4e6),CRLB1/sigma1, color = 'b', label= '4-Ton')
    plt.plot(SNR + 10*np.log10(4e6),CRLB2/sigma2, color = 'g', label= '8-Ton')
    plt.plot(SNR + 10*np.log10(4e6),CRLB3/sigma3, color = 'r', label= '16-Ton')
    plt.plot(SNR + 10*np.log10(4e6),CRLB4/(sigma4), color = 'c', label= '7-Ton')
    #plt.plot(SNR + 10*np.log10(4e6),CRLB5/sigma5, color = 'm', label= '1023-Ton')
    plt.ylabel(u'Schätzeffizienz')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylim(0,1.2)
    plt.xlim(66,140)
    plt.legend(loc=4)  
    fig1.show()
    
    
def Fehlerkurve_gemitteltePhasendifferenz(alpha):
    
    s1 = s.generate_Signal(4,1024,4e6,2e6) 
    s2 = s.generate_Signal(8,1024,4e6,2e6)
    s3 = s.generate_Signal(16,1024,4e6,2e6)
    s4 = s.generate_Signal(7,1022,4e6,2e6)
    
    
    randomPhase = np.linspace(0,np.pi,4)
    dt = np.linspace(0,1e-6,1001)
    
    t1_fehler = np.zeros(dt.size)
    t2_fehler = np.zeros(dt.size)
    t3_fehler = np.zeros(dt.size)
    t4_fehler = np.zeros(dt.size)
    
    t1 = np.zeros(1)
    t2 = np.zeros(1)
    t3 = np.zeros(1)
    t4 = np.zeros(1)
    
    
    
    for i,phi in enumerate(randomPhase):
        for j in np.arange(1001):
            for x in np.arange(1):
                
                r1 = ch.twowayChannel(0,dt[j],alpha,phi,1000,s1,4e6)
                r2 = ch.twowayChannel(0,dt[j],alpha,phi,1000,s2,4e6)
                r3 = ch.twowayChannel(0,dt[j],alpha,phi,1000,s3,4e6)
                r4 = ch.twowayChannel(0,dt[j],alpha,phi,1000,s4,4e6)
                
                t1[x] = est.gemittelte_Phasendifferenz(s1,r1,4e6,4)
                t2[x] = est.gemittelte_Phasendifferenz(s2,r2,4e6,8)
                t3[x] = est.gemittelte_Phasendifferenz(s3,r3,4e6,16)
                t4[x] = est.gemittelte_Phasendifferenz(s4,r4,4e6,7)
                
            t1_fehler[j] = np.mean(t1)
            t2_fehler[j] = np.mean(t2)
            t3_fehler[j] = np.mean(t3)
            t4_fehler[j] = np.mean(t4)
            
        dx = dt*3e8
        y1 = t1_fehler *3e8
        y2 = t2_fehler *3e8
        y3 = t3_fehler *3e8
        y4 = t4_fehler *3e8
        
        if(phi == np.pi):
            fig0 = plt.figure(9)
            plt.plot(dx,y1,label = '4-Ton')
            plt.plot(dx,y2,label = '8-Ton')
            plt.plot(dx,y3,label = '16-Ton')
            plt.plot(dx,y4,label = '7-Ton')
            plt.ylabel('rmse in m')
            plt.xlabel('Umweg in m')
            plt.ylim(-40,40)
            plt.legend()
            plt.title(u'Hüllkurve')
            fig0.show()
            
            
        
        fig1 = plt.figure(10)
        plt.plot(dx,y1, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 4-Ton')
        
        fig2 = plt.figure(11)
        plt.plot(dx,y2, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('Fehler in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 8-Ton')
        
        fig3 = plt.figure(12)
        plt.plot(dx,y3, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 16-Ton')
        
        fig4 = plt.figure(13)
        plt.plot(dx,y4, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 7-Ton')
        
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    
    return y1, y2, y3, y4
    
    
    
    
    
    
    
def Fehlerkurve_mseq_Phasendifferenz(alpha):
    
    s5 = s.generate_Signal(1023,1024,4e6,2e6)
    
    
    randomPhase = np.linspace(0,np.pi,4)
    dt = np.linspace(0,1e-6,1001)
    
    t5_fehler = np.zeros(dt.size)

    t5 = np.zeros(1)
    
    
    
    for i,phi in enumerate(randomPhase):
        for j in np.arange(1001):
            for x in np.arange(1):
                
                r5 = ch.twowayChannel(0,dt[j],alpha,phi,100,s5,4e6)
                
                
                t5[x] = est.gemittelte_Phasendifferenz(s5,r5,4e6,1023)

                
            t5_fehler[j] = np.mean(t5)

        dx = dt*3e8

        y5 = t5_fehler *3e8
        if(phi == 0):
            fig0 = plt.figure(9)
            plt.plot(dx,y5,label = '1023-Ton')
            plt.ylabel('rmse in m')
            plt.xlabel('Umweg in m')
            plt.ylim(-40,40)
            plt.legend()
            plt.title(u'Hüllkurve')
            fig0.show()
        
        fig = plt.figure(60)
        plt.plot(dx,y5, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 1023-Ton')
    fig.show()

    return y5  











      
 
def Fehlerkurve_LundR():
    
    s1 = s.generate_Signal(4,1024,4e6,2e6) 
    s2 = s.generate_Signal(8,1024,4e6,2e6)
    s3 = s.generate_Signal(16,1024,4e6,2e6)
    s4 = s.generate_Signal(7,1022,4e6,2e6)
    
    
    randomPhase = np.linspace(0,np.pi,4)
    dt = np.linspace(0,1e-6,1001)
    
    t1_fehler = np.zeros(dt.size)
    t2_fehler = np.zeros(dt.size)
    t3_fehler = np.zeros(dt.size)
    t4_fehler = np.zeros(dt.size)
    
    t1 = np.zeros(1)
    t2 = np.zeros(1)
    t3 = np.zeros(1)
    t4 = np.zeros(1)
    
    
    
    for i in np.arange(4):
        for j in np.arange(1001):
            for x in np.arange(1):
                
                r1 = ch.twowayChannel(0,dt[j],0.5,randomPhase[i],100,s1,4e6)
                r2 = ch.twowayChannel(0,dt[j],0.5,randomPhase[i],100,s2,4e6)
                r3 = ch.twowayChannel(0,dt[j],0.5,randomPhase[i],100,s3,4e6)
                r4 = ch.twowayChannel(0,dt[j],0.5,randomPhase[i],100,s4,4e6)
                
                t1[x] = est.LundR(s1,r1,4e6,4)
                t2[x] = est.LundR(s2,r2,4e6,8)
                t3[x] = est.LundR(s3,r3,4e6,16)
                t4[x] = est.LundR(s4,r4,4e6,7)
                
            t1_fehler[j] = np.mean(t1)
            t2_fehler[j] = np.mean(t2)
            t3_fehler[j] = np.mean(t3)
            t4_fehler[j] = np.mean(t4)
            
        dx = dt*3e8
        y1 = t1_fehler *3e8
        y2 = t2_fehler *3e8
        y3 = t3_fehler *3e8
        y4 = t4_fehler *3e8
        
        if(randomPhase[i] == 0):
            fig0 = plt.figure(19)
            plt.plot(dx,y1,label = '4-Ton')
            plt.plot(dx,y2,label = '8-Ton')
            plt.plot(dx,y3,label = '16-Ton')
            plt.plot(dx,y4,label = '7-Ton')
            plt.ylabel('rmse in m')
            plt.xlabel('Umweg in m')
            plt.ylim(-40,40)
            plt.legend()
            plt.title(u'Hüllkurve')
            fig0.show()
            
            
        
        fig1 = plt.figure(20)
        plt.plot(dx,y1, label = '$\phi= %i^{\circ}$' %(np.round(randomPhase[i]/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 4-Ton')
        
        fig2 = plt.figure(21)
        plt.plot(dx,y2, label = '$\phi= %i^{\circ}$' %(np.round(randomPhase[i]/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 8-Ton')
        
        fig3 = plt.figure(22)
        plt.plot(dx,y3, label = '$\phi= %i^{\circ}$' %(np.round(randomPhase[i]/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 16-Ton')
        
        fig4 = plt.figure(23)
        plt.plot(dx,y4, label = '$\phi= %i^{\circ}$' %(np.round(randomPhase[i]/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 7-Ton')
        
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    
    
    
    
def Fehlerkurve_mseq_LuR(alpha):
    
    #s_test = np.lib.pad(m_sequence(10),(2*1023,2*1023),'constant',constant_values=(0,0))
    
    s5 = s.generate_Signal(1023,1024,4e6,2e6)
    
    
    randomPhase = np.linspace(0,np.pi,4)
    dt = np.linspace(0,1e-6,301)
    
    t5_fehler = np.zeros(dt.size)

    t5 = np.zeros(1)
    
    
    
    for i,phi in enumerate(randomPhase):
        for j in np.arange(301):
            for x in np.arange(1):
                
                r5 = ch.twowayChannel(0,dt[j],alpha,phi,100,s5,4e6)
                
                
                t5[x] = est.LundR(s5,r5,4e6,1023)

                
            t5_fehler[j] = np.mean(t5)

        dx = dt*3e8

        y5 = t5_fehler *3e8
        
        if(phi == 0):
            fig0 = plt.figure(19)
            plt.plot(dx,y5,label = '1023-Ton')
            plt.ylabel('rmse in m')
            plt.xlabel('Umweg in m')
            plt.ylim(-40,40)
            plt.legend()
            plt.title(u'Hüllkurve')
            fig0.show()
        fig = plt.figure(61)
        plt.plot(dx,y5, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 1023-Ton')
    fig.show()

    return y5 
    
    







def Fehlerkurve_ESPRIT():
    s1 = s.generate_Signal(4,1024,4e6,2e6) 
    s2 = s.generate_Signal(8,1024,4e6,2e6)
    s3 = s.generate_Signal(16,1024,4e6,2e6)
    s4 = s.generate_Signal(7,1022,4e6,2e6)
    
    
    randomPhase = np.linspace(0,np.pi,4)
    dt = np.linspace(0,1e-6,1001)
    
    t1_fehler = np.zeros(dt.size)
    t2_fehler = np.zeros(dt.size)
    t3_fehler = np.zeros(dt.size)
    t4_fehler = np.zeros(dt.size)
    
    t1 = np.zeros(1)
    t2 = np.zeros(1)
    t3 = np.zeros(1)
    t4 = np.zeros(1)
    
    
    
    for i,phi in enumerate(randomPhase):
        for j in np.arange(1001):
            for x in np.arange(1):
                
                r1 = ch.twowayChannel(0,dt[j],0.5,phi,100,s1,4e6)
                r2 = ch.twowayChannel(0,dt[j],0.5,phi,100,s2,4e6)
                r3 = ch.twowayChannel(0,dt[j],0.5,phi,100,s3,4e6)
                r4 = ch.twowayChannel(0,dt[j],0.5,phi,100,s4,4e6)
                
                t1[x] = sub.esprit(s1,r1,4e6,4)
                t2[x] = sub.esprit(s2,r2,4e6,8)
                t3[x] = sub.esprit(s3,r3,4e6,16)
                t4[x] = sub.esprit(s4,r4,4e6,7)
                
            t1_fehler[j] = np.mean(t1)
            t2_fehler[j] = np.mean(t2)
            t3_fehler[j] = np.mean(t3)
            t4_fehler[j] = np.mean(t4)
            
        dx = dt*3e8
        y1 = t1_fehler *3e8
        y2 = t1_fehler *3e8
        y3 = t1_fehler *3e8
        y4 = t1_fehler *3e8
        
        if(phi == 0):
            fig0 = plt.figure(19)
            plt.plot(dx,y1,label = '4-Ton')
            plt.plot(dx,y2,label = '8-Ton')
            plt.plot(dx,y3,label = '16-Ton')
            plt.plot(dx,y4,label = '7-Ton')
            plt.ylabel('rmse in m')
            plt.xlabel('Umweg in m')
            plt.legend()
            plt.title('Mehrwegeausbreitung')
            fig0.show()
            
            
        
        fig1 = plt.figure(20)
        plt.plot(dx,y1, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 4-Ton')
        
        fig2 = plt.figure(21)
        plt.plot(dx,y2, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 8-Ton')
        
        fig3 = plt.figure(22)
        plt.plot(dx,y3, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 16-Ton')
        
        fig4 = plt.figure(23)
        plt.plot(dx,y4, label = '$\phi= %i^{\circ}$' %(np.round(phi/(2*np.pi)*360)))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m')
        plt.ylim(-40,40)
        plt.legend(ncol = 2, title = 'Phasenoffset')
        plt.title(u'Hüllkurve für einen 7-Ton')
        
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    

def gesamt_Kanalauswertung():
    
    s1 = s.generate_Signal(7,1022,4e6,2e6)
    
    umweg = np.linspace(0,0.5e-6,1001)
    SNR = np.linspace(0,30,4)
    
    t_fehler1 = np.zeros(umweg.size)
    t_fehler2 = np.zeros(umweg.size)
    temp1 = np.zeros(10)
    temp2 = np.zeros(10)
    
    for i in np.arange(4):
        print i
        for j in np.arange(1001):
            np.random.seed(2)
            for x in np.arange(10):
                r = ch.twowayChannel(0,umweg[j],0.5,0,SNR[i],s1,4e6)
                temp1[x] = sub.esprit(s1,r,4e6,7)
                temp2[x] = sub.rtmusic(s1,r,4e6,7)
            t_fehler1[j]=np.sqrt(np.mean(temp1**2))
            t_fehler2[j]=np.sqrt(np.mean(temp2**2))
        dx = umweg*3e8
        y1 = t_fehler1 *3e8
        y2 = t_fehler2 *3e8
        
        fig1 = plt.figure(100)
        plt.plot(dx,y1,label = '%i dB-Hz ' %((i*10)+66))   
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m') 
        plt.ylim(0,25)
        plt.legend(title = r'$\frac{C}{N_0}$')
        plt.title('ESPRIT')
        fig1.show()
        
        fig2 = plt.figure(101)
        plt.plot(dx,y2,label = '%i dB-Hz' %((i*10)+66))
        plt.ylabel('rmse in m')
        plt.xlabel('Umweg in m') 
        plt.ylim(0,25)
        plt.legend(title = r'$\frac{C}{N_0}$')
        plt.title('ESPRIT')
        fig2.show()



def Physikalische_Kanalauswertung():
    
   
    s1 = s.generate_Signal(7,1022,4e6,2e6)
    

    Umweg = np.linspace(0,0.5e-6,151)[::10]
    LOS = np.linspace(0,0.5e-6,151)[::10]

    
    t_fehler1 = np.zeros((LOS.size,Umweg.size))
    temp1 = np.zeros(100)
    
    for i in np.arange(LOS.size):
        print i
        for j in np.arange(Umweg.size):
            for x in np.arange(100):
                r = ch.twowayChannelfor3D(LOS[i],LOS[i]+Umweg[j], 0, 0, s1, 4e6)

                temp1[x] = sub.esprit(s1,r,4e6,7)
                

            t_fehler1[i,j] = np.sqrt(np.mean(temp1**2)) * 3e8
        
    return  t_fehler1
    

  
def gesamtAuwertung_Plot(Z):
    
    
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        X, Y = np.meshgrid(np.arange(Z.shape[1]), np.arange(Z.shape[0]))
        
        ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
        #ax.contourf(X, Y, Z, zdir='x', offset=16, cmap=cm.coolwarm)
        #ax.contourf(X, Y, Z, zdir='y', offset=16, cmap=cm.coolwarm)
        plt.show()
        
        
        
                
                

            
def music_mehrwegeauswertung(typ):
    
    
        signal = s.generate_Signal(typ,1024,4e6,2e6)
    
        randomPhase = np.linspace(0,np.pi,4)
        dt = np.linspace(0,1.4e-6,75)
        t = np.zeros(1)
        t_fehler = np.zeros(dt.size)
        
        for i,phi in enumerate(randomPhase):
            for j,tau in enumerate(dt):
                for x in np.arange(1):
                    r = ch.twowayChannel(0,tau,0.5,phi,0,signal,4e6)
                    t[x] = sub.esprit(signal,r,4e6,typ)
                    #t[x] = sub.rtmusic(signal,r,2e6,typ)
                t_fehler[j] = np.sqrt(np.mean(t**2))
                #t_fehler[j] = np.std(t)
                #t_fehler[j] = np.mean(t)
                
                
            dx = dt*3e8
            y = t_fehler * 3e8               
            
            plt.figure(2)
            plt.plot(dx,y)   
            plt.ylabel('Fehler in [m]')
            plt.xlabel('Umweg in [m]')  
            plt.text(dx[np.argmax(np.abs(y))], y[np.argmax(np.abs(y))]+2 , '$\phi= %i$' %(np.round(phi/(2*np.pi)*360)))
            plt.title('Error bei 10dB ueber Laufzeit')
                  
                  
                  
def music_SNR(typ):
    
    # Hadamard-sequenz mit 4 Subträgern und einer länge 2048 was einem burst von 512us entspricht
    signal = s.generate_Signal(typ,1024,2e6,1e6)
    #phis =  np.linespace(0,np.pi,7)
    dt = np.linspace(0,1.4e-6,15)
    t = np.zeros(100)
    sigma_t = np.zeros(dt.size)
    SNR = np.linspace(0,20,3)
    for i in SNR:
        for j in dt:
            
            for x in np.arange(100):
                r = ch.twowayChannel(0,j,0.5,0,i,signal,2e6)
                t[x] = sub.rtmusic(signal,r,2e6,typ)
            sigma_t[np.where(dt == j)] = np.mean(t)
                
        
        y = sigma_t * 3e8  
        x = dt * 3e8
        plt.figure(1)         
        plt.plot(x,y)
        plt.ylabel('Fehler in [m]')
        plt.xlabel('Umweg in [m]')



def plot_varianz(sigma1, sigma2, sigma3, sigma4, sigma5, estimator = 'Estimatortype'):
    
    x = np.linspace(0,100,101) + 10*np.log10(4e6)
    
        
    
    fig1 = plt.figure(51)
    plt.plot(x, sigma1, color='b', label = '4-Ton')
    plt.plot(x, sigma2, color='g', label = '8-Ton')
    plt.plot(x, sigma3, color='r', label = '16-Ton')
    plt.plot(x, sigma4*1.2, color='c', label = '7-Ton')
    #plt.plot(x, sigma5, color='m', label = '1023-Ton')
    plt.xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
    plt.ylabel('rmse in m')
    plt.ylim(0,3)
    plt.xlim(66,140)
    plt.legend()
    plt.title(estimator)        
    fig1.show()

def zoomed_plot(sigma1, sigma2, sigma3, sigma4, sigma5, estimator = 'Estimatortype'):    
        
        x = np.linspace(0,100,101) + 10*np.log10(4e6)
        
        fig, ax = plt.subplots()
        ax.plot(x, sigma1,color='b', label='4-Ton')
        ax.plot(x, sigma2,color='g', label='8-Ton')
        ax.plot(x, sigma3,color='r', label='16-Ton')
        ax.plot(x, sigma4,color='c', label='7-Ton')
        ax.plot(x, sigma5,color='m', label='1023-Ton')
        ax.set_xlabel(r'$\mathbf{\frac{C}{N_0}}$ in dB-Hz')
        ax.set_ylabel('rmse in m')
        ax.set_xlim(66,110)
        ax.set_title(estimator)
        plt.legend()
        
        #axins = zoomed_inset_axes(ax,7.5,loc=5)
        axins = inset_axes(ax,4,1,loc=7)
        axins.plot(x, sigma1,color='b', label='4-Ton')
        axins.plot(x, sigma2,color='g', label='8-Ton')
        axins.plot(x, sigma3,color='r', label='16-Ton')
        axins.plot(x, sigma4,color='c', label='7-Ton')
        axins.plot(x, sigma5,color='m', label='1023-Ton')
        x1,x2,y1,y2 = 66, 80 , 0, 2
        axins.set_xlim(x1,x2)
        axins.set_ylim(y1,y2)
        axins.set_xticks([70,80])
        axins.set_yticks([1,2])
        mark_inset(ax, axins, loc1=2,loc2=4,ec='0.5')
        fig.show()
        
    
    
    
    
def test_zoom():
    
    x = np.linspace(-2*np.pi,2*np.pi,101)  
    s = np.sin(x)
    fig, ax = plt.subplots()
    ax.plot(x,s)
    
    axins = zoomed_inset_axes(ax,3.5,loc=1)
    axins.plot(x,s)
    
    x1,x2,y1,y2 = 1 , 2, 0.8, 1
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    
    axins.set_yticks([])
    axins.set_xticks([])
    
    mark_inset(ax,axins, loc1 = 2, loc2 = 3, fc = 'none', ec = 'g')    
    
    fig.show()