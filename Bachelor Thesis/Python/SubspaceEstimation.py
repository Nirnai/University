# -*- coding: utf-8 -*-
"""
Created on Tue Mar 08 14:30:36 2016

@author: Nirnai
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig, eigh
from algos import totlstsq


def ChannelCovMatrix(s,r):
    
    r"""
    Calculates the Channel Covarience matrix with the recieved signal and the sent signal
    
    Parameters
    ----------
    s : Send signal
    r : recieved signal
    typ: typ of the spectrum (4-,8-,16-subcarrier or mseq)
    """
    
    

#    spec = spec[np.where(np.abs(spec_ref) > 1e-4)]
#    spec_ref = spec_ref[np.where(np.abs(spec_ref) > 1e-4)]
    
    S = np.fft.fftshift(np.fft.fft(s)[::5])
    R = np.fft.fftshift(np.fft.fft(r)[::5])
    
    R = R[np.where(np.abs(S) > 1e-4)]
    S = S[np.where(np.abs(S) > 1e-4)]
    
#    freq = np.fft.fftshift(np.fft.fftfreq(S.size,0.25))
        
#    # Bandbegrenzen des Signals
#    lower_bandlimit = np.where(freq>-1.001)[0][0]
#    upper_bandlimit = np.where(freq>0.999)[0][0]
#
#    S = S[lower_bandlimit:upper_bandlimit]
#    R = R[lower_bandlimit:upper_bandlimit]
#    S = S[np.where(np.abs(S) > 1e-4)]
#    R = R[np.where(np.abs(S) > 1e-4)]

    # Creates the Channel frequency response 
    Y = (R/S).reshape(-1,1)
    
    # Covariance Matrix of the Channel frequency response
    Cy = 1/len(Y)*Y.dot(np.conjugate(Y).T)    
    
    Cy = Cy[:-1,:-1] + Cy[1:,1:]
    
    
    
#    J = np.flipud(np.eye(len(Cy)))
#
#    Cy += J.dot(Cy.conj()).dot(J)
       
    return Cy
    
    
    

   
def esprit(s,r,fs,typ):

    # Creates the Channel covarience matrix
    Cy = ChannelCovMatrix(s,r)

    # Smoothing
    Cy = Cy[:-1,:-1] + Cy[1:,1:]
    
    # estimate signal subspace

    sigma, U = eigh(Cy)

    U_s = U[:, np.argsort(sigma)[-1:]]


    # estimate frequencies

    U_s1 = U_s[:-1]

    U_s2 = U_s[1:]

    Psi = totlstsq(U_s2, U_s1)

    omega = eig(Psi)[0]
    t = 2*typ* np.angle(omega)/(2*np.pi*fs)
    

    return t[np.argsort(np.abs(t))][0]
    
    
    
    
    
def music(s,r,fc):
    
    
    # Creates a Channel covariance matrix    
    Cy = ChannelCovMatrix(s,r,4)
    
    # Smoothing
    Cy = Cy[:-1,:-1] + Cy[1:,1:]
    
    S = np.fft.fftshift(np.fft.fft(s)[::5])
    R = np.fft.fftshift(np.fft.fft(r)[::5])
    f = np.fft.fftshift(np.fft.fftfreq(len(R), 1/2/2))
    f_Y = f[np.where(S>1e-4)]
    
    # Eigenwerte und Eigenvektoren aus der Hermitischen Covarianz Matrix berechnen
    w, v = np.linalg.eigh(Cy)
    
    # Sortiert die Eigenvektoren nach den größten Eigenwerten und schmeißt je nach anzahl der Signale die n letzten werte weg
    Vr = v[:,np.argsort(w)][:,0:-2]
    
    # normaler MUSIC:
    fs = np.linspace(-1*np.pi, 1*np.pi, num=360+1)
    S = np.empty(len(fs))
    for i, omega in enumerate(fs):
        
        a = np.exp(-1j * omega*f_Y[:-1])
        S[i] = 1/(a.conj().T.dot(Vr).dot(Vr.conj().T).dot(a))
        
    t = np.linspace(-1, 1, num=len(S))
    
    tau = t[np.argmax(S)]/fc
        
    plt.figure('MUSIC')
    plt.plot(np.linspace(-1, 1, num=len(S)), S/np.max(S))
    plt.ylim(0, 1.1)
        
    return tau
    
    
def rtmusic(s,r,fs,typ):
    
#    S = np.fft.fftshift(np.fft.fft(s)[::5])
#    S = S[np.where(np.abs(S) > 1e-4)]
#    S = np.abs(S)
        
    #A = np.diagflat(S)
    #A = A[:-1,:-1] + A[1:,1:]
    


    # Creates the Channel covarience matrix
    Cy = ChannelCovMatrix(s,r)
    
    
    # calculates the Eigenvalues w and their coresponding Eigenvektors v
    w, v = eigh(Cy)
    
    # Extracts all Eigenvektors wich belong to the Noisespace
    Vn = v[:,np.argsort(w)[:-1]]
    
    # Creates the Noise covarience matrix wich are the polinomial coefficients
    Cn = Vn.dot(Vn.T.conj())
    p = np.empty(shape = (2 * len(Cn) - 1), dtype=complex)
    
    for i in xrange(len(Cn)):  
        p[len(Cn) - 1 + i] = sum(np.diag(Cn, i))
        p[len(Cn) - 1 - i] = sum(np.diag(Cn, -i))
    
    
    # Estimate Frequencies    
    rt = np.roots(p)
    rt = np.select([abs(rt) < 1], [rt])
    omega = rt[np.argsort(abs(rt))[-1:]]
    
    
        
    t = 2*typ* np.angle(omega)/(2*np.pi*fs)
    
     
    return t[np.argsort(np.abs(t))][0]
    
    