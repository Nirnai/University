# -*- coding: utf-8 -*-
"""
 Created on Thu Jun 15 12:54:54 2017

@author: rnirn
"""
from os import listdir
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt

# majority vote
from collections import Counter

# %%
# Training set
img_names0 = listdir('./yaleBfaces/subset0')

T = np.zeros((50*50,len(img_names0)))
T_Label = np.zeros(len(img_names0))

for i,name in enumerate(img_names0):
    T[:,i] = sci.misc.imread('./yaleBfaces/subset0/' + name).reshape(-1)
    T_Label[i] = int(name[6:8])
  

# %%
# Testset
# Creating a dict S with each subset beeing one key
S = {}
S_Label = {}
for i in np.arange(4)+1:
    path = './yaleBfaces/subset%i'%i
    img_names = listdir(path)
    S['subset%i'%i] = np.zeros((2500,1))
    S_Label['subset%i'%i] = np.array([])
    temp = np.zeros((50*50,len(img_names)))
    for j,name in enumerate(img_names):
        temp[:,j] = sci.misc.imread(path+ '/' + name).reshape(-1)
        S_Label['subset%i'%i] = np.append(S_Label['subset%i'%i],int(name[6:8]))
    S['subset%i'%i] = np.concatenate((S['subset%i'%i],temp),axis=1)
    S['subset%i'%i] = S['subset%i'%i][:,1:]


# %% preprocessing data 
T_mean = np.repeat(np.mean(T, axis=0).reshape(1,70),2500,axis = 0)
T_std = np.repeat(np.std(T, axis=0).reshape(1,70),2500,axis = 0)
T = (T - T_mean)/T_std
    
T_generic = np.mean(T, axis = 1)
T = T - np.repeat(T_generic.reshape(2500,1), T.shape[1], axis=1)
#
for i in np.arange(4)+1:
    S_mean = np.repeat(np.mean(S['subset%i'%i], axis=0).reshape(1,S['subset%i'%i].shape[1]),2500,axis = 0)
    S_std = np.repeat(np.std(S['subset%i'%i], axis=0).reshape(1,S['subset%i'%i].shape[1]),2500,axis = 0)
    S['subset%i'%i] = (S['subset%i'%i]-S_mean)/S_std 
    S['subset%i'%i] = S['subset%i'%i] - np.repeat(T_generic.reshape(2500,1), S['subset%i'%i].shape[1], axis=1) 
#  
# 

# %%

def singluarVectors(T):
    U,s,V = np.linalg.svd(T)
    U = U[:,:20]
    return U

# %%
U = singluarVectors(T)
# %%
# error before preprocessing
def error_rate(T, S, T_Label, S_Label, U, k):
    
    # global variables to calculate error rate
    total_clasif = 0.0
    false_clasif = 0.0
    error_rate = np.zeros(4)
    
    # projector on a k-dimensional space
    U = U[:,:k]
    
    # projecting data on lower dimensional space
    T_proj = U.dot(U.T.dot(T)) 
    #T_proj = U.T.dot(T)
    S_proj = {}
    for i in np.arange(4)+1:
        S_proj['subset%i'%i] = U.dot(U.T.dot(S['subset%i'%i])) 
        #S_proj['subset%i'%i] = U.T.dot(S['subset%i'%i])
        for j in np.arange(S_proj['subset%i'%i].shape[1]):
            temp = np.repeat(S_proj['subset%i'%i][:,j:j+1],T_proj.shape[1], axis=1)
            distance = np.sqrt(np.sum((T_proj - temp)**2, axis = 0))
            # three k nearest neighbours
            k_nearest = distance.argsort()[:3]
            print k_nearest
            # choose label
            c = Counter(k_nearest)
            #if c.most_common()[0][1] > 1:
            predicted_label = T_Label[c.most_common()[0][0]]
            #else:
                #predicted_label = T_Label[k_nearest[0]]
            print predicted_label
            #weighted k nearest neighbors
            #predicted_label = np.round(np.sum(k_nearest * 1/distance[k_nearest]))
            
            if predicted_label != S_Label['subset%i'%i][j]:
                false_clasif = false_clasif + 1
            total_clasif = total_clasif + 1
        error_rate[i-1] = false_clasif/total_clasif
     
    return error_rate





# %%
# Outputs for data set

error = np.zeros((20,4))

for k in np.arange(20):
    error[k,:] = error_rate(T,S,T_Label, S_Label, U, k)



# %%
U2 = U[:,3:]
error2 = np.zeros((17,4))

for k in np.arange(17):
    error2[k,:] = error_rate(T,S,T_Label, S_Label, U2, k)


# %% subtask 1
plt.figure()
plt.subplot(231)
plt.title('first Eigenface')
plt.imshow(U[:,0].reshape(50,50),cmap='gray')

plt.subplot(232)
plt.title('second Eigenface')
plt.imshow(U[:,1].reshape(50,50),cmap='gray')

plt.subplot(233)
plt.title('third Eigenface')
plt.imshow(U[:,2].reshape(50,50),cmap='gray')

plt.subplot(234)
plt.title('fourth Eigenface')
plt.imshow(U[:,3].reshape(50,50),cmap='gray')

plt.subplot(235)
plt.title('fifth Eigenface')
plt.imshow(U[:,4].reshape(50,50),cmap='gray')

plt.subplot(236)
plt.title('sixt Eigenface')
plt.imshow(U[:,5].reshape(50,50),cmap='gray')
# %% subtask 2   

plt.figure()
plt.plot(error[:,0], label='subset1' )
plt.plot(error[:,1], label='subset2' )
plt.plot(error[:,2], label='subset3' )
plt.plot(error[:,3], label='subset4' )
plt.xlabel('k')
plt.ylabel('error rate in %')
plt.legend()




# %%  subtask 3   
plt.figure()

plt.plot(error2[:,0], label='subset1' )
plt.plot(error2[:,1], label='subset2' )
plt.plot(error2[:,2], label='subset3' )
plt.plot(error2[:,3], label='subset4' )
plt.xlabel('k')
plt.ylabel('error rate in %')
plt.legend()
