import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils.extmath import randomized_svd
from mpl_toolkits.mplot3d import Axes3D

m = 50
n = 200
X = np.random.randn(2, m)/10

for idx in range(n):
    t = np.random.randn(2, 2)
    tmp = t[:, 0:1]/np.linalg.norm(t[:, 0]) + t[:, 1]/np.linalg.norm(t[:, 1:2])/10
    X = np.hstack((X, tmp))
    
plt.figure()
plt.scatter(X[0, m+1:m+n], X[1, m+1:m+n])
plt.scatter(X[0, 1:m], X[1, 1:m], c='r')

plt.show()

#%% 
# Kernels 

sigma = 0.1
c = 0
d = 100
theta = 1

gauss_kernal = lambda x,y: np.exp(- np.linalg.norm(x-y)**2/(2*sigma**2))

poly_kernal = lambda x,y: (x.T.dot(y) + c)**d

wave_kernal = lambda x,y:  (np.sin(np.linalg.norm(x-y)/theta)*theta)/np.linalg.norm(x-y)

#%% generate gram matrix
def kgram(X, kappa):
    p, Nx = X.shape
    K = np.zeros((Nx,Nx))
    
    for i in np.arange(Nx):
        for j in np.arange(Nx):
            K[i,j] = kappa(X[:,i], X[:,j])
    return K


# %%
def kpca_v2(X, kappa, k):
    
    
    p, N_X = X.shape
    K = kgram(X, kappa)
    H = np.eye(N_X) - np.ones((N_X, N_X))/N_X
    K_centered = H.dot(K).dot(H)
    _, s, Vt = randomized_svd(K_centered, n_components=k)
    s = np.sqrt(s)
    S = np.dot(np.diag(s), Vt)
    
    return S
    
    
#%%
S = kpca_v2(X, gauss_kernal, 3)

#%%
plt.figure()
plt.scatter(S[0, m+1:m+n], S[1, m+1:m+n])
plt.scatter(S[0, 1:m], S[1, 1:m], c='r')
plt.show()

#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(S[0,:], S[1,:], S[2,:])
plt.show()