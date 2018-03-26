# kpca_demo
import numpy as np
import matplotlib.pyplot as plt

# data generation
    
# the following code generates 2 dimensional data. specifically, for each vector the first component is between 0 and alpha, while the second is between 0 and 1
n = 1000                       # number of data points
alpha = 2             # length/width ratio
s = np.array([alpha,1])
X = np.diag(s).dot(np.random.rand(2,n))          # uniformly distributed points on a rectangle

H = np.eye(n) - np.ones((n,n))/n           # create centering matrix

#%%
def custom_sdist(X):
    """
    Funktion that given a matrix X returns the squared pairwise distances 
    of the column vectors in matrix form
    """
    XX = np.dot(X.T, X)
    pdists = np.outer(np.diag(XX), np.ones(XX.shape[1]).T) + np.outer(np.ones(XX.shape[0]), np.diag(XX).T) - 2*XX
    return pdists

    
sigma = 1

def K(X):
    pdist = custom_sdist(X)
    K = np.exp(-(pdist)/(2*sigma**2))
    return K

k = 2                                   # number of eigenvectors

K = K(X)
K_centered = H.dot(K).dot(H)
_, s, Vt = np.linalg.svd(K_centered)
s = np.sqrt(s)

Y = np.dot(np.diag(s[:k]), Vt[:k, :])  #projection of the kernel matrix K onto the first two principal components

    
fig, axs = plt.subplots(1,2,figsize=(15,6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=.5, wspace=.1)
axs = axs.ravel()
for ip in range(k):
    tolerance = np.max(Y[ip,:])*0.05
    border = X[:,np.where((Y[ip,:]<tolerance) & (Y[ip,:]>-tolerance))]
    im = axs[ip].scatter(X[0,:], X[1,:], c=Y[ip,:])
    axs[ip].scatter(border[0,:], border[1,:], color = 'r')
    axs[ip].set_title('Color indicates value of PC {} at this point'.format(ip+1))
    
fig.colorbar(im)
plt.show()
print sigma    

# %%
alpha = np.arange(0,12)+1
sigma = np.array([1,1,2,4,6,9,12,17,21,25,31,36])

plt.figure()
plt.plot(alpha,sigma)
plt.xlabel('alpha')
plt.ylabel('sigma')
plt.show()