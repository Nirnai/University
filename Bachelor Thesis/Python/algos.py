from __future__ import division

from itertools import permutations

from numpy import array as ar, floor, hstack, inf, r_, sqrt, sum
from numpy.fft import fft2, ifft2
from numpy.linalg import lstsq, svd


def fftcorr2d(a, b):
    corrshape = r_[a.shape] + r_[b.shape] - 1
    return ifft2(fft2(a, corrshape) * fft2(b, corrshape), corrshape)


# Linear Interpolation Between Closest Ranks
# different from numpy.percentile!
def percentile(data, p):
    assert(len(data))
    assert(p >= 0 and p <= 100)

    data = sorted(data)
    idx = p * len(data) / 100 - 0.5
    if idx <= 0:
        return data[0]
    elif idx >= len(data) - 1:
        return data[-1]
    else:
        iidx = int(floor(idx))
        didx = idx - iidx
        return (1 - didx) * data[iidx] + didx * data[iidx + 1]


# Total Least Squares
def totlstsq(A, B):
    d = A.shape[1]
    _, _, V = svd(hstack((A, B)))
    V = V.T.conj()
    return lstsq(V[d:, d:].T, -V[:d, d:].T)[0].T


def bestrmsematch(lst1, lst2):
    minrmse = inf
    for perm in permutations(xrange(len(lst2))):
        perm = ar(perm)
        rmse = sqrt(sum((ar(lst1) - ar(lst2)[perm]) ** 2) / len(lst1))
        if rmse < minrmse:
            minrmse = rmse
            optperm = perm
    return minrmse, optperm


def wrap(data, lower, upper):
    return (data - lower) % (upper - lower) + lower
