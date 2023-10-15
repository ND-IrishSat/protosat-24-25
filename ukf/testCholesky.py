import scipy
import time
import numpy as np

def sigma(means, cov, n, scaling):
    '''
    sigma
        creates sigma point matrix based on formula that represents distribution of means and cov

    @params
        means: mean of gaussian of estimated states so far (maybe???). Also first column of sigma matrix (1 x n)
        cov: covariance matrix of state (n x n)
        n: dimensionality of model
        scaling: how far from mean we distribute our points, used in sigma point formula
    @returns
        sigmaMatrix: matrix of sigma points (2 * n + 1, n) 
    '''
    # intialize 2N + 1 sigma points to zeroes
    sigmaMatrix = np.zeros((2*n+1,n))
    temp = np.zeros((n, n))

    for i in range(len(cov)):
      for j in range(len(cov[i])):
        #sigma point formula from website
        temp[i][j] = cov[i][j]  * (n + scaling)

    temp = scipy.linalg.sqrtm(temp)
    # first column of sigma matrix is means
    sigmaMatrix[0] = means

    # traverse N (10) dimensions, calculating all sigma points
    for i in range(n):
        sigmaMatrix[i + 1] = np.add(means, temp[i])  # mean + covariance
        sigmaMatrix[i + 1 + n] = np.subtract(means, temp[i])  # mean - covariance

    # return the sigma matrix (21 columns)
    return sigmaMatrix




tot = 0
iterations = 1000
n = 10
scaling = 3-n
#100 iterations for original sigma: .00064
for i in range(iterations):
    start=time.time()
    means = np.random.rand(1, n)
    cov = np.random.rand(10,10)
    sigma(means, cov, n , scaling)
    tot += (time.time()-start)

print(tot/iterations)