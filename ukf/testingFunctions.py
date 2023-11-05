'''
testingFunctions.py
Authors: Andrew Gaylord and Alex Casillas
Last modified: 10/31/2023

Speed tests for arbitrary functions
'''

import time
import numpy as np
import scipy
import timeit
import run_UKF
import UKF_algorithm
import random



def speedTest (func, inputs, n):
    '''
    speedTest
        takes arbitrary function and finds average speed (in milliseconds) over n iterations

    @params
        func: function we want to test
        inputs: list of inputs for function
        n: number of iterations to run speed test on
    @returns
        average: average time to complete n trials (in milliseconds)
    '''
    start = time.time()

    for i in range(n):
        func(*inputs)
    
    total = time.time() - start

    average = (total / n) * 1000

    print("Average time over {} iterations: {} ms".format(n, average))
    return average


def alex(a, b, c):
    sum = 0
    for i in range(c):
        sum += a * b


# speedTest(alex, [3, 5, 10000], 10000)


n = 10
m = 9
scaling = 3-n
means = np.random.rand(n)
cov = np.random.rand(10,10)
r=np.zeros(n)
q=np.zeros(m)
for i in range(m):
    r[i]=random.random()
    q[i]=random.random() * .1
u_k = []
data = np.random.rand(6)

# speedTest(UKF_algorithm.sigma, [means, cov, n, scaling], 10000)
# about half a millisecond

speedTest(UKF_algorithm.UKF, [means, cov, r, q, u_k, data], 10000)
# about 7.9 ms





