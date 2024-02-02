'''
ideal_test_cases.py
Authors: Andrew Gaylord
Last modified 1/20/24

    tests two ideal test cases for UKF algorithm: nonmoving and constant angular velocity

'''

import numpy as np
import random
import pygame
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
import datetime as dt

from pyquaternion import Quaternion

import UKF_algorithm
from run_UKF import *
import hfunc

# from ukf.PySOL import spacecraft as sp
# from ukf.PySol.spacecraft import *
from PySOL.sol_sim import *
import PySOL.spacecraft as sp
import PySOL.orb_tools as ot

from filterpy.kalman import KalmanFilter
from filterpy.common import Q_discrete_white_noise



def run_basic_test():
    '''
    ideal test case for UKF: zeroed out, nonmoving cubesat
    0 reaction wheel speed and angular velocity, constant magnetic field
    no frame of reference transformation or gps implementation

    THINGS CHANGED FOR TEST CASE ONE:
        removed hfunc and passed u_k as magnetic field
        
        removed normalization

        removed zCov
        
        need to find optimal starting r, q for each case

        got rid of different dimensionality for m: instead, everything uses 7

        switched covariance calculation to option 1

        fixed crosscovariance calculation bug

        created diagonal starting covariance

        what is control input vector in EOMs?
    '''

    n = 7

    # quaternion and angular velocity should be zero
    start = np.array([0, 0, 1, 0, 0, 0, 0])
    # start = np.array([0, 1, 0, 0, 0, 0, 0])


    # start with diagonal covariance
    cov = np.identity(n) * 0.5
    # print("Starting cov: ", cov)


    # we want magnetomer reading to be constant, rest to be 0
    data = [0, 0, 1, 0, 0, 0, 0]
    # data = [0, 1, 0, 0, 0, 0]

    # data = [0, 1, 0, 0, 0, 0, 0]

    # gaussian noise
    noiseMagnitude = .005
    r = np.random.normal(0, noiseMagnitude, 1000)
    noiseMagnitude = .0025
    q = np.random.normal(0, noiseMagnitude, 1000)

    
    # functioning [0, 0, 1, 0, 0, 0, 0] noises
    r[:10] = [ 0.00072222,  0.00384547, -0.00737526,  0.00585633, -0.00058933,  0.00142955,
  0.0052347,   0.0036605,   0.00506041, -0.00103479]
    q[:10] = [-0.00375006, -0.00071075,  0.0004241,  -0.0014944,   0.00222291,  0.00088999,
  0.00338481,  0.0027458,   0.00346013,  0.00152124]
    # functioning [0, 1, 0, 0, 0, 0, 0] noises
#     r[:10] = [-0.00295791, -0.00030439, -0.00266072,  0.00867703, -0.0003597,  -0.0091533,
#  -0.0048222,  -0.00545613, -0.01036766, -0.00077585]
#     q[:10] = [ 0.00083288,  0.0012302,  -0.00149358, -0.00445953,  0.0005876,  -0.00232486,
#   0.00209268,  0.00337717,  0.00154955, -0.00357763]


    # edit code so that lat/long/hieght are not needed 
    # make u_k = magnetic field for this test only 
    # not even used in this test case bc hfunc is disabled
    u_k = np.zeros(3)

    # control input vector for eoms, zero for this test
    reaction_speeds = np.zeros(3)
    
    i = 0
    while(1):
        start, cov = UKF_algorithm.UKF(start, cov, r[i], q[i], u_k, reaction_speeds, data)

        game_visualize(np.array([start[:4]]), i)
        
        # uncomment for cube to move
        # if 100 > i > 25 and data[1] < 1:
        #     data[1] += .05
        # elif i > 100 and data[1] > 0:
        #     data[1] -= .05
        i += 1


def run_moving_test():
    ''''
    b field is the same but spinning

    need h func
    constant angular velocity and true b field
    same position (don't need lat/long again, j use true b field) and hence b field
    no reaciton wheel speed
    building data before we run:
        B true = [1, 0, 0] state. this can be u_k that is passed to hfunc. unchanging, omega
        B measurement = [changing accoridng to imu]. sen. calculate before running
        w true = [0, 0, 1]
    need true quaternion for every time step
    to find n number of b field measurement (sen) states: rotation (q) * B sen + noise
    W sen = w true + noise
    data = [0, B sen, W sen]

    '''

    # find n and use propogator to find quaternion
    n = 1000
    dt = 0.1
    dim = 7
    dimMes = dim - 1
    speed = 1
    # initQ = np.array([0, 0, 1, 0])
    # start = [0, 0, 1, 0, speed, 0, 0]
    initQ = np.array([1, 0, 0, 0])
    start = [1, 0, 0, 0, 0, speed, 0]
    w = np.array([0, speed, 0])
    cov = np.identity(dim) * 5e-10
    # constant B field
    B_true = np.array([0, 0, 1])
    reaction_speeds = np.zeros(3)

    # if model is less reliable: q > r
    # measurement noise (m x m i think)
    noiseMagnitude = 0.02
    # r = np.random.normal(0, noiseMagnitude, n)
    r = np.diag([noiseMagnitude] * dimMes)

    # process noise (n x n)
    noiseMagnitude = 0.005
    # q = np.random.normal(0, noiseMagnitude, n)

    '''
    q = Q_discrete_white_noise(dim=3, dt=dt, var=0.01**2, block_size=2)
    # https://filterpy.readthedocs.io/en/latest/common/common.html
    # https://github.com/rlabbe/filterpy/blob/master/filterpy/kalman/UKF.py

    q = np.zeros((dim, dim))
    q [0][:4] = [2.77777778e-12, 8.33333333e-11, 1.66666667e-09, 1.66666667e-08]
    q [1][:4] = [8.33333333e-11, 2.50000000e-09, 5.00000000e-08, 5.00000000e-07]
    q [2][:4] = [1.66666667e-09, 5.00000000e-08, 1.00000000e-06, 1.00000000e-05]
    q [3][:4] = [1.66666667e-08, 5.00000000e-07, 1.00000000e-05, 1.00000000e-05]
    q [4][4:] = [2.5e-09, 5.0e-08, 5.0e-07]
    q [5][4:] = [5.0e-08, 1.0e-06, 1.0e-05]
    q [6][4:] = [5.0e-07, 1.0e-05, 1.0e-05]

    q *= 100
        # redundant measurement noise covariance estimation

    # https://hackmd.io/@lancec/H1Ay6CizK

    
    # autocovariance least squares python
    https://quant.stackexchange.com/questions/8501/how-to-tune-kalman-filters-parameter
    https://www.sciencedirect.com/science/article/pii/S0959152407001631

    # robust/adaptive
    https://ieeexplore.ieee.org/abstract/document/6626597?casa_token=oroiTLL1ZggAAAAA:5xW15OXvGpCazSla1h_okSAkxqDG1Qd97YpvJShmNIPHWya9Xy44HEdF0boURmmVAItZsjN1XRc
    https://www.mdpi.com/1424-8220/18/3/808
    https://ieeexplore.ieee.org/document/7231764
    https://www.hindawi.com/journals/mpe/2015/218561/
    '''
    q = np.diag([noiseMagnitude] * dim)

    ''' NOISE SHOULD BE 2D MATRIX '''

    t0 = 0
    tf = 100
    i = 0

    # use attitude propagator to find actual quaternion for n steps
    propagator = AttitudePropagator(q_init=initQ, w_init=w)

    states = propagator.propagate_states(t0, tf, n)

    # calculate sensor b field for every time step
    # rotation matrix(q) * true B field + noise
    B_sens = np.array([np.matmul(hfunc.quaternion_rotation_matrix(states[0]), B_true)])
    # print(states[:10])
    # print("first calc: ", hfunc.quaternion_rotation_matrix(states[2]))
    # print(np.matmul(hfunc.quaternion_rotation_matrix(states[2]), B_true))
    
    for a in range(1, n):
        B_sens = np.append(B_sens, np.array([np.matmul(hfunc.quaternion_rotation_matrix(states[a]), B_true)]), axis=0)

    # need to add noise?

    
    while i < 1000:

        # create sensor data matrix of magnetomer reading and angular velocity
        data = [0] * dimMes
        data[0] = B_sens[i][0]
        data[1] = B_sens[i][1]
        data[2] = B_sens[i][2]
        data[3] = w[0]
        data[4] = w[1]
        data[5] = w[2]

        # run ukf algorithm with B_true as u_k instead of gps data instead
        # print(start, cov, r[i], q[i], list(B_true), reaction_speeds, data)
        # start, cov = UKF_algorithm.UKF(start, cov, r[i], q[i], list(B_true), reaction_speeds, data)
        start, cov = UKF_algorithm.UKF(start, cov, r, q, list(B_true), reaction_speeds, data)
        # start, cov = UKF_algorithm.UKF(start, cov, r, q, data[:3], reaction_speeds, data)


        print("Data: ", data)
        print("State: ", start[:4])
        print("Ideal: ", states[i])
        print("")

        game_visualize(np.array([start[:4]]), i)
        # game_visualize(np.array([states[i]]), i)
        i += 1






if __name__ == "__main__":

    # run_basic_test()
    run_moving_test()