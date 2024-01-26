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
    
    i = 1
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
    initQ = np.array([0, 0, 0, 1])
    start = [0, 0, 0, 1, 0, 0.25, 0]
    w = np.array([0, .25, 0])
    cov = np.identity(7) * 0.1
    # constant B field
    B_true = np.array([0, 0, 1])
    reaction_speeds = np.zeros(3)

    noiseMagnitude = .001
    r = np.random.normal(0, noiseMagnitude, n)
    noiseMagnitude = .001
    q = np.random.normal(0, noiseMagnitude, n)

    ''' SHOULD NOISE BE 2D MATRIX??? '''

    # noise for magnetomer + angular velocity readings?
    sensorNoise = np.random.normal(0, 0.005, 2000)
    print(sensorNoise[:10])

    # for a in range(20):
    #     print(r[a],", ", end=" ")
    # for b in range(20):
    #     print(q[b],", ", end=" ")

    t0 = 0
    tf = 100
    i = 0

    # use attitude propagator to find actual quaternion for n steps
    propagator = AttitudePropagator(q_init=initQ, w_init=w)

    states = propagator.propagate_states(t0, tf, n)


    # calculate sensor b field for every time step
    # rotation matrix(q) * true B field + noise
    B_sens = np.array([np.matmul(hfunc.quaternion_rotation_matrix(states[0]), B_true)])
    
    for a in range(1, n):
        B_sens = np.append(B_sens, np.array([np.matmul(hfunc.quaternion_rotation_matrix(states[a]), B_true)]), axis=0)

    # need to add noise
    # B_sens[:, 0] += sensorNoise[1000:]
    # B_sens[:, 1] += sensorNoise[::2]
    # B_sens[:, 2] += sensorNoise[:1000]

    
    while i < 1000:

        # create sensor data matrix of magnetomer reading and angular velocity
        data = [0] * 6
        data[0] = B_sens[i][0]
        data[1] = B_sens[i][1]
        data[2] = B_sens[i][2]
        data[3] = w[0]
        data[4] = w[1]
        data[5] = w[2]

        # run ukf algorithm with B_true as u_k instead of gps data instead
        # print(start, cov, r[i], q[i], list(B_true), reaction_speeds, data)
        start, cov = UKF_algorithm.UKF(start, cov, r[i], q[i], list(B_true), reaction_speeds, data)
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