'''
run_UKF.py
Authors: Andrew Gaylord, Claudia Kuczun, Micheal Paulucci, Alex Casillas, Anna Arnett
Last modified 10/7/23

Runs IrishSat UKF on generated or real-time data and simulates CubeSat using pygame

TODO:
    interface with gps sensor, find what frame it gives us (ECEF or ECI?)
    find correct value for zCov and noise (r, q)
    update EOMs with new inertia
    adding gps component/control input vector for EOMs?
    optimize for loops and numpy arrays
    test with different data sets
    remake sigma points?
'''

import numpy as np
import random
import pygame
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
import time

from pyquaternion import Quaternion

import UKF_algorithm

# from BNO055_magnetometer_basic import calibrate
# from BNO055_magnetometer_basic import get_data

##################################################################################################################
## AttitudePropagator - the important part, propagates states numerically using constant angular velocity model ##
##################################################################################################################
class AttitudePropagator:
    def __init__(self, q_init: np.ndarray, w_init: np.ndarray):
        ''' For generating a set of states for a constant angular velocity model, given an initial angular velocity vector

        Args:
            q_init (np.ndarray, (1x4)): initial quaternion
            init_w (np.ndarray, (1x3)): initial angular velocity

        '''
        self.q_init = q_init
        self.w = w_init

    def attitude_propagator(self, quaternion: np.ndarray, w_vect: np.ndarray, dt: float, normalize: bool):
        ''' Function to propagate quaternion/attitude given an angular velocity vector w. Equations of motion
        derived from https://ahrs.readthedocs.io/en/latest/filters/angular.html

        Form of propagtion equation is 
                    q_t+1 = K*q_t
        where K is
                ( cos(w*t/2)*I_4 + 1/w * sin(w*t/2) * Big_Omega(w) )

        See reference for more info

        Args:  
            quaternion (np.ndarray): quaternion before propagation step
            w_vect (np.ndarray): angular velocity array
            dt (float): timestep
            normalize (bool): boolean, where if True, normalize the new quaternion
        '''
        
        # Store norm of w_vect as separate variable, as w_magnitude. Also separate out components of w
        w_magnitude = np.linalg.norm(w_vect)
        w_x = w_vect[0]
        w_y = w_vect[1]
        w_z = w_vect[2]

        # Calculate factors for calculating matrix K
        cos_fact = np.cos(w_magnitude * dt/2)
        sin_fact = np.sin(w_magnitude * dt/2)
        big_omega =  np.array([[0, -w_x, -w_y, -w_z],
                               [w_x,  0,  w_z, -w_y],
                               [w_y, -w_z, 0,   w_x],
                               [w_z,  w_y, -w_x,  0]])

        # Putting together to calculate matrix K
        K = cos_fact * np.identity(4) + (1/w_magnitude) * sin_fact * big_omega

        # Propagate q_t+1
        new_quaternion = K @ quaternion

        # Normalize -- technically not needed but recommended to get rid off round off errors
        new_quaternion = new_quaternion / np.linalg.norm(new_quaternion)

        return new_quaternion

    def propagate_states(self, t_0: float, t_f: float, N: int):
        ''' Propagate attitude from its initial conditions from t = t_0 to t = t_f, using
        N number of steps

        Args:
            t_0 (float): starting time
            t_f (float): ending time
            N (int): number of steps between t_0 and t_f
        '''

        # Calculate dt from time interval and N number of steps
        dt = (t_f - t_0)/N 

        # Allocate state matrix and save initial quaternion to first row of state matrix
        quaternion_states = np.zeros((N+1, 4))
        quaternion_states[0] = self.q_init

        # Propagate N number of times
        for i in np.arange(1, N+1, 1):
            quaternion_states[i] = self.attitude_propagator(quaternion_states[i-1], self.w, dt, normalize = True)

        return quaternion_states

##############################################################################################################
## Display and 3D rendering stuff, not very important tbh, unless you're into it, jk jk haha... unless?   ####
##############################################################################################################
def rotate_points_by_quaternion(points, quaternion):
    ''' Points are n by 3
    
    '''
    v_prime = np.zeros(points.shape)
    
    quaternion = Quaternion(quaternion)
    
    for i in np.arange(len(points)):
        v = points[i]
        v_prime[i] = quaternion.rotate(v)
    
    return v_prime

def Draw(vertices, edges):
    
    glBegin(GL_LINES)

    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])

    glEnd()


def game_visualize(states, i):
    '''
    game_visualize 
        uses pygames and AttitudePropagator class to visualize simple cube with our data (written by Juwan)

    @params
        states: quaternion matrix to visualize (1 x 4)
        i: index of what step we are on (must start at 1 to properly initialize)
    ''' 
    vertices_cube = np.array([[1, -1, -1], [1, 1, -1],[-1, 1, -1],[-1, -1, -1],\
        [1, -1, 1],[1, 1, 1],[-1, -1, 1],[-1, 1, 1],[0,0,0],[0,-1,0],\
        [-0.8,0.8,1],[-0.8,0.6,1],[-0.6,0.6,1],[-0.6,0.8,1],\
        [-0.5, 0.6, 1], [-0.5,0.8,1], [-0.3,0.8,1],[-0.3,0.7,1],[-0.5,0.7,1], ])
    edges_cube = (
        (0,1),
        (0,3),
        (0,4),
        (2,1),
        (2,3),
        (2,7),
        (6,3),
        (6,4),
        (6,7),
        (5,1),
        (5,4),
        (5,7),
        (8,9),
        (0,6),
        (3,4),
        (10,11),
        (11,12),
        (12,13),
        (14,15),
        (15,16),
        (16,17),
        (17,18)
        )
    vertices_target_line = (
        (0,0,0),
        (0,-5,0)
    )
    edges_target_line = (
        (0,1),
        (0,0)
        )

    if(i == 1):
        pygame.init()
        display = (900, 700)
        pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
        gluPerspective(45, (display[0]/display[1]), 0.1, 50.0)
        glTranslatef(0.0, 0.0, -5)


    clock = pygame.time.Clock()
    Q_array = states[:, :4]  # array of quaternions (for each entry, take first four items -> array of [a,b,c,d])
    for i in range(0, len(Q_array)):
        Q_array[i][0] = 0
    i = 0

    num_states = states.shape[0]

    while True:
        clock.tick_busy_loop(15)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        vertices_cube_new = rotate_points_by_quaternion(vertices_cube, Q_array[i])
        Draw(vertices_cube_new, edges_cube)
        Draw(vertices_target_line, edges_target_line)
        pygame.display.flip()
        i += 1  # go to next index in quaternion

        if i == num_states:
            break


def check_zeros(data):
    '''
    check_zeros
        checks validity of data input

    @params
        data: input from sensor real-time (1 x 6)
    @returns
        true or false 
    '''
    if int(data[0]) == 0 and int(data[1]) == 0 and int(data[2]) == 0:
        return True


def run_ukf_textfile(start, cov, r, q, filename):
    '''
    run_ukf_textfile
        runs and visualizes UKF algorithm on input data file

    @params
        start: initialized state (1 x n)
        cov: initialized covariance matrix (n x n)
        r: noise vector for predictions (1 x n)
        q: noise vector for sensors (1 x m)
        filename: text file of cleaned sensor data to read from (any length)
    '''
    f = open(filename, "r")
    data = f.readline()
    splitData = data.split(",")
    splitData = [float(x) for x in splitData]
    # start[0] = splitData[0]
    # start[1] = splitData[1]
    # start[2] = splitData[2]
    # start[3] = splitData[3]
    i = 1
    u_k = []
    while(data):
        # get gps data and add time stamp
        u_k = get_gps_data()
        u_k = ecef_to_latlong(u_k[0], u_k[1], u_k[2])
        u_k.append(2023.8123)
        # run ukf and visualize output
        start, cov = UKF_algorithm.UKF(start, cov, r, q, u_k, splitData)
        game_visualize(np.array([start[:4]]), i)

        # continue to get data from file until empty
        data = f.readline()
        if(data == ''):
            break
        splitData = data.split(",")
        splitData = [float(x) for x in splitData]
        i+=1

    f.close()


def run_ukf_sensor(state, cov, r, q):
    '''
    run_ukf_sensor
        runs and visualizes UKF algorithm using real-time data from magnetometer/pi

    @params
        start: initialized state (1 x n)
        cov: initialized covariance matrix (n x n)
        r: noise vector for predictions (1 x n)
        q: noise vector for sensors (1 x m)
    '''

    # uncomment BNO055 imports to use

    # i = 1
    # u_k = []
    # calibrate()

    # while(1):
    #     time.sleep(0.5)
    #     data = get_data()
        # u_k = get_gps_data()
        # u_k = ecef_to_latlong(u_k[0], u_k[1], u_k[2])
        # u_k.append(2023.8123)
    #     # do not use data if B-field is all zeros 
    #     if check_zeros(data): continue 

    #     start, cov = UKF_algorithm.UKF(start, cov, r, q, u_k, data)
    #     game_visualize(np.array([start[:4]]), i)

    #     i += 1


if __name__ == "__main__":

    # Initialize noises and starting state/cov to random values
    n = 10
    m = 9

    r=np.zeros(n)
    q=np.zeros(m)
    for i in range(m):
        r[i]=random.random()
        q[i]=random.random() * .1

    start = np.random.rand(n)

    cov = np.zeros((n,n))
    for i in range(n):
        cov[i][i]=random.random()
    
    filename = "sensor_data_2.txt"

    # tests ukf with pre-generated and cleaned data file
    run_ukf_textfile(start, cov, r, q, filename)

    # must uncomment BNO055 imports to use in real-time with sensor
    # run_ukf_sensor(start, cov, r, q)
