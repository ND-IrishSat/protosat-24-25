'''
run_UKF.py
Authors: Andrew Gaylord, Claudia Kuczun, Micheal Paulucci, Alex Casillas, Anna Arnett
Last modified 9/26/23

Runs IrishSat UKF on generated or real-time data and simulates CubeSat using pygame
'''

# implement EOMS for real-time control in pygame???

import numpy as np
import random
import pygame
from pygame.locals import * # Bad coding practice >:( Change later
from OpenGL.GL import * # Bad coding practice >:( Change later
from OpenGL.GLU import * # Bad coding practice >:( Change later
import time

from pyquaternion import Quaternion

from UKF_algorithm import *
# change for readability: import __ as ___ and change names later on
# see above comments >:(

from BNO055_MAGNETOMETER_BASIC import calibrate
from BNO055_MAGNETOMETER_BASIC import get_data

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
    ''' Tryna implement PyGame and OpenGL to this bish 

    You don't stop
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
    # print("Q_ARRAY: ", Q_array)
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

# data: list
# typecasting for readability/a little functionality
def check_zeros(data):
    ''' Checks validity of data

        Args:
            data (???): data point

    '''
    if int(data[0]) == 0 and int(data[1]) == 0 and int(data[2]) == 0:
        return True

if __name__ == "__main__":
    states = []

    # Initialize
    r=np.zeros(10)
    q=np.zeros(9)

    for i in range(9):
        r[i]=random.random()
        q[i]=random.random() * .1
    
    start=np.zeros(10)
    for i in range(len(start)):
        start[i] = random.random()
    # start = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    cov = np.zeros((10,10))
    for i in range(10):
        cov[i][i]=random.random()

    # calibrate sensors before getting data
    calibrate()
    #get_data()

    cnt = 0
    for i in range(1, 150):
        
        '''
        EOMSData = EOMs(np.array(start))
        print(EOMSData)
        game_visualize(np.array([EOMSData[:4]]), i)
        start = EOMSData
        '''

        if cnt == 0: 
            data = get_data()
            cnt += 1
        time.sleep(0.5)
        data = get_data()
        if check_zeros(data): continue # do not use data if B-field is all zeros

        print(f"new data = {data}")
        start, cov = UKF(start, cov, r, q, data)
        game_visualize(np.array([start[:4]]), i)

       # print("STATE AFTER {} RUNTHROUGH: {}".format(i, start))
        #print("COV AFTER {} RUNTHROUGH: {}".format(i, cov))
        
