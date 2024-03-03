'''
run_UKF.py
Authors: Andrew Gaylord, Claudia Kuczun, Michael Paulucci, Alex Casillas, Anna Arnett
Last modified 10/7/23

Runs IrishSat UKF on generated or real-time data and simulates CubeSat using pygame

TODO:
    Hall sensor reading function: COMPLETED? Clean up + document motor_interface
    are motor scripts (simple something) and results of integration testing updated to github?
    find correct values of q and r noise (check out latex presentation/ask michael)
    find freshman who wants to learn UKF
    add ukf latex to github + clean up folder organization of repo
    big latex adcs document???

    add to latex documentation: testing cases/results, R and Q research/results, more explanation
    update EOMs to 4 reaction wheels (and gps data...?)
    interface with gps sensor, find what frame it gives us (ECEF or ECI?) and units?




    run_UKF out of ukf folder?
    create new file for simulator/visualizer?
'''

import numpy as np
import random
import pygame
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
import datetime as dt
import matplotlib.animation as animation

from pyquaternion import Quaternion

import UKF_algorithm
import gps_interface
import ideal_test_cases
import hfunc
# from happy_sensors import get_imu_data

# from ukf.PySOL import spacecraft as sp
# from ukf.PySol.spacecraft import *
from PySOL.sol_sim import *
import PySOL.spacecraft as sp
import PySOL.orb_tools as ot
# from ukf.PySOL import sol_sim
# from ukf.PySOL import orb_tools as ot

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
    ''' Rotates set of 3-vectors by the transformation described by the input quaternion

        Args:
            points (np.ndarray, (Nx3)): array of 3-vectors
            quaternion (np.ndarray, (1x4)): quaternion, using q0, q1, q2, q3 convention where q0 is the scalar component
    
    '''
    v_prime = np.zeros(points.shape)
    
    quaternion = Quaternion(quaternion)
    
    for i in np.arange(len(points)):
        v = points[i]
        v_prime[i] = quaternion.rotate(v)
    
    return v_prime

def Draw(vertices, edges):
    ''' Draw edges using defined vertices
    '''
    
    glBegin(GL_LINES)

    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])

    glEnd()

def game_visualize(states, a):
    ''' Visualization of time evolution of state using PyGame + OpenGL
    '''     

    # Initial vertices describe cube's position
    vertices_cube = np.array([[1, -1, -1], [1, 1, -1],[-1, 1, -1],[-1, -1, -1],\
        [1, -1, 1],[1, 1, 1],[-1, -1, 1],[-1, 1, 1]])
    
    # Vertices defining label z on one of the cube's face
    vertices_z_label = np.array([[-0.90, 0.90, 1], [-0.80, 0.90, 1], [-0.90, 0.80, 1], [-0.80, 0.80, 1]])
    edges_z_label = (
        (0, 1),
        (1, 2),
        (2, 3)
    )

    # Vertices defining label y on one of the cube's face
    vertices_y_label = np.array([[-0.85, 1, 0.85], [-0.70, 1, 0.85], [-0.775, 1, 0.75], [-0.775, 1, 0.70]])
    edges_y_label = (
        (0, 2),
        (1, 2),
        (2, 3)
    )

    # Vertices defining label x on one of the cube's face
    vertices_x_label = np.array([[1, -0.90, 0.90], [1, -0.80, 0.80], [1, -0.90, 0.80], [1, -0.80, 0.90]])
    edges_x_label = (
        (0, 1),
        (2, 3)
    )

    # List of tuples describing set of commands of how to draw edges. Ex: The first tuple, (0, 1), describes a command to draw a line from vertix 0 to vertix 1
    edges_cube = [
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
        ]

    # Define vertices of reaction wheels using lin space to define circle
    num = 100
    theta = np.linspace(0, 2*np.pi, num)
    vertices_RW_z = np.zeros((num, 3))
    vertices_RW_y = np.zeros((num, 3))
    vertices_RW_x = np.zeros((num, 3))

    sines = 0.5*np.cos(theta)
    cosines = 0.5*np.sin(theta)    

    vertices_RW_z[:, 0] = sines
    vertices_RW_z[:, 1] = cosines
    vertices_RW_z[:, 2] = 1.1*np.ones(num)

    vertices_RW_y[:, 0] = cosines
    vertices_RW_y[:, 1] = 1.1*np.ones(num)
    vertices_RW_y[:, 2] = sines

    vertices_RW_x[:, 0] = 1.1*np.ones(num)
    vertices_RW_x[:, 1] = sines
    vertices_RW_x[:, 2] = cosines

    # Define two rotation matrices to create fourth reaction wheel from the x-reaction wheel
    ang1 = np.pi/4
    ang2 = -np.pi/4

    c1 = np.cos(ang1)
    s1 = np.sin(ang1)

    Rz_mat = np.array([[c1, -s1, 0],
                       [s1, c1, 0],
                       [0, 0, 1]])
    
    c2 = np.cos(ang2)
    s2 = np.sin(ang2)
    
    Ry_mat = np.array([[c2, 0, s2],
                       [0, 1, 0],
                       [-s2, 0, c2]])
    
    vertices_RW_4 = np.zeros((num, 3))

    for i in np.arange(num):
        vertices_RW_4[i, :] = np.matmul(Rz_mat, np.matmul(Ry_mat, vertices_RW_x[i, :] + np.array([0.40, 0, 0])))

    # Define how to draw reaction wheel's edges
    edges_RW = []

    for i in np.arange(num-1):
        edges_RW = edges_RW + [(i, i+1)]

    edges_RW = edges_RW + [(num-1, 0)]

    # Vertices describing the lines showing the satellite's x, y, z axes / body frame
    vertices_x_axis = np.array([[0, 0, 0], [1.5, 0, 0]])
    vertices_y_axis = np.array([[0, 0, 0], [0, 1.5, 0]])
    vertices_z_axis = np.array([[0, 0, 0], [0, 0, 1.5]])

    edges_axis = (
        (0, 1),
        (0, 0)
    )

    # Initialize window, perspective, and "view" of camera
    if a == 0:
        pygame.init()
        display = (900, 700)
        pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
        gluPerspective(50, (display[0]/display[1]), 0.1, 50.0)
        glTranslatef(0.0, 0.0, -6)

    clock = pygame.time.Clock()
    Q_array = states[:, :4]  # array of quaternions (for each entry, take first four items -> array of [a,b,c,d])
    i = 0

    num_states = states.shape[0]

    # Draw the states iteratively in PyGame
    while True:
        clock.tick_busy_loop(25)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

        glEnable(GL_DEPTH_TEST)

        vertices_cube_new = rotate_points_by_quaternion(vertices_cube, Q_array[i])

        vertices_z_label_new = rotate_points_by_quaternion(vertices_z_label, Q_array[i])
        vertices_y_label_new = rotate_points_by_quaternion(vertices_y_label, Q_array[i])
        vertices_x_label_new = rotate_points_by_quaternion(vertices_x_label, Q_array[i])

        vertices_x_axis_new = rotate_points_by_quaternion(vertices_x_axis, Q_array[i])
        vertices_y_axis_new = rotate_points_by_quaternion(vertices_y_axis, Q_array[i])
        vertices_z_axis_new = rotate_points_by_quaternion(vertices_z_axis, Q_array[i])

        vertices_RW_z_new = rotate_points_by_quaternion(vertices_RW_z, Q_array[i])
        vertices_RW_y_new = rotate_points_by_quaternion(vertices_RW_y, Q_array[i])
        vertices_RW_x_new = rotate_points_by_quaternion(vertices_RW_x, Q_array[i])
        vertices_RW_4_new = rotate_points_by_quaternion(vertices_RW_4, Q_array[i])

        Draw(vertices_cube_new, edges_cube)

        Draw(vertices_z_label_new, edges_z_label)
        Draw(vertices_y_label_new, edges_y_label)
        Draw(vertices_x_label_new, edges_x_label)

        Draw(vertices_RW_z_new, edges_RW)
        Draw(vertices_RW_y_new, edges_RW)
        Draw(vertices_RW_x_new, edges_RW)
        Draw(vertices_RW_4_new, edges_RW)

        Draw(vertices_x_axis_new, edges_axis)
        Draw(vertices_y_axis_new, edges_axis)
        Draw(vertices_z_axis_new, edges_axis)

        pygame.display.flip()
        i += 1

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
    # startTime = 2022.321
    t0 = dt.datetime(2022, 3, 21, 0, 0, 0)
    # sim = sol_sim.Simulation(TIME = t0, mag_deg= 12)
    sim = Simulation(TIME = t0, mag_deg= 12)

    # how long we're simulating for
    duration = .02
    OE1 = ot.OE_array(f = 0, a = 6_800, e = 0.00068, i = 51, Om = 30, w = 30)
    sim.create_sc(OE_array= OE1, verbose = True, color = 'green', name = 'Low-Earth Orbit')


    DT = dt.timedelta(hours = duration)
    # resolution = timestep. Must match with rest of ukf
    sim.propogate(DT, resolution =  .1)
    orb_laln = sim.scs[0].state_mat.LALN
    orb_h = ot.calc_h(sim.scs[0].state_mat.R_ECEF)

    # print(sim.scs[0].state_mat.R_ECEF.shape)
    print(len(orb_laln))
    print(len(orb_h))

    # Get B field at current lat/long/altitude
    curr_date_time= np.array([2024.1066])
    lat = np.array([41.675])
    long = np.array([-86.252])
    alt = np.array([225.552]) # 740 feet (this is in meters)
    B_true = hfunc.bfield_calc(np.array([lat, long, alt, curr_date_time]))


    f = open(filename, "r")
    data = f.readline()
    splitData2 = np.array([float(x) for x in data.split(",")])
    # for test-still: accelerometer, gyro, magnetometer (microteslas)
    splitData = np.concatenate((splitData2[6:], splitData2[3:6]))

    reaction_speeds = np.zeros(3)
    i = 0
    while(data):
        # get gps data and add time stamp
        # gps_data = gps_interface.get_gps_data()
        # gps_data = gps_interface.ecef_to_latlong(gps_data[0], gps_data[1], gps_data[2]) # add time

        gps_data = np.array([np.array([orb_laln[i][0]]), np.array([orb_laln[i][1]]), np.array([orb_h[i]]), np.array([2022.257])])

        gps_data = B_true
        # run ukf and visualize output
        start, cov = UKF_algorithm.UKF(start, cov, q, r, gps_data, reaction_speeds, reaction_speeds, splitData)
        game_visualize(np.array([start[:4]]), i)

        # continue to get data from file until empty
        data = f.readline()
        if(data == ''):
            break
        splitData2 = [float(x) for x in data.split(",")]
        splitData = np.concatenate((splitData2[6:], splitData2[3:6]))

        i+=1

    f.close()


def run_ukf_sensor_iteration(state, cov, r, q, i):
    # data = get_imu_data()
    data = []
    gps_data = gps_interface.get_gps_data()
    gps_data = gps_interface.ecef_to_latlong(gps_data[0], gps_data[1], gps_data[2])
    gps_data.append(2023.8123)

    # do not use data if B-field is all zeros 
    if check_zeros(data): 
        return "Error" 

    state, cov = UKF_algorithm.UKF(state, cov, r, q, gps_data, data)
    # Visualize only if needed
    # game_visualize(np.array([state[:4]]), i) 

    return state, cov


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
    # gps_data = []
    # calibrate()

    # while(1):
    #     time.sleep(0.5)
    #     data = get_data()
        # gps_data = gps_interface.get_gps_data()
        # gps_data = gps_interface.ecef_to_latlong(gps_data[0], gps_data[1], gps_data[2])
        # gps_data.append(2023.8123)
    #     # do not use data if B-field is all zeros 
    #     if check_zeros(data): continue 

    #     start, cov = UKF_algorithm.UKF(start, cov, r, q, gps_data, data)
    #     game_visualize(np.array([start[:4]]), i)

    #     i += 1


if __name__ == "__main__":

    # Initialize noises and starting state/cov values
    n = 7
    m = n - 1

    start = np.array([1, 0, 0, 0, 0, 0, 0])

    cov = np.identity(n) * 5e-10


    # r: measurement noise (m x m)
    noiseMagnitude = 0.02
    r = np.diag([noiseMagnitude] * m)

    # q: process noise (n x n)
    noiseMagnitude = 0.005
    q = np.diag([noiseMagnitude] * n)

    
    # filename = "sensor_data_2.txt"
    filename = "test-still.txt"


    # tests ukf with pre-generated and cleaned data file
    run_ukf_textfile(start, cov, r, q, filename)
    # ideal_test_cases.run_moving_test()

    # must uncomment BNO055 imports to use in real-time with sensor
    # run_ukf_sensor(start, cov, r, q)
