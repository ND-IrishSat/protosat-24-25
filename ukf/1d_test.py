import numpy as np
import random
import pygame
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
import datetime as dt
from pyquaternion import Quaternion
from run_UKF import *
from PySOL.sol_sim import *
import PySOL.spacecraft as sp
import PySOL.orb_tools as ot
from adc.adc_pd_controller import pd_controller
import UKF_algorithm
import time
from happy_sensors import get_imu_data
import hfunc

MAX_PWM = 65535 # pwm val that gives max speed according to Tim


def main(target=[1,0,0,0]):

    '''NOTE
        Things that are left to do:
        - HALL SENSORS: method to convert Hall sensor outputs (number & time) to 4 reaction_speeds
        - VISUALIZATION: inputs to view it correctly
    '''

    dim = 7
    dim_mes = dim - 1
    state = 0
    cov = 0

    # TODO: Check if this noise is correcy for first test
    #       Based on ideal_test_cases.py format/setup
    # Maybe make them 0 for test?
    noise_mag = 0.02
    r = np.diag([noise_mag] * dim_mes)
    noise_mag = 0.005
    q = np.diag([noise_mag] * dim)
    u_k = 0

    # True magnetic field (constant in cage)
    B_true = np.array([0, 0, 1])

    # Magnetic field from sensors
    B_sens = np.array([np.matmul(hfunc.quaternion_rotation_matrix(state[0]), B_true)])
    for a in range(1, n): # TODO: get row of data from sensors!!
        B_sens = np.append(B_sens, np.array([np.matmul(hfunc.quaternion_rotation_matrix(state[a]), B_true)]), axis=0)

    
    # Infinite loop to run until you kill it
    while (1):
        # Starting time for loop 
        start_time = time.time()
        
        # TODO: Get reaction speeds from Hall sensors?
        #       Write function to convert frequency/time of Hall sensor script into RPM!
        #       (In whatever units needed for UKF)
        reaction_speeds = []
        reaction_speeds = get_hall_data() # NOTE: data will be in duty cycles (write function to convert from time & frequency to duty cycles/RPM?)

        # Get current imu data (accel*3, gyro*3, mag*3)
        imu_data = get_imu_data()
        angular_vel = imu_data[6:]

        # Data array to pass into
        data = [0] * dim_mes
        data[0] = B_sens[0]
        data[1] = B_sens[1]
        data[2] = B_sens[2]
        data[3] = angular_vel[0]
        data[4] = angular_vel[1]
        data[5] = angular_vel[2]

        # Run UKF
        # Again, same format as idael_test_cases.py
        state, cov = UKF_algorithm.UKF(state, cov, r, q, list(B_true), reaction_speeds, data)
        
        # Visualize if you want
        # game_visualize(np.array([starst[:4]]), i)
        
        # Run PD controller
        curr_state = state
        target = [1,0,0,0]

        # Sensor func to get currentangular velocity of cubesat
        omega = np.array([angular_vel[0], angular_vel[1], angular_vel[2]])
        kp = .05*MAX_PWM
        kd = .01*MAX_PWM
        
        # Run PD controller to generate output for reaction wheels
        pwm = pd_controller(curr_state, target, omega, kp, kd)

        # Run actuator (reaction wheels) script
        # TODO: Send commands to actuators

        # Ending time for loop
        end_time = time.time()

        # Output total runtime for single loop iteration
        print("Time to run single iteration of ADC loop: {:.2f}".format(end_time-start_time))


if __name__ == "__main__":
    target = [0]*4
    main(target)