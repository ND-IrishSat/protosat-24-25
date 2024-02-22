import numpy as np
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
from run_UKF import *
from PySOL.sol_sim import *
from adc_pd_controller import pd_controller
import UKF_algorithm
import time
from happy_sensors import get_imu_data, calibrate
from hfunc import *
import os

MAX_PWM = 65535 # pwm val that gives max speed according to Tim


def main(target=[1,0,0,0]):

    '''NOTE
        - HALL SENSORS: method to convert Hall sensor outputs (number & time) to 4 reaction_speeds
    '''

    dim = 7
    dim_mes = dim - 1

    # Calibrate the IMU
    success = calibrate()
    if not success:
        print("Something went wrong when calibrating & configuring!!!\n")

    # Read starting state from sensors
    # start_data = get_imu_data()
    state = np.array([1, 0, 0, 0, 0, 0, 0])
    cov = np.identity(dim) * 5e-10

    # Noises = 0 for 1D test
    noise_mag = 0
    r = np.diag([noise_mag] * dim_mes)
    noise_mag = 0
    q = np.diag([noise_mag] * dim)

    # Get B field at current lat/long/altitude
    curr_date_time= np.array([2024.1066])
    lat = np.array([41.675])
    long = np.array([-86.252])
    alt = np.array([225.552]) # 740 feet (this is in meters)
    B_true = bfield_calc(np.array([lat, long, alt, curr_date_time]))

    # Infinite loop to run until you kill it
    i = 0
    while (1):
        # Starting time for loop 
        start_time = time.time()
        
        # Get reaction speeds from Hall sensors? write function to convert frequency/time of Hall sensor script into RPM!
        reaction_speeds = get_hall_data() # TODO: data will be in duty cycles (write function to convert from time & frequency to duty cycles/RPM?)

        # Get current imu data (accel*3, gyro*3, mag*3)
        imu_data = get_imu_data()
        angular_vel = imu_data[:3]

        # Data array to pass into
        data = [0] * dim_mes
        data[0] = imu_data[6]
        data[1] = imu_data[7]
        data[2] = imu_data[8]
        data[3] = angular_vel[0]
        data[4] = angular_vel[1]
        data[5] = angular_vel[2]

        # Run UKF
        # Again, same format as idael_test_cases.py
        state, cov = UKF_algorithm.UKF(state, cov, q, r, B_true, reaction_speeds, data)
        
        # Visualize if you want
        game_visualize(np.array([state[:4]]), i)
        # Also, current quaternion of what we think it is
        print("Current state: ", state[:4])

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
        # Clear screen
        os.system('cls' if os.name == 'nt' else 'clear')

        # Increment iteration count
        i += 1


if __name__ == "__main__":
    target = [0]*4
    main(target)