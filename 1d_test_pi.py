'''
1d_test.py
Authors: Claudia Kuczun, Andrew Gaylord, Patrick Schwartz, Juwan Jeremy Jacob, Alex Casillas
Last updated: 2/25/24

Main script for 1D systems integration test
Collab test between GOAT lab and protoSat
Implements UKF, PD, hall sensors, IMU, visualizer, and actuators to rotate cubeSat on frictionless table on 1 axis 

'''

import sys, os
sys.path.extend([f'./{name}' for name in os.listdir(".") if os.path.isdir(name)])

import numpy as np
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
from ukf.PySOL.wmm import *
from ukf.simulator import *
from adc.adc_pd_controller_numpy import pd_controller
from ukf.UKF_algorithm import *
import time
from interface.happy_sensors import get_imu_data, calibrate
from ukf.hfunc import *
from interface.motors import *
from interface.hall import checkHall
from interface.init import initialize_setup

# MAX_DUTY and MAX_RPM are declared in motors.py


def main(target=[1,0,0,0]):
    '''NOTE
        - things to take out for UKF only test: motor initialization, reaction wheel update, pd section (127-146)
        - What is starting state guess? read from sensors?
        - Set proper target quaternion?
        - check with patrick that we want pwm[3] instead of pwm[2]
        - check that current_quaternion is correct (not entire state)
        - do we need to worry about hall sensors not reading if wheels not spinning?
    '''

    # Initialize setup for motors (I2C and GPIO)
    initialize_setup()
    print("initialized setup\n")

    # Initialize motor classes (for each of 3 reactions wheels) using global variables from motors.py
    x = Motor(pinX,dirX,hallX,default,default)
    y = Motor(pinY,dirY,hallY,default,default)
    z = Motor(pinZ,dirZ,hallZ,default,default)
    print("initialized motors\n")

    # Dimension of state and measurement space
    dim = 7
    dim_mes = dim - 1

    dt = 0.1

    # Calibrate the IMU
    success = calibrate()
    if not success:
        print("Something went wrong when calibrating & configuring!!!\n")

    # Read starting state from sensors?
    # start_data = get_imu_data()
    
    # initialize starting state and covariance
    state = np.array([1, 0, 0, 0, 0, 0, 0]) #[q0, q1, q2, q3, omega_x, omega_y, omega_z]
    cov = np.identity(dim) * 5e-10

    # r: measurement noise (m x m)
    noise_mag = .1
    # r = np.diag([noise_mag] * dim_mes)
    r = np.array([noise_mag, 0, 0, 0, 0, 0],
                 [0, noise_mag, 0, 0, 0, 0],
                 [0, 0, noise_mag, 0, 0, 0],
                 [0, 0, 0, noise_mag, 0, 0],
                 [0, 0, 0, 0, noise_mag, 0],
                 [0, 0, 0, 0, 0, noise_mag])

    # q: process noise (n x n)
    # Should depend on dt
    noise_mag = .5
    # q = np.diag([noise_mag] * dim)
    q = np.array([[dt, 3*dt/4, dt/2, dt/4, 0, 0, 0],
                [3*dt/4, dt, 3*dt/4, dt/2, 0, 0, 0],
                [dt/2, 3*dt/4, dt, 3*dt/4, 0, 0, 0],
                [dt/4, dt/2, 3*dt/4, dt, 0, 0, 0],
                [0, 0, 0, 0, dt, 2*dt/3, dt/3],
                [0, 0, 0, 0, 2*dt/3, dt, 2*dt/3],
                [0, 0, 0, 0, dt/3, 2*dt/3, dt]
    ])
    q = q * noise_mag
    
    # Current lat/long/altitude, which doesn't change for this test
    curr_date_time= np.array([2024.1266])
    lat = np.array([41.675])
    long = np.array([-86.252])
    alt = np.array([.225552]) # 740 feet (km)

    # calculate true B field at stenson-remick
    B_true = bfield_calc(np.array([lat, long, alt, curr_date_time]))

    # Set negative of last element to match magnetometer
    # Convert to North East Down to North East Up, matching X, Y, Z reference frame of magnetometer
    B_true[2] *= -1

    # Initialize current step and last step reaction wheel speeds
    # For this test they're 1x3: x, y, skew (z) wheels
    old_reaction_speeds = np.array([0,0,0])
    reaction_speeds = np.array([0,0,0])

    # Initialize pwm speeds
    pwm = np.array([0,0,0,0])

    # Infinite loop to run until you kill it
    i = 0
    while (1):
        # Starting time for loop 
        start_time = time.time()

        # Store reaction wheel speeds of last iteration
        old_reaction_speeds = reaction_speeds

        # Get reaction speeds (in duty cycles) from Hall sensors
        #NOTE: wheel must be spinning to check Hall
        x_speed = checkHall(x.hallList[0]) 
        y_speed = checkHall(y.hallList[1])
        z_speed = checkHall(z.hallList[2])
        reaction_speeds = [*x_speed, *y_speed, *z_speed]

        # Get current imu data (mag*3, gyro*3)
        imu_data = get_imu_data()

        # Angular velocity comes from gyro
        angular_vel = imu_data[3:]

        # Data array to pass into
        data = [0] * dim_mes
        data[0] = imu_data[0]
        data[1] = imu_data[1]
        data[2] = imu_data[2]
        data[3] = angular_vel[0]
        data[4] = angular_vel[1]
        data[5] = angular_vel[2]

        # Run UKF
        state, cov = UKF(state, cov, q, r, B_true, reaction_speeds, old_reaction_speeds, data)
        
        # Visualize if you want
        game_visualize(np.array([state[:4]]), i) # not working on Pi 0, but working on  Pi 4
        # print current quaternion estimate
        print("Current state: ", state[:4])

        # Run PD controller
        curr_quaternion = state[:4]
        target = [1,0,0,0]

        # get current angular velocity of cubesat
        omega = np.array([angular_vel[0], angular_vel[1], angular_vel[2]])
        # PD gains parameters (dependant upon max pwm/duty cycles)
        kp = .05*MAX_DUTY
        kd = .01*MAX_DUTY
        
        # Run PD controller to generate output for reaction wheels
        pwm = pd_controller(curr_quaternion, target, omega, kp, kd)

        # Get the pwm signals
        x.target = pwm(0)
        y.target = pwm(1)
        z.target = pwm(3)

        # Check directions & alter speeds
        x.setSpeed(x.target)
        y.setSpeed(y.target)
        z.setSpeed(z.target)

        # Ending time for loop
        end_time = time.time()

        # Output total runtime for single loop iteration
        #print("Time to run single iteration of ADC loop: {:.2f}".format(end_time-start_time))
        # Clear screen
        #os.system('cls' if os.name == 'nt' else 'clear')

        # Increment iteration count
        i += 1

    # Bring the wheels to a stop
    speed = 0
    x.target = speed
    x.setSpeed(x.target)
    y.target = speed
    y.setSpeed(y.target)
    z.target = speed
    z.setSpeed(z.target)

    # Confirm it is done
    print("done with script! ending...")

    GPIO.cleanup()


if __name__ == "__main__":
    target = [1, 0, 0, 0]
    main(target)
