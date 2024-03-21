'''
1d_test.py
Authors: Claudia Kuczun, Andrew Gaylord, Patrick Schwartz, Juwan Jeremy Jacob, Alex Casillas
Last updated: 2/25/24

Main script for 1D systems integration test
Collab test between GOAT lab and protoSat
Implements UKF, PD, hall sensors, IMU, visualizer, and actuators to rotate cubeSat on frictionless table on 1 axis 

'''
import numpy as np
from pygame.locals import * 
from OpenGL.GL import *
from OpenGL.GLU import *
from ukf.PySOL.wmm import *
from ukf.simulator import *
# from adc.adc_pd_controller import pd_controller_numpy
from ukf.UKF_algorithm import *
import time
from interface.happy_sensors import get_imu_data, calibrate
# from interface.motor_interface import Motor
from ukf.hfunc import *

MAX_PWM = 65535 # pwm val that gives max speed according to Tim
MAX_RPM = 8100 # according to Tim


def main(target=[1,0,0,0]):

    '''NOTE
        - HALL SENSORS: need method to convert Hall sensor outputs (number & time) to 4 reaction_speeds
        - things to take out for UKF only test: motor initialization, reaction wheel update, pd section (127-146)
        - What is starting state guess?
        - Set proper target quaternion?
        - check with patrick that we want pwm[3] instead of pwm[2]
        - What is max RPM (8100 or 8800?)
    '''

    dim = 7
    dim_mes = dim - 1

    # Calibrate the IMU
    success = calibrate()
    if not success:
        print("Something went wrong when calibrating & configuring!!!\n")

    # Read starting state from sensors
    # start_data = get_imu_data()
    state = np.array([1, 0, 0, 0, 0, 0, 0]) #[q0, q1, q2, q3, omega_x, omega_y, omega_z]
    cov = np.identity(dim) * 5e-10

    # Noises = 0 for 1D test. when we test with hall sensors, ask micheal about which values to test
    noise_mag = 0.85
    #r = np.diag([noise_mag] * dim_mes)
    mag = 0.1 
    omg = [1,3,1]
    
    r = [[mag,0,0,0,0,0],
	 [0,mag,0,0,0,0],
	 [0,0,mag,0,0,0],
	 [0,0,0,omg[0],0,0],
	 [0,0,0,0,omg[1],0],
	 [0,0,0,0,0,omg[2]]]
    
    q = [[1e-1,0,0,0,0,0,0],
	 [0,1e-1,0,0,0,0,0],
	 [0,0,1e-1,0,0,0,0],
	 [0,0,0,1e-1,0,0,0],
	 [0,0,0,0,5e4,0,0],
	 [0,0,0,0,0,5e4,0],
	 [0,0,0,0,0,0,5e4]]
    
    # current lat/long/altitude, which doesn't change for this test
    curr_date_time= np.array([2024.1266])
    lat = np.array([41.675])
    long = np.array([-86.252])
    alt = np.array([.225552]) # 740 feet (this is in km!!)
    # calculate true B field at stenson-remick
    B_true = bfield_calc(np.array([lat, long, alt, curr_date_time]))

    # set negative of last element to match magnetometer
    # Convert to North East Down to North East Up, matching X, Y, Z reference frame of magnetometer
    B_true[2] *= -1

    # inialize current step and last step reaction wheel speeds
    # for this test they're 1x3: x, y, skew wheels
    old_reaction_speeds = np.array([0,0,0])
    reaction_speeds = np.array([0,0,0])

    # initialize pwm speeds
    pwm = np.array([0,0,0,0])

    # Infinite loop to run until you kill it
    i = 0
    while (1):
        # Starting time for loop 
        start_time = time.time()

        # set pins to control reaction wheel motors
        # xMotorPin = 10
        # yMotorPin = 11
        # skewMotorPin = 12

        # # Initialize Motor class from Tim's script to read and write to motors
        # x = Motor(xMotorPin,24,[1,2,3],0,0)
        # y = Motor(yMotorPin,24,[1,2,3],0,0)
        # skew = Motor(skewMotorPin,24,[1,2,3],0,0) # pwm(2) is for z
        
        #old_reaction_speeds = reaction_speeds
        old_reaction_speeds = np.array([0, 0, 0])
        # Get reaction speeds from Hall sensors? write function to convert frequency/time of Hall sensor script into RPM!
        # TODO: data will be in duty cycles (write function to convert from time & frequency to duty cycles/RPM?)
        #reaction_speeds = get_hall_data() 

        # Temporary solution instead of function: get reaction wheel speeds from pwm input
        # we basically know how fast they'll be spinning based on what we told it last cycle
        # Scale from pwm to rad/s
        # scale = MAX_RPM/MAX_PWM*2*np.pi/60
        # TODO: is this pwm[3] or pwm[2]???
        # reaction_speeds = scale*np.array([pwm(0),pwm(1),pwm(3)])

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
        state, cov = UKF(state, cov, q, r, B_true, reaction_speeds, old_reaction_speeds, data)
        print('state: ',state)
        print('')
        
        # Visualize if you want
        game_visualize(np.array([state[:4]]), i) # not working on Pi 0, but working on  Pi 4
        # print current quaternion estimate
        # print("Current state: ", state[:4])

        # Run PD controller
        # curr_state = state
        # target = [1,0,0,0]

        # # Sensor func to get currentangular velocity of cubesat
        # omega = np.array([angular_vel[0], angular_vel[1], angular_vel[2]])
        # # PD gains parameters
        # kp = .05*MAX_PWM
        # kd = .01*MAX_PWM
        
        # # Run PD controller to generate output for reaction wheels
        # pwm = pd_controller(curr_state, target, omega, kp, kd)

        # # Set speed of each of the reaction wheels
        # x.target = pwm(0)
        # y.target = pwm(1)
        # skew.target = pwm(3)
        # x.setSpeed(x.target)
        # y.setSpeed(y.target)
        # skew.setSpeed(skew.target)

        # Ending time for loop
        end_time = time.time()

        # Output total runtime for single loop iteration
        #print("Time to run single iteration of ADC loop: {:.2f}".format(end_time-start_time))
        # Clear screen
        #os.system('cls' if os.name == 'nt' else 'clear')

        # Increment iteration count
        i += 1


if __name__ == "__main__":
    target = [0]*4
    main(target)
