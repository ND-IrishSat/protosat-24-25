'''
bertha.py
Authors: 

Combined script for interfacing with our BERTHA test bed.
WIP
First step, connect to VN100 IMU and read data

'''

# add to path variable so that subdirectory modules can be imported
import sys, os
sys.path.extend([f'./{name}' for name in os.listdir(".") if os.path.isdir(name)])

import os
import numpy as np

from vnpy import *

import interface.vn100.vn100_interface as vn

import ukf.simulator as simulator

# import ahrs

import time

def visualize_data(i, quaternion):
    '''
    use our cubesat simulator to show quaternion (orientation)

    @params:
        quaternion: orientation from one data read of VN100 IMU (1 x 4)
    '''

    simulator.game_visualize(np.array([quaternion]), i)


'''
read VN100 Software Documentation file, VN100 Quick Start Guide (both in VN100 References Folder in Drive)
    also open index.html from file explorer

get consistant data from VN100 IMU and generate data sets

'''
    
if __name__ == "__main__":

    # ==============================================================================
    # connect to VN100 IMU. run setup.py if needed, check and print sensor info, etc

	# declare sensor object
    # vn.connect()
    # vn.read_gps()
    ez = EzAsyncData.connect('COM5', 115200)


    # vn.print_mag_calibration()
    # ==============================================================================
    # once we're connected to IMU, set up a loop to read a stream of data

    # keep track of our iteration count
    i = 0
    while i < 10:

        print(ez.sensor.read_imu_measurements())
        print(ez.sensor.read_ins_solution_lla())

        # this returns a compositeData object
        data = ez.current_data

        print(dir(data))

        # print(data.position_gps_ecef)

        print(data.position_estimated_ecef)

        print(data.angular_rate)



    #   vn.read_ang_rates();
    #     visualize_data(i, quat)
    #     #  time.sleep(5);
        i += 1
        # optional: save to text file in form of magnetometer (magnetic field), angular velocity (gyroscope), and acceleration (accelerometer)


    # ==============================================================================
