'''
bertha.py
Authors: 

Combined script for interfacing with our BERTHA test bed.
WIP
First step, connect to VN100 IMU and read datapip

'''

# add to path variable so that subdirectory modules can be imported
import sys, os
sys.path.extend([f'./{name}' for name in os.listdir(".") if os.path.isdir(name)])

import os
import numpy as np

from vnpy import *

import ukf.simulator as simulator
# provides functions for us to transform state space (quaternion) to measurement space
import ukf.hfunc as hfunc

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

    s = VnSensor()
print(type(s))

ez = EzAsyncData.connect('/dev/ttyUSB0', 115200)
print(s.is_connected)
print(s.port)




    # ==============================================================================
    # connect to VN100 IMU. run setup.py if needed, check and print sensor info, etc

	# declare sensor object
    s = VnSensor()
    s.connect('COM5', 115200)
    print("CONNECTED")

    # ==============================================================================
    # once we're connected to IMU, set up a loop to read a stream of data

    # keep track of our iteration count
    i = 0

    while True:
        #print(s.read_yaw_pitch_roll())
        #print(s.read_magnetic_measurements())
        #  print(s.read_quaternion_magnetic_acceleration_and_angular_rates());
        #print((s.read_quaternion_magnetic_acceleration_and_angular_rates().mag))
        #print("mag")
        #print((s.read_quaternion_magnetic_acceleration_and_angular_rates().gyro))
        #print("gyro")
        #print((s.read_quaternion_magnetic_acceleration_and_angular_rates().accel))
        #print("accel")

        quat = s.read_quaternion_magnetic_acceleration_and_angular_rates().quat
        quat = [quat.w, quat.x, quat.y, quat.z]
        #print("quat\n")

        visualize_data(i, quat)

        #  time.sleep(5);
        i += 1


        # optional: save to text file in form of magnetometer (magnetic field), angular velocity (gyroscope), and acceleration (accelerometer)


    # ==============================================================================
