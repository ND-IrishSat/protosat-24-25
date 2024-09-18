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

# from interface.vn100.vn100_library.vnpy import *
#import Serial

import ukf.simulator as simulator
# provides functions for us to transform state space (quaternion) to measurement space
import ukf.hfunc as hfunc

# TODO: read this explanation, VN100 Software Documentation file, and VN100 Quick Start Guide (both in VN100 References Folder in Drive)
# explaination from same company as VN100:
    # https://www.vectornav.com/resources/inertial-navigation-primer/theory-of-operation/theory-ahrs

# TODO: pottentially research ahrs attitude estimation library
# https://ahrs.readthedocs.io/en/latest/
# provides methods to estimate orientation (quaternion) from IMU data
    # https://ahrs.readthedocs.io/en/latest/filters.html
	# look at Madwick and Algebraic Quaternion Algorithm (AQUA)


def visualize_data(data):
    '''
    use our cubesat simulator to show n sets of magnetometer data
    need to convert to quaternion first so our simulator can use it

    @params:
        data: magnetometer data from VN100 IMU (n x 3)
    '''
    i = 0
    while (i < len(data)):
        # convert to quaternion
        quaternion = [0] + data[i]
        simulator.game_visualize(np.array([quaternion]), i)

        i += 1

if __name__ == "__main__":
	
	# declare sensor object
	# s = VnSensor()
	# print(type(s))
	
	# attempt connection: baud rate is bugged
	# ez = EzAsyncData.connect('/dev/ttyUSB0', 115200)
	
	# alternate connection
	#s.connect('/dev/ttyUSB0', 9600)
	
	# info from sensor object
	# print(s.is_connected)
	# print(s.port)
	
	
	# s.change_baudrate(9600)
	#print(s.read_serial_number())
	
	#print(s.port)
	# print(s.read_model_number())


    # once you can populate a data array with magnetometer data, uncomment this
    visualize_data(data)





# example of reading from a data text file

# filename = "interface/sample_data/test-still.txt"
# # find absolute path to text file
# script_path = os.path.abspath(__file__)
# script_dir = os.path.split(script_path)[0]
# abs_file_path = os.path.join(script_dir, filename)
# f = open(abs_file_path, "r")

# data = f.readline()
# splitData2 = np.array([float(x) for x in data.split(",")])
# # for test-still: accelerometer, gyro, magnetometer (microteslas)
# # data is in the form of magnetic field (bx, by, zy) and angular velocity (wx, wy, wz)
# splitData = np.concatenate((splitData2[6:], splitData2[3:6]))
