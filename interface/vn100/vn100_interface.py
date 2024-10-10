'''
vn100_interface.py
Authors: Andrew Gaylord

test script for connecting to and reading data from vn100 imu

'''
from vnpy import *

s = VnSensor() 

def get_intance(): #returns the singleton instance of the VnSensor 
	return s

def connect():
	'''
	  connect(): connecting to the VnSensor version using specified COM port and printing connected to show 
	'''
	s.connect('COM5', 115200)
	print("CONNECTED")

def read_quat():#returnes a 1 by 4 array (w, x, y, z) of quaternion magnetic_acceleration_and_angular_rates
	'''
	readquat(): represnting the orientation of the object in 3d space through providing data about magentic acceleration and angular rates.

	@return 
	it returns a 1 by 4 array (w,x,y,z). w is the scalar part and x, y, z returns to vector part 
	'''
	
	quat = s.read_quaternion_magnetic_acceleration_and_angular_rates().quat
	quat = [quat.w, quat.x, quat.y, quat.z]
	return quat

def read_mag():
	mag_reading = s.read_magnetic_measurements()
	return [mag_reading.x, mag_reading.y, mag_reading.z]

def read_gyro():
	gyro_reading = s.read_angular_rate_measurements()
	return [gyro_reading.x, gyro_reading.y, gyro_reading.z]

def read_accel():
	accel_reading = s.read_acceleration_measurements()
	return [accel_reading.x, accel_reading.y, accel_reading.z]

def read_all():
	allData = s.read_magnetic_acceleration_and_angular_rates()
	mag = allData.mag
	gyro = allData.gyro
	return [mag, gyro]

def print_mag_calibration():
	mag_cal = s.read_calculated_magnetometer_calibration()
	print("b: ", mag_cal.b) # b and c are objects
	print("c: ", mag_cal.c)