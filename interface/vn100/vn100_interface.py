'''
vn100_interface.py
Authors: Andrew Gaylord

test script for connecting to and reading data from vn100 imu

'''
from vnpy import *

s = VnSensor()

def get_intance():
	return s

def connect():
	s.connect('COM5', 115200)
	if(s.verify_sensor_connectivity()):
		print("CONNECTED")

def disconnect():
	s.disconnect()
	print("DISCONNECTED")

def read_quat():
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

def read_gps():
	#gps_reading = s.read_gps_solution_lla()
	gps_reading = s.read_gps_solution_ecef()
	print("GPS: ", gps_reading);