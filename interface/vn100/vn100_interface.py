'''
vn100_interface.py
Authors: Andrew Gaylord

test script for connecting to and reading data from vn100 imu

'''
from vnpy import *

s = VnSensor()

'''
Alternate connection route:

ez = EzAsyncData.connect('COM5', 115200)

# call the sensor object through the object and read the relevant data
ez.sensor.read_imu_measurements()

# this returns a compositeData object
# lots of helper functions exist to extract the data from this object
data = ez.current_data

print("angular: ", data.angular_rate)

print("quat: ", data.any_attitude.quat)

'''

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
	quat = s.read_attitude_quaternion()
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

def read_temp():
	# read_imu_measurements returns mag, accel, gyro, temp, and pressure
	return s.read_imu_measurements().temp

def read_pressure():
	return s.read_imu_measurements().pressure