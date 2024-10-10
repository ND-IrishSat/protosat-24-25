'''
vn100_interface.py
Authors: Andrew Gaylord

test script for connecting to and reading data from vn100 imu

Script for changing connectivity to VN100 IMU sensor, and reading and/or getting data from the sensor.
Data such as quaternion, magnetic, angle, acceleration, and temperature.


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


def get_intance(): #returns the singleton instance of the VnSensor 
	return s

def connect():
	'''
	  connect(): connecting to the VnSensor version using specified COM port and printing connected to show 
	'''
	s.connect('COM5', 115200)
	if(s.verify_sensor_connectivity()):
		print("CONNECTED")

def disconnect():
	s.disconnect()
	print("DISCONNECTED")

def read_quat():#returnes a 1 by 4 array (w, x, y, z) of quaternion magnetic_acceleration_and_angular_rates
	'''
	readquat(): represnting the orientation of the object in 3d space through providing data about magentic acceleration and angular rates.

	@return 
	it returns a 1 by 4 array (w,x,y,z). w is the scalar part and x, y, z returns to vector part 
	'''
	
	quat = s.read_attitude_quaternion()
	quat = [quat.w, quat.x, quat.y, quat.z]
	return quat

def read_mag():
	'''
	readmag(): represnting the current magnetic fields reading from the sensor as three values corresponding to x,y,z dimensions 

	@return 
	it returns a truple containing 3 float values in microtesla (x,y,z). x, y, z returns to vector part. 
	'''
	mag_reading = s.read_magnetic_measurements()
	return [mag_reading.x, mag_reading.y, mag_reading.z]

def read_gyro():
	'''
	read_gyro(): detects the devaiation of the object and is used to measure angular rates from the sensor as x, y, z three values. 

	@return 
	it returns a 1 by 3 array (x, y, z), reading the angular rate in three different dimensions. 
	'''
	gyro_reading = s.read_angular_rate_measurements()
	return [gyro_reading.x, gyro_reading.y, gyro_reading.z]
	

def read_accel():
	'''
	read_accel(): detecs the acceleration of the object and displays 3 values 

	@return 
	it returns a 1 by 3 array (x,y,z), reading the acceleration in three different dimensions. 
	'''
	accel_reading = s.read_acceleration_measurements()
	return [accel_reading.x, accel_reading.y, accel_reading.z]

def read_all():
	'''
	Use
		A get function for magnetic acceleration and angular orentation
	@return
		Returns magnetic acceleration and angle orentation as 1 by 6 array, first the x,y,z of magnetic then x,y,z of angular
	'''
	allData = s.read_magnetic_acceleration_and_angular_rates()
	mag = allData.mag
	gyro = allData.gyro
	return [mag, gyro]

def print_mag_gyro_quat():
	# returns string with mag, gyro, and quat data (in that order)
	# to use when creating sample data txt files
	mag = [ str(a) for a in read_mag()]
	gyro = [ str(b) for b in read_gyro()]
	quat = [ str(c) for c in read_quat()]
	all_three = mag + gyro + quat # create list with all data
	all_three_string = ", ".join(all_three)
	return all_three_string

def print_data_to_file(count, file_name):
	i = 0 # keep track of our iteration count
	f = open(file_name, "a+")
	while i < count:
		i += 1
		# save to text file in form of magnetometer (magnetic field), angular velocity (gyroscope), and acceleration (accelerometer)
		f.write(print_mag_gyro_quat()) # put mag, gyro, quat data into text file
		if (i < count): f.write("\n") # add newline to separate data sets

	f.close()
	#source = f'./{file_name}'
	#destination = './new_sensor_tests'
	#os.rename(source, destination)
	return

def print_mag_calibration():
	'''
	Note: Unfinished
	
	Use
		A method that prints out the calculated calibration for both B and C
	'''
	mag_cal = s.read_calculated_magnetometer_calibration()
	print("b: ", mag_cal.b) # b and c are objects
	print("c: ", mag_cal.c)

#Returns temp from the sensor as 
def read_temp():
	'''
	read_temp(): measures the temperature in Celcius 
	@return
	returs the temperature reading of the sensor as a float 
	'''
	return s.read_imu_measurements().temp

#Returns the pressure reading of the sensor as a float 
def read_pressure():
	'''
	read_pressure(): measures the pressure in hectopascals (hPa)
	@return
	returs the pressure reading of the sensor as a float 
	'''
	return s.read_imu_measurements().pressure