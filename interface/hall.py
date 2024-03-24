'''
hall.py
Authors: Tim Roberts

interfaces with hall sensors, allowing us to read motor speed of reaction wheels

'''

import RPi.GPIO as GPIO
import numpy as np
from sklearn.linear_model import LinearRegression
import time

# Constants
c = 9100
# max duty cycles
k = 65535
# pin numbers for hall sensors [x,y,z]
hallList = [13,8,16]


# Check motor speed (in duty cycles)
def checkHall(sensor): 
	count = 0
	data = []
	while count < 10:
		if GPIO.event_detected(sensor):
			data.append(time.perf_counter())
			count += 1
	t = np.array(data)
	c = np.arange(10).reshape(-1,1)
	model = LinearRegression().fit(c,t)
	frequency = 1/model.coef_
	speed += ((frequency * 15) / c) * k 
	# Return the speed according to Hall sensor (duty cycles)
	return speed

# example call
# res = checkHall(hallList[0])
