import RPi.GPIO as GPIO
import numpy as np
from sklearn.linear_model import LinearRegression
import time


# Constants
k = 65535
c = 9100
hallList = [13,8,16] # pin numbers for hall sensors [x,y,z]


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
	return speed

# example call
# res = checkHall(hallList[0])
