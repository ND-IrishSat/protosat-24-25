'''
motors.py
Authors: Tim Roberts

motor interface class 
checks hall sensors readings and sets reaction wheel speeds
initialization must be done through init.py

'''
from board import SCL, SDA
import busio
from adafruit_pca9685 import PCA9685
import RPi.GPIO as GPIO
import time
import numpy as np
from sklearn.linear_model import LinearRegression
from hall import checkHall
import random

# Constants
MAX_RPM = 9100
# max duty cycles
MAX_DUTY = 65535
default = 0
enable = 11
pinX = 10
pinY = 9
pinZ = 8
dirX = 0
dirY = 5
dirZ = 6
hallX = [13,19,26]
hallY = [8,7,1]
hallZ = [16,20,21]

# motor class!
class Motor():
    # IO setup
    def __init__(self, pin, dir, hallList, current, target):
        self.pin = pin
        self.dirPin = dir
        self.hallList = hallList
        self.hData = [[],[],[]]
        self.current = current
        self.target = target
        GPIO.setup(self.dirPin,GPIO.OUT)
        GPIO.setup(self.hallList[0],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.hallList[1],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.hallList[2],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.add_event_detect(self.hallList[0], GPIO.FALLING)
        GPIO.add_event_detect(self.hallList[1], GPIO.FALLING)
        GPIO.add_event_detect(self.hallList[2], GPIO.FALLING)


    # Check Hall sensor (motor) speed
    def checkSpeed(self): 
        # use hall.py to check speed (duty cycles) for first hall sensor
        return checkHall(self.hallList[0])
        # count = 0
        # data = []
        # while count < 10:
        #     if GPIO.event_detected(self.hallList[0]):
        #         data.append(time.perf_counter())
        #         count += 1
        # t = np.array(data)
        # c = np.arange(10).reshape(-1,1)
        # model = LinearRegression().fit(c,t)
        # frequency = 1/model.coef_
        # speed += ((frequency * 15) / c) * MAX_DUTY
        # # Return the speed according to Hall sensor (duty cycles)
        # return speed
    

    # Set the speed
    def setSpeed(self, newVal):
        self.target = newVal
        pca.channels[self.pin].duty_cycle = abs(self.target)


    # Set the direction
    def setDir(self, val):
        GPIO.output(self.dir, val)


    def checkDir(self):
        # check if our target is in a direction that's opposite of our current direction
        # if so, bring it down to zero before switching directions
        if self.current != self.target:
            if self.target < 0 and self.current > 0:
                self.setSpeed(0)
                self.setDir(False)
                self.setSpeed(self.target())
            elif self.target > 0 and self.current < 0:
                self.setSpeed(0)
                self.setDir(True)
                self.setSpeed(self.target())
        else:
            self.setSpeed(self.target())


    # convert duty cycle to RPM
    def convertToRPM(self): 
        return (self.current / MAX_DUTY) * MAX_RPM

    # convert RPM to duty cycle
    def convertToDuty(RPM):
        if RPM > MAX_RPM:
            RPM = MAX_RPM
        return (RPM / MAX_RPM) * MAX_DUTY