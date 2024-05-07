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
from motorsandmags import mag
import random


# # GPIO
# GPIO.setmode(GPIO.BCM)
# GPIO.setup(enable,GPIO.OUT)
# GPIO.output(enable,True)
# # note: must also run GPIO.cleanup() at end of script

# I2C
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500


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
    def __init__(self, pin, direction, hallList, lastSpeed, target):
        self.pin = pin
        self.dirPin = direction
        self.hallList = hallList
        self.hData = [[],[],[]]
        self.lastSpeed = lastSpeed
        self.target = target
        self.count = 0
        self.M_switch = 0
        self.rate = 0
        GPIO.setup(self.dirPin,GPIO.OUT)
        GPIO.setup(self.hallList[0],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.hallList[1],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.hallList[2],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        GPIO.add_event_detect(self.hallList[0], GPIO.FALLING)
        GPIO.add_event_detect(self.hallList[1], GPIO.FALLING)
        GPIO.add_event_detect(self.hallList[2], GPIO.FALLING)
    
    
    # Set the speed
    def setSpeed(self):
        # self.target = newVal
        pca.channels[self.pin].duty_cycle = abs(self.target)


    # Set the direction
    def setDir(self, val):
        GPIO.output(self.dirPin, val)


    # Send output signal to actuators (& checks direction)
    def changeSpeed(self):
        temp_target = self.target
        if self.target < 0 and self.lastSpeed > 0:
            self.target = 0
            self.setSpeed()
            time.sleep(.75)

            # reverses the dirPin output, causing negative spin
            self.setDir(False)
            self.target = temp_target
            self.setSpeed()

        elif self.target > 0 and self.lastSpeed < 0:
            self.target = 0
            self.setSpeed()
            time.sleep(.75)

            # defines dirPin output = True as positive spin
            self.setDir(True)
            self.target = temp_target
            self.setSpeed()
        else:
            if abs(self.target) > .87 * MAX_DUTY:
                self.count+=1
                if self.count >= 10:
                    #self.M_switch = 1
                    mag.magOn(50,0, pca)
                    mag.magOn(50,1, pca)
            else:
                self.count = 0
                if abs(self.target) < .87 * MAX_DUTY: 
                    #self.M_switch = 0
                    mag.magOff(0, pca)
                    mag.magOff(1, pca)
            self.target -= 1 * self.rate # NOTE: M_switch doesn't exist?
            self.setSpeed()
    
	

    # convert duty cycle to RPM
    def convertToRPM(self): 
        return (self.target / MAX_DUTY) * MAX_RPM

    # convert RPM to duty cycle
    def convertToDuty(RPM):
        if RPM > MAX_RPM:
            RPM = MAX_RPM
        return (RPM / MAX_RPM) * MAX_DUTY
