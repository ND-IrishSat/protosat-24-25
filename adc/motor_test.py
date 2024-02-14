# Imports
"""
from board import SCL, SDA
import busio
from adafruit_pca9685 import PCA9685
import RPi.GPIO as GPIO
import time
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500
"""
import numpy as np
import time


# Constants
maxSpeed = 8100     # Max speed in rpm
currentSpeed = 2500 # Current speed in rpm
k = 65535
dirX = 0
dirY = 0
dirZ = 0


# Motor class for actuators
class Motor():
    # Create class instance
    def __init__(self, pin, dir, hallList, current, target):
        self.pin = pin
        self.dirPin = dir
        self.hallList = hallList
        self.current = current
        self.target = target
        #GPIO.setup(23,GPIO.IN)

    # Check speed
    def setSpeed(self, tar):
        if abs(self.current - tar) >= 10000:
            for i in np.linspace(self.current,tar,num=10):
                #pca.channels[self.pin].duty_cycle = abs(i)
                #add some delay, make thread? multiprocess?
                print(i)
        else:
            #pca.channels[self.pin].duty_cycle = abs(tar)
            pass

    # Set direction 
    def setDir(self, val):
        GPIO.output(self.dir, val)

    # Check reaction wheel direction for change
    def checkDir(self):
        # No change needed
        if self.current == 0 and self.target == 0:
            pass
        # Change needed
        elif self.current != self.target:
            # Switch negative to positive direction (T -> F) speed
            if self.target < 0 and self.current > 0:
                self.setSpeed(0)
                self.setDir(False)
                self.setSpeed(self.target())
                # Switch positive to negative direction (F -> T) speed
            elif self.target > 0 and self.current < 0:
                self.setSpeed(0)
                self.setDir(True)
                self.setSpeed(self.target())
        # Same speed
        else:
            self.setSpeed(self.target())

    # Read hall sensor data
    def hallRead(self):
        pass


# Convert RPM to duty cycle
def convert(RPM):
    if RPM > maxSpeed:
        RPM = maxSpeed
    return (RPM / maxSpeed) * k


# Test class and methods
target = convert(int(input("ENTER RPM: ")))
x = Motor(10,24,[1,2,3],0,target)                                                                 
x.setSpeed()
