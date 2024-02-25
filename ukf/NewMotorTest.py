from board import SCL, SDA
import busio
from adafruit_pca9685 import PCA9685
import RPi.GPIO as GPIO
import time
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500
import numpy as np
import time

MAX_SPEED = 8100 # max speed in rpm
CURRENT_SPEED = 2500 # current speed in rpm
K = 65535 # max duty cycle

class Motor():
    def __init__(self, pin, dir, hallList, current, target):
        self.pin = pin
        self.dirPin = dir
        self.hallList = hallList
        self.current = current
        self.target = target
        #GPIO.setup(23,GPIO.IN)

    def setSpeed(self, target=None): #checks speed
        if not target:
            target = self.target
        if abs(self.current - target) >= 10000:
            for i in np.linspace(self.current,target,num=10):
                pca.channels[self.pin].duty_cycle = abs(i)
                # add some delay, make thread? multiprocess?
        else:
            pca.channels[self.pin].duty_cycle = abs(target)
            pass

    def setDir(self, val): #sets direction
        GPIO.output(self.dir, val)

    def checkDir(self):
        if self.current == 0 and self.target == 0:
            return
        elif self.current != self.target:
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

    def hallRead(self): #read hall data
        pass

def convert(RPM): #converts RPM to duty cycle
    if RPM > MAX_SPEED:
        RPM = MAX_SPEED
    return (RPM / MAX_SPEED) * K

# target = convert(int(input("ENTER RPM: ")))
# x = Motor(10,24,[1,2,3],0,target)
# x.setSpeed()
# dirX = 0
# dirY = 0
# dirZ = 0
