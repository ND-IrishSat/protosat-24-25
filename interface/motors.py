from board import SCL, SDA
import busio
from adafruit_pca9685 import PCA9685
import RPi.GPIO as GPIO
import time
import numpy as np
from sklearn.linear_model import LinearRegression
import random

i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500

k = 65535
c = 9100
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
maxSpeed = 65535

#GPIO setup
GPIO.setmode(GPIO.BCM)
GPIO.setup(enable,GPIO.OUT)
GPIO.output(enable,True)


class Motor():
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

    def setSpeed(self,newVal): #sets speed
        self.target = newVal
        pca.channels[self.pin].duty_cycle = abs(self.target)

    def checkHall(self): #check motor speed
        count = 0
        data = []
        while count < 10:
            if GPIO.event_detected(self.hallList[0]):
                data.append(time.perf_counter())
                count += 1
        t = np.array(data)
        c = np.arange(10).reshape(-1,1)
        model = LinearRegression().fit(c,t)
        frequency = 1/model.coef_
        finAveSpeed += ((frequency * 15) / c) * k
        return finAveSpeed

    def setDir(self, val): #sets direction
        GPIO.output(self.dir, val)

    def checkDir(self): #checks direction
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

    # def convertToRPM(self): #converts duty cycle to RPM
    #     return (self.current / k) * maxSpeed

x = Motor(pinX,dirX,hallX,default,default)
y = Motor(pinY,dirY,hallY,default,default)
z = Motor(pinZ,dirZ,hallZ,default,default)

GPIO.cleanup()
