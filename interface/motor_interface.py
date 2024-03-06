'''
motor_interface.py
Authors: Tim Roberts
Last modified: 3/3/24

Creates motor object that can read data from hall sensors and send signals to reaction wheels

'''

"""
from board import SCL, SDA
import busio
from adafruit_pca9685 import PCA9685
import RPi.GPIO as GPIO
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500
"""
import numpy as np
from sklearn.linear_model import LinearRegression
import random
import time

maxSpeed = 9100 #max speed in rpm
k = 65535 #max duty cycles
c = 6553.5
default = 5000
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

#GPIO setup
GPIO.setmode(GPIO.BCM)
GPIO.setup(enable,GPIO.OUT,True)


class Motor():
    def __init__(self, pin, dir, hallList, current, target):
        self.pin = pin
        self.dirPin = dir
        self.hallList = hallList
        self.hData = []
        self.current = current
        self.target = target
        #GPIO.setup(self.dirPin,GPIO.OUT,False)
        #GPIO.setup(self.hallList[0],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        #GPIO.setup(self.hallList[1],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)
        #GPIO.setup(self.hallList[2],GPIO.IN,pull_up_down=GPIO.PUD_DOWN)

    def setTarget(self,newVal):
        # update target speed
        self.target = newVal

    def setSpeed(self):
        # set reaction wheel speed speed to self.target
        if abs(self.current - self.target) >= 10000:
            for i in np.linspace(self.current,self.target,num=10):
                pca.channels[self.pin].duty_cycle = abs(i)
                # add some delay, make thread? multiprocess?
        else:
            pca.channels[self.pin].duty_cycle = abs(self.target)
            pass

    def process(self,channel):
        self.hData.append(time.perfCounter())

    def setCurrent(self): 
        #read hall sensor data, find current speed, and update self.current
        temp = []
        finAveSpeed = 0
        for i in self.hallList:
            for l in range(0,100):
                event = GPIO.add_event_detect(i,GPIO.RISING,callback=self.process())
            t = np.array(self.hData.copy())
            c = np.arange(100).reshape(-1,1)
            model = LinearRegression().fit(c,t)
            frequency = 1/model.coef_
            finAveSpeed += ((frequency * 15) / c) * k
        self.current = finAveSpeed
        self.hData = []

    def setDir(self, val): 
        # sets direction of wheel
        #GPIO.output(self.dir, val)
        pass

    def checkDir(self):
        # check if our target is in a direction that's opposite of our current direction
        # if so, bring it down to zero before switching directions
        if self.current == 0 and self.target == 0:
            pass
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

    def convertToRPM(self): 
        #converts duty cycle to RPM
        return (self.current / k) * maxSpeed

    def convertToDuty(RPM):
        # Convert RPM to duty cycle
        if RPM > maxSpeed:
            RPM = maxSpeed
        return (RPM / maxSpeed) * k


# Test class and methods
target = convert(int(input("ENTER RPM: ")))
x = Motor(10,24,[1,2,3],0,target)                                                                 
x.setSpeed()


x = Motor(pinX,dirX,hallX,default,default)
y = Motor(pinY,dirY,hallY,default,default)
z = Motor(pinZ,dirZ,hallZ,default,default)

target = random.randrange(10000,20000)
print(target)
dif = 100
temp = 10000
while dif > 10:
    while dif > 100:
        temp += 100
        x.setTarget(temp)
        x.setSpeed()
        time.sleep(0.5)
        x.setCurrent()
        dif = abs(target - x.current)
    temp += 10
    x.setTarget(temp)
    x.setSpeed()
    time.sleep(0.5)
    x.setCurrent()
    dif = abs(target - x.current)
#GPIO.cleanup()
