'''
motor_test.py
Author: Andrew Gaylord

tests that we can read sensor data from the hall sensors 
and run our motors at desired rpm
'''

from motors import *
from hall import checkHall

# Initialize motor classes (for each of 3 reactions wheels) using global variables from motors.py
x = Motor(pinX,dirX,hallX,default,default)
# y = Motor(pinY,dirY,hallY,default,default)
# z = Motor(pinZ,dirZ,hallZ,default,default)

# i2c initialization
i2c_bus = busio.I2C(SCL, SDA)
pca = PCA9685(i2c_bus)
pca.frequency = 1500

# GPIO setup
GPIO.setmode(GPIO.BCM)
GPIO.setup(enable,GPIO.OUT)
GPIO.output(enable,True)

count = 0
while(count < 100000):
    x_speed = checkHall(Motor.hallList[0]) 
    # y_speed = checkHall(Motor.hallList[1])
    # z_speed = checkHall(Motor.hallList[2])

    print("hall sensor reading: ", x_speed)

    speed = 0

    if count > 1000:
        speed = 100

    # Get the pwm signals
    x.target = 100
    # y.target = pwm(1)
    # z.target = pwm(3)

    # Check directions & alter speeds
    x.changeSpeed()
    # y.changeSpeed()
    # z.changeSpeed()

GPIO.cleanup()
