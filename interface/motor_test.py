'''
motor_test.py
Author: Andrew Gaylord

tests that we can read sensor data from the hall sensors 
and run our motors at desired rpm
'''

from motors import *
from hall import checkHall
from init import initialize_setup


# Initialize for execution
initialize_setup()

# Initialize motor classes (for each of 3 reactions wheels) using global variables from motors.py
x = Motor(pinX,dirX,hallX,default,default)
# y = Motor(pinY,dirY,hallY,default,default)
# z = Motor(pinZ,dirZ,hallZ,default,default)

count = 0

while (count < 100):    
    # Intialize the speed to 10000
    speed = 10000
    
    # Get the pwm signals
    x.target = speed
    # y.target = pwm(1)
    # z.target = pwm(3)

    # Alter speed
    x.setSpeed(x.target)
    # y.setSpeed(y.target)
    # z.setSpeed(z.target)

    # Check the current speed
    x_speed = checkHall(x.hallList[0]) 
    # y_speed = checkHall(y.hallList[1])
    # z_speed = checkHall(z.hallList[2])

    print(f"hall sensor reading: {x_speed}")
    count += 1


# Bring the wheels to a stop
speed = 0
x.target = speed
x.setSpeed(x.target)
print("done with script! ending...")

GPIO.cleanup()
