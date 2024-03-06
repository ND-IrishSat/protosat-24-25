'''
BNOO55_magnetometer_basic.py
Authors: Claudia Kuczun
Last modified: 11/5/2023

Interfaces with and calibrates BN0055 magnetometer to get data for IrishSat Kalman Filter
'''

# SPDX-FileCopyrightText: 2021 ladyada for Adafruit Industries
# SPDX-License-Identifier: MIT

import time
import board
import busio
import serial
from Adafruit_BNO055 import BNO055

def calibrate():
    # code to get data from imu :)
    i2c = board.I2C()
    sensor = BNO055.BNO055_I2C(i2c)
    if sensor.calibrated: print("CALIBRATED!")
    print("calibration status:", sensor.calibration_status)

    print("calibrating the sensors")
    while True:
        if sensor.calibration_status[1] < 2:
            print("need to calibrate gyro! hold the sensor still")
            time.sleep(3)
        elif sensor.calibration_status[3] < 2:
            print("need to calibrate magnetometer! move the sensor in figure 8 config")
            time.sleep(4)

        if sensor.calibration_status[1] >= 2 and sensor.calibration_status[3] >= 1:
            break

        print(sensor.calibration_status, "\n")


def get_data():
    i2c = board.I2C()
    sensor = BNO055.BNO055_I2C(i2c)
    res = sensor.magnetic + sensor.gyro
    print(f"res = {res}")
    return res

if __name__ == "__main__":
    calibrate()
    print(get_data())
