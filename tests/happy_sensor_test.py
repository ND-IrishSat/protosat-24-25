

import time
from mpu9250_jmdev.registers import *
from mpu9250_jmdev.mpu_9250 import MPU9250


def sensor_test():
    mpu = MPU9250(
        address_ak = AK8963_ADDRESS,
        address_mpu_master = MPU9050_ADDRESS_68,
        address_mpu_slave = None,
        bus = 1,
        gfs = GFS_1000,
        afs = AFS_8G,
        mfs = AK8963_BIT_16,
        mode = AK8963_MODE_C100HZ)
    time.sleep(5)
    print("Calibrating gyroscope and accelerometer...")
    # mpu.calibrate()
    # mpu.configure()

    # print("Calibrating magnetometer...")
    # mpu.calibrateAK8963()
    # mpu.configure() # always configure after calibration

    # while True:
        # print("\nMPU 9250 running...")
        # print("accelerometer: ", mpu.readAccelerometerMaster())
        # print("gyroscope: ", mpu.readGyroscopeMaster())
        # print("magnetometer: ", mpu.readMagnetometerMaster())

        # time.sleep(1)

    # TODO
    # 1. create test files named with certain movements (holding still, turning slowly around each axis)    #    csv files
    # 2. output array of 6 values => [magnetometer, gyroscope]

    count = 0 
    line = ['Gyroscope']
    line2 = ['Magnetometer']

    while count < 2000:
        line1, line2 = [],[]
        f = open("testfile-still.csv", "w")
        print(mpu.readGyroscopeMaster())
        print(mpu.readGyroscopeMaster())
        print("appedning data now")
        line.append(mpu.readGyroscopeMaster())
        #f.write(str(line) +"\n")
        line.append(mpu.readMagnetometerMaster())
        print(line)
        #line2.append(mpu.readMagnetometerMaster())
        #f.write(str(line2))
        count += 1
    f.close()
    
if __name__ == "__main__":
    sensor_test()

