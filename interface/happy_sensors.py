import time
from mpu9250_jmdev.registers import *
from mpu9250_jmdev.mpu_9250 import MPU9250

mpu = MPU9250(
        address_ak = AK8963_ADDRESS,
        address_mpu_master = MPU9050_ADDRESS_68,
        address_mpu_slave = None,
        bus = 1,
        gfs = GFS_1000,
        afs = AFS_8G,
        mfs = AK8963_BIT_16,
        mode = AK8963_MODE_C100HZ)


def calibrate():
    # Calibrate and configure the IMU
    try:
        mpu.configure()
        mpu.calibrate()
        mpu.configure()
        print("SENSOR: Calibrated & configured...\n\n\nBeginning data read...\n\n")
        return True
    
    # Something went wrong
    except:
        return False
    


def get_imu_data():
    data = [*mpu.readMagnetometerMaster(),
            *mpu.readGyroscopeMaster()]
    print(data)
    
    return data
