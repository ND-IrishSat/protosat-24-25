from mpu9250_jmdev.registers import *
from mpu9250_jmdev.mpu_9250 import MPU9250


def test_on_pi():
    
        mpu = MPU9250(
                address_ak = AK8963_ADDRESS,
                address_mpu_master = MPU9050_ADDRESS_68,
                address_mpu_slave = None,
                bus = 1,
                gfs = GFS_1000,
                afs = AFS_8G,
                mfs = AK8963_BIT_16,
                mode = AK8963_MODE_C100HZ)

        mpu.configure()

        print("JERE")

        while True:
                print("MPU 9250 running...")
                print("accelerometer: ", mpu.readAccelerometerMaster())
                print("gyroscope: ", mpu.readGyroscopeMaster())
                print("magnetometer: ", mpu.readMagnetometerMaster())

                time.sleep(1)



def get_imu_data():
        # Create IMU instance
        mpu = MPU9250(
                address_ak = AK8963_ADDRESS,
                address_mpu_master = MPU9050_ADDRESS_68,
                address_mpu_slave = None,
                bus = 1,
                gfs = GFS_1000,
                afs = AFS_8G,
                mfs = AK8963_BIT_16,
                mode = AK8963_MODE_C100HZ)
        
        # Get the data
        res = []
        res.append(mpu.readAccelerometerMaster())
        res.append(mpu.readGyroscopeMaster())
        res.append(mpu.readMagnetometerMaster())


        return res

if __name__ == "__main__":
        print(get_imu_data())
