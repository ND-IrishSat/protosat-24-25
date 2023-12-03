from mpu9250_i2c import *
import time


t0 = time.time()
start_bool = False # boolean for connection
while (time.time()-t0)<5: # wait for 5-sec to connect to IMU
    try:
        print("TRYING")
        from mpu9250_i2c import *
        start_bool = True # True for forthcoming loop
        break 
    except:
        continue
        
        
imu_devs   = ["ACCELEROMETER","GYROSCOPE","MAGNETOMETER"]
imu_labels = ["x-dir","y-dir","z-dir"]
imu_units  = ["g","g","g","dps","dps","dps","uT","uT","uT"]

print("PAST TRYING") 
f = open("testfile-long-random.csv", "w")
count = 0
while count < 7000:
    count += 1
    print(count)
    if start_bool==False: # make sure the IMU was started
        print("IMU not Started, Check Wiring") # check wiring if error
        break
    ##################################
    # Reading and Printing IMU values
    ##################################
    #
    try:
        ax,ay,az,wx,wy,wz = mpu6050_conv() # read and convert mpu6050 data
        mx,my,mz = AK8963_conv() # read and convert AK8963 magnetometer data
        
        line = str([ax,ay,az,wx,wy,wz,mx,my,mz])
        line += "\n"
        print(line)
        f.write(line)
    except:
        continue 
    #
    ##################################
    # Reading and Printing IMU values
    ##################################
    #
    # print(50*"-")
    # for imu_ii,imu_val in enumerate([ax,ay,az,wx,wy,wz,mx,my,mz]):
        # if imu_ii%3==0:
            # print(20*"_"+"\n"+imu_devs[int(imu_ii/3)]) # print sensor header
        # #
        # ###############
        # # Print Data
        # ###############
        # #
        # print("{0}: {1:3.2f} {2}".format(imu_labels[imu_ii%3],imu_val,imu_units[imu_ii]))
        
    #time.sleep() # wait between prints
    count += 1
f.close()
