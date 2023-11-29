#!/usr/bin/env


# Imports
from ukf.UKF_algorithm import *
from ukf.run_UKF import *
from happy_sensors import get_imu_data


def initialize():
    n = 10
    m = 9

    r=np.zeros(n)
    q=np.zeros(m)
    for i in range(m):
        r[i]=random.random()
        q[i]=random.random() * .1

    start = np.random.rand(n)

    cov = np.zeros((n,n))
    for i in range(n):
        cov[i][i]=random.random()
    
    return r, q, start, cov


"""
Main method:
Runs the ADCS loop
    1. Get sensor & GPS data -> IMU (accelerometer, gyroscope, magentometer)
    2. Run UKF -> calculate current state/orientation
        a. Get EoMs input
        b. Get sensor input
    3. Run PID controller
    4. Generate RPM
    5. Calculate RPM limiting
    6. Etc... -> hardware items
"""
def main():
    # Get start state for UKF
    r, q, state, cov = initialize()
    i = 1

    # For now, run infinite loop
    while True:

        # ADCS run from textfile data?
        # imu_data = get_imu_data() # Sensors
        # filename = "./ukf/sensor_data_2.txt" # Textfile
        # run_ukf_textfile(start, cov, r, q, filename)

        ### ADCS with live sensor data
        ## 1/2. Run UKF (get sensor data within UKF) with live sensors: [3 accelerometer, 3 gyroscope, 3 magnetometer]
        state, cov = run_ukf_sensor_iteration(state, cov, r, q, i)

        ## 3. Feed output of UKF into PID script
        #pid_out = pid_run()

        ## 4. Output of PID script to affect actuators (reaction wheels/magnetorquers)
        # TODO  
       
        i += 1   

    return


if __name__ == "__main__":
    # mode = "0" # Textfule
    # mode = "1" # Sensor
    main()
