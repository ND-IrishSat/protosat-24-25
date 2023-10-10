'''
UKF_algorithm.py
Authors: Andrew Gaylord, Claudia Kuczun, Micheal Paulucci, Alex Casillas, Anna Arnett
Last modified 10/7/23

UKF algorithm for IrishSat based on following resource:
https://towardsdatascience.com/the-unscented-kalman-filter-anything-ekf-can-do-i-can-do-it-better-ce7c773cf88d

Variables needed throughout UKF process:
  n = dimensionality of model (10)
  m = dimension of measurement space that excludes first 0 of quaternion (9)
  r = noise vector for predictions (we choose this & make it) (1 x n)
  q = noise vector for sensors (provided on data sheet for sensors) (1 x m)
  scaling = parameter for sigma point generation (3 - n)
  means = mean of gaussian of estimated states so far (maybe???) (1 x n)
  cov = covariance matrix (n x n)
  predMeans = matrix of predicted means (1 x n)
  predCovid = matrix of predicted covariance (n x n)
  g = matrix of predicted sigma points (state space using EOMs) (2*n+1 x n)
  h = matrix of transformed sigma points (in the measurement space) (2*n+1 x m)
  q_wmm = B field represented as quaternion (1 x 4)
  meanInMes = means in the measurement space (1 x m)
  covidInMes = covariance matrix of points in measurement space (m x m)
  z = sensor data (1 x n except switched to 1 x m for some calculations)
  kalman = kalman gain for each step (n x m)
'''

import numpy as np
import math
import bigEOMS
import random
import scipy
import scipy.linalg
#from statsmodels import *
#from statsmodels.stats import correlation_tools


def sigma(means, cov, n, scaling):
    '''
    sigma
        creates sigma point matrix based on formula that represents distribution of means and cov

    @params
        means: mean of gaussian of estimated states so far (maybe???). Also first column of sigma matrix (1 x n)
        cov: covariance matrix of state (n x n)
        n: dimensionality of model
        scaling: how far from mean we distribute our points, used in sigma point formula
    @returns
        sigmaMatrix: matrix of sigma points (2 * n + 1, n) 
    '''
    # intialize 2N + 1 sigma points to zeroes
    sigmaMatrix = np.zeros((2*n+1,n))
    temp = np.zeros((n, n))

    for i in range(len(cov)):
      for j in range(len(cov[i])):
        #sigma point formula from website
        temp[i][j] = cov[i][j]  * (n + scaling)

    temp = scipy.linalg.sqrtm(temp)
    # first column of sigma matrix is means
    sigmaMatrix[0] = means

    # traverse N (10) dimensions, calculating all sigma points
    for i in range(n):
        sigmaMatrix[i + 1] = np.add(means, temp[i])  # mean + covariance
        sigmaMatrix[i + 1 + n] = np.subtract(means, temp[i])  # mean - covariance

    # return the sigma matrix (21 columns)
    return sigmaMatrix


def EOMs(state, u_k):
    '''
    EOMs
        uses the current state of the system to output the new predicted state of the system based on physics Equations of Motions
        change dt within to control time step
        u_k: control input vector for magnetorquers?

    @params
        state: column of sigma point matrix to propogate (1 x n)
            [a b c d w_x w_y w_z theta_dot_RW1 theta_dot_RW2 theta_dot_RW3]
    @returns
        x_predicted: next step based on Euler's method (1 x n)
    '''

    # func: instance of EoMs object, defined in eoms.py
    func = bigEOMS.bigEOMS()

    # u_k: Control input vector (Currently, magnetorquers are not being used, all set to 0)
                # [t_motor1, t_motor2, t_motor3, M_mt1, M_mt2, M_mt3]
    u_k = np.zeros(6)

    # I_body_tensor: Moment of inertia tensor of the cubesat
                # [[I_XX  I_XY  I_XZ]
                #  [I_YX  I_YY  I_YZ]
                #  [I_ZX  I_ZY  I_ZZ]]
    I_body_tensor = [[1728.7579, -60.6901, -8.7583],
                     [-60.6901, 1745.997, 53.4338],
                     [-8.7583, 53.4338, 1858.2584]]
    
    # I_RW: Moment of inertias of the three reaction wheels
                # [I_RW1 I_RW2 I_RW3]
    I_RW = [578.5944, 578.5944, 578.5944]

    # dt: The timestep between the current state and predicted state
    dt = 0.1

    # Initialize prediction
    x_predicted = np.zeros(len(state))

    # Grab components of state vector
    a = state[0]
    b = state[1]
    c = state[2]
    d = state[3]
    w_x = state[4]
    w_y = state[5]
    w_z = state[6]
    theta_dot_RW1 = state[7]
    theta_dot_RW2 = state[8]
    theta_dot_RW3 = state[9]

    # Grab moment of inertias
    I_xx = I_body_tensor[0][0]
    I_xy = I_body_tensor[0][1]
    I_xz = I_body_tensor[0][2]
    I_yx = I_body_tensor[1][0]
    I_yy = I_body_tensor[1][1]
    I_yz = I_body_tensor[1][2]
    I_zx = I_body_tensor[2][0]
    I_zy = I_body_tensor[2][1]
    I_zz = I_body_tensor[2][2]
    I_RW1_XX = I_RW[0]
    I_RW2_YY = I_RW[1]
    I_RW3_ZZ = I_RW[2]

    # Grab components of control input
    t_motor1 = u_k[0]
    t_motor2 = u_k[1]
    t_motor3 = u_k[2]
    M_mt1 = u_k[3]
    M_mt2 = u_k[4]
    M_mt3 = u_k[5]

    # Do Euler's method to get next state
    x_predicted[0] = state[0] + dt * func.adot(a, b, c, d, w_x, w_y, w_z)
    x_predicted[1] = state[1] + dt * func.bdot(a, b, c, d, w_x, w_y, w_z)
    x_predicted[2] = state[2] + dt * func.cdot(a, b, c, d, w_x, w_y, w_z)
    x_predicted[3] = state[3] + dt * func.ddot(a, b, c, d, w_x, w_y, w_z)
    x_predicted[4] = state[4] + dt * func.w_dot_x(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    x_predicted[5] = state[5] + dt * func.w_dot_y(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    x_predicted[6] = state[6] + dt * func.w_dot_z(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    x_predicted[7] = state[7] + dt * func.theta_ddot_RW1(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    x_predicted[8] = state[8] + dt * func.theta_ddot_RW2(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    x_predicted[9] = state[9] + dt * func.theta_ddot_RW3(
        M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
        I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)

    return x_predicted


"""
Functions of getting predicted measurement from magnetometer (quaternion of B-field in body frame) from state (quaternion of body frame with reference to ECI)
"""
def observe_b_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm):
    return a_kf*(a_kf*b_wmm - c_kf*d_wmm + c_wmm*d_kf) - c_kf*(a_kf*d_wmm - b_kf*c_wmm + b_wmm*c_kf) + b_kf*(b_kf*b_wmm + c_kf*c_wmm + d_kf*d_wmm) + d_kf*(a_kf*c_wmm + b_kf*d_wmm - b_wmm*d_kf)


def observe_c_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm):
    return a_kf*(a_kf*c_wmm + b_kf*d_wmm - b_wmm*d_kf) + b_kf*(a_kf*d_wmm - b_kf*c_wmm + b_wmm*c_kf) + c_kf*(b_kf*b_wmm + c_kf*c_wmm + d_kf*d_wmm) - d_kf*(a_kf*b_wmm - c_kf*d_wmm + c_wmm*d_kf)


def observe_d_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm):
    return a_kf*(a_kf*d_wmm - b_kf*c_wmm + b_wmm*c_kf) - b_kf*(a_kf*c_wmm + b_kf*d_wmm - b_wmm*d_kf) + c_kf*(a_kf*b_wmm - c_kf*d_wmm + c_wmm*d_kf) + d_kf*(b_kf*b_wmm + c_kf*c_wmm + d_kf*d_wmm)


def H_func(state, q_wmm): 
    '''
    H_func
        transforms sigma points from state space to measurement space by running transformation and removing unnecessary first element of quaternion

    @params
        state: estimate in state space (1 x n)
        q_wmm: B field of ECI frame (measurement space) represented as quaternion (1 x 4, first element 0)
            [0, bx, by, bz] normalized by dividing by magnitude
    @returns
        za: state estimate in measurement space (1 x m)
    '''

    # Grab components from vectors
    a_kf = state[0]
    b_kf = state[1]
    c_kf = state[2]
    d_kf = state[3]
    w_x = state[4]
    w_y = state[5]
    w_z = state[6]
    theta_dot_RW1 = state[7]
    theta_dot_RW2 = state[8]
    theta_dot_RW3 = state[9]

    a_wmm = q_wmm[0]
    b_wmm = q_wmm[1]
    c_wmm = q_wmm[2]
    d_wmm = q_wmm[3]

    # Perform observation function, only needed for quaternion components. The rest have 1 to 1 mapping
    a_B_BF = 0  # For now, assume this will always be zero. Trust P. Wensing!!!
    b_B_BF = observe_b_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm)
    c_B_BF = observe_c_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm)
    d_B_BF = observe_d_B_BF(a_kf, b_kf, c_kf, d_kf, a_wmm, b_wmm, c_wmm, d_wmm)

    # Observation vector
    za = np.array([
        b_B_BF, c_B_BF, d_B_BF, w_x, w_y, w_z, theta_dot_RW1,
        theta_dot_RW2, theta_dot_RW3
    ])

    
    #normalized version if we need it idk
    #should already be normal, but isn't occasually due to rounding errors. Could normalize sometimes
    '''
    normalize=1/math.sqrt(abs(a_B_BF*a_B_BF + b_B_BF*b_B_BF + c_B_BF*c_B_BF + d_B_BF*d_B_BF))
    Nza = np.array([
        a_B_BF * normalize, b_B_BF * normalize, c_B_BF * normalize, 
        d_B_BF * normalize, w_x, w_y, w_z, theta_dot_RW1,
        theta_dot_RW2, theta_dot_RW3
    ])
    '''
    
    return za


def quaternionMultiply(a, b):
    '''
    quaternionMultiply
        custom function to perform quaternion multiply on two passed-in matrices

    @params
        a, b: quaternion matrices (1 x 4)
    @returns
        multiplied quaternion matrix
    '''
    return [[a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]],
            [a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]],
            [a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]],
            [a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]]]


def UKF(passedMeans, passedCov, r, q, u_k, data):
    '''
    UKF
        estimates state at for time step based on sensor data, noise, and equations of motion

    @params
        passedMeans: means of previous states (1 x n)
        passedCov: covariance matrix of state (n x n)
        r: noise vector of predictions (1 x n)
        q: noise vector of sensors (1 x m)
        u_k: control input vector for EOMs step (gps data)
        data: magnetometer (magnetic field) and gyroscope (angular velocity) data reading from sensor (1 x 6)
    @returns
        means: calcuated state estimate at current time (1 x n)
        cov: covariance matrix (n x n)
    '''

    # initialize vars (top of file for descriptions)
    n = 10
    m = 9
    cov = passedCov
    predCovid = np.zeros((n,n))
    meanInMes = np.zeros(m)
    covidInMes = np.zeros(m)
    h = np.zeros((2 * n + 1,m))
    g = np.zeros((n * 2 + 1, n))
    q_wmm = []
    q_wmm.append(0)
    q_wmm.append(data[0])
    q_wmm.append(data[1])
    q_wmm.append(data[2])
    # z = passedMeans
    z=[]
    z.append(0)
    z.append(data[0])
    z.append(data[1])
    z.append(data[2])
    z.append(data[3])
    z.append(data[4])
    z.append(data[5])
    z.append(passedMeans[7])
    z.append(passedMeans[8])
    z.append(passedMeans[9])

    scaling = 3-n
    w1 = scaling / (n + scaling) # weight for first value
    w2 = 1 / (2 * ( n + scaling)) # weight for all other values
    
    # track the average of the estimated K states so far
    means = passedMeans


    """
    Calculate mean of Gaussian (populates the global predicted means matrix)
    1. Store temporary sigma points
    2. Apply the EOMs to the temporary (stored) sigma points
    2. Calculate the means of sigma points without weights
    4. Calculate the new predicted means by applying predetermined weights
    Let the sigma matrix = the starting sigma point matrix
    """
    # initialize the means array to zeroes
    predMeans = np.zeros(n)

    sigTemp = sigma(means, cov, n, scaling)  # temporary sigma points
    # oldSig = sigTemp

    for i in range(1, n * 2 + 1):  # generate 2N+1 sigma points
        x = EOMs(sigTemp[i], u_k)  # use state estimation equations
        g[i] = x  # add the next entry to g matrix
        predMeans = np.add(predMeans,x)  # calculate means of sigma points w/out weights

    # apply weights to predicted means
    predMeans *= w2   # w2 for later weights
    x = EOMs(sigTemp[0], u_k)  # calculate EoMs for first sigma point
    g[0] = x  # add first sigma point to first index in g(x)
    predMeans = np.add(predMeans, x*w1)  # w1 for first weight
    

    """
    Calculate predicted covariance of Gaussian
    """
    # for all sigma points
    for i in range(1, n * 2 + 1):
        arr = np.subtract(g[i], predMeans)[np.newaxis]
        # subtract the predicted mean from the transformed sigma points

        arr = np.matmul(arr.transpose(), arr)
        # matrix multiplication: multiply the matrix by itself transposed!
        predCovid = np.add(predCovid, arr)

    arr = np.subtract(g[0], predMeans)[np.newaxis]
    predCovid*=w2
    d = np.matmul(arr.transpose(), arr)*w1  # d: separates out first value

    # add d back to predicted covariance matrix
    predCovid = np.add(predCovid, d)
  

    """ 
    Mean to measurement
    """
    # create arbitrary covariance for sensors
    zCov = [[.1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, .1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, .1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, .1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, .1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, .1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, .1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, .1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, .1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, .1]
    ]
    for i in range(n):
        zCov[i][i] = .2
    # zCov = cov
    # z = means
    # create temporary sigma points
    sigTemp = sigma(z, zCov, n, scaling)
    
    # pass the sigma point to the h function
    for i in range(1, n * 2 + 1):
        x = H_func(sigTemp[i], q_wmm)
        # x = sigTemp[i] 
        '''works'''
        # transforms sigma points into measurement space
        h[i] = x  # store the calculated point in the h matrix
        meanInMes = np.add(meanInMes, x)  # update mean in measurement mean

    meanInMes *= w2  # weight for later value

    x = H_func(sigTemp[0], q_wmm)  # get first mapped point
    # x = sigTemp[0]
    '''works'''
    h[0] = x  # set the first element in h matrix

    # adjust the means in measurement space for first value
    meanInMes = np.add(meanInMes, [i * w1 for i in x])


    """
    Creates covariance matrix in measurement space
    """
    for i in range(1, n * 2 + 1):
        arr = np.subtract(h[i], meanInMes)[np.newaxis]
        arr = np.matmul(arr.transpose(), arr)
        covidInMes = np.add(covidInMes, arr)
    
    arr = np.subtract(h[0], meanInMes)[np.newaxis]
    d = np.matmul(arr.transpose(), arr)  #ordering?

    for i in range(m):
        for j in range(m):
            covidInMes[i][j] *= w2
            d[i][j] *= w1
    covidInMes = np.add(covidInMes, d)
    '''remove/add sensor noise here '''
    # covidInMes=np.add(covidInMes,q) 


    '''
    Cross covariance matrix (t) between state space and predicted space

    Remake sigma points here now that we have new data up to the group
    '''
    # sig = sigma(means, cov, n, scaling)
    sig = sigTemp
    crossCo = np.zeros((n,m))

    for i in range(1, n * 2 + 1):
        arr1 = np.subtract(sig[i], predMeans)[np.newaxis]
        arr2 = np.subtract(h[i], meanInMes)[np.newaxis]
        arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        crossCo = np.add(crossCo, arr1)
        # arr1 = np.subtract(h[i], meanInMes)[np.newaxis]
        # arr2 = np.subtract(sig[i], predMeans)[np.newaxis]
        # arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        # crossCo = np.add(crossCo, arr1)
    '''switch ordering??'''

    # arr1 = np.subtract(h[-1], meanInMes)[np.newaxis]
    arr1 = np.subtract(sig[0], predMeans)[np.newaxis]
    arr2 = np.subtract(h[0], meanInMes)[np.newaxis]

    d = np.matmul(arr1.transpose(), arr2)

    for i in range(n):
        for j in range(m):
            crossCo[i][j] *= w2
            d[i][j] *= w1

    crossCo = np.add(crossCo, d)


    """
    Kalman gain and final update
    """
    # calculate kalman gain by multiplying cross covariance matrix and transposed predicted covariance
    # nxm
    kalman = np.matmul(crossCo, np.linalg.inv(covidInMes))

    z = z[1:]
    # updated final mean = predicted + kalman(measurement - predicted in measurement space)
    means = np.add(predMeans, np.matmul(kalman, np.subtract(z, meanInMes)))

    # 2 options for updated cov:
    
    # updated covariance = predicted covariance * (n identity matrix - kalman * cross covariance)
    # cov = np.matmul(np.subtract(np.identity(m), np.matmul(kalman, crossCo)), predCovid)
    # this one doesn't work with different n and m for some reason

    # updated covariance = predicted covariance - (kalman * covariance in measurement * transposed kalman)
    cov = np.subtract(predCovid, np.matmul(np.matmul(kalman, covidInMes), kalman.transpose()))

    return [means, cov]
