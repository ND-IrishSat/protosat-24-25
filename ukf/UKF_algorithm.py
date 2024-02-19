'''
UKF_algorithm.py
Authors: Andrew Gaylord, Claudia Kuczun, Micheal Paulucci, Alex Casillas, Anna Arnett
Last modified 2/13/24

Unscented Kalman Filter algorithm for IrishSat based on following resource:
The Unscented Kalman Filter for Nonlinear Estimation 
Eric A. Wan and Rudolph van der Merwe 
Oregon Graduate Institute of Science & Technology 

Variables needed throughout UKF process:
  n = dimensionality of model (7)
  m = dimension of measurement space. Frame of reference of the satellite/sensors (6)
  q: process noise covariance matrix (n x n)
  r: measurement noise covariance matrix (m x m)
  scaling = parameter for sigma point generation (3 - n)
  means = estimated current state (1 x n)
  cov = covariance matrix of current state (n x n)
  predMeans = matrix of predicted means based on physics EOMs (1 x n)
  predCov = matrix of predicted covariance (n x n)
  g = matrix of predicted sigma points (state space using EOMs) (2*n+1 x n)
  h = matrix of transformed sigma points (in the measurement space using hfunc) (2*n+1 x m)
  mesMeans = means in the measurement space after nonlinear transformation (hfunc) (1 x m)
  mesCov = covariance matrix of points in measurement space (m x m)
  data: magnetometer (magnetic field) and gyroscope (angular velocity) data reading from sensor at each step (1 x m)
  kalman = kalman gain for each step, measures how much we should change based on our trust of sensors vs our model (n x m)
'''

import numpy as np
import math
import bigEOMS
import scipy
import scipy.linalg
from hfunc import *
from typing import Optional


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

    # 1) sigma point generation
    # first column of sigma matrix is means
    sigmaMatrix[0] = means

    # take the square root of the inside
    temp = scipy.linalg.sqrtm(np.multiply(cov, (n + scaling)))

    # traverse n dimensions, calculating all other sigma points
    # means + sqrt for 1 to n
    sigmaMatrix[1:(n+1)] = np.add(means, temp)
    # means - sqrt for n + 1 to 2*n
    sigmaMatrix[(n+1):(2*n+1)] = np.subtract(means, temp)

    # return the sigma matrix (2 * n + 1 columns)
    return sigmaMatrix


def EOMs(state, reaction_speeds):
    '''
    EOMs
        uses the current state of the system to output the new predicted state of the system based on physics Equations of Motions
        change dt within to control time step
        reaction_speeds: reaction wheel velocities (1 x 3)
        gps_data: Control input vector (Currently, magnetorquers are not being used, all set to 0)
                [t_motor1, t_motor2, t_motor3, M_mt1, M_mt2, M_mt3]

    @params
        state: column of sigma point matrix to propogate (1 x n)
            [a b c d w_x w_y w_z theta_dot_RW1 theta_dot_RW2 theta_dot_RW3]
    @returns
        x_predicted: next step based on Euler's method (1 x n)
    '''

    # func: instance of EoMs object, defined in eoms.py
    func = bigEOMS.bigEOMS()

    # gps_data: Control input vector (Currently, magnetorquers are not being used, all set to 0)
                # [t_motor1, t_motor2, t_motor3, M_mt1, M_mt2, M_mt3]
    gps_data = np.zeros(6)
    # gps_data[0] = reaction_speeds[0]
    # gps_data[1] = reaction_speeds[1]
    # gps_data[2] = reaction_speeds[2]


    # I_body_tensor: Moment of inertia tensor of the cubesat
                # [[I_XX  I_XY  I_XZ]
                #  [I_YX  I_YY  I_YZ]
                #  [I_ZX  I_ZY  I_ZZ]]
    #I_body_tensor = [[1728.7579, -60.6901, -8.7583],
    #                 [-60.6901, 1745.997, 53.4338],
    #                 [-8.7583, 53.4338, 1858.2584]]
    I_body_tensor = np.identity(3)

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
    theta_dot_RW1 = reaction_speeds[0]
    theta_dot_RW2 = reaction_speeds[1]
    theta_dot_RW3 = reaction_speeds[2]

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
    t_motor1 = gps_data[0]
    t_motor2 = gps_data[1]
    t_motor3 = gps_data[2]
    M_mt1 = gps_data[3]
    M_mt2 = gps_data[4]
    M_mt3 = gps_data[5]

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
    # x_predicted[7] = state[7] + dt * func.theta_ddot_RW1(
    #     M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
    #     I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
    #     I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    # x_predicted[8] = state[8] + dt * func.theta_ddot_RW2(
    #     M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
    #     I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
    #     I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)
    # x_predicted[9] = state[9] + dt * func.theta_ddot_RW3(
    #     M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx,
    #     I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY,
    #     I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3)

    return x_predicted


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


def generateMeans(func, controlVector, sigmaPoints, w0, w1, n, dimensionality):
    '''
    generateMeans
        generate mean after passing sigma point distribution through a transformation function using equation 3b) and 4b)
        also stores and returns all transformed sigma points
            
    @params
        func: transformation/predictive function we are passing sigma points through (H_func or EOMs)
        controlVector: additional input needed for func (gps_data or q_wmm)
        sigmaPoints: sigma point matrix (2xn+1 x n)
        w0, w1: weight for first and all other sigma points, respectively
        n: dimensionality of model 
        dimensionality: dimensionality of what state we are generating for (n or m)

    @returns
        means: mean of distribution in state or measurement space (1 x n or 1 x m)
        transformedSigma: sigma matrix of transformed points (n*2+1 x n or n*2+1 x m)
    '''
    # initialize means and new sigma matrix with correct dimensionality
    means = np.zeros(dimensionality)
    transformedSigma = np.zeros((2 * n + 1, dimensionality))

    # pass all sigma points to the transformation function
    for i in range(1, n * 2 + 1):
        # 3a) and 4a)
        x = func(sigmaPoints[i], controlVector)
        # store calculated sigma point in transformed sigma matrix
        transformedSigma[i] = x
        # update mean with point
        means = np.add(means, x)

    # apply weight to mean without first point
    means *= w1

    # pass first sigma point through transformation function
    x = func(sigmaPoints[0], controlVector) 
    
    # store new point as first element in transformed sigma matrix
    transformedSigma[0] = x

    # adjust the means for first value and multiply by correct weight
    means = np.add(means, x*w0)

    return means, transformedSigma


def generateCov(means, transformedSigma, w0, w1, n, noise):
    '''
    generateCov
        generates covariance matrix from equation 3c) and 4c) based on means and sigma points
        uncoment noise addition once r and q are figured out
        
    @params
        means: means in state or measurement space (1 x n or 1 x m)
        transformedSigma: stored result of passing sigma points through the EOMs or H_func (n*2+1 x m or n*2+1 x n)
        w0, w1: weight for first and all other sigma points, respectively
        n: dimensionality of model 
        noise: noise value array to apply to our cov matrix (r or q)
        
    @returns
        cov: covariance matrix in state or measurement space (n x n or m x m)
    '''
    # find dimension of cov by looking at size of sigma point array
    # prediction points will have n columns, measurement points will have m columns
    covDimension = transformedSigma.shape[1]
    
    # initialize cov with proper dimensionality
    cov = np.zeros((covDimension, covDimension))

    # for all transformed sigma points, apply covariance formula
    for i in range(1, n * 2 + 1):
        # subtract mean from sigma point and multiply by itself transposed
        arr = np.subtract(transformedSigma[i], means)[np.newaxis]
        arr = np.matmul(arr.transpose(), arr)
        cov = np.add(cov, arr)
    
    # separate out first value and update with correct weight
    arr = np.subtract(transformedSigma[0], means)[np.newaxis]
    d = np.matmul(arr.transpose(), arr) * w0

    # use other weight for remaining values
    cov *= w1

    # add back first element
    cov = np.add(cov, d)

    # add noise to covariance matrix
    cov = np.add(cov, noise)

    return cov


def generateCrossCov(predMeans, mesMeans, g, h, w0, w1, n):
    '''
    generateCrossCov
        use equation 5a) to generate cross covariance between our means and sigma points in our state and measurement space

    @params
        predMeans: predicted means based on EOMs (1 x n)
        mesMeans: predicted means in measurement space (1 x m)
        g: sigma point matrix that has passed through the EOMs (n*2+1 x n)
        h: sigma point matrix propogated through non-linear transformation h func (n*2+1 x m)
        w0: weight for first value
        w1: weight for other values
        n: dimensionality of model
    
    @returns
        crossCov: represents uncertainty between our state and measurement space estimates (n x m)
    '''
    m = len(mesMeans)
    crossCov = np.zeros((n,m))

    for i in range(1, n * 2 + 1):
        arr1 = np.subtract(g[i], predMeans)[np.newaxis]
        arr2 = np.subtract(h[i], mesMeans)[np.newaxis]
        arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        crossCov = np.add(crossCov, arr1)

    arr1 = np.subtract(g[0], predMeans)[np.newaxis]
    arr2 = np.subtract(h[0], mesMeans)[np.newaxis]

    # seperate out first element
    d = np.matmul(arr1.transpose(), arr2)

    # multiply by weights for first and other values
    crossCov = np.multiply(crossCov, w1)
    d = np.multiply(d, w0)

    # add first value back into cross covariance
    crossCov = np.add(crossCov, d)

    return crossCov


def UKF(means, cov, q, r, gps_data, reaction_speeds, data):
    '''
    UKF
        estimates state at time step based on sensor data, noise, and equations of motion

    @params
        means: means of previous states (1 x n)
        cov: covariance matrix of state (n x n)
        q: process noise covariance matrix (n x n)
        r: measurement noise covariance matrix (m x m)
        gps_data: control input vector for hfunc (gps data: longitude, latitude, height, time)
        reaction_speeds: control input for reaction wheel speeds (1 x 3)
        data: magnetometer (magnetic field) and gyroscope (angular velocity) data reading from sensor (1 x m)

    @returns
        means: calculated state estimate at current time (1 x n)
        cov: covariance matrix (n x n)
    '''

    # dimensionality of state space = dimension of means
    n = len(means)
    # dimensionality of measurement space = dimension of measurement noise
    m = len(r)
    
    # scaling parameters
    # alpha and k scale points around the mean. To capture the kurtosis of a gaussian distribution, a=1 and k=3-n should be used
    #   If a decrease in the spread of the SPs is desired, use κ = 0 and alpha < 1
    #   If an increase in the spread of the SPs is desired, use κ > 0 and alpha = 1
    alpha = 0.001
    k = 0
    # beta minimizes higher order errors in covariance estimation
    beta = 2
    # equation ??
    scaling = alpha * alpha * (n + k) - n

    # 1) sigma point generation
    sigmaPoints = sigma(means, cov, n, scaling)
    

    # 2) weights calculation
    w0_m = scaling / (n + scaling) # weight for first value for means
    w0_c = scaling / (n + scaling) + (1 - alpha * alpha + beta) # weight for first value for covariance
    w1 = 1 / (2 * (n + scaling)) # weight for all other values


    # 3) predictive step
    # 3a) and 3b): pass sigma points through EOMs (g) and generate mean in state space
    predMeans, g = generateMeans(EOMs, reaction_speeds, sigmaPoints, w0_m, w1, n, n)
    
    # print("PREDICTED MEANS: ", predMeans)
    
    # 3c) generate predicted covariance + process noise q
    predCov = generateCov(predMeans, g, w0_c, w1, n, q)

    # print("PRED COVID: ", predCov)


    # finds true B field based on gps data
    if len(gps_data) == 4:
        Bfield = bfield_calc(gps_data)
    else: 
        # for ideal tests only, use gps_data as b field vector and skip calculating it from the gps data
        Bfield = gps_data


    # print("BFIELD: ", Bfield)

    # 4) non linear transformation
    # 4a) and 4b): non linear transformation of predicted sigma points g into measurement space (h), and mean generation
    mesMeans, h = generateMeans(hfunc, Bfield, g, w0_m, w1, n, m)

    # print("MEAN IN MEASUREMENT: ", mesMeans)

    # 4c) measurement covariance + measurement noise r
    mesCov = generateCov(mesMeans, h, w0_c, w1, n, r)


    # 5) measurement updates
    # 5a) cross covariance: compare our different sets of sigma points and our predicted/measurement means
    crossCov = generateCrossCov(predMeans, mesMeans, g, h, w0_c, w1, n)

    # print("covariance in measurement: ", mesCov)
    # print("cross covariance: ", crossCov)

    # 5b) calculate kalman gain (n x m) by multiplying cross covariance matrix and transposed predicted covariance
    kalman = np.matmul(crossCov, np.linalg.inv(mesCov))

    # print("KALMAN: ", kalman)

    # 5c) updated final mean = predicted + kalman(measurement data - predicted means in measurement space)
    means = np.add(predMeans, np.matmul(kalman, np.subtract(data, mesMeans)))

    # normalize the quaternion to reduce small calculation errors over time
    means[0:4] = means[0:4]/np.linalg.norm(means[0:4])

    # updated covariance = predicted covariance * (n identity matrix - kalman * cross covariance)
    # cov = np.matmul(np.subtract(np.identity(n), np.matmul(kalman, np.transpose(crossCov))), predCov)

    # 5d) updated covariance = predicted covariance - kalman * measurement cov * transposed kalman
    cov = np.subtract(predCov, np.matmul(np.matmul(kalman, mesCov), kalman.transpose()))


    # print("MEANS AT END: ", means)
    # print("COV AT END: ", cov)
    return [means, cov]
