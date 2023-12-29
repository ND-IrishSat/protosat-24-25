'''
UKF_algorithm.py
Authors: Andrew Gaylord, Claudia Kuczun, Micheal Paulucci, Alex Casillas, Anna Arnett
Last modified 10/7/23

Unscented Kalman Filter algorithm for IrishSat based on following resource:
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
  meanInMes = means in the measurement space (1 x m)
  covidInMes = covariance matrix of points in measurement space (m x m)
  z = sensor data (1 x n except switched to 1 x m for some calculations)
  kalman = kalman gain for each step (n x m)
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

    # for i in range(len(cov)):
    #   for j in range(len(cov[i])):
    #     temp[i][j] = cov[i][j]  * (n + scaling)

    # sigma point formula from website
    temp = np.multiply(cov, (n + scaling))

    # take the square root of the matrix
    temp = scipy.linalg.sqrtm(temp)

    # first column of sigma matrix is means
    sigmaMatrix[0] = means

    # traverse N (10) dimensions, calculating all sigma points
    # for i in range(n):
    #     sigmaMatrix[i + 1] = np.add(means, temp[i])  # mean + covariance
    #     sigmaMatrix[i + 1 + n] = np.subtract(means, temp[i])  # mean - covariance
    sigmaMatrix[1:(n+1)] = np.add(means, temp)
    sigmaMatrix[(n+1):(2*n+1)] = np.subtract(means, temp)

    # return the sigma matrix (21 columns)
    return sigmaMatrix


def EOMs(state, reaction_speeds):
    '''
    EOMs
        uses the current state of the system to output the new predicted state of the system based on physics Equations of Motions
        change dt within to control time step
        reaction_speeds: reaction wheel velocities (1 x 3)
        u_k: Control input vector (Currently, magnetorquers are not being used, all set to 0)
                [t_motor1, t_motor2, t_motor3, M_mt1, M_mt2, M_mt3]

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
    u_k[0] = reaction_speeds[0]
    u_k[1] = reaction_speeds[1]
    u_k[2] = reaction_speeds[2]


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


def generateMeans(func, controlVector, sigmaPoints, w1, w2, n, dimensionality):
    # Bfield: Optional[np.array]
    '''
    generateMeans
        generate mean after passing sigma point distribution through a transformation function
        also stores and returns all transformed sigma points
            
    @params
        func: transformation/predictive function we are passing sigma points through (H_func or EOMs)
        controlVector: additional input needed for func (u_k or q_wmm)
        sigmaPoints: sigma point matrix (2xn+1 x n)
        w1, w2: weight for first and all other sigma points, respectively
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
        x = func(sigmaPoints[i], controlVector)
        # store calculated sigma point in transformed sigma matrix
        transformedSigma[i] = x
        # update mean with point
        means = np.add(means, x)

    # apply weight to mean without first point
    means *= w2

    # pass first sigma point through transformation function
    x = func(sigmaPoints[0], controlVector) 
    # store new point as first element in transformed sigma matrix
    transformedSigma[0] = x

    # adjust the means for first value and multiply by correct weight
    means = np.add(means, x*w1)

    return means, transformedSigma


def generateCov(means, transformedSigma, w1, w2, n, noise):
    '''
    generateCov
        generates covariance matrix from website equation based on means and sigma points
        uncoment noise addition once r and q are figured out
        
    @params
        means: means in state or measurement space (1 x n or 1 x m)
        transformedSigma: stored result of passing sigma points through the EOMs or H_func (n*2+1 x m or n*2+1 x n)
        w1, w2: weight for first and all other sigma points, respectively
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
    d = np.matmul(arr.transpose(), arr) * w1

    # use other weight for remaining values
    cov *= w2

    # add back first element
    cov = np.add(cov, d)

    # add noise to covariance matrix
    cov = np.add(cov, noise)

    return cov



def UKF(passedMeans, passedCov, r, q, u_k, reaction_speeds, data):
    '''
    UKF
        estimates state at time step based on sensor data, noise, and equations of motion

    @params
        passedMeans: means of previous states (1 x n)
        passedCov: covariance matrix of state (n x n)
        r: noise vector of predictions (1 x n)
        q: noise vector of sensors (1 x m)
        u_k: control input vector for hfunc (gps data: longitude, latitude, height, time)
        reaction_speeds: control input for reaction wheel speeds (1 x 3)
        data: magnetometer (magnetic field) and gyroscope (angular velocity) data reading from sensor (1 x 6)
    @returns
        means: calculated state estimate at current time (1 x n)
        cov: covariance matrix (n x n)
    '''

    # initialize vars (top of file for descriptions)
    n = 7
    m = n
    cov = passedCov
    predMeans = np.zeros(n)
    predCovid = np.zeros((n,n))
    meanInMes = np.zeros(m)
    covidInMes = np.zeros((m, m))
    crossCo = np.zeros((n,m))
    g = np.zeros((n * 2 + 1, n))
    h = np.zeros((n * 2 + 1, m))
    
    # z = np.array([0.0] + data)
    z = data
    

    scaling = 3-n
    w1 = scaling / (n + scaling) # weight for first value
    w2 = 1 / (2 * (n + scaling)) # weight for all other values
    
    # track the average of the estimated K states so far
    means = passedMeans


    """
    Calculate mean of Gaussian (populates the global predicted means matrix)
    1. Store temporary sigma points
    2. Apply the EOMs to the temporary (stored) sigma points
    2. Calculate the means of sigma points without weights
    4. Calculate the new predicted means by applying predetermined weights
    Initialize sigma matrix to the starting sigma point matrix
    """

    sigTemp = sigma(means, cov, n, scaling)  # temporary sigma points
    # print("SIGMA POINTS", sigTemp)

    predMeans, g = generateMeans(EOMs, reaction_speeds, sigTemp, w1, w2, n, n)
    
    # print("PREDICTED MEANS: ", predMeans)

    """
    Calculate predicted covariance of Gaussian
    """
    predCovid = generateCov(predMeans, g, w1, w2, n, r)


    """ 
    Mean to measurement
    """
    # create arbitrary covariance for sensors
    # zCov = [[.1, 0, 0, 0, 0, 0, 0],
    #         [0, .1, 0, 0, 0, 0, 0],
    #         [0, 0, .1, 0, 0, 0, 0],
    #         [0, 0, 0, .1, 0, 0, 0],
    #         [0, 0, 0, 0, .1, 0, 0],
    #         [0, 0, 0, 0, 0, .1, 0],
    #         [0, 0, 0, 0, 0, 0, .1]
    # ]
    # for i in range(n):
    #     zCov[i][i] = .005

    # zCov = cov
  
    # create temporary sigma points
    # sigTemp = sigma(z, zCov, n, scaling)

    # Bfield = bfield_calc(u_k)
    # Bfield = u_k[:3]

    # print("BFIELD: ", Bfield)

    # update with gps control vector instead of q_wmm
    meanInMes, h = generateMeans(hfunc, u_k, sigTemp, w1, w2, n, m)

    # print("MEAN IN MEASUREMENT: ", meanInMes)

    # print("Transformed sigma points: ", h)


    """
    Creates covariance matrix in measurement space
    """
    covidInMes = generateCov(meanInMes, h, w1, w2, n, q)


    '''
    Cross covariance matrix (t) between state space and predicted space

    Remake sigma points here now that we have new data up to the group?
    '''
    # sig = sigma(means, cov, n, scaling)
    sig = sigTemp

    # use formula from website to compare our different sets of sigma points and our predicted/measurement means
    for i in range(1, n * 2 + 1):
        arr1 = np.subtract(sig[i], predMeans)[np.newaxis]
        arr2 = np.subtract(h[i], meanInMes)[np.newaxis]
        arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        # print("     crossCo calc: ", arr1)
        crossCo = np.add(crossCo, arr1)
        # arr1 = np.subtract(h[i], meanInMes)[np.newaxis]
        # arr2 = np.subtract(sig[i], predMeans)[np.newaxis]
        # arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        # crossCo = np.add(crossCo, arr1)
    '''switch ordering??'''

    arr1 = np.subtract(sig[0], predMeans)[np.newaxis]
    arr2 = np.subtract(h[0], meanInMes)[np.newaxis]
    # arr1 = np.subtract(h[i], meanInMes)[np.newaxis]
    # arr2 = np.subtract(sig[i], predMeans)[np.newaxis]

    # seperate out first element
    d = np.matmul(arr1.transpose(), arr2)

    # multiply by weights for first and other values
    crossCo = np.multiply(crossCo, w2)
    d = np.multiply(d, w1)

    # add first value back into cross covariance
    crossCo = np.add(crossCo, d)

    """
    Kalman gain and final update
    """
    # calculate kalman gain by multiplying cross covariance matrix and transposed predicted covariance
    # n x m
    # print("covariance in measurement: ", covidInMes)
    # print("cross covariance: ", crossCo)
    kalman = np.matmul(crossCo, np.linalg.inv(covidInMes))

    # print("KALMAN: ", kalman)

    # z = z[1:]
  
    # updated final mean = predicted + kalman(measurement - predicted in measurement space)
    means = np.add(predMeans, np.matmul(kalman, np.subtract(z, meanInMes)))
    # normalize the quaternion?
    # normal = np.linalg.norm(means[0:4])
    # means[0:4] = means[0:4]/normal

    # 2 options for updated cov:
    
    # option 1:
    # updated covariance = predicted covariance * (n identity matrix - kalman * cross covariance)
    cov = np.matmul(np.subtract(np.identity(m), np.matmul(kalman, crossCo)), predCovid)
    # this one doesn't work with different n and m for some reason

    # option 2:
    # updated covariance = predicted covariance - (kalman * covariance in measurement * transposed kalman)
    # cov = np.subtract(predCovid, np.matmul(np.matmul(kalman, covidInMes), kalman.transpose()))

    # print("MEANS AT END: ", means)
    # print("COV AT END: ", cov)
    return [means, cov]
