'''
1D_UKF_algorithm.py
edited by Juwan

Unscented Kalman Filter algorithm for IrishSat based on following resource:
The Unscented Kalman Filter for Nonlinear Estimation 
Eric A. Wan and Rudolph van der Merwe 
Oregon Graduate Institute of Science & Technology 

Variables needed throughout UKF process:
  n = dimensionality of model (7)
  m = dimension of measurement space. Frame of reference of the satellite/sensors (6)
  q: process noise covariance matrix (n x n)
  r: measurement noise covariance matrix (m x m)
  scaling = parameter for sigma point generation (equal to alpha^2 * (n + k) - n)
  means = estimated current state (1 x n)
  cov = covariance matrix of current state (n x n)
  predMeans = matrix of predicted means based on physics EOMs (1 x n)
  predCov = matrix of predicted covariance (n x n)
  f = matrix of predicted sigma points (state space using EOMs) (2*n+1 x n)
  h = matrix of transformed sigma points (in the measurement space using hfunc) (2*n+1 x m)
  mesMeans = means in the measurement space after nonlinear transformation (hfunc) (1 x m)
  mesCov = covariance matrix of points in measurement space (m x m)
  data: magnetometer (magnetic field) and gyroscope (angular velocity) data reading from sensor at each step (1 x m)
  kalman = kalman gain for each step, measures how much we should change based on our trust of sensors vs our model (n x m)

NOTE: This implementation uses simplified rotational dynamics, i.e. 2D dynamics
Simplifications: quaternion (4x1) -> psi (angle, +psi_dot points out of the page)
                 w_sat (3x1) -> psi_dot
'''

import numpy as np
import scipy
import scipy.linalg
from hfunc_1D import *
from typing import Optional


def sigma(means, cov, n, scaling):
    '''
    sigma
        creates sigma point matrix that is a representative sampling of the mean and covariance of the system (eq 5-7)
        this allows for efficient transferal of information through our nonlinear function as we estimate our measurement space

    @params
        means: mean of estimated states so far. Also first column of sigma matrix (eq 5) (1 x n)
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


class DEMO_1D_EOMS():
    def __init__(self, I_body: np.ndarray, I_w_spin: float, I_w_trans: float):
        
        # Initially store moment of inertia tensor w/o reaction wheel inertias!
        self.I_body = I_body

        # Store principal moment of inertia for reaction wheels about spin axis and about axis transverse to spin axis respectively
        self.I_w_spin = I_w_spin
        self.I_w_trans = I_w_trans

        # changes for test
        #alpha = 1/np.sqrt(3)
        #beta = 1/np.sqrt(3)
        #gamma = 1/np.sqrt(3) 

        # self.rw_config = np.array([[1, 0, alpha], [0, 1, beta], [0, 0, gamma]])
        # Double check to make sure rw_config for dynamics for 1D test is correct!
        # make sure x axis aligns with wheel
        self.rw_config = np.identity(3)

        # Update: self.rw_config is almost certainly wrong here, need to rewrite! Looks like it's a rotation about Z
        # by theta_1D = 135 degrees
        #theta_1D = 135*np.pi/180.0
        #self.rw_config = np.array([[np.cos(theta_1D), np.sin(theta_1D), 0], [-np.sin(theta_1D), np.cos(theta_1D), 0], [1/np.sqrt(2), 0, 1/np.sqrt(2)]])
        
        # Calculate contributions of reaction wheel to moment of inertia tensor due to principal moment transverse to the spin axis
        for i in np.arange(self.rw_config.shape[1]):
            self.I_body = self.I_body + I_w_trans*(np.identity(3) - np.matmul(self.rw_config[:, i], np.transpose(self.rw_config[:, i]))) 

    def eoms(self, psi: float, psi_dot: float, alpha_rw: float, dt: float):
        
        # Grab moment of inertia about z axis
        Izz = self.I_body[2, 2]

        # Calculate psi_ddot, i.e. angular acceleration
        psi_ddot = -(self.I_w_spin/Izz) * alpha_rw

        # Propagate
        new_psi = psi + psi_dot * dt
        new_psi_dot = psi_dot + psi_ddot * dt

        #print(new_psi)
        #print(new_psi_dot)
        
        new_state = np.array([new_psi, new_psi_dot])

        return new_state
    
    
def generatePredMeans(eomsClass, sigmaPoints, w0, w1, reaction_speed, old_reaction_speed, n):
    '''
    generatePredMeans
        generate mean (eq 9) after passing sigma point distribution through a transformation function (eq 8)
        also stores and returns all transformed sigma points
            
    @params
        eomsClass: EOMs class to pass our sigma points through
        sigmaPoints: sigma point matrix (2xn+1 x n)
        w0, w1: weight for first and all other sigma points, respectively
        reaction_speeds/old_reaction_speeds: reaction wheel speeds for current and last time step (float for constrained demo)
        n: dimensionality of state space

    @returns
        means: mean of distribution in state or measurement space (1 x n or 1 x m)
        transformedSigma: sigma matrix of transformed points (n*2+1 x n or n*2+1 x m)
    '''
    # initialize means and new sigma matrix with correct dimensionality
    means = np.zeros(n)
    transformedSigma = np.zeros((2 * n + 1, n))

    # CHANGE TO PASS TO FUNCITON THROUGHOUT??
    dt = 0.1
    # calculate angular acceleration using old and current reaction wheel speeds
    alpha = (reaction_speed - old_reaction_speed) / dt

    # pass all sigma points to the transformation function
    for i in range(1, n * 2 + 1):
        # 3a) and 4a)
        x = eomsClass.eoms(sigmaPoints[i][0], sigmaPoints[i][1], alpha, dt)

        # store calculated sigma point in transformed sigma matrix
        transformedSigma[i] = x
        # update mean with point
        means = np.add(means, x)

    # apply weight to mean without first point
    means *= w1

    # pass first sigma point through transformation function
    x = eomsClass.eoms(sigmaPoints[i][0], sigmaPoints[i][1], alpha, dt)
    
    # store new point as first element in transformed sigma matrix
    transformedSigma[0] = x

    # adjust the means for first value and multiply by correct weight
    means = np.add(means, x*w0)

    return means, transformedSigma


def generateMesMeans(func, controlVector, sigmaPoints, w0, w1, n, dimensionality):
    '''
    generateMesMeans
        generate mean (eq 12) after passing sigma point distribution through non-linear transformation function (eq 11)
        also stores and returns all transformed sigma points
            
    @params
        func: transformation function we are passing sigma points through (H_func)
        controlVector: additional input needed for func: true magnetic field (1 x 3)
        sigmaPoints: sigma point matrix (2xn+1 x n)
        w0, w1: weight for first and all other sigma points, respectively
        n: dimensionality of model 
        dimensionality: dimensionality of what state we are generating for (measurement space: m)

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
        generates covariance matrix from eq 10 and 13 based on means and sigma points
        
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


def generateCrossCov(predMeans, mesMeans, f, h, w0, w1, n):
    '''
    generateCrossCov
        use equation 14 to generate cross covariance between our means and sigma points in our state and measurement space

    @params
        predMeans: predicted means based on EOMs (1 x n)
        mesMeans: predicted means in measurement space (1 x m)
        f: sigma point matrix that has passed through the EOMs (n*2+1 x n)
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
        arr1 = np.subtract(f[i], predMeans)[np.newaxis]
        arr2 = np.subtract(h[i], mesMeans)[np.newaxis]
        arr1 = np.matmul(arr1.transpose(), arr2)  # ordering?
        crossCov = np.add(crossCov, arr1)

    arr1 = np.subtract(f[0], predMeans)[np.newaxis]
    arr2 = np.subtract(h[0], mesMeans)[np.newaxis]

    # seperate out first element
    d = np.matmul(arr1.transpose(), arr2)

    # multiply by weights for first and other values
    crossCov = np.multiply(crossCov, w1)
    d = np.multiply(d, w0)

    # add first value back into cross covariance
    crossCov = np.add(crossCov, d)

    return crossCov


def UKF(means, cov, q, r, gps_data, reaction_speed, old_reaction_speed, data):
    '''
    UKF
        estimates state at time step based on sensor data, noise, and equations of motion

    @params
        means: means of previous states (1 x n)
        cov: covariance matrix of state (n x n)
        q: process noise covariance matrix (n x n)
        r: measurement noise covariance matrix (m x m)
        gps_data: given B field
        reaction_speeds: control input for EOMs (1 x 4)
        old_reaction_speeds: speeds for past step, used to find angular acceleration (1 x 4)
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

    # eq 1: scaling factor lambda
    scaling = alpha * alpha * (n + k) - n

    # eq 2-4: weights calculation
    w0_m = scaling / (n + scaling) # weight for first value for means
    w0_c = scaling / (n + scaling) + (1 - alpha * alpha + beta) # weight for first value for covariance
    w1 = 1 / (2 * (n + scaling)) # weight for all other values


    # eq 5-7: sigma point generation
    sigmaPoints = sigma(means, cov, n, scaling)


    # prediction step
    # intertia constants from juwan
    I_body = np.array([[2337899.19, -14882.35, 38212.04],
              [-14882.35, 5112345.28,19754.53],
              [38212.04, 19754.53, 5387496.72]])
    I_body = I_body * 1e-7

    # Reorder moment of inertia matrix. This is as the moment of inertia
    # was calculated in CAD where y is pointing up vertically, while our body frame
    # uses z as pointing up vertically. Transform from CAD to body frame is rotation
    # about X by 90 degrees
    Ixx = I_body[0,0]
    Ixy = I_body[0,1]
    Ixz = I_body[0,2]
    Iyx = I_body[1,0]
    Iyy = I_body[1,1]
    Iyz = I_body[1,2]
    Izx = I_body[2,0]
    Izy = I_body[2,1]
    Izz = I_body[2,2]

    # Define rotated I (given by [R][I][R]^T)
    I_body = np.array([[Ixx, -Ixz, Ixy],
                       [-Izx, Izz, -Izy],
                       [Iyx, -Iyz, Iyy]])

    # Define I_spin and I_trans of RW
    I_spin = 0.319 * 1.82899783e-5 # kg (m^2)
    I_trans = 0

    # intialize 1D EOMs using intertia measurements of cubeSat
    EOMS = DEMO_1D_EOMS(I_body, I_spin, I_trans)
    
    # eq 8-9: pass sigma points through EOMs (f) and generate mean in state space
    predMeans, f = generatePredMeans(EOMS, sigmaPoints, w0_m, w1, reaction_speed, old_reaction_speed, n)
    
    # print("PREDICTED MEANS: ", predMeans)
    
    # eq 10: generate predicted covariance + process noise q
    predCov = generateCov(predMeans, f, w0_c, w1, n, q)

    # print("PRED COVID: ", predCov)


    if len(gps_data) == 4:
        # finds true B field based on gps data
        Bfield = bfield_calc(gps_data)
    else: 
        # for ideal tests only, use gps_data as b field vector and skip calculating it from the gps data
        Bfield = gps_data


    # print("BFIELD: ", Bfield)

    # non linear transformation
    # eq 11-12: non linear transformation of predicted sigma points f into measurement space (h), and mean generation
    mesMeans, h = generateMesMeans(hfunc, Bfield, f, w0_m, w1, n, m)

    # print("MEAN IN MEASUREMENT: ", mesMeans)

    # eq 13: measurement covariance + measurement noise r
    mesCov = generateCov(mesMeans, h, w0_c, w1, n, r)


    # measurement updates
    # eq 14: cross covariance. compare our different sets of sigma points and our predicted/measurement means
    crossCov = generateCrossCov(predMeans, mesMeans, f, h, w0_c, w1, n)

    # print("covariance in measurement: ", mesCov)
    # print("cross covariance: ", crossCov)

    # eq 15: calculate kalman gain (n x m) by multiplying cross covariance matrix and transposed predicted covariance
    kalman = np.matmul(crossCov, np.linalg.inv(mesCov))

    # print("KALMAN: ", kalman)

    # eq 16: updated final mean = predicted + kalman(measurement data - predicted means in measurement space)
    means = np.add(predMeans, np.matmul(kalman, np.subtract(data, mesMeans)))

    # normalize the quaternion to reduce small calculation errors over time

    # eq 17: updated covariance = predicted covariance - kalman * measurement cov * transposed kalman
    cov = np.subtract(predCov, np.matmul(np.matmul(kalman, mesCov), kalman.transpose()))


    # print("MEANS AT END: ", means)
    # print("COV AT END: ", cov)
    return [means, cov]
