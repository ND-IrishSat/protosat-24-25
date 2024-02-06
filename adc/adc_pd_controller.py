import numpy as np
import matplotlib.pyplot as plt


MAX_SPEED = 8100


def quaternionMultiply(a, b):
    '''
    quaternionMultiply
        custom function to perform quaternion multiply on two passed-in matrices

    @params
        a, b: quaternion matrices (4 x 1)
    @returns
        multiplied quaternion matrix
    '''
    return [[a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]],
            [a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]],
            [a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]],
            [a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]]]

def delta_q(state,target):
    '''
    delta_q
        returns error quaternion by taking quaternion product (x) of current 
        quaternion and inverse of goal quaternion

    @params
        state, target: quaternion matrices (4 x 1) [q0 ; q1:3]
    @returns
        error quaternion
    '''
    # Takes the inverse of the quaternion
    quat_inv_target = [[target[0]],
                       [-target[1]],
                       [-target[2]],
                       [-target[3]]] / ((np.sqrt(target[0]**2 + target[1]**2 + target[2]**2 + target[3]**2) )**2)
    return quaternionMultiply(state,quat_inv_target)

def pd_controller(state,target, omega, kp, kd):
    '''
    
    pd_controller
        proportional derivative controller as described in Fundamentals book pgs 289 -293

    @params
        state: current quaternion of cubesat (4 x 1) [q0 ; q1:3]
        target: goal quaternion of cubesat (4 x 1) [q0 ; q1:3]
        omega: current angular velocity of cubesat
        kp: proportional gain
        kd: derivative gain
    
    @returns
        4 pwm inputs for each of the 4 motors

    '''
    # pwm output for 3 reaction wheels
    delta_q_out = delta_q(state, target)
    pwm = -kp * np.sign(delta_q_out[0]) * delta_q_out[1:3] - kd * omega

    # transformation matrix for NASA configuration
    W =  [[1, 0, 0, 1/np.sqrt(3)],
             [0, 1, 0, 1/np.sqrt(3)],
             [0, 0, 1, 1/np.sqrt(3)]]
    
    # convert output for 3 rx wheels to 4
    pwm = W*pwm

    # Ensure pwm is always within limits
    if pwm > 1:
        pwm = 1
    elif pwm < -1:
        pwm = -1

    # MAYBE: Ensure pwm is always within limits of RPMs
    if pwm > MAX_SPEED:
        pwm = MAX_SPEED
    elif pwm < -MAX_SPEED:
        pwm = -MAX_SPEED

    return pwm