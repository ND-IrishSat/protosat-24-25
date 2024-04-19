'''
1D_hfunc.py

Simplified h function relating magnetometer to angle in XY plane 
'''

import numpy as np
import matplotlib.pyplot as plt


def hfunc(state, Bx_true, By_true):
    '''
    hfunc
        transformation from state space to measurement space using magnetic field
    '''

    # find angle in xy plane
    psi = state[0] # rad
    rotationMatrix = np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]])
    
    # calculate measured B-field
    Btrue = np.array([Bx_true, By_true])
    Bmeas = np.matmul(rotationMatrix, Btrue)
    
    return np.append(Bmeas, state[1])

def angle2quat(psi):
    ''' Takes in angle psi (Euler angle about Z axis) and returns quaternion
    '''

    # Build rotation matrix (Z rotation matrix)
    Rmat = np.array([[np.cos(psi), -np.sin(psi), 0],
                     [np.sin(psi), np.cos(psi), 0],
                     [0, 0, 1]])
    
    # Grab components of rotation matrices
    m00 = Rmat[0,0]
    m01 = Rmat[0,1]
    m02 = Rmat[0,2]
    m10 = Rmat[1,0]
    m11 = Rmat[1,1]
    m12 = Rmat[1,2]
    m20 = Rmat[2,0]
    m21 = Rmat[2,1]
    m22 = Rmat[2,2]

    if m00 < -m11:
        t = 1 - m00 - m11 + m22
        quat = np.array([m20 + m02, m12 + m21, t, m01 - m10])
    else:
        t = 1 + m00 + m11 + m22
        quat = np.array([m12 - m21, m20 - m02, m01 - m10, t])
    
    quat = (0.5/np.sqrt(t)) * quat

    return quat