import numpy as np
import matplotlib.pyplot as plt
from wmm import WMM

from mpl_toolkits.mplot3d import Axes3D

def hfunc(state, controls):

    quaternion = state[:4]
    rotationMatrix = quaternion_rotation_matrix(quaternion)

    #need to convert this from xyz gps (jiwan has the code already on his person)
    lat = controls[0] 
    long = controls[1]
    height = controls[2] 

    time = controls[3] #formatted as 2023.percentage of the year in month type stuff

    wmm_model = WMM(12, 'WMMcoef.csv')
    wmm_model.calc_gcc_components(lat, long, height, time, degrees=True)
    Bfield1 = wmm_model.get_Bfield()
    return np.concatenate( np.matmul(rotationMatrix,Bfield1), state[4:])



def quaternion_rotation_matrix(Q):
    """
    Covert a quaternion into a full three-dimensional rotation matrix.
 
    Input
    :param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 
 
    Output
    :return: A 3x3 element matrix representing the full 3D rotation matrix. 
             This rotation matrix converts a point in the local reference 
             frame to a point in the global reference frame.
    """
    # Extract the values from Q
    q0 = Q[0]
    q1 = Q[1]
    q2 = Q[2]
    q3 = Q[3]
     
    # First row of the rotation matrix
    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)
     
    # Second row of the rotation matrix
    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)
     
    # Third row of the rotation matrix
    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1
     
    # 3x3 rotation matrix
    rot_matrix = np.array([[r00, r01, r02],
                           [r10, r11, r12],
                           [r20, r21, r22]])
                            
    return rot_matrix



n = 10
state = np.random.rand(n)
transformed = np.array()
# need some kind of generator for random long/lat coordinates, height, and time
#   arrays needed for controls: 
    # lat_gd (np.array): array holding the geodesic latitude associated with a state
    # lon (np.array): array holding the longtitude associated with a state
    # h_ellp (np.array): array holding the estimated heights above the ellipsoid in m
    # t (np.array): array of times associated with an array of states, given in decimal years
controls = [[]]

# Perform observation function, only needed for quaternion components. The rest have 1 to 1 mapping
transformed.append(transformed, np.array(hfunc(state, controls)))

transformed.append(transformed, np.array(state[4:]))

# transformed = np.array([np.array(hfunc(state, q_wmm)), np.array(state[4:])])




if __name__ == '__main__':

    q = np.array([1,0,1,1])
    val = np.linalg.norm(q)
    q = q/val
    rotationMatrix = quaternion_rotation_matrix(q)

    print('quaternion: ',q,'\nrotation matrix: ', rotationMatrix)

    original = [1,0,0]
    rotated = np.matmul(rotationMatrix,original)

    #PLOTTING (DOESNT MATTER)
    original = np.concatenate(([0, 0, 0], original))
    rotated = np.concatenate(([0, 0, 0], rotated))

    print(original)

    soa = np.array([original, rotated])

    X, Y, Z, U, V, W = zip(*soa)

    print(X, Y, Z, U, V, W)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(X, Y, Z, U, V, W)
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    ax.set_zlim([-2, 2])
    plt.show()