''' 
bigEOMS.py
This module contains the class definition for the improved equations of motion for IrishSat
    
    IrishSat
    Juwan Jeremy Jacobe
    University of Notre Dame 
    14 April 2022
'''
# c g/cm^2

class bigEOMS:
    ''' Equations of motions for the state of the system of ProtoSat, where the state is made up of quaternion components
    a, b, c, d, angular velocities of the body w_x, w_y, w_z, and angular velocities of the reaction wheels theta_dot_RW1, theta_dot_RW2,
    and theta_dot_RW3
    
    '''
    
    def __init__(self):
        pass
        
    def adot(self, a, b, c, d, w_x, w_y, w_z):
        return - (w_x*b)/2 - (w_y*c)/2 - (w_z*d)/2
        
        # The below return statements for adot, bdot, cdot, and ddot were the EoMs according to Wensing
        # with extra terms for "normalizing" the quaternion back to a unit quaternion. However, implementing them 
        # as is leads to the propagator diverging. For now, the terms without those added for v2 will not be used,
        # and normalizing will be done outside the propagator
        #return -a**3 - a*b**2 - a*c**2 - a*d**2 + a - (w_x*b)/2 - (w_y*c)/2 - (w_z*d)/2
        
    def bdot(self, a, b, c, d, w_x, w_y, w_z):  
        return (w_x*a)/2 + (w_z*c)/2 - (w_y*d)/2
        
        #return -a**2*b + (w_x*a)/2 - b**3 - b*c**2 - b*d**2 + b + (w_z*c)/2 - (w_y*d)/2 - 1
        
    def cdot(self, a, b, c, d, w_x, w_y, w_z):  
        return (w_y*a)/2 - (w_z*b)/2 + (w_x*d)/2
        
        #return -a**2*c + (w_y*a)/2 - b**2*c - (w_z*b)/2 - c**3 - c*d**2 + c + (w_x*d)/2 - 1
        
    def ddot(self, a, b, c, d, w_x, w_y, w_z): 
        
        return (w_z*a)/2 + (w_y*b)/2- (w_x*c)/2
        
        #return -a**2*d + (w_z*a)/2 - b**2*d + (w_y*b)/2 - c**2*d - (w_x*c)/2 - d**3 + d - 1
        
    def w_dot_x(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return (I_xy**2*I_yz*w_y**2 - I_xz**2*I_zy*w_z**2 + I_yz*I_zy**2*w_y**2 - I_yz**2*I_zy*w_z**2 + I_RW2_YY*I_RW3_ZZ*M_mt1 + I_RW2_YY*I_xz*M_mt3 + I_RW3_ZZ*I_xy*M_mt2 - I_RW3_ZZ*I_yy*M_mt1 - 
        I_RW2_YY*I_zz*M_mt1 + I_xy*I_yz*M_mt3 - I_xz*I_yy*M_mt3 - I_xy*I_zz*M_mt2 + I_xz*I_zy*M_mt2 + I_yy*I_zz*M_mt1 - I_yz*I_zy*M_mt1 - I_RW2_YY*I_RW3_ZZ*t_motor1 - I_RW2_YY*I_xz*t_motor3 - 
        I_RW3_ZZ*I_xy*t_motor2 + I_RW3_ZZ*I_yy*t_motor1 + I_RW2_YY*I_zz*t_motor1 - I_xy*I_yz*t_motor3 + I_xz*I_yy*t_motor3 + I_xy*I_zz*t_motor2 - I_xz*I_zy*t_motor2 - I_yy*I_zz*t_motor1 + 
        I_yz*I_zy*t_motor1 + I_RW2_YY*I_RW3_ZZ*I_yz*w_z**2 - I_RW2_YY*I_RW3_ZZ*I_zy*w_y**2 + I_RW2_YY*I_xy*I_xz*w_y**2 - I_RW3_ZZ*I_xy*I_xz*w_z**2 - I_RW2_YY*I_xz*I_yx*w_x**2 + I_RW3_ZZ*I_xy*I_zx*w_x**2 - 
        I_RW3_ZZ*I_yy*I_yz*w_z**2 - I_RW2_YY*I_yz*I_zz*w_z**2 + I_RW3_ZZ*I_yy*I_zy*w_y**2 + I_RW2_YY*I_zy*I_zz*w_y**2 - I_xy*I_xz*I_yy*w_y**2 - I_xy*I_yx*I_yz*w_x**2 + I_xz*I_yx*I_yy*w_x**2 + I_xy*I_xz*I_zz*w_z**2 - 
        I_xy*I_zx*I_zz*w_x**2 + I_xz*I_zx*I_zy*w_x**2 + I_yy*I_yz*I_zz*w_z**2 - I_yy*I_zy*I_zz*w_y**2 + I_RW2_YY*I_xz**2*w_y*w_z - I_RW3_ZZ*I_xy**2*w_y*w_z - I_RW3_ZZ*I_yy**2*w_y*w_z + 
        I_RW2_YY*I_zz**2*w_y*w_z + I_xz*I_yy**2*w_x*w_y - I_xy*I_yz**2*w_x*w_z - I_xz**2*I_yy*w_y*w_z + I_xz*I_zy**2*w_x*w_y - I_xy*I_zz**2*w_x*w_z + I_xy**2*I_zz*w_y*w_z - 
        I_yy*I_zz**2*w_y*w_z + I_yy**2*I_zz*w_y*w_z - I_RW2_YY*I_RW3_ZZ**2*w_y*theta_dot_RW3 + I_RW2_YY**2*I_RW3_ZZ*w_z*theta_dot_RW2 - I_RW2_YY**2*I_xz*w_x*theta_dot_RW2 + I_RW3_ZZ**2*I_xy*w_x*theta_dot_RW3 + 
        I_RW3_ZZ**2*I_yy*w_y*theta_dot_RW3 - I_RW2_YY**2*I_zz*w_z*theta_dot_RW2 + I_RW2_YY*I_RW3_ZZ*I_yx*w_x*w_z + I_RW2_YY*I_RW3_ZZ*I_yy*w_y*w_z - I_RW2_YY*I_RW3_ZZ*I_zx*w_x*w_y - 
        I_RW2_YY*I_RW3_ZZ*I_zz*w_y*w_z + I_RW2_YY*I_xx*I_xz*w_x*w_y - I_RW3_ZZ*I_xx*I_xy*w_x*w_z - I_RW2_YY*I_xz*I_yy*w_x*w_y - I_RW2_YY*I_xz*I_yz*w_x*w_z + I_RW3_ZZ*I_xy*I_zy*w_x*w_y - 
        I_RW3_ZZ*I_yx*I_yy*w_x*w_z + I_RW3_ZZ*I_xy*I_zz*w_x*w_z - I_RW2_YY*I_yx*I_zz*w_x*w_z + I_RW3_ZZ*I_yy*I_zx*w_x*w_y - I_RW2_YY*I_yy*I_zz*w_y*w_z + I_RW3_ZZ*I_yy*I_zz*w_y*w_z + 
        I_RW2_YY*I_zx*I_zz*w_x*w_y + I_xx*I_xy*I_yz*w_x*w_y - I_xx*I_xz*I_yy*w_x*w_y + I_xy*I_xz*I_yz*w_y*w_z + I_xx*I_xy*I_zz*w_x*w_z - I_xx*I_xz*I_zy*w_x*w_z - I_xy*I_yy*I_yz*w_x*w_y - 
        I_xy*I_xz*I_zy*w_y*w_z + I_xz*I_yy*I_yz*w_x*w_z - I_xy*I_zy*I_zz*w_x*w_y + I_yx*I_yy*I_zz*w_x*w_z - I_yx*I_yz*I_zy*w_x*w_z + I_xz*I_zy*I_zz*w_x*w_z - I_yy*I_yz*I_zy*w_y*w_z - 
        I_yy*I_zx*I_zz*w_x*w_y + I_yz*I_zx*I_zy*w_x*w_y + I_yz*I_zy*I_zz*w_y*w_z + I_RW1_XX*I_RW2_YY*I_xz*w_y*theta_dot_RW1 - I_RW1_XX*I_RW3_ZZ*I_xy*w_z*theta_dot_RW1 - I_RW2_YY*I_RW3_ZZ*I_yy*w_z*theta_dot_RW2 + 
        I_RW2_YY*I_RW3_ZZ*I_zz*w_y*theta_dot_RW3 + I_RW1_XX*I_xy*I_yz*w_y*theta_dot_RW1 - I_RW1_XX*I_xz*I_yy*w_y*theta_dot_RW1 - I_RW2_YY*I_xy*I_yz*w_x*theta_dot_RW2 + I_RW2_YY*I_xz*I_yy*w_x*theta_dot_RW2 + 
        I_RW1_XX*I_xy*I_zz*w_z*theta_dot_RW1 - I_RW1_XX*I_xz*I_zy*w_z*theta_dot_RW1 - I_RW3_ZZ*I_xy*I_zz*w_x*theta_dot_RW3 + I_RW3_ZZ*I_xz*I_zy*w_x*theta_dot_RW3 + I_RW2_YY*I_yy*I_zz*w_z*theta_dot_RW2 - 
        I_RW2_YY*I_yz*I_zy*w_z*theta_dot_RW2 - I_RW3_ZZ*I_yy*I_zz*w_y*theta_dot_RW3 + I_RW3_ZZ*I_yz*I_zy*w_y*theta_dot_RW3)/(I_RW2_YY*I_RW3_ZZ*I_xx - I_RW1_XX*I_RW2_YY*I_RW3_ZZ + I_RW1_XX*I_RW3_ZZ*I_yy + 
        I_RW1_XX*I_RW2_YY*I_zz - I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + I_xx*I_yy*I_zz - I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + I_xy*I_yz*I_zx + 
        I_xz*I_yx*I_zy - I_xz*I_yy*I_zx)
        
    def w_dot_y(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return -(I_xz*I_yx**2*w_x**2 + I_xz*I_zx**2*w_x**2 - I_xz**2*I_zx*w_z**2 - I_yz**2*I_zx*w_z**2 - I_RW1_XX*I_RW3_ZZ*M_mt2 + I_RW3_ZZ*I_xx*M_mt2 - I_RW1_XX*I_yz*M_mt3 - I_RW3_ZZ*I_yx*M_mt1 + I_RW1_XX*I_zz*M_mt2 + 
        I_xx*I_yz*M_mt3 - I_xz*I_yx*M_mt3 - I_xx*I_zz*M_mt2 + I_xz*I_zx*M_mt2 + I_yx*I_zz*M_mt1 - I_yz*I_zx*M_mt1 + I_RW1_XX*I_RW3_ZZ*t_motor2 - I_RW3_ZZ*I_xx*t_motor2 + I_RW1_XX*I_yz*t_motor3 + 
        I_RW3_ZZ*I_yx*t_motor1 - I_RW1_XX*I_zz*t_motor2 - I_xx*I_yz*t_motor3 + I_xz*I_yx*t_motor3 + I_xx*I_zz*t_motor2 - I_xz*I_zx*t_motor2 - I_yx*I_zz*t_motor1 + I_yz*I_zx*t_motor1 + I_RW1_XX*I_RW3_ZZ*I_xz*w_z**2 - 
        I_RW1_XX*I_RW3_ZZ*I_zx*w_x**2 - I_RW3_ZZ*I_xx*I_xz*w_z**2 - I_RW1_XX*I_xy*I_yz*w_y**2 + I_RW1_XX*I_yx*I_yz*w_x**2 - I_RW1_XX*I_xz*I_zz*w_z**2 + I_RW3_ZZ*I_xx*I_zx*w_x**2 - I_RW3_ZZ*I_yx*I_yz*w_z**2 + 
        I_RW3_ZZ*I_yx*I_zy*w_y**2 + I_RW1_XX*I_zx*I_zz*w_x**2 + I_xx*I_xy*I_yz*w_y**2 - I_xy*I_xz*I_yx*w_y**2 - I_xx*I_yx*I_yz*w_x**2 + I_xx*I_xz*I_zz*w_z**2 - I_xx*I_zx*I_zz*w_x**2 + I_yx*I_yz*I_zz*w_z**2 - 
        I_yx*I_zy*I_zz*w_y**2 + I_yz*I_zx*I_zy*w_y**2 - I_RW3_ZZ*I_xx**2*w_x*w_z + I_RW1_XX*I_yz**2*w_x*w_z - I_RW3_ZZ*I_yx**2*w_x*w_z + I_RW1_XX*I_zz**2*w_x*w_z + I_xx**2*I_yz*w_x*w_y - 
        I_xx*I_yz**2*w_x*w_z - I_xz**2*I_yx*w_y*w_z - I_xx*I_zz**2*w_x*w_z + I_xx**2*I_zz*w_x*w_z + I_yz*I_zx**2*w_x*w_y + I_yx**2*I_zz*w_x*w_z - I_yx*I_zz**2*w_y*w_z - I_RW1_XX*I_RW3_ZZ**2*w_x*theta_dot_RW3 
        + I_RW1_XX**2*I_RW3_ZZ*w_z*theta_dot_RW1 + I_RW3_ZZ**2*I_xx*w_x*theta_dot_RW3 - I_RW1_XX**2*I_yz*w_y*theta_dot_RW1 + I_RW3_ZZ**2*I_yx*w_y*theta_dot_RW3 - I_RW1_XX**2*I_zz*w_z*theta_dot_RW1 + I_RW1_XX*I_RW3_ZZ*I_xx*w_x*w_z + 
        I_RW1_XX*I_RW3_ZZ*I_xy*w_y*w_z - I_RW1_XX*I_RW3_ZZ*I_zy*w_x*w_y - I_RW1_XX*I_RW3_ZZ*I_zz*w_x*w_z - I_RW3_ZZ*I_xx*I_xy*w_y*w_z - I_RW1_XX*I_xx*I_yz*w_x*w_y - I_RW1_XX*I_xz*I_yz*w_y*w_z - 
        I_RW1_XX*I_xx*I_zz*w_x*w_z + I_RW1_XX*I_yy*I_yz*w_x*w_y - I_RW1_XX*I_xy*I_zz*w_y*w_z + I_RW3_ZZ*I_xx*I_zy*w_x*w_y + I_RW3_ZZ*I_xx*I_zz*w_x*w_z - I_RW3_ZZ*I_yx*I_yy*w_y*w_z + 
        I_RW3_ZZ*I_yx*I_zx*w_x*w_y + I_RW3_ZZ*I_yx*I_zz*w_y*w_z + I_RW1_XX*I_zy*I_zz*w_x*w_y - I_xx*I_xz*I_yx*w_x*w_y + I_xx*I_xz*I_yz*w_y*w_z - I_xx*I_xz*I_zx*w_x*w_z - I_xx*I_yy*I_yz*w_x*w_y + 
        I_xz*I_yx*I_yy*w_x*w_y + I_xx*I_xy*I_zz*w_y*w_z - I_xy*I_xz*I_zx*w_y*w_z + I_xz*I_yx*I_yz*w_x*w_z - I_xx*I_zy*I_zz*w_x*w_y + I_xz*I_zx*I_zy*w_x*w_y - I_yx*I_yz*I_zx*w_x*w_z + 
        I_xz*I_zx*I_zz*w_x*w_z + I_yx*I_yy*I_zz*w_y*w_z - I_yy*I_yz*I_zx*w_y*w_z - I_yx*I_zx*I_zz*w_x*w_y + I_yz*I_zx*I_zz*w_y*w_z - I_RW1_XX*I_RW3_ZZ*I_xx*w_z*theta_dot_RW1 + I_RW1_XX*I_RW2_YY*I_yz*w_x*theta_dot_RW2 
        - I_RW2_YY*I_RW3_ZZ*I_yx*w_z*theta_dot_RW2 + I_RW1_XX*I_RW3_ZZ*I_zz*w_x*theta_dot_RW3 + I_RW1_XX*I_xx*I_yz*w_y*theta_dot_RW1 - I_RW1_XX*I_xz*I_yx*w_y*theta_dot_RW1 - I_RW2_YY*I_xx*I_yz*w_x*theta_dot_RW2 + 
        I_RW2_YY*I_xz*I_yx*w_x*theta_dot_RW2 + I_RW1_XX*I_xx*I_zz*w_z*theta_dot_RW1 - I_RW1_XX*I_xz*I_zx*w_z*theta_dot_RW1 - I_RW3_ZZ*I_xx*I_zz*w_x*theta_dot_RW3 + I_RW3_ZZ*I_xz*I_zx*w_x*theta_dot_RW3 + 
        I_RW2_YY*I_yx*I_zz*w_z*theta_dot_RW2 - I_RW2_YY*I_yz*I_zx*w_z*theta_dot_RW2 - I_RW3_ZZ*I_yx*I_zz*w_y*theta_dot_RW3 + I_RW3_ZZ*I_yz*I_zx*w_y*theta_dot_RW3)/(I_RW2_YY*I_RW3_ZZ*I_xx - I_RW1_XX*I_RW2_YY*I_RW3_ZZ + 
        I_RW1_XX*I_RW3_ZZ*I_yy + I_RW1_XX*I_RW2_YY*I_zz - I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + I_xx*I_yy*I_zz - I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + 
        I_xy*I_yz*I_zx + I_xz*I_yx*I_zy - I_xz*I_yy*I_zx)
        
    def w_dot_z(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return (I_xy*I_yx**2*w_x**2 - I_xy**2*I_yx*w_y**2 + I_xy*I_zx**2*w_x**2 - I_yx*I_zy**2*w_y**2 + I_RW1_XX*I_RW2_YY*M_mt3 - I_RW2_YY*I_xx*M_mt3 - I_RW1_XX*I_yy*M_mt3 + I_RW1_XX*I_zy*M_mt2 + I_RW2_YY*I_zx*M_mt1 + 
        I_xx*I_yy*M_mt3 - I_xy*I_yx*M_mt3 - I_xx*I_zy*M_mt2 + I_xy*I_zx*M_mt2 + I_yx*I_zy*M_mt1 - I_yy*I_zx*M_mt1 - I_RW1_XX*I_RW2_YY*t_motor3 + I_RW2_YY*I_xx*t_motor3 + I_RW1_XX*I_yy*t_motor3 - I_RW1_XX*I_zy*t_motor2 - 
        I_RW2_YY*I_zx*t_motor1 - I_xx*I_yy*t_motor3 + I_xy*I_yx*t_motor3 + I_xx*I_zy*t_motor2 - I_xy*I_zx*t_motor2 - I_yx*I_zy*t_motor1 + I_yy*I_zx*t_motor1 + I_RW1_XX*I_RW2_YY*I_xy*w_y**2 - I_RW1_XX*I_RW2_YY*I_yx*w_x**2 - 
        I_RW2_YY*I_xx*I_xy*w_y**2 - I_RW1_XX*I_xy*I_yy*w_y**2 + I_RW2_YY*I_xx*I_yx*w_x**2 + I_RW1_XX*I_yx*I_yy*w_x**2 - I_RW1_XX*I_xz*I_zy*w_z**2 + I_RW2_YY*I_yz*I_zx*w_z**2 + I_RW1_XX*I_zx*I_zy*w_x**2 - I_RW2_YY*I_zx*I_zy*w_y**2 + 
        I_xx*I_xy*I_yy*w_y**2 - I_xx*I_yx*I_yy*w_x**2 + I_xx*I_xz*I_zy*w_z**2 - I_xy*I_xz*I_zx*w_z**2 - I_xx*I_zx*I_zy*w_x**2 + I_yx*I_yz*I_zy*w_z**2 - I_yy*I_yz*I_zx*w_z**2 + I_yy*I_zx*I_zy*w_y**2 - I_RW2_YY*I_xx**2*w_x*w_y + 
        I_RW1_XX*I_yy**2*w_x*w_y + I_RW1_XX*I_zy**2*w_x*w_y - I_RW2_YY*I_zx**2*w_x*w_y - I_xx*I_yy**2*w_x*w_y + I_xx**2*I_yy*w_x*w_y - I_xx*I_zy**2*w_x*w_y + I_xx**2*I_zy*w_x*w_z - I_xy**2*I_zx*w_y*w_z + 
        I_yy*I_zx**2*w_x*w_y + I_yx**2*I_zy*w_x*w_z - I_yy**2*I_zx*w_y*w_z - I_RW1_XX*I_RW2_YY**2*w_x*theta_dot_RW2 + I_RW1_XX**2*I_RW2_YY*w_y*theta_dot_RW1 + I_RW2_YY**2*I_xx*w_x*theta_dot_RW2 - 
        I_RW1_XX**2*I_yy*w_y*theta_dot_RW1 - I_RW1_XX**2*I_zy*w_z*theta_dot_RW1 + I_RW2_YY**2*I_zx*w_z*theta_dot_RW2 + I_RW1_XX*I_RW2_YY*I_xx*w_x*w_y + I_RW1_XX*I_RW2_YY*I_xz*w_y*w_z - I_RW1_XX*I_RW2_YY*I_yy*w_x*w_y - 
        I_RW1_XX*I_RW2_YY*I_yz*w_x*w_z - I_RW2_YY*I_xx*I_xz*w_y*w_z - I_RW1_XX*I_xx*I_yy*w_x*w_y + I_RW2_YY*I_xx*I_yy*w_x*w_y - I_RW1_XX*I_xz*I_yy*w_y*w_z + I_RW2_YY*I_xx*I_yz*w_x*w_z - 
        I_RW1_XX*I_xx*I_zy*w_x*w_z - I_RW1_XX*I_xy*I_zy*w_y*w_z + I_RW1_XX*I_yy*I_yz*w_x*w_z + I_RW2_YY*I_yx*I_zx*w_x*w_z + I_RW2_YY*I_yy*I_zx*w_y*w_z + I_RW1_XX*I_zy*I_zz*w_x*w_z - 
        I_RW2_YY*I_zx*I_zz*w_y*w_z - I_xx*I_xy*I_yx*w_x*w_y + I_xx*I_xz*I_yy*w_y*w_z - I_xy*I_xz*I_yx*w_y*w_z - I_xx*I_xy*I_zx*w_x*w_z + I_xy*I_yx*I_yy*w_x*w_y + I_xx*I_xy*I_zy*w_y*w_z - 
        I_xx*I_yy*I_yz*w_x*w_z + I_xy*I_yx*I_yz*w_x*w_z + I_xy*I_zx*I_zy*w_x*w_y - I_yx*I_yy*I_zx*w_x*w_z - I_xx*I_zy*I_zz*w_x*w_z + I_xy*I_zx*I_zz*w_x*w_z + I_yx*I_yy*I_zy*w_y*w_z - 
        I_yx*I_zx*I_zy*w_x*w_y - I_yx*I_zy*I_zz*w_y*w_z + I_yy*I_zx*I_zz*w_y*w_z - I_RW1_XX*I_RW2_YY*I_xx*w_y*theta_dot_RW1 + I_RW1_XX*I_RW2_YY*I_yy*w_x*theta_dot_RW2 + I_RW1_XX*I_RW3_ZZ*I_zy*w_x*theta_dot_RW3 - 
        I_RW2_YY*I_RW3_ZZ*I_zx*w_y*theta_dot_RW3 + I_RW1_XX*I_xx*I_yy*w_y*theta_dot_RW1 - I_RW1_XX*I_xy*I_yx*w_y*theta_dot_RW1 - I_RW2_YY*I_xx*I_yy*w_x*theta_dot_RW2 + I_RW2_YY*I_xy*I_yx*w_x*theta_dot_RW2 + 
        I_RW1_XX*I_xx*I_zy*w_z*theta_dot_RW1 - I_RW1_XX*I_xy*I_zx*w_z*theta_dot_RW1 - I_RW3_ZZ*I_xx*I_zy*w_x*theta_dot_RW3 + I_RW3_ZZ*I_xy*I_zx*w_x*theta_dot_RW3 + I_RW2_YY*I_yx*I_zy*w_z*theta_dot_RW2 - 
        I_RW2_YY*I_yy*I_zx*w_z*theta_dot_RW2 - I_RW3_ZZ*I_yx*I_zy*w_y*theta_dot_RW3 + I_RW3_ZZ*I_yy*I_zx*w_y*theta_dot_RW3)/(I_RW2_YY*I_RW3_ZZ*I_xx - I_RW1_XX*I_RW2_YY*I_RW3_ZZ + I_RW1_XX*I_RW3_ZZ*I_yy + I_RW1_XX*I_RW2_YY*I_zz - 
        I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + I_xx*I_yy*I_zz - I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + I_xy*I_yz*I_zx + I_xz*I_yx*I_zy - I_xz*I_yy*I_zx)
        
    def theta_ddot_RW1(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return (I_RW1_XX*I_xz**2*I_zy*w_z**2 - I_RW1_XX*I_xy**2*I_yz*w_y**2 - I_RW1_XX*I_yz*I_zy**2*w_y**2 + I_RW1_XX*I_yz**2*I_zy*w_z**2 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*M_mt1 - I_RW1_XX*I_RW2_YY*I_xz*M_mt3 - I_RW1_XX*I_RW3_ZZ*I_xy*M_mt2 + 
        I_RW1_XX*I_RW3_ZZ*I_yy*M_mt1 + I_RW1_XX*I_RW2_YY*I_zz*M_mt1 - I_RW1_XX*I_xy*I_yz*M_mt3 + I_RW1_XX*I_xz*I_yy*M_mt3 + I_RW1_XX*I_xy*I_zz*M_mt2 - I_RW1_XX*I_xz*I_zy*M_mt2 - I_RW1_XX*I_yy*I_zz*M_mt1 + 
        I_RW1_XX*I_yz*I_zy*M_mt1 + I_RW1_XX*I_RW2_YY*I_xz*t_motor3 + I_RW1_XX*I_RW3_ZZ*I_xy*t_motor2 + I_RW2_YY*I_RW3_ZZ*I_xx*t_motor1 + I_RW1_XX*I_xy*I_yz*t_motor3 - I_RW1_XX*I_xz*I_yy*t_motor3 - I_RW3_ZZ*I_xx*I_yy*t_motor1 + 
        I_RW3_ZZ*I_xy*I_yx*t_motor1 - I_RW1_XX*I_xy*I_zz*t_motor2 + I_RW1_XX*I_xz*I_zy*t_motor2 - I_RW2_YY*I_xx*I_zz*t_motor1 + I_RW2_YY*I_xz*I_zx*t_motor1 + I_xx*I_yy*I_zz*t_motor1 - I_xx*I_yz*I_zy*t_motor1 - I_xy*I_yx*I_zz*t_motor1 + 
        I_xy*I_yz*I_zx*t_motor1 + I_xz*I_yx*I_zy*t_motor1 - I_xz*I_yy*I_zx*t_motor1 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yz*w_z**2 + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zy*w_y**2 - I_RW1_XX*I_RW2_YY*I_xy*I_xz*w_y**2 + I_RW1_XX*I_RW3_ZZ*I_xy*I_xz*w_z**2 + 
        I_RW1_XX*I_RW2_YY*I_xz*I_yx*w_x**2 - I_RW1_XX*I_RW3_ZZ*I_xy*I_zx*w_x**2 + I_RW1_XX*I_RW3_ZZ*I_yy*I_yz*w_z**2 + I_RW1_XX*I_RW2_YY*I_yz*I_zz*w_z**2 - I_RW1_XX*I_RW3_ZZ*I_yy*I_zy*w_y**2 - I_RW1_XX*I_RW2_YY*I_zy*I_zz*w_y**2 + 
        I_RW1_XX*I_xy*I_xz*I_yy*w_y**2 + I_RW1_XX*I_xy*I_yx*I_yz*w_x**2 - I_RW1_XX*I_xz*I_yx*I_yy*w_x**2 - I_RW1_XX*I_xy*I_xz*I_zz*w_z**2 + I_RW1_XX*I_xy*I_zx*I_zz*w_x**2 - I_RW1_XX*I_xz*I_zx*I_zy*w_x**2 - 
        I_RW1_XX*I_yy*I_yz*I_zz*w_z**2 + I_RW1_XX*I_yy*I_zy*I_zz*w_y**2 - I_RW1_XX*I_RW2_YY*I_xz**2*w_y*w_z + I_RW1_XX*I_RW3_ZZ*I_xy**2*w_y*w_z + I_RW1_XX*I_RW3_ZZ*I_yy**2*w_y*w_z - I_RW1_XX*I_RW2_YY*I_zz**2*w_y*w_z - 
        I_RW1_XX*I_xz*I_yy**2*w_x*w_y + I_RW1_XX*I_xy*I_yz**2*w_x*w_z + I_RW1_XX*I_xz**2*I_yy*w_y*w_z - I_RW1_XX*I_xz*I_zy**2*w_x*w_y + I_RW1_XX*I_xy*I_zz**2*w_x*w_z - I_RW1_XX*I_xy**2*I_zz*w_y*w_z + 
        I_RW1_XX*I_yy*I_zz**2*w_y*w_z - I_RW1_XX*I_yy**2*I_zz*w_y*w_z + I_RW1_XX*I_RW2_YY*I_RW3_ZZ**2*w_y*theta_dot_RW3 - I_RW1_XX*I_RW2_YY**2*I_RW3_ZZ*w_z*theta_dot_RW2 + I_RW1_XX*I_RW2_YY**2*I_xz*w_x*theta_dot_RW2 - 
        I_RW1_XX**2*I_RW2_YY*I_xz*w_y*theta_dot_RW1 - I_RW1_XX*I_RW3_ZZ**2*I_xy*w_x*theta_dot_RW3 + I_RW1_XX**2*I_RW3_ZZ*I_xy*w_z*theta_dot_RW1 - I_RW1_XX*I_RW3_ZZ**2*I_yy*w_y*theta_dot_RW3 + I_RW1_XX*I_RW2_YY**2*I_zz*w_z*theta_dot_RW2 - 
        I_RW1_XX**2*I_xy*I_yz*w_y*theta_dot_RW1 + I_RW1_XX**2*I_xz*I_yy*w_y*theta_dot_RW1 - I_RW1_XX**2*I_xy*I_zz*w_z*theta_dot_RW1 + I_RW1_XX**2*I_xz*I_zy*w_z*theta_dot_RW1 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yx*w_x*w_z - 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yy*w_y*w_z + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zx*w_x*w_y + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zz*w_y*w_z - I_RW1_XX*I_RW2_YY*I_xx*I_xz*w_x*w_y + I_RW1_XX*I_RW3_ZZ*I_xx*I_xy*w_x*w_z + 
        I_RW1_XX*I_RW2_YY*I_xz*I_yy*w_x*w_y + I_RW1_XX*I_RW2_YY*I_xz*I_yz*w_x*w_z - I_RW1_XX*I_RW3_ZZ*I_xy*I_zy*w_x*w_y + I_RW1_XX*I_RW3_ZZ*I_yx*I_yy*w_x*w_z - I_RW1_XX*I_RW3_ZZ*I_xy*I_zz*w_x*w_z + 
        I_RW1_XX*I_RW2_YY*I_yx*I_zz*w_x*w_z - I_RW1_XX*I_RW3_ZZ*I_yy*I_zx*w_x*w_y + I_RW1_XX*I_RW2_YY*I_yy*I_zz*w_y*w_z - I_RW1_XX*I_RW3_ZZ*I_yy*I_zz*w_y*w_z - I_RW1_XX*I_RW2_YY*I_zx*I_zz*w_x*w_y - 
        I_RW1_XX*I_xx*I_xy*I_yz*w_x*w_y + I_RW1_XX*I_xx*I_xz*I_yy*w_x*w_y - I_RW1_XX*I_xy*I_xz*I_yz*w_y*w_z - I_RW1_XX*I_xx*I_xy*I_zz*w_x*w_z + I_RW1_XX*I_xx*I_xz*I_zy*w_x*w_z + I_RW1_XX*I_xy*I_yy*I_yz*w_x*w_y + 
        I_RW1_XX*I_xy*I_xz*I_zy*w_y*w_z - I_RW1_XX*I_xz*I_yy*I_yz*w_x*w_z + I_RW1_XX*I_xy*I_zy*I_zz*w_x*w_y - I_RW1_XX*I_yx*I_yy*I_zz*w_x*w_z + I_RW1_XX*I_yx*I_yz*I_zy*w_x*w_z - I_RW1_XX*I_xz*I_zy*I_zz*w_x*w_z + 
        I_RW1_XX*I_yy*I_yz*I_zy*w_y*w_z + I_RW1_XX*I_yy*I_zx*I_zz*w_x*w_y - I_RW1_XX*I_yz*I_zx*I_zy*w_x*w_y - I_RW1_XX*I_yz*I_zy*I_zz*w_y*w_z + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yy*w_z*theta_dot_RW2 - 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zz*w_y*theta_dot_RW3 + I_RW1_XX*I_RW2_YY*I_xy*I_yz*w_x*theta_dot_RW2 - I_RW1_XX*I_RW2_YY*I_xz*I_yy*w_x*theta_dot_RW2 + I_RW1_XX*I_RW3_ZZ*I_xy*I_zz*w_x*theta_dot_RW3 - I_RW1_XX*I_RW3_ZZ*I_xz*I_zy*w_x*theta_dot_RW3 - 
        I_RW1_XX*I_RW2_YY*I_yy*I_zz*w_z*theta_dot_RW2 + I_RW1_XX*I_RW2_YY*I_yz*I_zy*w_z*theta_dot_RW2 + I_RW1_XX*I_RW3_ZZ*I_yy*I_zz*w_y*theta_dot_RW3 - I_RW1_XX*I_RW3_ZZ*I_yz*I_zy*w_y*theta_dot_RW3)/(I_RW1_XX*(I_RW2_YY*I_RW3_ZZ*I_xx - 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ + I_RW1_XX*I_RW3_ZZ*I_yy + I_RW1_XX*I_RW2_YY*I_zz - I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + I_xx*I_yy*I_zz - 
        I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + I_xy*I_yz*I_zx + I_xz*I_yx*I_zy - I_xz*I_yy*I_zx))
        
    def theta_ddot_RW2(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return (I_RW2_YY*I_xz*I_yx**2*w_x**2 + I_RW2_YY*I_xz*I_zx**2*w_x**2 - I_RW2_YY*I_xz**2*I_zx*w_z**2 - I_RW2_YY*I_yz**2*I_zx*w_z**2 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*M_mt2 + I_RW2_YY*I_RW3_ZZ*I_xx*M_mt2 - I_RW1_XX*I_RW2_YY*I_yz*M_mt3 - 
        I_RW2_YY*I_RW3_ZZ*I_yx*M_mt1 + I_RW1_XX*I_RW2_YY*I_zz*M_mt2 + I_RW2_YY*I_xx*I_yz*M_mt3 - I_RW2_YY*I_xz*I_yx*M_mt3 - I_RW2_YY*I_xx*I_zz*M_mt2 + I_RW2_YY*I_xz*I_zx*M_mt2 + I_RW2_YY*I_yx*I_zz*M_mt1 - 
        I_RW2_YY*I_yz*I_zx*M_mt1 + I_RW1_XX*I_RW2_YY*I_yz*t_motor3 + I_RW1_XX*I_RW3_ZZ*I_yy*t_motor2 + I_RW2_YY*I_RW3_ZZ*I_yx*t_motor1 - I_RW2_YY*I_xx*I_yz*t_motor3 + I_RW2_YY*I_xz*I_yx*t_motor3 - I_RW3_ZZ*I_xx*I_yy*t_motor2 + 
        I_RW3_ZZ*I_xy*I_yx*t_motor2 - I_RW1_XX*I_yy*I_zz*t_motor2 + I_RW1_XX*I_yz*I_zy*t_motor2 - I_RW2_YY*I_yx*I_zz*t_motor1 + I_RW2_YY*I_yz*I_zx*t_motor1 + I_xx*I_yy*I_zz*t_motor2 - I_xx*I_yz*I_zy*t_motor2 - I_xy*I_yx*I_zz*t_motor2 + 
        I_xy*I_yz*I_zx*t_motor2 + I_xz*I_yx*I_zy*t_motor2 - I_xz*I_yy*I_zx*t_motor2 + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xz*w_z**2 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zx*w_x**2 - I_RW2_YY*I_RW3_ZZ*I_xx*I_xz*w_z**2 - I_RW1_XX*I_RW2_YY*I_xy*I_yz*w_y**2 + 
        I_RW1_XX*I_RW2_YY*I_yx*I_yz*w_x**2 - I_RW1_XX*I_RW2_YY*I_xz*I_zz*w_z**2 + I_RW2_YY*I_RW3_ZZ*I_xx*I_zx*w_x**2 - I_RW2_YY*I_RW3_ZZ*I_yx*I_yz*w_z**2 + I_RW2_YY*I_RW3_ZZ*I_yx*I_zy*w_y**2 + I_RW1_XX*I_RW2_YY*I_zx*I_zz*w_x**2 + 
        I_RW2_YY*I_xx*I_xy*I_yz*w_y**2 - I_RW2_YY*I_xy*I_xz*I_yx*w_y**2 - I_RW2_YY*I_xx*I_yx*I_yz*w_x**2 + I_RW2_YY*I_xx*I_xz*I_zz*w_z**2 - I_RW2_YY*I_xx*I_zx*I_zz*w_x**2 + I_RW2_YY*I_yx*I_yz*I_zz*w_z**2 - I_RW2_YY*I_yx*I_zy*I_zz*w_y**2 + 
        I_RW2_YY*I_yz*I_zx*I_zy*w_y**2 - I_RW2_YY*I_RW3_ZZ*I_xx**2*w_x*w_z + I_RW1_XX*I_RW2_YY*I_yz**2*w_x*w_z - I_RW2_YY*I_RW3_ZZ*I_yx**2*w_x*w_z + I_RW1_XX*I_RW2_YY*I_zz**2*w_x*w_z + I_RW2_YY*I_xx**2*I_yz*w_x*w_y - 
        I_RW2_YY*I_xx*I_yz**2*w_x*w_z - I_RW2_YY*I_xz**2*I_yx*w_y*w_z - I_RW2_YY*I_xx*I_zz**2*w_x*w_z + I_RW2_YY*I_xx**2*I_zz*w_x*w_z + I_RW2_YY*I_yz*I_zx**2*w_x*w_y + I_RW2_YY*I_yx**2*I_zz*w_x*w_z - 
        I_RW2_YY*I_yx*I_zz**2*w_y*w_z - I_RW1_XX*I_RW2_YY*I_RW3_ZZ**2*w_x*theta_dot_RW3 + I_RW1_XX**2*I_RW2_YY*I_RW3_ZZ*w_z*theta_dot_RW1 + I_RW2_YY*I_RW3_ZZ**2*I_xx*w_x*theta_dot_RW3 + I_RW1_XX*I_RW2_YY**2*I_yz*w_x*theta_dot_RW2 - 
        I_RW1_XX**2*I_RW2_YY*I_yz*w_y*theta_dot_RW1 + I_RW2_YY*I_RW3_ZZ**2*I_yx*w_y*theta_dot_RW3 - I_RW2_YY**2*I_RW3_ZZ*I_yx*w_z*theta_dot_RW2 - I_RW1_XX**2*I_RW2_YY*I_zz*w_z*theta_dot_RW1 - I_RW2_YY**2*I_xx*I_yz*w_x*theta_dot_RW2 + 
        I_RW2_YY**2*I_xz*I_yx*w_x*theta_dot_RW2 + I_RW2_YY**2*I_yx*I_zz*w_z*theta_dot_RW2 - I_RW2_YY**2*I_yz*I_zx*w_z*theta_dot_RW2 + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xx*w_x*w_z + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xy*w_y*w_z - 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zy*w_x*w_y - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zz*w_x*w_z - I_RW2_YY*I_RW3_ZZ*I_xx*I_xy*w_y*w_z - I_RW1_XX*I_RW2_YY*I_xx*I_yz*w_x*w_y - I_RW1_XX*I_RW2_YY*I_xz*I_yz*w_y*w_z - 
        I_RW1_XX*I_RW2_YY*I_xx*I_zz*w_x*w_z + I_RW1_XX*I_RW2_YY*I_yy*I_yz*w_x*w_y - I_RW1_XX*I_RW2_YY*I_xy*I_zz*w_y*w_z + I_RW2_YY*I_RW3_ZZ*I_xx*I_zy*w_x*w_y + I_RW2_YY*I_RW3_ZZ*I_xx*I_zz*w_x*w_z - 
        I_RW2_YY*I_RW3_ZZ*I_yx*I_yy*w_y*w_z + I_RW2_YY*I_RW3_ZZ*I_yx*I_zx*w_x*w_y + I_RW2_YY*I_RW3_ZZ*I_yx*I_zz*w_y*w_z + I_RW1_XX*I_RW2_YY*I_zy*I_zz*w_x*w_y - I_RW2_YY*I_xx*I_xz*I_yx*w_x*w_y + 
        I_RW2_YY*I_xx*I_xz*I_yz*w_y*w_z - I_RW2_YY*I_xx*I_xz*I_zx*w_x*w_z - I_RW2_YY*I_xx*I_yy*I_yz*w_x*w_y + I_RW2_YY*I_xz*I_yx*I_yy*w_x*w_y + I_RW2_YY*I_xx*I_xy*I_zz*w_y*w_z - I_RW2_YY*I_xy*I_xz*I_zx*w_y*w_z + 
        I_RW2_YY*I_xz*I_yx*I_yz*w_x*w_z - I_RW2_YY*I_xx*I_zy*I_zz*w_x*w_y + I_RW2_YY*I_xz*I_zx*I_zy*w_x*w_y - I_RW2_YY*I_yx*I_yz*I_zx*w_x*w_z + I_RW2_YY*I_xz*I_zx*I_zz*w_x*w_z + I_RW2_YY*I_yx*I_yy*I_zz*w_y*w_z - 
        I_RW2_YY*I_yy*I_yz*I_zx*w_y*w_z - I_RW2_YY*I_yx*I_zx*I_zz*w_x*w_y + I_RW2_YY*I_yz*I_zx*I_zz*w_y*w_z - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xx*w_z*theta_dot_RW1 + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_zz*w_x*theta_dot_RW3 + 
        I_RW1_XX*I_RW2_YY*I_xx*I_yz*w_y*theta_dot_RW1 - I_RW1_XX*I_RW2_YY*I_xz*I_yx*w_y*theta_dot_RW1 + I_RW1_XX*I_RW2_YY*I_xx*I_zz*w_z*theta_dot_RW1 - I_RW1_XX*I_RW2_YY*I_xz*I_zx*w_z*theta_dot_RW1 - 
        I_RW2_YY*I_RW3_ZZ*I_xx*I_zz*w_x*theta_dot_RW3 + I_RW2_YY*I_RW3_ZZ*I_xz*I_zx*w_x*theta_dot_RW3 - I_RW2_YY*I_RW3_ZZ*I_yx*I_zz*w_y*theta_dot_RW3 + I_RW2_YY*I_RW3_ZZ*I_yz*I_zx*w_y*theta_dot_RW3)/(I_RW2_YY*(I_RW2_YY*I_RW3_ZZ*I_xx - 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ + I_RW1_XX*I_RW3_ZZ*I_yy + I_RW1_XX*I_RW2_YY*I_zz - I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + 
        I_xx*I_yy*I_zz - I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + I_xy*I_yz*I_zx + I_xz*I_yx*I_zy - I_xz*I_yy*I_zx))
        
    def theta_ddot_RW3(self, M_mt1, M_mt2, M_mt3, t_motor1, t_motor2, t_motor3, w_x, w_y, w_z, I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz, I_RW1_XX, I_RW2_YY, 
        I_RW3_ZZ, theta_dot_RW1, theta_dot_RW2, theta_dot_RW3):
        
        return (I_RW3_ZZ*I_xy**2*I_yx*w_y**2 - I_RW3_ZZ*I_xy*I_yx**2*w_x**2 - I_RW3_ZZ*I_xy*I_zx**2*w_x**2 + I_RW3_ZZ*I_yx*I_zy**2*w_y**2 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*M_mt3 + I_RW2_YY*I_RW3_ZZ*I_xx*M_mt3 + 
        I_RW1_XX*I_RW3_ZZ*I_yy*M_mt3 - I_RW1_XX*I_RW3_ZZ*I_zy*M_mt2 - I_RW2_YY*I_RW3_ZZ*I_zx*M_mt1 - I_RW3_ZZ*I_xx*I_yy*M_mt3 + I_RW3_ZZ*I_xy*I_yx*M_mt3 + I_RW3_ZZ*I_xx*I_zy*M_mt2 - I_RW3_ZZ*I_xy*I_zx*M_mt2 - 
        I_RW3_ZZ*I_yx*I_zy*M_mt1 + I_RW3_ZZ*I_yy*I_zx*M_mt1 + I_RW1_XX*I_RW2_YY*I_zz*t_motor3 + I_RW1_XX*I_RW3_ZZ*I_zy*t_motor2 + I_RW2_YY*I_RW3_ZZ*I_zx*t_motor1 - I_RW2_YY*I_xx*I_zz*t_motor3 + I_RW2_YY*I_xz*I_zx*t_motor3 - 
        I_RW3_ZZ*I_xx*I_zy*t_motor2 + I_RW3_ZZ*I_xy*I_zx*t_motor2 - I_RW1_XX*I_yy*I_zz*t_motor3 + I_RW1_XX*I_yz*I_zy*t_motor3 + I_RW3_ZZ*I_yx*I_zy*t_motor1 - I_RW3_ZZ*I_yy*I_zx*t_motor1 + I_xx*I_yy*I_zz*t_motor3 - 
        I_xx*I_yz*I_zy*t_motor3 - I_xy*I_yx*I_zz*t_motor3 + I_xy*I_yz*I_zx*t_motor3 + I_xz*I_yx*I_zy*t_motor3 - I_xz*I_yy*I_zx*t_motor3 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xy*w_y**2 + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yx*w_x**2 + 
        I_RW2_YY*I_RW3_ZZ*I_xx*I_xy*w_y**2 + I_RW1_XX*I_RW3_ZZ*I_xy*I_yy*w_y**2 - I_RW2_YY*I_RW3_ZZ*I_xx*I_yx*w_x**2 - I_RW1_XX*I_RW3_ZZ*I_yx*I_yy*w_x**2 + I_RW1_XX*I_RW3_ZZ*I_xz*I_zy*w_z**2 - I_RW2_YY*I_RW3_ZZ*I_yz*I_zx*w_z**2 - 
        I_RW1_XX*I_RW3_ZZ*I_zx*I_zy*w_x**2 + I_RW2_YY*I_RW3_ZZ*I_zx*I_zy*w_y**2 - I_RW3_ZZ*I_xx*I_xy*I_yy*w_y**2 + I_RW3_ZZ*I_xx*I_yx*I_yy*w_x**2 - I_RW3_ZZ*I_xx*I_xz*I_zy*w_z**2 + I_RW3_ZZ*I_xy*I_xz*I_zx*w_z**2 + 
        I_RW3_ZZ*I_xx*I_zx*I_zy*w_x**2 - I_RW3_ZZ*I_yx*I_yz*I_zy*w_z**2 + I_RW3_ZZ*I_yy*I_yz*I_zx*w_z**2 - I_RW3_ZZ*I_yy*I_zx*I_zy*w_y**2 + I_RW2_YY*I_RW3_ZZ*I_xx**2*w_x*w_y - I_RW1_XX*I_RW3_ZZ*I_yy**2*w_x*w_y - 
        I_RW1_XX*I_RW3_ZZ*I_zy**2*w_x*w_y + I_RW2_YY*I_RW3_ZZ*I_zx**2*w_x*w_y + I_RW3_ZZ*I_xx*I_yy**2*w_x*w_y - I_RW3_ZZ*I_xx**2*I_yy*w_x*w_y + I_RW3_ZZ*I_xx*I_zy**2*w_x*w_y - I_RW3_ZZ*I_xx**2*I_zy*w_x*w_z + 
        I_RW3_ZZ*I_xy**2*I_zx*w_y*w_z - I_RW3_ZZ*I_yy*I_zx**2*w_x*w_y - I_RW3_ZZ*I_yx**2*I_zy*w_x*w_z + I_RW3_ZZ*I_yy**2*I_zx*w_y*w_z + I_RW1_XX*I_RW2_YY**2*I_RW3_ZZ*w_x*theta_dot_RW2 - 
        I_RW1_XX**2*I_RW2_YY*I_RW3_ZZ*w_y*theta_dot_RW1 - I_RW2_YY**2*I_RW3_ZZ*I_xx*w_x*theta_dot_RW2 + I_RW1_XX**2*I_RW3_ZZ*I_yy*w_y*theta_dot_RW1 - I_RW1_XX*I_RW3_ZZ**2*I_zy*w_x*theta_dot_RW3 + 
        I_RW1_XX**2*I_RW3_ZZ*I_zy*w_z*theta_dot_RW1 + I_RW2_YY*I_RW3_ZZ**2*I_zx*w_y*theta_dot_RW3 - I_RW2_YY**2*I_RW3_ZZ*I_zx*w_z*theta_dot_RW2 + I_RW3_ZZ**2*I_xx*I_zy*w_x*theta_dot_RW3 - I_RW3_ZZ**2*I_xy*I_zx*w_x*theta_dot_RW3 + 
        I_RW3_ZZ**2*I_yx*I_zy*w_y*theta_dot_RW3 - I_RW3_ZZ**2*I_yy*I_zx*w_y*theta_dot_RW3 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xx*w_x*w_y - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xz*w_y*w_z + 
        I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yy*w_x*w_y + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yz*w_x*w_z + I_RW2_YY*I_RW3_ZZ*I_xx*I_xz*w_y*w_z + I_RW1_XX*I_RW3_ZZ*I_xx*I_yy*w_x*w_y - 
        I_RW2_YY*I_RW3_ZZ*I_xx*I_yy*w_x*w_y + I_RW1_XX*I_RW3_ZZ*I_xz*I_yy*w_y*w_z - I_RW2_YY*I_RW3_ZZ*I_xx*I_yz*w_x*w_z + I_RW1_XX*I_RW3_ZZ*I_xx*I_zy*w_x*w_z + I_RW1_XX*I_RW3_ZZ*I_xy*I_zy*w_y*w_z - 
        I_RW1_XX*I_RW3_ZZ*I_yy*I_yz*w_x*w_z - I_RW2_YY*I_RW3_ZZ*I_yx*I_zx*w_x*w_z - I_RW2_YY*I_RW3_ZZ*I_yy*I_zx*w_y*w_z - I_RW1_XX*I_RW3_ZZ*I_zy*I_zz*w_x*w_z + I_RW2_YY*I_RW3_ZZ*I_zx*I_zz*w_y*w_z + 
        I_RW3_ZZ*I_xx*I_xy*I_yx*w_x*w_y - I_RW3_ZZ*I_xx*I_xz*I_yy*w_y*w_z + I_RW3_ZZ*I_xy*I_xz*I_yx*w_y*w_z + I_RW3_ZZ*I_xx*I_xy*I_zx*w_x*w_z - I_RW3_ZZ*I_xy*I_yx*I_yy*w_x*w_y - 
        I_RW3_ZZ*I_xx*I_xy*I_zy*w_y*w_z + I_RW3_ZZ*I_xx*I_yy*I_yz*w_x*w_z - I_RW3_ZZ*I_xy*I_yx*I_yz*w_x*w_z - I_RW3_ZZ*I_xy*I_zx*I_zy*w_x*w_y + I_RW3_ZZ*I_yx*I_yy*I_zx*w_x*w_z + 
        I_RW3_ZZ*I_xx*I_zy*I_zz*w_x*w_z - I_RW3_ZZ*I_xy*I_zx*I_zz*w_x*w_z - I_RW3_ZZ*I_yx*I_yy*I_zy*w_y*w_z + I_RW3_ZZ*I_yx*I_zx*I_zy*w_x*w_y + I_RW3_ZZ*I_yx*I_zy*I_zz*w_y*w_z - 
        I_RW3_ZZ*I_yy*I_zx*I_zz*w_y*w_z + I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_xx*w_y*theta_dot_RW1 - I_RW1_XX*I_RW2_YY*I_RW3_ZZ*I_yy*w_x*theta_dot_RW2 - I_RW1_XX*I_RW3_ZZ*I_xx*I_yy*w_y*theta_dot_RW1 + 
        I_RW1_XX*I_RW3_ZZ*I_xy*I_yx*w_y*theta_dot_RW1 + I_RW2_YY*I_RW3_ZZ*I_xx*I_yy*w_x*theta_dot_RW2 - I_RW2_YY*I_RW3_ZZ*I_xy*I_yx*w_x*theta_dot_RW2 - I_RW1_XX*I_RW3_ZZ*I_xx*I_zy*w_z*theta_dot_RW1 + 
        I_RW1_XX*I_RW3_ZZ*I_xy*I_zx*w_z*theta_dot_RW1 - I_RW2_YY*I_RW3_ZZ*I_yx*I_zy*w_z*theta_dot_RW2 + I_RW2_YY*I_RW3_ZZ*I_yy*I_zx*w_z*theta_dot_RW2)/(I_RW3_ZZ*(I_RW2_YY*I_RW3_ZZ*I_xx - I_RW1_XX*I_RW2_YY*I_RW3_ZZ + 
        I_RW1_XX*I_RW3_ZZ*I_yy + I_RW1_XX*I_RW2_YY*I_zz - I_RW3_ZZ*I_xx*I_yy + I_RW3_ZZ*I_xy*I_yx - I_RW2_YY*I_xx*I_zz + I_RW2_YY*I_xz*I_zx - I_RW1_XX*I_yy*I_zz + I_RW1_XX*I_yz*I_zy + I_xx*I_yy*I_zz - I_xx*I_yz*I_zy - I_xy*I_yx*I_zz + 
        I_xy*I_yz*I_zx + I_xz*I_yx*I_zy - I_xz*I_yy*I_zx))
