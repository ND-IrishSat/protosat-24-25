'''
params.py

Reference file that contains all important variables for our system
Here, all factors that we would want to change are easily accessible
'''

import numpy as np


QUAT_INITIAL = np.array([1.0, 0.0, 0.0, 0.0])
RW_INITIAL = np.array([0.0, 0.0, 0.0, 0.0])
RW_EXTERNAL_TORQUE = np.array([0.0, 0.0, 0.0, 0.0])

# Pulse Width Modulation (PWM) signal that generates the max speed in our motors
MAX_PWM = 65535
# (TODO) max torque that our wheels can handle (Nm)
MAX_RW_TORQUE = 0.02

# =======  OPTIONS  =======================================

IDEAL_KNOWN = True
DATA_FILE = "data.txt"

# should be true if not doing controls
RUN_STATISTICAL_TESTS = False

# 0 = create pdf report, 1 = show 3D animation visualization, 2 = both, 3 = none
RESULT = 2

OUTPUT_FILE = "output.pdf"

# =======  CONTROLS  =======================================

# target orientation for if we're simulating controls
TARGET = np.array([1.0, 0.0, 0.5, 0.0])

# gains for our PID controller
KP = MAX_PWM * 4.0e-8       # Proportional gain
KI = MAX_PWM * 1e-9         # Integral gain
KD = MAX_PWM * 9e-9         # Derivative gain

TF = 10
DT = .02

# =======  UKF  =================================================

COVARIANCE_INITIAL_MAG = 5e-7

# Filter process noise Q represents uncertainty in our state transition model
# Higher values mean we trust our model less and measurements more
PROCESS_NOISE_MAG = 0.00001
PROCESS_NOISE_K = 10

# Filter measurement noise R
# Higher values mean we trust our measurements less
MEASUREMENT_MAGNETOMETER_NOISE = 0.001
MEASUREMENT_GYROSCOPE_NOISE = 0.01

# scaling parameters used to calculate scaling factor and sigma point weights
# alpha and k scale points around the mean. To capture the kurtosis of a gaussian distribution, a=1 and k=3-n should be used
#   If a decrease in the spread of the SPs is desired, use κ = 0 and alpha < 1
#   If an increase in the spread of the SPs is desired, use κ > 0 and alpha = 1
ALPHA = 0.001
K = 0
# beta minimizes higher order errors in covariance estimation
BETA = 2

# =======  SENSORS  ==================================================

# noise sd = noise density * sqrt(sampling rate)
# vn100 imu sampling rate from user manual = 200 Hz

# mag noise density from vn100 website = 140 uGauss /sqrt(Hz)
SENSOR_MAGNETOMETER_SD = (140 * 10e-6) * np.sqrt(200)

# gyro noise density from vn100 website = 0.0035 /s /sqrt(Hz)
SENSOR_GYROSCOPE_SD = 0.0035 * np.sqrt(200)

# =======  PHYSICS  ================================================

# principal moment of inertia for reaction wheels about spin axis and about axis transverse to spin axis respectively
SPIN_AXIS_INERTIA = 5.1e-7
# SPIN_AXIS_INERTIA = 1e-7
# TODO: what does this mean?
TRANSVERSE_AXIS_INERTIA = 0.0

# moment of inertia tensor of 2U CubeSat (w/o reaction wheel inertias) (kg m^2)
CUBESAT_BODY_INERTIA = (1e-7) * np.array([[46535.388, 257.834, 536.12],
                                          [257.834, 47934.771, -710.058],
                                          [536.12, -710.058, 23138.181]])

# Moments of Inertia of reaction wheels [g cm^2] - measured
Iw1 = (1/2)*38*1.8**2 # I_disc = 1/2 * M * R^2
Iw2 = Iw1
Iw3 = Iw1
Iw4 = Iw1
# Moment of inertia tensor of rxn wheels [kg m^2]
# this gets multiplied by SPIN_AXIS_INERTIA during EOMs calculation
RW_CONFIG_INERTIA = np.array([[Iw1, 0, 0, 0],
                              [0, Iw2, 0, 0],
                              [0, 0, Iw3, 0],
                              [0, 0, 0, Iw4]])

# Transformation matrix for NASA config given in Fundamentals pg 153-154
TRANSFORMATION = np.array([[1, 0, 0, 1/np.sqrt(3)],
                           [0, 1, 0, 1/np.sqrt(3)],
                           [0, 0, 1, 1/np.sqrt(3)]])

# motor model parameters (Maxon DCX 8 M (9 volts)) used for controls sim
# TODO: find accurate numbers for these
MAX_CURRENT = 1
MIN_CURRENT = -1
THERMAL_RESISTANCE = 0.01  # °C per A^2 (or Kelvin per A^2). how much the current flowing through the system causes heat generation
COOLING_CONSTANT = 0.1     # 1/s (rate of cooling). how quickly the temperature difference between the system and its surroundings dissipates
WHEEL_COUPLING_FACTOR = 0.5  # coupling between ambient and reaction wheel temperature

RWA = 3.54      # Ohms, winding resistance at ambient temperature
LW = 0.424e-3  # Henry
KT = 8.82e-3   # Torque constant Nm/A
KV = KT    # Voltage constant V*s/rad
JM = 5.1*(1e-7)   # Kg m^2
BM = 3.61e-6   # [N·m·s/rad] Viscous friction
ALPHA_CU = 0.00393 # copper's temperature coefficient [1/K]
Rha = 16.5      # K/W
Rwh = 2.66     # K/W
Cwa = 2.31/Rwh     # Thermal Capacitance
Cha = 162/Rha      # Thermal Capacitance



# ======= QUATERNION TOLERANCES ============================
# Quaternion error tolerances define how close we need to be to our target orientation
QUAT_ERROR_TOLERANCE = 0.01  # Maximum acceptable quaternion error magnitude
                            # 0.01 ≈ 1.15 degrees of rotation error
                            # sqrt(1 - cos(theta/2)) for small angles
ANGULAR_RATE_TOLERANCE = 0.001  # rad/s, maximum acceptable angular rate when "settled"

# ======= ENVIRONMENTAL DISTURBANCES ======================
# Typical disturbance torques for a 2U CubeSat in LEO
GRAVITY_GRADIENT_MAX = 1e-7  # Nm, maximum gravity gradient torque
                            # Varies with orbit altitude and satellite orientation
                            # Typically ~10^-7 Nm for 2U in 400km orbit

SOLAR_PRESSURE_MAX = 1e-8   # Nm, maximum solar radiation pressure torque
                           # Depends on surface area, reflectivity, sun angle
                           # Typically ~10^-8 Nm for 2U

AERO_DRAG_MAX = 2e-7       # Nm, maximum aerodynamic drag torque
                          # Varies with altitude (atmospheric density)
                          # Typically ~10^-7 Nm at 400km

MAGNETIC_RESIDUAL = 1e-6   # Am^2, residual magnetic dipole of spacecraft
                          # Creates torque when interacting with Earth's field
                          # Typical value for CubeSat with basic magnetic cleanliness

# Combined disturbance for simulation
TOTAL_DISTURBANCE_MAX = (GRAVITY_GRADIENT_MAX + 
                        SOLAR_PRESSURE_MAX + 
                        AERO_DRAG_MAX +
                        MAGNETIC_RESIDUAL * np.linalg.norm([19e-6, 1.7e-6, 49e-6]))  # Total maximum disturbance torque


    # # Add environmental disturbances
    # disturbance_direction = np.random.rand(3) - 0.5  # Random direction
    # disturbance_direction = disturbance_direction / np.linalg.norm(disturbance_direction)
    # disturbance_magnitude = np.random.uniform(0, TOTAL_DISTURBANCE_MAX)
    # disturbance_torque = disturbance_magnitude * disturbance_direction
    
    # # Add disturbance to dynamics
    # angular_acceleration += np.dot(self.J_B_inv, disturbance_torque)