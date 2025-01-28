import math
import numpy as np
import matplotlib.pyplot as plt

"""
Practice EKF 

Authors: Jacobo Wiesner, Martin Castellanos

Practice EKF calculating the altitude of a falling body using three variables:
altitude, velocity and ballistic coefficient.

params:
"""

def get_system_parameters():
  """
  Returns the system parameters and constants.
  """
  rho0 = 0.0034       # Air density at sea level
  g = 32.2            # Acceleration due to gravity (ft/s^2)
  k = 22000           # Scale height constant for air density
  
  # Get covariance matrices
  Q, R = get_covariance_matrices()
  
  params = {
    'rho0': rho0,
    'g': g,
    'k': k,
    'Q': Q,
    'R': R
  }
  return params

def state_transition(x, t, params):
  """
  Computes the derivatives of the state variables.

  Parameters:
  - x: State vector [x1, x2, x3]
  - t: Current time (unused as system is time-invariant)
  - params: Dictionary of system parameters

  Returns:
  - dxdt: Derivative of state vector
  """
  x1, x2, x3 = x
  rho0 = params['rho0']
  g = params['g']
  k = params['k']
  
  exp_term = np.exp(-x1 / k)
  
  dx1_dt = x2
  dx2_dt = (rho0 * exp_term * x2**2) / (2 * x3) - g
  dx3_dt = 0  # No process noise in state equations for this step
  
  dxdt = np.array([dx1_dt, dx2_dt, dx3_dt])
  return dxdt

def measurement_function(x):
  """
  Computes the measurement based on the current state.

  Parameters:
  - x: State vector [x1, x2, x3]

  Returns:
  - y: Measurement vector
  """
  x1 = x[0]
  y = np.array([x1])
  return y

def get_initial_conditions():
  """
  Returns the initial true state, initial estimate, and initial covariance.

  Returns:
  - x0: Initial true state vector
  - x_hat0: Initial estimated state vector
  - P0: Initial covariance matrix
  """
  x0 = np.array([100000, -6000, 1/2000])           # True initial state
  x_hat0 = np.array([100010, -6100, 1/2500])      # Initial estimate
  P0 = np.array([
    [500,      0,          0],
    [0,   20000,          0],
    [0,        0, 1/250000]
  ])                                             # Initial covariance
  return x0, x_hat0, P0

def get_simulation_settings():
  """
  Returns the simulation settings.

  Returns:
  - dt: Integration step size (seconds)
  - measurement_interval: Time between measurements (seconds)
  - num_simulations: Number of Monte Carlo simulations
  """
  dt = 0.0004            # 0.4 milliseconds
  measurement_interval = 0.5  # 0.5 seconds
  num_simulations = 1000
  return dt, measurement_interval, num_simulations

def compute_jacobian_A(x, params):
  """
  Computes the Jacobian matrix A = ∂f/∂x at state x.

  Parameters:
  - x: State vector [x1, x2, x3]
  - params: Dictionary of system parameters

  Returns:
  - A: Jacobian matrix (3x3)
  """
  x1, x2, x3 = x
  rho0 = params['rho0']
  k = params['k']
  
  exp_term = np.exp(-x1 / k)
  x2_squared = x2 ** 2
  x3_squared = x3 ** 2

  # Compute A21, A22, A23
  A21 = - (rho0 * exp_term * x2_squared) / (2 * k * x3)
  A22 = (rho0 * exp_term * x2) / x3
  A23 = - (rho0 * exp_term * x2_squared) / (2 * x3_squared)

  # Assemble the Jacobian matrix A
  A = np.array([
    [0,    1,    0],
    [A21, A22, A23],
    [0,    0,    0]
  ])
  return A

def compute_jacobian_H():
  """
  Computes the Jacobian matrix H = ∂h/∂x.

  Returns:
  - H: Measurement Jacobian matrix (1x3)
  """
  H = np.array([[1, 0, 0]])  # Since h(x) = x1 + v
  return H

def compute_jacobians(x, params):
  """
  Computes both the state transition Jacobian A and measurement Jacobian H.

  Parameters:
  - x: State vector [x1, x2, x3]
  - params: Dictionary of system parameters

  Returns:
  - A: State transition Jacobian matrix (3x3)
  - H: Measurement Jacobian matrix (1x3)
  """
  A = compute_jacobian_A(x, params)
  H = compute_jacobian_H()
  return A, H

def get_covariance_matrices():
  """
  ---This can be adapted to make covariance change based on system occurences
  Returns the process noise covariance matrix Q and measurement noise covariance matrix R.
  
  Process Noise:
    E[w_i^2(t)] = 0 for i = 1, 2, 3
    Hence, Q is a zero matrix.
  
  Measurement Noise:
    E[v^2(t)] = 100
    Hence, R is a scalar or a 1x1 matrix with value 100.
  
  Returns:
    Q (numpy.ndarray): Process noise covariance matrix (3x3)
    R (numpy.ndarray): Measurement noise covariance matrix (1x1)
  """
  # Process Noise Covariance Matrix Q
  Q = np.zeros((3, 3))  # Since E[w_i^2(t)] = 0 for all i
  
  # Measurement Noise Covariance Matrix R
  R = np.array([[100]])  # Measurement noise variance
  
  return Q, R

def ekf_predict(x_hat, P, params, dt):
  """
  Performs the EKF prediction step.

  Parameters:
  - x_hat: Current state estimate vector
  - P: Current covariance matrix
  - params: Dictionary of system parameters (includes Q)
  - dt: Time step for integration

  Returns:
  - x_hat_pred: Predicted state estimate vector
  - P_pred: Predicted covariance matrix
  """
  # Extract Q from parameters
  Q = params['Q']
  
  # Compute the Jacobian A at current x_hat
  A = compute_jacobian_A(x_hat, params)
  
  # Predict the next state using Euler integration
  x_dot = state_transition(x_hat, 0, params)  # Assuming time-invariant
  x_hat_pred = x_hat + x_dot * dt
  
  # Predict the next covariance
  P_pred = P + (A @ P + P @ A.T + Q) * dt
  
  return x_hat_pred, P_pred

def simulate_true_state(x0, params, dt, total_time):
  """
  Simulates the true state of the system over time using Euler's method.
  
  Parameters:
  - x0 (numpy.ndarray): Initial true state vector [x1, x2, x3]
  - params (dict): System parameters
  - dt (float): Time step size (seconds)
  - total_time (float): Total simulation time (seconds)
  
  Returns:
  - true_states (list of numpy.ndarray): True state at each time step
  """
  true_states = [x0.copy()]
  x = x0.copy()
  num_steps = int(total_time / dt)
  
  for _ in range(num_steps):
    # Compute state derivatives
    dxdt = state_transition(x, 0, params)  # t is unused in state_transition
    
    # Update state using Euler's method
    x = x + dxdt * dt
    
    true_states.append(x.copy())

  return true_states

def ekf_update(x_hat_pred, P_pred, y, params):
  """
  Performs the EKF update step.

  Parameters:
  - x_hat_pred: Predicted state estimate vector
  - P_pred: Predicted covariance matrix
  - y: Measurement vector
  - params: Dictionary of system parameters (includes R)

  Returns:
  - x_hat_updated: Updated state estimate vector
  - P_updated: Updated covariance matrix
  """
  # Extract R from parameters
  R = params['R']
  
  # Compute the Jacobian H
  H = compute_jacobian_H()
  
  # Compute the innovation (measurement residual)
  y_pred = measurement_function(x_hat_pred)
  innovation = y - y_pred
  
  # Compute the innovation covariance
  S = H @ P_pred @ H.T + R
  
  # Compute the Kalman Gain
  K = P_pred @ H.T @ np.linalg.inv(S)
  
  # Update the state estimate
  x_hat_updated = x_hat_pred + K @ innovation
  
  # Update the covariance estimate
  P_updated = (np.eye(len(x_hat_pred)) - K @ H) @ P_pred
  
  return x_hat_updated, P_updated

def run_ekf(total_time):
  """
  Runs the hybrid Extended Kalman Filter over the simulation period.
  
  Returns:
  - state_estimates: List of state estimates at each measurement time
  - covariances: List of covariance matrices at each measurement time
  """
  # Get initial conditions and parameters
  x0, x_hat0, P0 = get_initial_conditions()
  params = get_system_parameters()
  Q, R = get_covariance_matrices()
  params['Q'] = Q
  params['R'] = R
  
  # Simulation settings
  dt, measurement_interval, num_simulations = get_simulation_settings() # num simpulations would be used if we wanted to run simulation many times
  #total_time = 100  # Define total simulation time as needed
  time_steps = np.arange(0, total_time, dt)
  print(f'Time Steps: {time_steps}')
  #measurement_times = np.arange(0, total_time, measurement_interval) - unnecessary because np.close is used to check for measurment interval
  
  # Initialize state and covariance estimates
  x_hat = x_hat0.copy()
  P = P0.copy()
  x_true = x0.copy()
  
  # Lists to store estimates
  measurement_times = []
  true_altitudes = []
  true_velocities = []
  true_ballistics = []
    
  estimated_altitudes = []
  estimated_velocities = []
  estimated_ballistics = []
  
  for t in time_steps:
    # Prediction step
    x_hat, P = ekf_predict(x_hat, P, params, dt)
    
    # Update true state
    dxdt_true = state_transition(x_true, t, params)
    x_true = x_true + dxdt_true * dt
    
    # Measurement update at specified intervals
    if np.isclose(t % measurement_interval, 0, atol=dt/10):
      # Generate noisy measurement
      y = measurement_function(x_true) + np.random.normal(0, np.sqrt(R[0, 0]))
      
      # Update step
      x_hat, P = ekf_update(x_hat, P, y, params)
      
      # Store values
      measurement_times.append(t)
            
      true_altitudes.append(x_true[0])
      true_velocities.append(x_true[1])
      true_ballistics.append(1.0 / x_true[2])
      
      estimated_altitudes.append(x_hat[0])
      estimated_velocities.append(x_hat[1])
      estimated_ballistics.append(1.0 / x_hat[2])
  
  return measurement_times, true_altitudes, true_velocities, true_ballistics, estimated_altitudes, estimated_velocities, estimated_ballistics

def main():
  # Define total simulation time (seconds)
  total_time = 16  # Adjust as needed
  
  # Run a single EKF simulation
  measurement_times, true_altitudes, true_velocities, true_ballistics, estimated_altitudes, estimated_velocities, estimated_ballistics = run_ekf(total_time)
  
  fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 15))
    
  # 1. Altitude vs Time
  axes[0].plot(measurement_times, true_altitudes, label='True Altitude', marker='o')
  axes[0].plot(measurement_times, estimated_altitudes, label='Estimated Altitude', marker='x')
  axes[0].set_xlabel('Time (s)')
  axes[0].set_ylabel('Altitude')
  axes[0].set_title('True vs. Estimated Altitude Over Time')
  axes[0].grid(True)
  axes[0].legend()
  
  # 2. Velocity vs Time
  axes[1].plot(measurement_times, true_velocities, label='True Velocity', marker='o')
  axes[1].plot(measurement_times, estimated_velocities, label='Estimated Velocity', marker='x')
  axes[1].set_xlabel('Time (s)')
  axes[1].set_ylabel('Velocity (ft/s)')
  axes[1].set_title('True vs. Estimated Velocity Over Time')
  axes[1].grid(True)
  axes[1].legend()
  
  # 3. Ballistic Coefficient vs Time
  axes[2].plot(measurement_times, true_ballistics, label='True Ballistic Coefficient', marker='o')
  axes[2].plot(measurement_times, estimated_ballistics, label='Estimated Ballistic Coefficient', marker='x')
  axes[2].set_xlabel('Time (s)')
  axes[2].set_ylabel('Ballistic Coefficient')
  axes[2].set_title('True vs. Estimated Ballistic Coefficient Over Time')
  axes[2].grid(True)
  axes[2].legend()
  
  # Adjust layout so plots don’t overlap
  plt.tight_layout()
  plt.show()

if __name__ == "__main__":
  main()
