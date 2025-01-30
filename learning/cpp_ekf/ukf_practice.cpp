#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "matplotlibcpp.h"
#include <numpy/arrayobject.h>

using namespace std;

// Constants for the system
struct SystemParameters {
    double g = 9.81;              // Gravity (m/s^2)
    double rho = 1.225;           // Air density (kg/m^3)
    double scaleHeight = 8000;    
    double processNoise = 1e-5;   // Process noise covariance
    double measurementNoise = 25; // Measurement noise covariance
};

// State transition function
vector<double> stateTransition(const vector<double>& x, double dt, const SystemParameters& params) {
    double v = x[1]; // Velocity
    double alt = x[0]; // Altitude
    double drag = params.rho * exp(-alt / params.scaleHeight) * v * fabs(v);

    vector<double> x_next(2);
    x_next[0] = alt + v * dt;
    x_next[1] = v - params.g * dt - drag * dt;

    return x_next;
}

// Measurement function
double measurementFunction(const vector<double>& x) {
    return x[0]; // Measure altitude only
}

// Jacobian of state transition function
vector<vector<double>> computeJacobianA(const vector<double>& x, double dt, const SystemParameters& params) {
    double v = x[1];
    double alt = x[0];
    double drag = params.rho * exp(-alt / params.scaleHeight);

    vector<vector<double>> A(2, vector<double>(2));
    A[0][0] = 1.0;
    A[0][1] = dt;
    A[1][0] = drag * dt / params.scaleHeight * v * fabs(v);
    A[1][1] = 1.0 - 2 * drag * dt;

    return A;
}

// Jacobian of measurement function
vector<vector<double>> computeJacobianH() {
    return {{1.0, 0.0}}; // Measurement matrix
}

// EKF Prediction Step
void ekfPredict(vector<double>& x, vector<vector<double>>& P, const SystemParameters& params, double dt) {
    x = stateTransition(x, dt, params);

    vector<vector<double>> A = computeJacobianA(x, dt, params);
    vector<vector<double>> Q = {{params.processNoise, 0.0}, {0.0, params.processNoise}};

    vector<vector<double>> AP(P.size(), vector<double>(P[0].size()));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            AP[i][j] = 0.0;
            for (size_t k = 0; k < A.size(); ++k) {
                AP[i][j] += A[i][k] * P[k][j];
            }
        }
    }

    vector<vector<double>> ATP(P.size(), vector<double>(P[0].size()));
    for (size_t i = 0; i < AP.size(); ++i) {
        for (size_t j = 0; j < AP[i].size(); ++j) {
            ATP[i][j] = 0.0;
            for (size_t k = 0; k < A.size(); ++k) {
                ATP[i][j] += AP[i][k] * A[j][k];
            }
        }
    }

    for (size_t i = 0; i < P.size(); ++i) {
        for (size_t j = 0; j < P[i].size(); ++j) {
            P[i][j] = ATP[i][j] + Q[i][j];
        }
    }
}

// EKF Update Step
void ekfUpdate(vector<double>& x, vector<vector<double>>& P, double z, const SystemParameters& params) {
    auto H = computeJacobianH();
    double R = params.measurementNoise;

    // Compute innovation
    double y = z - measurementFunction(x);

    // Compute Kalman gain
    double S = P[0][0] + R;
    vector<double> K = {P[0][0] / S, P[1][0] / S};

    // Update state
    x[0] += K[0] * y;
    x[1] += K[1] * y;

    // Update covariance
    vector<vector<double>> I = {{1.0, 0.0}, {0.0, 1.0}};
    for (size_t i = 0; i < P.size(); ++i) {
        for (size_t j = 0; j < P[i].size(); ++j) {
            P[i][j] = (I[i][j] - K[i] * H[0][j]) * P[i][j];
        }
    }
}


namespace plt = matplotlibcpp;

// Add your EKF code here...

int main() {
    SystemParameters params;

    double dt = 0.1; // Time step
    int steps = 100;

    vector<double> trueAltitudes, estimatedAltitudes, time;
    vector<double> trueVelocities, estimatedVelocities;

    // Initialize true state, estimated state, and covariance
    vector<double> trueState = {10000, 0}; // Altitude and velocity
    vector<double> estimatedState = {9500, 10}; // Initial guess
    vector<vector<double>> P = {{1000, 0}, {0, 100}}; // Initial covariance

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> noise(0, sqrt(params.measurementNoise));

    for (int t = 0; t < steps; ++t) {
        // Simulate true state
        trueState = stateTransition(trueState, dt, params);
        double measurement = measurementFunction(trueState) + noise(gen);

        // EKF steps
        ekfPredict(estimatedState, P, params, dt);
        ekfUpdate(estimatedState, P, measurement, params);

        // Store results for plotting
        time.push_back(t * dt);
        trueAltitudes.push_back(trueState[0]);
        estimatedAltitudes.push_back(estimatedState[0]);
        trueVelocities.push_back(trueState[1]);
        estimatedVelocities.push_back(estimatedState[1]);
    }

    // Plotting Altitude
    plt::figure();
    plt::named_plot("True Altitude", time, trueAltitudes, "r-");
    plt::named_plot("Estimated Altitude", time, estimatedAltitudes, "b--");
    plt::xlabel("Time (s)");
    plt::ylabel("Altitude (m)");
    plt::title("EKF Altitude Estimation");
    plt::legend();
    plt::show();

    // Plot Velocity
    plt::figure();
    plt::named_plot("True Velocity", time, trueVelocities, "g-");
    plt::named_plot("Estimated Velocity", time, estimatedVelocities, "m--");
    plt::xlabel("Time (s)");
    plt::ylabel("Velocity (m/s)");
    plt::title("EKF Velocity Estimation");
    plt::legend();
    plt::show();

}
