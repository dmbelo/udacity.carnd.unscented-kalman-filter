#include "ukf.h"
#include <math.h>

// Assuming only radar measurements

#define _USE_MATH_DEFINES
#define NX 5 // Number of states
#define NZ 3 // Number of measurements
#define NA 7 // Number of augmented states
#define NS 15 // Number of sigma points

UKF::UKF(double lambda_, double sigma_v_dot_, double sigma_psi_dot2_) {

    cout << "UKF Constructor" << endl;
    lambda = lambda_;
    sigma_v_dot = sigma_v_dot_; 
    sigma_psi_dot2 = sigma_psi_dot2_;
    is_initialized = false;
    dt = 0.0;

    x = VectorXd::Zero(NX);
    P = MatrixXd::Zero(NX, NX);
    xa = VectorXd::Zero(NA);
    xs = MatrixXd::Zero(NX, NS);
    xsa = MatrixXd::Zero(NA, NS);
    Pa = MatrixXd::Zero(NA, NA);
    // A = MatrixXd::Zero(NA, NA);

}

UKF::~UKF() {

    cout << "UKF Destructor" << endl;

}

void UKF::ProcessMeasurement(MeasurementPackage measurement) {

    if (is_initialized == false) {
        // Initialize

        if (measurement.sensor_type_ == MeasurementPackage::RADAR) {
            float rho = measurement.raw_measurements_(0);
	        float phi = measurement.raw_measurements_(1);
            float rho_dot = measurement.raw_measurements_(2);
            float px = rho * cos(phi);
            float py = rho * sin(phi);
            float v = rho_dot;
            float psi = phi;
            float psi_dot = 0.0;
            x << px, py, v, psi, psi_dot;

        }

        previous_timestamp = measurement.timestamp_;    
        is_initialized = true;
        return;

    }

    // Calculate timestep
    dt = (measurement.timestamp_ - previous_timestamp) / 1000000.0;
    previous_timestamp = measurement.timestamp_;

    // Generate sigma 
    // xsa - augmented state sigma points [7x15] at k
    GenerateAugmentedSigmaPoints();

    // Predict sigma points
    // xs - state sigma points [5x15] at k+1
    PredictSigmaPoints(dt);

    // Predict mean/covariance of predicted state
    // x - predicted state mean vector [5x1]
    // P - predicted state covariance [5x5]
    // PredictStateMeanAndCovariance();

    // Predict measurement
    // Use the (predicted?) xs sigma points
    // Will have to deal with different measurement types here {RADAR/LIDAR}
    // zs - predicted measurement sigma points [3x15]
    // PredictMeasurement();

    // Predict mean/covariance of predicted measurements
    // zp - predicted measurement mean vector [3x1]
    // S - predicted measurement covariance [3x3]
    // PredictMeasurementMeanAndCovariance();
    
    // Update state
    // z - measurement [3x1]
    // x - updated state mean vector [5x1]
    // P - updated state covariance matrix [5x5] 
    // UpdateState();

}

void UKF::GenerateAugmentedSigmaPoints() {

    xa.head(NX) = x;
    Pa.topLeftCorner(NX, NX) = P;
    Pa(5, 5) = sigma_v_dot * sigma_v_dot;
    Pa(6, 6) = sigma_psi_dot2 * sigma_psi_dot2;

    MatrixXd A = (Pa.llt().matrixL());
    A *= sqrt(lambda + NA); //square root matrix

    xsa.col(0) = xa;
    xsa.col(1) = xa + A.col(0);
    xsa.col(2) = xa + A.col(1);
    xsa.col(3) = xa + A.col(2);
    xsa.col(4) = xa + A.col(3);
    xsa.col(5) = xa + A.col(4);
    xsa.col(6) = xa + A.col(5);
    xsa.col(7) = xa + A.col(6);
    xsa.col(8) = xa - A.col(0);
    xsa.col(9) = xa - A.col(1);
    xsa.col(10) = xa - A.col(2);
    xsa.col(11) = xa - A.col(3);
    xsa.col(12) = xa - A.col(4);
    xsa.col(13) = xa - A.col(5);
    xsa.col(14) = xa - A.col(6);

}

void UKF::PredictSigmaPoints(double dt) {

    for (int i = 0; i < NA; i++) 
    {
        CTRVProcessModel(xs.col(i), xsa.col(i).head(NX), xsa.col(i).tail(NA-NX), dt);
    }

}

void UKF::PredictStateMeanAndCovariance() {

}

void UKF::PredictMeasurement() {

}

void UKF::PredictMeasurementMeanAndCovariance() {

}

void UKF::UpdateState() {

}

void UKF::CTRVProcessModel(Ref<VectorXd> xp, Ref<VectorXd> x, Ref<VectorXd> nu, double dt) {

    // Constant Turn Rate and Velocity (CTRV) Model

    // States
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double psi = x(3);
    double psi_dot = x(4);

    // Noise
    double nu_accel = nu(0);
    double nu_yaw_accel = nu(1);

    // Minimize repeated calcs
    double ppd = psi + psi_dot * dt;
    double sp = sin(psi);
    double cp = cos(psi);
    double dt2 = dt * dt;

    if (abs(psi_dot) > 0.001) 
    {
        double vpd = v/psi_dot;
        px = px + vpd * (sin(ppd) - sp) + 0.5 * dt2 * cp * nu_accel;
        py = py + vpd * (-cos(ppd) + cp) + 0.5 * dt2 * sp * nu_accel;
        v = v + dt * nu_accel;
        psi = psi + psi_dot * dt + 0.5 * dt2 * nu_yaw_accel;
        psi_dot = psi_dot + dt * nu_yaw_accel;
    }
    else
    {
        px = px + v * cp * dt + 0.5 * dt2 * cp * nu_accel;
        py = py + v * sp * dt + 0.5 * dt2 * sp * nu_accel;
        v = v + dt * nu_accel;
        psi = psi + 0.5 * dt2 * nu_yaw_accel;
        psi_dot = psi_dot + dt * nu_yaw_accel;
    }

    xp(0) = px;
    xp(1) = py;
    xp(2) = v;
    xp(3) = psi;
    xp(4) = psi_dot;

}


void UKF::RadarMeasurementModel() {

}

void UKF::LidarMeasurementModel() {

}