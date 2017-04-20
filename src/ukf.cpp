#include "ukf.h"
#include <math.h>

// Assuming only radar measurements

#define _USE_MATH_DEFINES
#define NX 5       // Number of states
#define NZ_RADAR 3 // Dimensionality of radar measurements
#define NZ_LIDAR 2 // Dimensionality of lidar measurement
#define NA 7       // Number of augmented states
#define NS 15      // Number of sigma points (2 * NA + 1)

UKF::UKF(double lambda_, double sigma_v_dot_, double sigma_psi_dot2_,
         double std_laspx, double std_laspy, double std_radr, 
         double std_radphi, double std_radrd) {

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
    zp = VectorXd::Zero(NZ);
    zs = MatrixXd::Zero(NZ, NS);
    S = MatrixXd::Zero(NZ, NZ);
    Tc = MatrixXd::Zero(NX, NZ);
    K = MatrixXd::Zero(NX, NZ);

    // Calculate weights
    weights = VectorXd::Zero(NS);    
    weights(0) = lambda / (lambda + NA);  
    for (int i = 1; i < NS; i++) {

        weights(i) = 0.5 / (lambda + NA);

    }

    // Assemble measurement covariance matrix
    R_radar = MatrixXd(NZ, NZ);
    R_radar << std_radr*std_radr, 0, 0,
         0, std_radphi*std_radphi, 0,
         0, 0, std_radrd*std_radrd;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage measurement) {

    if (is_initialized == false) {
        // Initialize

        if (measurement.sensor_type_ == MeasurementPackage::RADAR) {

            float rho = measurement.raw_measurements_(0);
	        float phi = measurement.raw_measurements_(1);
            float rho_dot = measurement.raw_measurements_(2);
            float px = rho * cos(phi);
            float py = rho * sin(phi);
            float v = 0.0;
            float psi = 0.0;
            float psi_dot = 0.0;
            x << px, py, v, psi, psi_dot;
            P = MatrixXd::Identity(NX, NX);

        }

        previous_timestamp = measurement.timestamp_;    
        is_initialized = true;
        return;

    }
    else if (measurement.sensor_type_ == MeasurementPackage::LASER) {
        return;
    }

    // Calculate timestep
    dt = (measurement.timestamp_ - previous_timestamp) / 1000000.0;
    previous_timestamp = measurement.timestamp_;

    // cout << "***************************************************" << endl;
    // cout << "State Mean Vector" << endl;
    // cout << x << endl;
    // cout << "State Covariance Vector" << endl;
    // cout << P << endl;

    // Generate sigma 
    // xsa - augmented state sigma points [7x15] at k
    GenerateAugmentedSigmaPoints();
    // cout << "Sigma Points" << endl;
    // cout << xsa << endl;

    // Predict sigma points
    // xs - state sigma points [5x15] at k+1
    PredictStateSigmaPoints(dt);
    // cout << "Predicted State Sigma Points" << endl;
    // cout << xs << endl;

    // Predict mean/covariance of predicted state
    // x - predicted state mean vector [5x1]
    // P - predicted state covariance [5x5]
    CalculateStateMeanAndCovariance();
    // cout << "Predicted State Mean Vector" << endl;
    // cout << x << endl;
    // cout << "Predicted State Covariance Matrix" << endl;
    // cout << P << endl;

    
    if (sensor_type == MeasurementPackage::RADAR) {

        // Predict measurement
        // Use the (predicted?) xs sigma points
        // Will have to deal with different measurement types here {RADAR/LIDAR}
        // zs - predicted measurement sigma points [3x15]
        PredictRadarMeasurementSigmaPoints();
        // Predict mean/covariance of predicted measurements
        // zp - predicted measurement mean vector [3x1]
        // S - predicted measurement covariance [3x3]
        CalculateRadarMeasurementMeanAndCovariance();
        // Update state
        // z - measurement [3x1]
        // x - updated state mean vector [5x1]
        // P - updated state covariance matrix [5x5] 
        UpdateStateRadar();

    }
    else if (sensor_type == MeasurementPackage::LASER) {

        PredictLidarMeasurementSigmaPoints();
        CalculateRadarMeasurementMeanAndCovariance();
        UpdateStateLidar();

    }
    
    // cout << "Updated State Mean Vector" << endl;
    cout << x.transpose() << endl;
    // cout << "Updated State Covariance Matrix" << endl;
    // cout << P << endl << endl;

}

void UKF::GenerateAugmentedSigmaPoints() {

    xa.head(NX) = x;
    Pa.topLeftCorner(NX, NX) = P;
    Pa(5, 5) = sigma_v_dot * sigma_v_dot;
    Pa(6, 6) = sigma_psi_dot2 * sigma_psi_dot2;

    MatrixXd A = Pa.llt().matrixL();
    A *= sqrt(lambda + NA); //square root matrix

    xsa.col(0) = xa;

    for (int i = 0; i < NA; i++) {

        xsa.col(i + 1) = xa + A.col(i);
        xsa.col(i + 1 + NA) = xa - A.col(i);

    }

}

void UKF::PredictStateSigmaPoints(double dt) {

    for (int i = 0; i < NS; i++) 
    {
        CTRVProcessModel(xs.col(i), xsa.col(i).head(NX), xsa.col(i).tail(NA-NX), dt);
    }

}

void UKF::CalculateStateMeanAndCovariance() { 

    x.fill(0.0); // Reset to state mean to zero
    for (int i = 0; i < NS; i++) {

        x = x + weights(i) * xs.col(i);

    }

    VectorXd x_diff = VectorXd(NX);
    P.fill(0.0); // Reset state covariance to zero
    for (int i = 0; i < NS; i++) {

        // state difference
        x_diff = xs.col(i) - x;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
        while (x_diff(3) <- M_PI) x_diff(3) += 2.0 * M_PI;
        P = P + weights(i) * x_diff * x_diff.transpose();

    }

}

void UKF::PredictRadarMeasurementSigmaPoints() {


    for (int i = 0; i < NS; i++) 
    {
        RadarMeasurementModel(zs_radar.col(i), xsa.col(i).head(NX));
    }

}

void UKF::PredictLidarMeasurementSigmaPoints() {

    for (int i = 0; i < NS; i++) 
    {
        LidarMeasurementModel(zs_lidar.col(i), xsa.col(i).head(NX));
    }

}

void UKF::CalculateRadarMeasurementMeanAndCovariance() {

    // Predicted measurement mean
    zp_radar.fill(0.0);
    for (int i = 0; i < NS; i++) {

        zp_radar = zp_radar + weights(i) * zs_radar.col(i);
        
    }

    // Predicted measurement covariance
    VectorXd z_diff = VectorXd(NZ_RADAR);
    S.fill(0.0);
    for (int i = 0; i < NS; i++) {

        z_diff = zs_radar.col(i) - zp_radar;

        // angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

        S_radar = S_radar + weights(i) * z_diff * z_diff.transpose();

    }

    // Add measurement noise covariance matrix
    S_radar = S_radar + R_radar;

}

void UKF::CalculateLidarMeasurementMeanAndCovariance() {

    // Predicted measurement mean
    zp_lidar.fill(0.0);
    for (int i = 0; i < NS; i++) {

        zp_lidar = zp_lidar + weights(i) * zs_lidar.col(i);
        
    }

    // Predicted measurement covariance
    VectorXd z_diff = VectorXd(NZ_LIDAR);
    S.fill(0.0);
    for (int i = 0; i < NS; i++) {

        z_diff = zs_lidar.col(i) - zp_lidar;
        S_lidar = S_lidar + weights(i) * z_diff * z_diff.transpose();

    }

    // Add measurement noise covariance matrix
    S_lidar = S_lidar + R_lidar;

}

void UKF::UpdateStateRadar(VectorXd z) {

    VectorXd z_diff = VectorXd(NZ_RADAR);
    VectorXd x_diff = VectorXd(NX);
    
    Tc_radar.fill(0.0);
    for (int i = 0; i < NS; i++) {

        //residual
        z_diff = zs_radar.col(i) - zp_radar;
        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1) <- M_PI) z_diff(1) += 2.*M_PI;

        // state difference
        x_diff = xs.col(i) - x;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3) <- M_PI) x_diff(3) += 2.*M_PI;

        Tc_radar = Tc_radar + weights(i) * x_diff * z_diff.transpose();

    }

    K_radar = Tc_radar * S_radar.inverse(); // Kalman gain
    z_diff = z - zp_radar;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) <- M_PI) z_diff(1) += 2.0 * M_PI;

    // Update state mean and covariance matrix
    x = x + K_radar * z_diff;
    P = P - K_radar * S_radar * K_radar.transpose();

}

void UKF::UpdateStateLidar(VectorXd z) {

    VectorXd z_diff = VectorXd(NZ_LIDAR);
    VectorXd x_diff = VectorXd(NX);
    
    Tc_lidar.fill(0.0);
    for (int i = 0; i < NS; i++) {

        //residual
        z_diff = zs_lidar.col(i) - zp_lidar;

        // state difference
        x_diff = xs.col(i) - x;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3) <- M_PI) x_diff(3) += 2.*M_PI;

        Tc_lidar = Tc_lidar + weights(i) * x_diff * z_diff.transpose();

    }

    K_lidar = Tc_lidar * S_lidar.inverse(); // Kalman gain
    z_diff = z - zp_lidar;

    // Update state mean and covariance matrix
    x = x + K_lidar * z_diff;
    P = P - K_lidar * S_lidar * K_lidar.transpose();

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

    if (fabs(psi_dot) > 0.001) 
    {
        double vpd = v / psi_dot;
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

void UKF::RadarMeasurementModel(Ref<VectorXd> zp, Ref<VectorXd> x) {

    double px = x(0);
    double py = x(1);
    double v  = x(2);
    double psi = x(3);

    double vx = cos(psi) * v;
    double vy = sin(psi) * v;

    double rho = sqrt(px * px + py * py);

    zp(0) = rho; // rho
    zp(1) = atan2(py, px); // phi
    if (fabs(rho) > 0.0001) {
        zp(2) = (px * vx + py * vy ) / rho; // rho_dot
    }
    else {
        zp(2) = 0.0;
    }
    
}

void UKF::LidarMeasurementModel(Ref<VectorXd> zp, Ref<VectorXd> x) {

    zp(0) = x(0); // px
    zp(1) = x(1); // py

}