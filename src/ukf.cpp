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
    
    // Radar data 
    zp_radar = VectorXd::Zero(NZ_RADAR);
    zs_radar = MatrixXd::Zero(NZ_RADAR, NS);
    S_radar = MatrixXd::Zero(NZ_RADAR, NZ_RADAR); 
    // Assemble measurement covariance matrix
    R_radar = MatrixXd(NZ_RADAR, NZ_RADAR);
    R_radar << std_radr*std_radr, 0,                     0,
               0,                 std_radphi*std_radphi, 0,
               0,                 0,                     std_radrd*std_radrd;
    Tc_radar = MatrixXd::Zero(NX, NZ_RADAR);   
    K_radar = MatrixXd::Zero(NX, NZ_RADAR);
    
    // Lidar data
    zp_lidar = VectorXd::Zero(NZ_LIDAR);
    zs_lidar = MatrixXd::Zero(NZ_LIDAR, NS);
    S_lidar = MatrixXd::Zero(NZ_LIDAR, NZ_LIDAR);
    R_lidar = MatrixXd(NZ_LIDAR, NZ_LIDAR);
    R_lidar << std_laspx*std_laspx, 0,
               0,                   std_laspy*std_laspy;
    Tc_lidar = MatrixXd::Zero(NX, NZ_LIDAR);
    K_lidar = MatrixXd::Zero(NX, NZ_LIDAR);

    // Calculate weights
    weights = VectorXd::Zero(NS);    
    weights(0) = lambda / (lambda + NA);  
    for (int i = 1; i < NS; i++) {

        weights(i) = 0.5 / (lambda + NA);

    }

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
            x << px, py, 0.0, 0.0, 0.0;

        }
        else if (measurement.sensor_type_ == MeasurementPackage::LASER) {

            x << measurement.raw_measurements_, 0.0, 0.0, 0.0;

        }

        // P = MatrixXd::Identity(NX, NX);
        P << 0.15, 0, 0, 0, 0,
             0, 0.15, 0, 0, 0,
             0, 0, 1, 0, 0, 
             0, 0, 0, 0.01, 0,
             0, 0, 0, 0, 0.01;

        previous_timestamp = measurement.timestamp_;    
        is_initialized = true;
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

    
    if (measurement.sensor_type_ == MeasurementPackage::RADAR) {

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
        UpdateStateRadar(measurement.raw_measurements_);

    }
    else if (measurement.sensor_type_ == MeasurementPackage::LASER) {

        PredictLidarMeasurementSigmaPoints();
        CalculateLidarMeasurementMeanAndCovariance();
        UpdateStateLidar(measurement.raw_measurements_);

    }
    
    // cout << "Updated State Mean Vector" << endl;
    // cout << x.transpose() << endl;
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
    S_radar.fill(0.0);
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
    S_lidar.fill(0.0);
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
    
    zp(1) = 0.0;
    if (fabs(py) > 0.001 || fabs(px) > 0.001) {
        zp(1) = atan2(py, px); // phi
    }
    
    
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