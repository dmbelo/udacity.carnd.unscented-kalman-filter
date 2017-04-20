#ifndef UKF_H
#define UKF_H

#include <vector>
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include "measurement_package.h"

using namespace std;
using namespace Eigen;

class UKF {
public:

  double lambda; // Spreading parameters for augmentation
  double sigma_v_dot; // Process noise standard deviation long. accel
  double sigma_psi_dot2; // Process noise standard deviation yaw accel 
  double dt; // Detla time (s)
  
  long long previous_timestamp; // Last time stamp
  
  bool is_initialized;

  VectorXd x; // State mean vector [NX x 1]
  MatrixXd P; // State covariance matrix [NX x NX]
  MatrixXd Pa; // Augmented state covariance matrix [NA x NA]
  VectorXd xa; // Augmented vector [NX+2 x 1]
  MatrixXd xs; // State sigma points [NX x NS]
  MatrixXd xsa; // Augmented state sigma points [NA x NS]
  VectorXd zp; // Predicted measurement mean vector [NZ x 1]
  MatrixXd S; // Predicted measurement covariance matrix [NZ x NZ]
  MatrixXd R; // Measurement covariance matrix [NZ x NZ]
  MatrixXd zs; // Measurement sigma points [NZ x NS]
  VectorXd weights; // Weights for mean/covariance calc during augmentation 
  MatrixXd Tc; // Cross correlation matrix
  MatrixXd K; // Kalman gain


  // Construcor
	UKF(double lambda_, 
      double sigma_v_dot_,
      double sigma_psi_dot2_,
      double std_laspx,
      double std_laspy,
      double std_radr,
      double std_radphi,
      double std_radrd);

	virtual ~UKF(); // Destructor
  void ProcessMeasurement(MeasurementPackage measurement);
  void GenerateAugmentedSigmaPoints();
  void PredictStateSigmaPoints(double dt);
  void CalculateStateMeanAndCovariance();
  void PredictRadarMeasurementSigmaPoints();
  void PredictLidarMeasurementSigmaPoints();
  void CalculateRadarMeasurementMeanAndCovariance();
  void CalculateLidarMeasurementMeanAndCovariance();
  void UpdateStateRadar(VectorXd z);
  void UpdateStateLidar(VectorXd z);
  void CTRVProcessModel(Ref<VectorXd> xp, Ref<VectorXd> x, Ref<VectorXd> nu, double dt);
  void RadarMeasurementModel(Ref<VectorXd> zp, Ref<VectorXd> x);
  void LidarMeasurementModel(Ref<VectorXd> zp, Ref<VectorXd> x);
};

#endif /* UKF_H */
