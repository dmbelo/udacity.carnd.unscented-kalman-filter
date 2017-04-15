#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;

class UKF {
public:
  double lambda; // Spreading parameters for augmentation
  double sigma_v_dot; // Process noise standard deviation long. accel
  double sigma_psi_dot2; // Process noise standard deviation yaw accel 

  VectorXd x;
  MatrixXd P;
  VectorXd xa;
  VectorXd xs;
  MatrixXd xsa;
  MatrixXd Pa;
  MatrixXd A;

	UKF(double lambda_, double sigma_v_dot_, double sigma_psi_dot2_); // Construcor
	virtual ~UKF(); // Destructor
	void Initialize();
  void ProcessMeasurement(double dt);
  void GenerateAugmentedSigmaPoints();
  void PredictSigmaPoints(double dt);
  void PredictStateMeanAndCovariance();
  void PredictMeasurement();
  void PredictMeasurementMeanAndCovariance();
  void UpdateState();
  void CTRVProcessModel(VectorXd* x, VectorXd* nu, double dt);
  void RadarMeasurementModel();
  void LidarMeasurementModel();

  // void GenerateSigmaPoints(MatrixXd* Xsig_out);
  // void AugmentedSigmaPoints(MatrixXd* Xsig_out);
  // void SigmaPointPrediction(MatrixXd* Xsig_out);
  // void PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred);
  // void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out);
  // void UpdateState(VectorXd* x_out, MatrixXd* P_out);
  // void PredictionFunction(VectorXd* x, VectorXd* x_pred, double dt);

};

#endif /* UKF_H */
