#include <iostream>
#include "ukf.h"

UKF::UKF() {
  //TODO Auto-generated constructor stub
  Init();
}

UKF::~UKF() {
  //TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions: 
*******************************************************************************/

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  //set state dimension
  int n_x = 5;

  //define spreading parameter
  double lambda = 3 - n_x;

  //set example state
  VectorXd x = VectorXd(n_x);
  x <<   5.7441,
         1.3800,
         2.2049,
         0.5015,
         0.3528;

  //set example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  /*******************************************************************************
   * Student part begin
  ******************************************************************************/

  //set first column of sigma point matrix
  Xsig.col(0)  = x;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x + sqrt(lambda+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x - sqrt(lambda+n_x) * A.col(i);
  }

  /*******************************************************************************
   * Student part end
  ******************************************************************************/

  //print result
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  //write result
  *Xsig_out = Xsig;

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set example state
  VectorXd x = VectorXd(n_x);
  x <<   5.7441,
         1.3800,
         2.2049,
         0.5015,
         0.3528;

  //create example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/
 
  //create augmented mean state (last two elements are zero bc noise is zero mean)
  x_aug.head(5) = x;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P;
  P_aug.bottomRightCorner(2, 2) << std_a * std_a, 0, 0, std_yawdd * std_yawdd;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd tmp = sqrt(lambda + n_aug) * A;

  Xsig_aug.col(0) = x_aug;
  Xsig_aug.col(1) = x_aug + tmp.col(0);
  Xsig_aug.col(2) = x_aug + tmp.col(1);
  Xsig_aug.col(3) = x_aug + tmp.col(2);
  Xsig_aug.col(4) = x_aug + tmp.col(3);
  Xsig_aug.col(5) = x_aug + tmp.col(4);
  Xsig_aug.col(6) = x_aug + tmp.col(5);
  Xsig_aug.col(7) = x_aug + tmp.col(6);
  Xsig_aug.col(8) = x_aug - tmp.col(0);
  Xsig_aug.col(9) = x_aug - tmp.col(1);
  Xsig_aug.col(10) = x_aug - tmp.col(2);
  Xsig_aug.col(11) = x_aug - tmp.col(3);
  Xsig_aug.col(12) = x_aug - tmp.col(4);
  Xsig_aug.col(13) = x_aug - tmp.col(5);
  Xsig_aug.col(14) = x_aug - tmp.col(6);

  
  /*******************************************************************************
   * Student part end
   ******************************************************************************/ 
  
  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;

  /* expected result:
   Xsig_aug =
  5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
    1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
  2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
  0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
  0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
       0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
       0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641
  */

}

void UKF::PredictionFunction(VectorXd* x, VectorXd* x_pred, double dt) {
  /*
   * x_pred = f(x, dt)
   * Effect of process noise on position is approximated as a linear 
   * superposition (this holds well as long as yaw rate acceleation is not 
   * significant) 
   */

  double px = (*x)(0);
  double py = (*x)(1);
  double v = (*x)(2);
  double psi = (*x)(3);
  double psi_dot = (*x)(4);
  double nu_accel = (*x)(5);
  double nu_yaw_accel = (*x)(6);

  double vpd = v/psi_dot;
  double ppd = psi + psi_dot * dt;
  double sp = sin(psi);
  double cp = cos(psi);
  double dt2 = dt * dt;

  if (psi_dot > 0.001) 
  {
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

  (*x_pred)(0) = px;
  (*x_pred)(1) = py;
  (*x_pred)(2) = v;
  (*x_pred)(3) = psi;
  (*x_pred)(4) = psi_dot;

}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  double delta_t = 0.1; //time diff in sec
  /*******************************************************************************
   * Student part begin
  ******************************************************************************/

  //predict sigma points
  for (int i = 0; i < 15; i++) 
  {
    VectorXd x_ = Xsig_aug.col(i);
    VectorXd x_pred_ = Xsig_pred.col(i);
    PredictionFunction(&x_, &x_pred_, delta_t);
    Xsig_pred.col(i) = x_pred_;
  };

  //avoid division by zero

  //write predicted sigma points into right column
  

  /*******************************************************************************
   * Student part end
  ******************************************************************************/

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}

