#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // Initialize
  VectorXd rmse(4);
  rmse << 0.0, 0.0, 0.0, 0.0;

  /* Check that estimations and groud_truth are not zero and are equal to each
  other in size */
  if(estimations.size() != ground_truth.size() || estimations.size() == 0 || ground_truth.size() == 0) {
    
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;

  }

  // Accumulate square of residuals
  for(int i = 0; i < estimations.size(); i++) {

    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;

  }

  // Calculate mean
  rmse = rmse / estimations.size();

  // Calculate squared root
  rmse = rmse.array().sqrt();

  return rmse;

}
