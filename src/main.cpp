#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main() {

	//Create a UKF instance
	UKF ukf(5, 0.2, 0.2);
	cout << ukf.lambda << endl;

/*******************************************************************************
* Programming assignment calls
*******************************************************************************/

    // MatrixXd Xsig_pred = MatrixXd(15, 5);
    // ukf.SigmaPointPrediction(&Xsig_pred);

	return 0;
}