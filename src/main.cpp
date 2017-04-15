#include <fstream>
#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
// using Eigen::MatrixXd;
// using Eigen::VectorXd;
// using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {
	
	check_arguments(argc, argv);

	string in_file_name_ = argv[1];
	ifstream in_file_(in_file_name_.c_str(), ifstream::in);

	string out_file_name_ = argv[2];
	ofstream out_file_(out_file_name_.c_str(), ofstream::out);

	check_files(in_file_, in_file_name_, out_file_, out_file_name_);

	/**********************************************
	*  Set Measurements                          *
	**********************************************/

	vector<MeasurementPackage> measurement_pack_list;
	vector<GroundTruthPackage> gt_pack_list;

	string line;

	// prep the measurement packages (each line represents a measurement at a
	// timestamp)


	return 0;
}