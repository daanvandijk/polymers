#include "spring.h"
#include <iostream>
#include <assert.h>

using namespace std;

double spring::potential(MatrixXd &R, double k, double l) {
  assert(k > 0);
  assert(l > 0);
  double V = 0.0;
  MatrixXd L = MatrixXd::Zero(R.rows(), R.cols()-1);
  // loop through particles
  for (int ip = 0; ip < R.cols()-1; ip++) {
    L.col(ip) = (R.col(ip)-R.col(ip+1))/polymer::d(R,ip,ip+1);
  } 
  for (int ip = 0; ip < R.cols()-1; ip++) {
    V += 0.5*k*(R.col(ip)-R.col(ip+1)-l*L.col(ip)).transpose() 
      * (R.col(ip)-R.col(ip+1)-l*L.col(ip));
  } 
  return V;
}

void spring::force(MatrixXd &R, double k, double l, Eigen::MatrixXd &F) { 
  assert(k > 0);
  assert(R.cols() == F.cols());
  assert(R.rows() == F.rows());

  // loop through particles
  for (int ip = 0; ip < R.cols(); ip++) {
    if (ip > 0) {
      F.col(ip) = F.col(ip) - k*(R.col(ip)-R.col(ip-1))*(1-l/polymer::d(R,ip,ip-1));
    }
    if (ip < R.cols() - 1) {
      F.col(ip) = F.col(ip) - k*(R.col(ip)-R.col(ip+1))*(1-l/polymer::d(R,ip,ip+1));
    }
  }
}

void spring::force_numerical(MatrixXd &R, double k, double l, MatrixXd &F) { 
  assert(k > 0);
  assert(R.cols() == F.cols());
  assert(R.rows() == F.rows());
  MatrixXd R2 = R;

  double delta = 1e-9;
  // loop through dimensions, typically 2 or 3 for a polymer
  for (int d = 0; d < R.rows(); d++) {
    // loop through particles
    for (int i = 0; i < R.cols(); i++) {
      R2(d,i) += delta;
      F(d, i) -= (potential(R2,k,l)-potential(R,k,l))/delta;
      R2(d,i) = R(d,i);
    }
  }
}

void spring::test() {
  cout << "Spring unit test\n";

  // test spring force
  // compare numerical differentation of the potential with the
  // analytical force
  {
    MatrixXd R(3, 10);
    MatrixXd F1 = MatrixXd::Zero(3, 10);
    MatrixXd F2 = MatrixXd::Zero(3, 10);

    polymer::generate(1.0, 0.1, 0.0, 0.3, R);
    force(R, 1.0, 1.0, F1);
    force_numerical(R, 1.0, 1.0, F2);
    assert(polymer::compare(F1, F2) < 1e-2);
  }
}
