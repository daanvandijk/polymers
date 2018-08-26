#include "rousetransform.h"
#include "polymer.h"

using namespace Eigen;

double phi(double n, double p, double N) {
  return cos(p * M_PI * (n+0.5) / N);
}

void rouse::transform(MatrixXd &R, MatrixXd &X) {
  assert(R.cols() == X.cols());
  assert(R.rows() == X.rows());
  int N = R.cols();
  int dims = R.rows();

  for (int id = 0; id < dims; id++) {
    for (int p = 0; p < N; p++) {
      X(id, p) = 0.0;
      for (int n = 0; n < N; n++) {
        X(id, p) += phi(n, p, N) * R(id, n) / ((double) N);
      }
    }
  }
}

Eigen::MatrixXd rouse::transform(Eigen::MatrixXd &R) {
  MatrixXd X = MatrixXd::Zero(R.rows(), R.cols());
  transform(R, X);
  return X;
}

void rouse::invtransform(MatrixXd &X, MatrixXd &R) {
  assert(R.cols() == X.cols());
  assert(R.rows() == X.rows());
  int N = R.cols();
  int dims = R.rows();

  for (int id = 0; id < dims; id++) {
    for (int n = 0; n < N; n++) {
      R(id, n) = phi(n, 0, N) * X(id, 0);
      for (int p = 1; p < N; p++) {
        R(id, n) += 2 * phi(n, p, N) * X(id, p);
      }
    }
  }
}

Eigen::MatrixXd rouse::invtransform(Eigen::MatrixXd &X) {
  MatrixXd R = MatrixXd::Zero(X.rows(), X.cols());
  invtransform(X, R);
  return R;
}

void rouse::test() {
  printf("Rouse transform unit test\n");

  MatrixXd R = MatrixXd::Zero(3, 10);
  MatrixXd X = MatrixXd::Zero(3, 10);
  MatrixXd Xinv = MatrixXd::Zero(3, 10);

  // init R
  polymer::generate(1.0, 0.1, 0.0, 0.5, R);

  transform(R, X);
  invtransform(X, Xinv);
  assert(polymer::compare(R, Xinv) < 1e-9);
}
