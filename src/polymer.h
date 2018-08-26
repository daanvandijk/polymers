#include <Eigen/Dense>

#ifndef _POLYMER_
#define _POLYMER_
namespace polymer { 
  // for working with matrices of polymer with different time states
  void save_state(int index_datapoint, Eigen::MatrixXd &R, Eigen::MatrixXd &Rcurrent);
  // returns mse (mean squared error)
  double compare(Eigen::MatrixXd s1, Eigen::MatrixXd s2); 
  // calculate dot product of row a and b of matrix s
  double dot(Eigen::MatrixXd &s, int a, int b); 
  void generate(double mu_l, double sigma_l, double mu_theta, double sigma_theta,Eigen::MatrixXd &s); 
  Eigen::MatrixXd generate(double mu_l, double sigma_l, double mu_theta, double sigma_theta, int Np); 
  // Euclidean metric/distance function
  // Calculates distance between two particles a and b
  double d(Eigen::MatrixXd &m, int a, int b); 
  // Distance between particles from a time integration matrix
  // The dimensions of this matrix are:
  // (3 \ctimes datapoints) \times Np
  double d(Eigen::MatrixXd &m, int a, int b, int time_offset);

  double angle(Eigen::MatrixXd &m, int ip);
  double angle(Eigen::MatrixXd &m, int ip, int time_offset);
  void test(); 
}
#endif
