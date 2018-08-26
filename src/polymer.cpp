#include <assert.h>
#include <iostream>
#include <math.h>
#include <random>
#include <chrono>
#include "polymer.h"

using namespace std;
using namespace Eigen;

void polymer::save_state(int index_datapoint, MatrixXd &R, Eigen::MatrixXd &Rcurrent) {
  assert(index_datapoint >= 0 && 3*index_datapoint < R.rows());
  assert(Rcurrent.rows() == 3);
  assert(Rcurrent.cols() == R.cols());
  assert(R.rows() % 3 == 0);

  // loop through dimensions
  for (int id = 0; id < 3; id++)
    R.row(id+3*index_datapoint) = Rcurrent.row(id);
}

// returns mse (mean squared error)
double polymer::compare(MatrixXd s1, MatrixXd s2) {
  assert(s1.rows() == s2.rows() && s1.cols() == s2.cols());
  double flag = 0.0;
  // loop through particles
  for (int i = 0; i < s1.cols(); i++) {
    // scale with number of particles
    flag += sqrt((s1.col(i)-s2.col(i)).squaredNorm()); 
  }
  return flag / (double)s1.cols();
}

// calculate dot product of row a and b of matrix s
double polymer::dot(MatrixXd &s, int a, int b) {
  assert(0 <= a && a < s.cols() && 0 <= b && b < s.cols());
  return s.col(a).dot(s.col(b));
}

MatrixXd polymer::generate(double mu_l, double sigma_l, double mu_theta, double sigma_theta, int Np) {
  MatrixXd pol = MatrixXd::Zero(3, Np);
  generate(mu_l, sigma_l, mu_theta, sigma_theta, pol);
  return pol;
}

void polymer::generate(double mu_l, double sigma_l, double mu_theta, double sigma_theta, MatrixXd &s) {
  assert(s.cols() > 2);
  assert(s.rows() == 3);

  // init random variables
  std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<double> distr_l(mu_l, sigma_l);
  std::normal_distribution<double> distr_theta(mu_theta, sigma_theta);

  s(0, 1) = distr_l(generator);
  double t = distr_theta(generator);
  for (int i = 1; i < s.cols() - 1; i++) {
    // Rotation matrix M:
    MatrixXd M = MatrixXd::Zero(2,2);
    M << 
      cos(t), -sin(t),
      sin(t), cos(t);

    double l = distr_l(generator);
    s.block(0,i+1,2,1) = 
      s.block(0,i,2,1)+ 
      (l/d(s,i,i-1))*M*(s.block(0,i,2,1)-s.block(0,i-1,2,1));
    t += distr_theta(generator);
  }
}

// Euclidean metric/distance function
// Calculates distance between two particles a and b
double polymer::d(MatrixXd &m, int a, int b) {
  assert(0 <= a && a < m.cols() && 0 <= b && b < m.cols());
  return sqrt((m.col(a) - m.col(b)).squaredNorm()); 
}

// Euclidean metric/distance function
// Calculates distance between two particles a and b
double polymer::d(MatrixXd &m, int a, int b, int time_offset) {
  assert(0 <= a && a < m.cols() && 0 <= b && b < m.cols());
  assert(m.rows() % 3 == 0);
  assert(0 <= time_offset && 3*time_offset+3 <= m.rows());
  return sqrt(
      (m.block(time_offset*3, a, 3, 1) - m.block(time_offset*3, b, 3, 1))
      .squaredNorm());
}


double polymer::angle(MatrixXd &m, int ip) {
  return polymer::angle(m, ip, 0);
}

double polymer::angle(MatrixXd &m, int ip, int time_offset) {
  assert(m.rows() >= (time_offset+1)*3);
  assert(m.cols() >= 3 + ip);
  Vector3d va, vb;
  va = m.block(3*time_offset, ip+1, 3, 1)-m.block(3*time_offset, ip, 3, 1);
  vb = m.block(3*time_offset, ip+2, 3, 1)-m.block(3*time_offset, ip+1, 3, 1);
  return atan2(sqrt(va.cross(vb).squaredNorm()), va.dot(vb));
}

void polymer::test() {
  cout << "Polymer unit test\n";

  // test distance function
  {
    MatrixXd state(3,3);
    generate(1.0, 0.0, 0.0, 0.3, state);
    assert(d(state, 0, 0) == 0.0);
    assert(d(state, 1, 0) == d(state, 0, 1));
    assert(d(state, 0, 2) <= d(state, 1, 0) + d(state, 1, 2));

    for (int i = 0; i < 3; i++) {
      state(0, i) = i;
      state(1, i) = i;
      state(2, i) = i;
    }
    assert(d(state, 0, 1) == sqrt(3));
    assert(d(state, 1, 2) == sqrt(3));
    assert(d(state, 0, 2) == sqrt(12));
  }

  // test 2nd distance function
  {
    MatrixXd state(6,3);
    for (int d = 0; d < 6; d++)
      for (int i = 0; i < 3; i++)
        state(d, i) = d * i;

    assert(d(state, 0, 1, 0) + d(state, 1, 2, 0) >= d(state, 0, 2, 0));
    assert(d(state, 0, 1, 1) + d(state, 1, 2, 1) >= d(state, 0, 2, 1));
  }

  // test dot product
  {
    MatrixXd s(3, 2);
    for (int i = 0; i < 3; i++) {
      s(i,0) = 1.0;
      s(i,1) = 1.0;
    }
    assert(dot(s,0,1) == 3);
    assert(dot(s,0,0) == 3);
    assert(dot(s,1,1) == 3);
    s(0,1) = 2;
    assert(dot(s,0,1) == 4);
    assert(dot(s,1,1) == 6);
  }

  // test angle function
  for (double theta = 0; theta < 3.0; theta += 0.01) {
    MatrixXd R(3,3);
    R << 
      -1.0, 0.0, cos(theta),
      0.0, 0.0, sin(theta),
      0.0, 0.0, 0.0;
    assert(abs(polymer::angle(R, 0) - theta) < 1.0e-8);
  }

  // test 2nd angle function
  for (double theta = 0; theta < 3.0; theta += 0.01) {
    MatrixXd R(6,3);
    R << 
      // time index 0
      -1.0, 0.0, cos(theta),
      0.0, 0.0, sin(theta),
      0.0, 0.0, 0.0,
      // time index 1
      -1.0, 0.0, cos(theta+0.01),
      0.0, 0.0, sin(theta+0.01),
      0.0, 0.0, 0.0;
    assert(abs(polymer::angle(R, 0, 0) - theta) < 1.0e-8);
    assert(abs(polymer::angle(R, 0, 1) - (theta+0.01)) < 1.0e-8);
  }

  // test polymer generate function
  {
    vector<double> l_vec = {1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
    int Np = 1000;
    MatrixXd pol = MatrixXd::Zero(3, Np);

    for (double l : l_vec) {
      double mu_l = l;
      double sigma_l = l / 5;
      double mu_theta = 0.0;
      double sigma_theta = 0.3;
      polymer::generate(mu_l, sigma_l, mu_theta, sigma_theta, pol);

      double mu_dR_estimate = 0.0;
      double sigma_dR_estimate = 0.0;
      for (int ip = 0; ip < Np-1; ip++)
        mu_dR_estimate += polymer::d(pol, ip, ip+1)/((double)Np-1);
      for (int ip = 0; ip < Np-1; ip++)
        sigma_dR_estimate += pow(polymer::d(pol, ip, ip+1)-mu_l,2.0);
      sigma_dR_estimate = sqrt(sigma_dR_estimate/((double)Np-2.0));

      assert(abs(mu_l-mu_dR_estimate)<2*sigma_l/sqrt(Np-1));
      assert(abs(mu_l-mu_dR_estimate)<2*sigma_dR_estimate/sqrt(Np-1));
    }

    //// The angles should be folded normal distributed
    //// Using the recursive method described here:
    //// https://en.wikipedia.org/wiki/Folded_normal_distribution#Parameter_estimation
    //// sigma^2 = {\sum_{i=1}^n x_i^2 \over n} - \mu^2
    //double mu_theta_estimate = 0.0;
    //double sigma_theta_estimate = 2.0; // guess
    //double sigma_theta_estimate2 = 0.0; // guess
    //do {
    //double sum = 0.0;
    //for (int ip = 0; ip < Np-2; ip++) {
    //sum += pow(polymer::angle(pol, ip), 2.0)/((double)Np-2);
    //}
    //mu_theta_estimate = sqrt(sum-pow(sigma_theta_estimate,2.0));
    //sigma_theta_estimate2 = sum - pow(mu_theta_estimate,2.0);
    //}
    //while (abs(sigma_theta_estimate - sigma_theta_estimate2) < 1.0e-9);

    //cout << mu_theta << endl;
    //cout << mu_theta_estimate << endl;
  }
}
