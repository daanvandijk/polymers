#include "analysis.h"
#include "polymer.h"
#include <iostream>

using namespace std;
using namespace Eigen;

void analysis::cm(MatrixXd &R, MatrixXd &cm) {
  assert(R.rows() % 3 == 0);
  int Nt = R.rows() / 3;
  int Np = R.cols();
  cm = MatrixXd::Zero(3, Nt);

  // loop through time
  for (int it = 0; it < Nt; it++) {
    // loop through dimensions
    for (int dim = 0; dim < 3; dim++) {
      cm(dim, it) = R.row(dim + it * 3).mean();
    }
  }
}

void analysis::dR(MatrixXd &R, MatrixXd &dR) {
  assert(R.rows() % 3 == 0);
  int Nt = R.rows() / 3;
  int Np = R.cols();
  dR = MatrixXd::Zero(Np-1, Nt);

  // loop through time
  for (int it = 0; it < Nt; it++) {
    // loop through particles
    for (int ip = 0; ip < Np-1; ip++) {
      dR(ip, it) = polymer::d(R,ip,ip+1,it); 
    }
  }
}

void analysis::theta(MatrixXd &R, MatrixXd &theta) {
  assert(R.rows() % 3 == 0);
  assert(R.cols() >= 3);
  int Nt = R.rows() / 3;
  int Np = R.cols();
  theta = MatrixXd::Zero(Np-2, Nt);

  // loop through time
  for (int it = 0; it < Nt; it++) {
    // loop thourgh particles
    for (int ip = 0; ip < Np-2; ip++) {
      theta(ip, it) = polymer::angle(R, ip, it);
    }
  }
}

void analysis::tau(MatrixXd &R, MatrixXd &dR, MatrixXd &tau) {
  assert(R.rows() % 3 == 0);
  int Nt = R.rows() / 3;
  int Np = R.cols();

  tau = MatrixXd(1, Nt);
  // loop through time
  for (int it = 0; it < Nt; it++)
    tau(0, it) = polymer::d(R, 0, Np-1, it) / dR.col(it).sum(); 
}

analysis::ensemble::ensemble(int N) {
  ensemble::N = N;
  ensemble::i = 0;
}

void analysis::monomer(MatrixXd &R, MatrixXd &m, int ip) {
  // R dimensions: (3*Nt) x Np
  // m dimensions: 3 x Nt
  int Np = R.cols();
  int Nt = R.rows() / 3;
  assert(0 <= ip && ip < Np);

  m = MatrixXd(3, Nt);

  // loop through time
  for (int it = 0; it < Nt; it++) {
    m(0, it) = R(3*it, ip);
    m(1, it) = R(3*it+1, ip);
    m(2, it) = R(3*it+2, ip);
  } 
}

void analysis::ensemble::add_trial(MatrixXd &R, MatrixXd &time) {
  assert(R.rows() % 3 == 0);
  // If this is the first trial, 
  // we have to figure out the dimensions
  if (i == 0) {
    ensemble::Nt = R.rows() / 3;
    ensemble::Np = R.cols();

    // Copy time matrix
    ensemble::time = time;

    // Init matrices
    cm_mean = MatrixXd::Zero(3, Nt);
    cm_sigma = MatrixXd::Zero(1, Nt);
    dR_mean = MatrixXd::Zero(Np-1, Nt);
    dR_sigma = MatrixXd::Zero(1, Nt);
    theta_mean = MatrixXd::Zero(Np-2, Nt);
    theta_sigma = MatrixXd::Zero(1, Nt);
    tau_mean = MatrixXd::Zero(1, Nt);
    tau_sigma = MatrixXd::Zero(1, Nt);
    mm_mean = MatrixXd::Zero(3, Nt);
    mm_sigma = MatrixXd::Zero(1, Nt);
  } 
  assert(R.rows() == 3 * Nt);
  assert(R.cols() == Np);
  assert(0 <= i && i < 2*N);

  // Do work...
  {
    MatrixXd cm, dR, theta, tau, mm;
    analysis::cm(R, cm);
    analysis::dR(R, dR);
    analysis::theta(R, theta);
    analysis::tau(R, dR, tau);
    analysis::monomer(R, mm, Np / 2);
    if (0 <= i && i < N) {
      cm_mean += cm / (double) N;
      dR_mean += dR / (double) N;
      theta_mean += theta / (double) N;
      tau_mean += tau / (double) N;
      mm_mean += mm / (double) N;
    }
    else if (N <= i && i < 2*N) {
      // sigma_cm^2 = \langle (X-mu_X)^2 \rangle
      // with X-mu_X = R_cm(t) - R_cm(0)
      for (int k = 0; k < Nt; k++) {
        cm_sigma(0, k) = 
          (cm_mean.col(0) - cm.col(k)).array().pow(2.0) // dims: 3 x 1
          .sum(); // dims: 1 x 1
      }
      cm_sigma += (cm_mean - cm).array().pow(2.0) // dims: 3 x Nt
        .matrix().colwise().sum(); // dims: 1 x Nt 
      dR_sigma += (dR_mean - dR).array().pow(2.0) // dims: (Np-1) x Nt
        .matrix().colwise().sum(); // dims: 1 x Nt 
      theta_sigma += (theta_mean - theta).array().pow(2.0) // dims: (Np-2) x Nt
        .matrix().colwise().sum(); // dims: 1 x Nt 
      tau_sigma += (tau_mean - tau).array().pow(2.0).matrix(); // dims: 1 x Nt
      mm_sigma += (mm_mean - mm).array().pow(2.0) // dims: 3 x Nt
        .matrix().colwise().sum(); // dims: 1 x Nt
    }
  }

  // last timestep
  if (i == 2*N-1) {
    cm_sigma = cm_sigma.array().sqrt().matrix() / sqrt((double)N-1.0);
    dR_sigma = cm_sigma.array().sqrt().matrix() / sqrt((double)N-1.0);
    theta_sigma = cm_sigma.array().sqrt().matrix() / sqrt((double)N-1.0);
    tau_sigma = cm_sigma.array().sqrt().matrix() / sqrt((double)N-1.0);
    mm_sigma = mm_sigma.array().sqrt().matrix() / sqrt((double)N-1.0);

    // do the fitting
    fit(time);
  }

  // Increment index trial
  i++;
}

void analysis::ensemble::fit(MatrixXd &time) {
  // fitting of center of mass
  // polynomial of order 2, i.e.
  // y = c_1 x^2 + c_2 x + c_3
  // c = (X^T X)^(-1) X^T y
  // X = [1 x_1 x_1^2 ...
  //      1 x_2 x_2^2
  //      .
  //      .
  //      1 x_n x_n^2 ...]
  // see https://en.wikipedia.org/wiki/Polynomial_regression#Matrix_form_and_calculation_of_estimates  
  MatrixXd X = MatrixXd::Ones(Nt, 3);
  for (int i = 0; i < Nt; i++) {
    for (int j = 1; j < 3; j++) {
      X(i,j) = pow(time(i), j);
    }
  }

  //cout << "X: " << X.rows() << " x " << X.cols() << endl;
  //cout << "cm_sigma: " << cm_sigma.rows() << " x " << cm_sigma.cols() << endl;
  cm_sigma_fit = (X.transpose() * X).inverse() * X.transpose() * 
    cm_sigma.transpose();
  //printf("cm_sigma = %e + %e t + %e t^2\n", 
      //cm_sigma_fit(0), 
      //cm_sigma_fit(1), 
      //cm_sigma_fit(2));
}

void analysis::test() {
  cout << "Analysis test\n";
  // test cm and dR
  {
    MatrixXd R(6, 4);
    R << 
      -2, -1, 1, 2,
      -2, -1, 1, 2,
      -2, -1, 1, 2,
      -1, 0, 1, 2,
      -1, 0, 1, 2,
      -1, 0, 1, 2;
    MatrixXd cm, dR;
    analysis::cm(R, cm);
    analysis::dR(R, dR);
    assert(cm(0, 0) == 0 && cm(1, 0) == 0 && cm(2, 0) == 0);
    assert(cm(0, 1) == 0.5 && cm(1, 1) == 0.5 && cm(2, 1) == 0.5);
    assert(dR(0, 0) == sqrt(3) && dR(0,1) == sqrt(3));
  }

  // test theta
  for (double theta = 0.0; theta < 3.0; theta += 0.01) {
    MatrixXd R(6, 3);
    R << 
      -1, 0, cos(theta),
      0, 0, sin(theta),
      0, 0, 0,  
      -1, 0, cos(theta+0.01),
      0, 0, sin(theta+0.01),
      0, 0, 0; 
    MatrixXd t;
    analysis::theta(R, t);
    assert(abs(t(0,0) - theta) < 0.001 && abs(t(0,1) - theta-0.01) < 0.001);
  }

  // test tau
  {
    MatrixXd R(6, 3);
    R <<
      0, 1, 2,
      0, 0, 0,
      0, 0, 0,
      0, 1, 2,
      0, 0.5, -0.5,
      0, 0, 0;
    MatrixXd dR, tau;
    analysis::dR(R, dR);
    analysis::tau(R, dR, tau);
    assert(tau(0,0) == 1.0 && 0 <= tau(0,1) && tau(0, 1) < 1.0);
  }

  // test the ensemble class
  {
    int N = 1000;
    ensemble ensemble(N);

    // 2 timesteps, 3 particles
    MatrixXd R = MatrixXd::Zero(3 * 2, 3);
    MatrixXd time = MatrixXd::Zero(1, 2);
    time(0, 0) = 0;
    time(0, 1) = 1;
    double theta = 1.0;
    double delta = 0.05;
    double sigma = 1.0;
    R << 
      // timestep 1:
      -1.0, 0.0, cos(theta),
      0.0, 0.0, sin(theta),
      0.0, 0.0, 0.0,
      // timestep 2:
      -1.0, 0.0, cos(theta),
      0.0, 0.0, sin(theta),
      0.0, 0.0, 0.0;
    R.block(3, 0, 3, 3).array() += delta;

    // init random variables
    std::random_device rd;
    //std::mt19937 generator;
    std::normal_distribution<double> distr_cm(0.0, sigma);

    for (int i = 0; i < 2*N; i++) {
      MatrixXd R_ = R.array() + distr_cm(rd);
      ensemble.add_trial(R_, time);
    }

    // test for the mean cm-values
    // sqrt(n) law 95%-confidence interval
    // loop through time
    for (int it = 0; it < 2; it++) {
      // loop through dimensions
      for (int id = 0; id < 3; id++) {
        double mu = R.row(id+3*it).mean();
        assert(abs(mu-ensemble.cm_mean(id,it)) < 2*sigma/sqrt(N));
      }
    }
  }

  // test polymer generate function
  {
    double mu_l = 1e-7;
    double sigma_l = 3e-8;
    double  mu_theta = 0.0;
    double sigma_theta = 0.01;

    int Np = 10;
    int N = 100;
    ensemble ensemble(N);
    for (int i = 0; i < 2*N; i++) {
      MatrixXd pol = MatrixXd::Zero(3, Np);
      MatrixXd time = MatrixXd::Zero(1,1);
      polymer::generate(mu_l, sigma_l, mu_theta, sigma_theta, pol);
      ensemble.add_trial(pol, time);
    }

    for (int i = 0; i < Np-1; i++)
      assert(abs(mu_l-ensemble.dR_mean(i,0)) <= 2*sqrt(Np-1)*sigma_l/sqrt(N));
    for (int i = 0; i < Np-1; i++)
      assert(abs(mu_l-ensemble.dR_mean(i,0)) <= 2*sqrt(Np-1)*ensemble.dR_sigma(0,0)/sqrt(N));
  }
}

