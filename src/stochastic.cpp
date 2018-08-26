#include "stochastic.h"
#include <iostream>
#include <assert.h>
#include <random>

using namespace std;
using namespace Eigen;

void stochastic::generate_thermal(double k_BT, MatrixXd &xi) {
  assert(xi.rows() == 3);
  int Np = xi.cols();

  std::random_device rd;
  std::normal_distribution<double> thermal(0.0, 3*k_BT);

  // loop through dimensions
  for (int id = 0; id < 3; id++) {
    // loop through particles
    for (int ip = 0; ip < Np; ip++) {
      xi(id, ip) = thermal(rd);
    }
  }
}

void stochastic::generate_active_init(double C, 
    double dt, 
    int Nt,
    int Np,
    std::vector<normrnd_endless> &gens) {
  // E[xi_A(t) \cdot xi_A(t')] = 3C \exp(-|t-t'|)
  // The factor 3 is from the number of dimensions
  VectorXd v(Nt);
  for (int it = 0; it < Nt; it++)
    v(it) = C*exp(-it*dt);
  MatrixXd Sigma;
  stochastic::covariance_matrix(v, Sigma); 
  for (int ig = 0; ig < 3*Np; ig++) {
    gens.push_back(stochastic::normrnd_endless{ Sigma });
  }
}

void stochastic::generate_active(MatrixXd &xi_A, 
    std::vector<normrnd_endless> &gens) {
  assert(xi_A.rows() == 3);
  assert(gens.size() == xi_A.cols() * 3);
  for (int ip = 0; ip < xi_A.cols(); ip++) {
    for (int id = 0; id < 3; id++) {
      xi_A(id, ip) = gens[ip + xi_A.cols() * id]();
    }
  }
}

void stochastic::conditional_normal(
    VectorXd const &mu_X,
    VectorXd const &mu_Y,
    MatrixXd const &Sigma_XX,
    MatrixXd const &Sigma_XY,
    MatrixXd const &Sigma_YY,
    VectorXd const &x,
    VectorXd &mu_C,
    MatrixXd &Sigma_C) {

  mu_C = mu_Y + 
    Sigma_XY.transpose() * Sigma_XX.inverse() * (x - mu_X);
  Sigma_C = Sigma_YY -
    Sigma_XY.transpose() * Sigma_XX.inverse() * Sigma_XY;

}

stochastic::normrnd::normrnd(VectorXd const& mu, 
    MatrixXd const& Sigma) : mu(mu) {
  assert(mu.size() == Sigma.rows() && mu.size() == Sigma.cols());
  normal_transform = MatrixXd(mu.size(), mu.size());
  LLT<Eigen::MatrixXd> chol_solver(Sigma);

  if (chol_solver.info() == Success) {
    normal_transform = chol_solver.matrixL();
  }
  else {
    SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(Sigma);
    normal_transform = eigen_solver.eigenvectors() *
      eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
  }
}

VectorXd stochastic::normrnd::operator()() const {
  static std::mt19937 gen{ std::random_device{}() };
  static std::normal_distribution<double> dist;

  return mu + normal_transform * 
    VectorXd{ mu.size() }.
    unaryExpr([&](double x) { return dist(gen); });
}

stochastic::normrnd_endless::normrnd_endless(VectorXd const& mu, 
    MatrixXd const& Sigma) : mu(mu) {
  assert(mu.size() == Sigma.rows() && mu.size() == Sigma.cols());
  assert(mu.size() % 2 == 0);
  normrnd_endless::mu = mu;
  normrnd_endless::Sigma = Sigma;
  i = 0;
}

double stochastic::normrnd_endless::operator()() {
  i++;
  if (previous_realisation.size() == 0) { 
    // first realisation
    i = 0;
    normrnd N(mu.head( mu.size() / 2 ), 
        Sigma.topLeftCorner( mu.size() / 2, mu.size() / 2));
    previous_realisation = N();
  }
  else if (previous_realisation.size() == i) {
    // we need to generate more datapoints
    i = 0;
    VectorXd mu_C;
    MatrixXd Sigma_C;
    conditional_normal(
        mu.head(mu.size()/2),
        mu.tail(mu.size()/2),
        Sigma.topLeftCorner(mu.size()/2,mu.size()/2),
        Sigma.topRightCorner(mu.size()/2,mu.size()/2),
        Sigma.bottomRightCorner(mu.size()/2,mu.size()/2),
        previous_realisation, 
        mu_C,
        Sigma_C);
    normrnd N(mu_C, Sigma_C); 
    previous_realisation = N();
  }
  return previous_realisation(i);
}

void stochastic::covariance_matrix(VectorXd &v, 
    MatrixXd &C) {
  // Init C vector
  C = MatrixXd::Zero(v.size(), v.size()); 

  for (int offset = 0; offset < v.size(); offset++)
    for (int i = 0; i < v.size()-offset; i++) 
    {
      C(i+offset, i) = v(offset);
      C(i, i+offset) = v(offset);
    }
}

void stochastic::estimate_conditional_norm(MatrixXd &m, 
    VectorXd &mu, 
    MatrixXd &Sigma) {
  int L = m.cols();
  int N = m.rows();
  mu = VectorXd(L);
  Sigma = MatrixXd::Zero(L, L);

  // calculate mu
  for (int i = 0; i < L; i++) {
    mu(i) = m.col(i).mean();
  }

  // calculate Sigma
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      for (int n = 0; n < N; n++) {
        Sigma(i, j) += m(n, i) * m(n, j) / ((double)N);
      }
    }
  }
}

void stochastic::test() {
  cout << "Stochastic forces test\n";

  // testing standard normal random variables
  // X = N(0,1) 
  // P(|X| < 2) \approx 5%
  // tsting normal random variables
  // X = N(mu, sigma)
  // P(|X-mu| < 2*sigma) \approx 5%
  // test random forces generation
  {
    int N = 1000;
    double k_BT = 0.5; // (1/3) sigma
    MatrixXd xi(3, N);
    stochastic::generate_thermal(k_BT, xi); // N(0, 3*k_BT) for each element
    double fraction = 0.0;
    for (int d = 0; d < 3; d++) {
      for (int i = 0; i < N; i++) {
        fraction += (abs(xi(d,i)) < 2*3*k_BT) * 1.0/((double)N);
      }
    }
    assert(0.93 < fraction/3.0 && fraction/3.0 < 0.97); 
  }

  // test covariance matrix generation
  {
    for (int N = 1; N < 15; N++) {
      VectorXd v(N);
      for (int i = 0; i < N; i++)
        v(i) = i;
      MatrixXd C;
      covariance_matrix(v, C);
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
          assert(C(i, j) == abs(i-j));
    }
  }

  // test conditional normal
  {
    MatrixXd covar(2,2), Sigma_C;
    VectorXd mu(2), x(1), mu_C;
    covar << 1, .5, .5, 1;
    mu << 0, 0;
    x << 0;
    conditional_normal(
        mu.head(1),
        mu.tail(1),
        covar.topLeftCorner(1, 1),
        covar.bottomLeftCorner(1, 1),
        covar.bottomRightCorner(1, 1),
        x,
        mu_C,
        Sigma_C);
    assert(abs(Sigma_C(0, 0) - 0.75) < 0.001);
    assert(abs(mu_C(0)) < 0.001);
  }

  // test normrnd
  {
    MatrixXd covar;
    VectorXd v(2); 
    v << 1, 0.5;
    covariance_matrix(v, covar);
    normrnd sample { covar };

    int N = 1000;
    MatrixXd m(N, 2);
    for (int i = 0; i < N; i++)
      m.row(i) = sample();
    VectorXd mu;
    MatrixXd Sigma;
    estimate_conditional_norm(m, mu, Sigma);
    assert(abs(mu(0)) < 2/sqrt(N));
    assert(abs(mu(1)) < 2/sqrt(N));
    assert(abs(Sigma(0, 1) - covar(0, 1)) < 4/sqrt(N));
  }

  // test normrnd_endless
  {
    int Nt = 1e5;
    int N = 100;
    double C = 1.0;
    double dt = 1e-6;

    VectorXd v(N);
    for (int it = 0; it < N; it++)
      v(it) = C * exp(-it * dt);
    MatrixXd Sigma;
    covariance_matrix(v, Sigma);

    normrnd_endless endless(Sigma);

    for (int it = 0; it < Nt; it++) {
      assert(abs(endless()) < 1e3);
    }
  }

  // test normrnd_endless
  {
    int M = 1000;
    int N = 4;
    double C = 1.0;
    double dt = 0.1;

    VectorXd v(N);
    for (int it = 0; it < N; it++)
      v(it) = C * exp(-it * dt);
    MatrixXd Sigma;
    covariance_matrix(v, Sigma);

    MatrixXd realisations(M, N);
    for (int m = 0; m < M; m++) {
      normrnd_endless endless(Sigma);

      for (int it = 0; it < N; it++) {
        realisations(m, it) = endless();
        assert(abs(realisations(m, it)) < 1e3);
      }
    }

    VectorXd mu_estimate;
    MatrixXd Sigma_estimate;
    estimate_conditional_norm(realisations, mu_estimate, Sigma_estimate);
    // todo: implement assert tests
    // I'm to lazy now, and my bloodpressure is too high
    //cout << "mu_estimate: \n" << mu_estimate << endl;
    //cout << "Sigma_estimate: \n" << Sigma_estimate << endl;
  }

}
