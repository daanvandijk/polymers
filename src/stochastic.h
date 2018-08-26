#include "polymer.h"
#include <vector>

#ifndef _STOCHASTIC_
#define _STOCHASTIC_
namespace stochastic {
  using namespace Eigen;
  void generate_thermal(double k_BT, MatrixXd &xi); 
  // CovarianceMatrix Creates a covariance matrix from a vector v
  // For a given vector v = [v_0, ..., v_{n-1}] this function returns a
  // symmetric matrix C 
  // [v_0, v_1, v_2, v_3, ... 
  //  v_1, v_0, v_1, v_2, ...
  //  v_2, v_1, v_0, v_1, ...
  //  .
  //  .
  //  .
  //  v_n, v_{n-1}, ...]
  void covariance_matrix(VectorXd &v, MatrixXd &C);
  // Returns the conditional distribution of Y, i.e. mu_C, Sigma_C,
  // given a realisation of X
  void conditional_normal(
      VectorXd const &mu_X,
      VectorXd const &mu_Y,
      MatrixXd const &Sigma_XX,
      MatrixXd const &Sigma_XY,
      MatrixXd const &Sigma_YY,
      VectorXd const &x,
      VectorXd &mu_C,
      MatrixXd &Sigma_C);
  void test();
  // Generate multivariate normall distribution realisations
  struct normrnd {
    normrnd(VectorXd const &mu, Eigen::MatrixXd const &Sigma); 
    normrnd(MatrixXd const &Sigma) :
      normrnd(VectorXd::Zero(Sigma.rows()), Sigma) { }
    VectorXd mu;
    MatrixXd normal_transform; 
    VectorXd operator()() const;
  };
  // Generate an endless stream of multivariate distributed variables
  // This is done by calculating the conditional normal distribution
  // between two consecutive realisation. 
  // Note: the size of mu should be devisable by 2.
  class normrnd_endless {
    private:
      int i;
      VectorXd previous_realisation;
      VectorXd mu;
      MatrixXd Sigma;

    public:
      normrnd_endless(
          VectorXd const &mu, 
          MatrixXd const &Sigma); 
      normrnd_endless(MatrixXd const &Sigma) :
        normrnd_endless(VectorXd::Zero(Sigma.rows()), Sigma) { }
      MatrixXd normal_transform; 
      double operator()();
  };
  // Estimates the mean mu and the covariance matrix Sigma from a collection of samples
  // m = [v_1; 
  //      v_2;
  //      .
  //      .
  //      v_n]
  void estimate_conditional_norm(MatrixXd &m, 
      VectorXd &mu, 
      MatrixXd &Sigma);

  void generate_active_init(double C, 
      double dt, 
      int Nt,
      int Np,
      std::vector<normrnd_endless> &gens);
  void generate_active(MatrixXd &xi_A, 
      std::vector<normrnd_endless> &gens);
}
#endif 
