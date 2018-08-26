#include <Eigen/Dense>
#include <assert.h>
#include <random>

#ifndef _ANALYSIS_
#define _ANALYSIS_
namespace analysis {
  // The following function are for one polymer trial
  void cm(Eigen::MatrixXd &R,Eigen::MatrixXd &cm); 
  void dR(Eigen::MatrixXd &R,Eigen::MatrixXd &dR); 
  void theta(Eigen::MatrixXd &R,Eigen::MatrixXd &theta); 
  void tau(Eigen::MatrixXd &R,Eigen::MatrixXd &dR, Eigen::MatrixXd &tau); 
  //void end_to_end(Eigen::MatrixXd &R, Eigen::MatrixXd &end_to_end);
  void monomer(Eigen::MatrixXd &R,Eigen::MatrixXd &m, int ip);

  class ensemble {
  private:
    // Total number of trials in ensemble N
    int N; 
    // Index of current trial
    int i;
    // Number of particles in trial
    int Np;
    // Number of timesteps
    int Nt;
    // polyfitting the data
    void fit(Eigen::MatrixXd &time);

  public:
    ensemble(int N); 
    void add_trial(Eigen::MatrixXd &R,Eigen::MatrixXd &time);
    Eigen::MatrixXd 
      cm_mean,     // dims: 3 x Nt
      dR_mean,     // dims: (Np-1) x Nt
      theta_mean,  // dims: (Np-2) x Nt
      mm_mean,     // dims: 3 x Nt
      dR_sigma,    // length: 1 x Nt
      cm_sigma,
      theta_sigma,
      tau_mean,
      tau_sigma,
      mm_sigma,
      time;
    Eigen::Vector3d cm_sigma_fit; // dims: 3
  };

  void test();
};
#endif
