#ifndef _EXPERIMENT_
#define _EXPERIMENT_
#include <Eigen/Dense> 
#include <cereal/cereal.hpp>

struct polymer_trial {
  Eigen::MatrixXd R;
  Eigen::MatrixXd time;

  template<class Archive>
    void serialize(Archive &archive) {
      archive(R, time);
    }
};

struct experiment_parameters {
  // title will also be used as filename, 
  // so be careful with [\/.] 
  std::string title;
  int Ntrials;
  int Np;
  int Nt;
  double a;
  double k;
  double l;
  double k_BT;
  double C;
  int datapoints;
  double dt;
  double mu_l;
  double sigma_l;
  double mu_theta;
  double sigma_theta;
  std::string space; // linear or log
  std::string method; // ef, rk4
  std::string tex_title; // doesn't get archive'd 

  bool check() const {
    return (Nt >= 0) && 
      (a >= 0) && 
      (k >= 0) && 
      (l >= 0) && 
      (k_BT >= 0) && 
      (C >= 0) &&
      (datapoints > 0) &&
      (space == "linear" || space == "log") && 
      (method == "ef" || method == "rk4");
  }

  template<class Archive>
    void serialize(Archive &archive) {
      archive(
          CEREAL_NVP(title),
          CEREAL_NVP(Ntrials),
          CEREAL_NVP(Np),
          CEREAL_NVP(Nt),
          CEREAL_NVP(a),
          CEREAL_NVP(k),
          CEREAL_NVP(l),
          CEREAL_NVP(k_BT),
          CEREAL_NVP(C),
          CEREAL_NVP(datapoints),
          CEREAL_NVP(dt),
          CEREAL_NVP(mu_l),
          CEREAL_NVP(sigma_l),
          CEREAL_NVP(mu_theta),
          CEREAL_NVP(sigma_theta),
          CEREAL_NVP(space),
          CEREAL_NVP(method));
    }

  // equality comparison. doesn't modify object. therefore const.
  bool operator==(const experiment_parameters& x) const {
    return 
      (title == x.title) && 
      (Ntrials == x.Ntrials) &&
      (Nt == x.Nt) && 
      (Np == x.Np) && 
      (a == x.a) && 
      (k == x.k) && 
      (l == x.l) && 
      (k_BT == x.k_BT) && 
      (C == x.C) && 
      (datapoints == x.datapoints) && 
      (dt == x.dt) && 
      (mu_l == x.mu_l) && 
      (sigma_l == x.sigma_l) && 
      (mu_theta == x.mu_theta) && 
      (sigma_theta == x.sigma_theta) &&
      (space == x.space) && 
      (method == x.method);
  }
};

// [ data_header | size = s ] [ polymer_trial #1 ] ... [ polymer_trial #s ]
struct data_header {
  experiment_parameters p;
  Eigen::MatrixXd Ri;
  int size; // num of trials behind header 

  template<class Archive>
    void serialize(Archive &archive) {
      archive(p, Ri, size);
    }
};

#endif
