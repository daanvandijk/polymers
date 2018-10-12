#ifndef _MAIN_
#define _MAIN_
#include <vector>
#include <string>
#include <algorithm>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include "cereal_matrix.h"
#include "polymer.h"
#include "time_integration.h"
#include "analysis.h"

// display help information 
void help();

struct polymer_trial {
  MatrixXd R;
  MatrixXd time;

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
  std::string tex_title; // doesn't get archive'd 

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
        CEREAL_NVP(space));
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
      (space == x.space);
  }
};

// [ data_header | size = s ] [ polymer_trial #1 ] ... [ polymer_trial #s ]
struct data_header {
  experiment_parameters p;
  MatrixXd Ri;
  int size; // num of trials behind header 

  template<class Archive>
  void serialize(Archive &archive) {
    archive(p, Ri, size);
  }
};

std::ostream& operator << (std::ostream &os, const data_header &s);

// the interesting functions
void perform_experiments(const char *experiments_path);
void list_experiments(const char *experiments_path);
void perform_experiment(data_header const &exp, polymer_trial &trial, int i);
analysis::ensemble* read_experiment(const char *path);
polymer_trial* read_experiment_raw(const char* path);
void tex(const char *experiments_path);

std::vector<data_header> get_experiments(const char *experiments_path);
#endif
