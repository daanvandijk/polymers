#include "polymer.h"
#include "spring.h"
#include "angle.h"
#include <mutex>

#ifndef _TIME_INTEGRATION_
#define _TIME_INTEGRATION_
using namespace Eigen;

namespace time_integration {
  static std::mutex mtx;
  void one_polymer(int Nt, 
      double dt, 
      double a,
      double k,
      double l,
      double k_BT,
      double C,
      MatrixXd const &Ri, 
      int Ndata, 
      MatrixXd &R, 
      MatrixXd &time,
      std::string space,
      bool gui_update); 
  void one_polymer(int Nt, 
      double dt, 
      double a,
      double k,
      double l,
      double k_BT,
      double C,
      MatrixXd const &Ri, 
      int Ndata, 
      MatrixXd &R, 
      std::string space,
      MatrixXd&time); 
  void test(); 
  class datapoints {
    protected:
      Eigen::ArrayXi indices;
      int next_it; // next time index
      int next_id; // next index datapoint
      void linspace(double a, double b, Eigen::ArrayXd &space);
      void logspace(double a, double b, Eigen::ArrayXd &space);
    public:
      datapoints(int Nt, int Npoints, double dt, std::string type);
      bool add_point(int it) const { return next_it == it; }
      void propogate(int it);
      int get_id() { return next_id; }
  };
};

#endif
