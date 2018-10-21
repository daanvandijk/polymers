#include "polymer.h"
#include "spring.h"
#include "angle.h"
#include "experiment.h"
#include <mutex>

#ifndef _TIME_INTEGRATION_
#define _TIME_INTEGRATION_
namespace time_integration {
  static std::mutex mtx;
  void one_polymer(
      const experiment_parameters &p,
      MatrixXd const &Ri, 
      MatrixXd &R, 
      Eigen::MatrixXd&time,
      bool gui_update); 
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
