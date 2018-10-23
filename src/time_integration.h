#include "polymer.h"
#include "spring.h"
#include "angle.h"
#include "experiment.h"
#include "stochastic.h"
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

  class integration_method {
    protected:
      experiment_parameters p;
      bool spring_forces, angle_forces, thermal_forces, active_forces;
      Eigen::MatrixXd xi_T, xi_A;
      // active force generators
      std::vector<stochastic::normrnd_endless> gens;

    public:
      integration_method(const experiment_parameters &p); 
      virtual void step(Eigen::MatrixXd &Rold, Eigen::MatrixXd &Rnew)=0;
  };

  class milstein : public integration_method {
    private:

    public:
      milstein(const experiment_parameters &p);
      void step(Eigen::MatrixXd &Rold, Eigen::MatrixXd &Rnew);
  };

  class ef : public integration_method { 
    private:
      Eigen::MatrixXd F;

    public:
      ef(const experiment_parameters &p);
      void step(Eigen::MatrixXd &Rold, Eigen::MatrixXd &Rnew);
  };

  class rk4 : public integration_method {
    private:
      Eigen::MatrixXd k1, k2, k3, k4, buf;
    public:
      rk4(const experiment_parameters &p);
      void step(Eigen::MatrixXd &Rold, Eigen::MatrixXd &Rnew);
      ~rk4();
  };

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
