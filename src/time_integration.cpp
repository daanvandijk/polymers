#include "time_integration.h"
#include "stochastic.h"
#include <iostream>
#include <chrono>
#include <random>
#include <stdlib.h>
#include <math.h>

using namespace std;
using namespace Eigen;

time_integration::integration_method::integration_method(const experiment_parameters &p) {
  assert(p.check() == true);
  this->p = p;
  
  // Check if this simulation uses thermal forces
  if (p.k_BT == 0) thermal_forces = false;
  else thermal_forces = true; 
  // Check if this simulation uses active forces
  if (p.C == 0) active_forces = false;
  else active_forces = true; 
  // Check if this simulation uses spring forces
  if (p.k == 0) spring_forces = false;
  else spring_forces = true;
  // Check if this simulation uses angle forces
  if (p.a == 0) angle_forces = false;
  else angle_forces = true;

  if (thermal_forces) {
    xi_T = MatrixXd::Zero(3, p.Np);
  }
  
  // I don't want to use all my mem for active forces :P
  if (active_forces) {
    // todo: M could maybe be smaller
    int M = 100;
    stochastic::generate_active_init(p.C, p.dt, M, p.Np, gens);
    xi_A = MatrixXd::Zero(3, p.Np);
  }
}

void time_integration::datapoints::linspace(double a, double b, Eigen::ArrayXd &space) {
  assert(space.size() > 1);
  int n = space.size();
  double delta = (b-a)/((double)n-1.0);
  for (int i = 0; i < n; i++) {
    space[i] = a + delta*i;
  }
}

void time_integration::datapoints::logspace(double a, double b, Eigen::ArrayXd &space) {
  assert(space.size() > 1);
  assert(a > 0.0 && b > 0.0);
  int n = (int) space.size();
  double factora = log(a);
  double factorb = log(b);
  double delta = (factorb-factora)/((double)n-1.0);
  for (int i = 0; i < n; i++) {
    space[i] = exp(factora + delta*i);
  }
}

void time_integration::one_polymer(const experiment_parameters &p, MatrixXd const &Ri, MatrixXd &R, MatrixXd &time, bool gui_update) {
  std::chrono::time_point<std::chrono::system_clock> start, end; 
  MatrixXd Rold, Rnew;
  double current_time;
  std::chrono::duration<double> elapsed;

  // Check if input is consistent
  assert(p.check() == true);
  assert(Ri.cols() > 0 && Ri.rows() == 3);
  assert(R.rows() == 3 * p.datapoints && R.cols() == Ri.cols());

  // clock for reporting progress
  start = std::chrono::system_clock::now();
  // Do time integration
  Rold = Ri; 
  Rnew = MatrixXd::Zero(3, p.Np); 

  // the index for the datapoints
  time_integration::datapoints points(p.Nt, p.datapoints, p.dt, p.space);

  integration_method *method;
  if (p.method == "rk4") 
    method = new rk4(p);
  else if (p.method == "ef") 
    method = new ef(p);
  else if (p.method == "milstein") 
    method = new milstein(p);
  else {
    fprintf(stderr, "Error: unknown method '%s'\n", p.method.c_str());
    throw;
  }

  // loop through time
  for (int it = 0; it < p.Nt; it++) {
    method->step(Rold, Rnew);

    // check if we need to add another datapoint
    if (points.add_point(it)) {
      // calculate current time
      current_time = it * p.dt;
      // save current time
      time(0, points.get_id()) = current_time;
      // save current state
      polymer::save_state(points.get_id(), R, Rold);

      for (int i = 0; i < Rold.rows(); i++) {
        for (int j = 0; j < Rold.cols(); j++) {
          if (isnan(Rold(i,j)) || isinf(Rold(i,j))) throw std::exception();
        }
      }

      // report progress 
      // (only when the previous report was more than 3 seconds ago)
      end = std::chrono::system_clock::now();
      elapsed = end - start;
      if (elapsed.count() > 0.5) {
        if (gui_update == true) {
          // avoid data race
          mtx.lock();
          printf("\r%.2f%%", round(100 * (double)it/(double)p.Nt)); 
          fflush(stdout);
          mtx.unlock();
        }
        // reset timer
        start = std::chrono::system_clock::now();
      }
      points.propogate(it);
    } 
    Rold = Rnew;
  }
  if (gui_update) printf("\n");

  free(method);
}

time_integration::ef::ef(const experiment_parameters &p) : time_integration::integration_method::integration_method(p) {
  F = MatrixXd::Zero(3, p.Np);
}

void time_integration::ef::step(MatrixXd &Rold, MatrixXd &Rnew) {
  if (spring_forces) spring::force(Rold, p.k, p.l, Rnew);
  if (angle_forces) angle::force(Rold, p.a, Rnew);

  Rnew = Rold + p.dt * F;
  if (thermal_forces == true) {
    stochastic::generate_thermal(p.k_BT, xi_T);
    Rnew += sqrt(p.dt) * xi_T;
  }
  if (active_forces == true) {
    stochastic::generate_active(xi_A, gens);
    Rnew += sqrt(p.dt) * xi_A;
  }
}

time_integration::milstein::milstein(const experiment_parameters &p) : time_integration::integration_method(p) {

}

void time_integration::milstein::step(MatrixXd &Rold, MatrixXd &Rnew) {
  assert(false);
}

time_integration::rk4::rk4(const experiment_parameters &p) : time_integration::integration_method::integration_method(p) {
  k1 = MatrixXd::Zero(3, p.Np); 
  k2 = MatrixXd::Zero(3, p.Np); 
  k3 = MatrixXd::Zero(3, p.Np); 
  k4 = MatrixXd::Zero(3, p.Np); 
}

void time_integration::rk4::step(MatrixXd &Rold, MatrixXd &Rnew) {
  if (spring_forces) spring::force(Rold, p.k, p.l, k1);
  if (angle_forces) angle::force(Rold, p.a, k1);
  k1 = p.dt * k1;
  buf = Rold + 0.5*k1;
  if (spring_forces) spring::force(buf, p.k, p.l, k2);
  if (angle_forces) angle::force(buf, p.a, k2);
  k2 = p.dt * k2;
  buf = Rold + 0.5*k2;
  if (spring_forces) spring::force(buf, p.k, p.l, k3);
  if (angle_forces) angle::force(buf, p.a, k3);
  k3 = p.dt * k3;
  buf = Rold + k3;
  if (spring_forces) spring::force(buf, p.k, p.l, k4);
  if (angle_forces) angle::force(buf, p.a, k4);
  k4 = p.dt * k4;
  Rnew = Rold + (k1 + 2*k2 + 2*k3 + k4)/6; 
  if (thermal_forces == true) {
    stochastic::generate_thermal(p.k_BT, xi_T);
    Rnew += sqrt(p.dt) * xi_T;
  }
  if (active_forces == true) {
    stochastic::generate_active(xi_A, gens);
    Rnew += sqrt(p.dt) * xi_A;
  }

  buf = MatrixXd::Zero(3, p.Np);
}

time_integration::rk4::~rk4() {
  k1.setZero();
  k2.setZero();
  k3.setZero();
  k4.setZero();
}

time_integration::datapoints::datapoints(int Nt, int datapoints, double dt, std::string type) {
  assert(Nt > datapoints);
  ArrayXd space(datapoints);
  if (type == "linear") {
    linspace(0, dt*(Nt-1), space);
  }
  else if (type == "log") {
    logspace(dt, dt*(Nt-1), space);
  }
  else {
    cout << "Error: time_integration::datapoints constructor; undefined type" << endl;
    throw;
  }

  indices = ArrayXi(datapoints);
  for (int i = 0; i < datapoints; i++) {
    indices[i] = round(space[i]/dt);
  }

  indices[0] = 0;

  // make sure that there are no double indices
  for (int i = 0; i < datapoints-1; i++) {
    if (indices[i] == indices[i+1]) {
      for (int j = i+1; indices[i] == indices[j] && j < datapoints; j++) {
        indices[j]++;
      }
    }
  }

  // set next index
  next_it = indices[0];
  next_id = 0;
}

void time_integration::datapoints::propogate(int it) {
  if (add_point(it)) {
    // find next index
    if (next_id < indices.size()-1) {
      next_id++;
    }
    next_it = indices[next_id];
  }
}

void time_integration::test() {
  cout << "Time integration test\n";

  // check if datapoints class is working correctly
  {
    std::string spaces[] = {"linear", "log"};

    for (auto space : spaces) {
      datapoints points(1e3, 10, 0.01, space);
      for (int i = 0; i < 1e3; i++) {
        bool added = points.add_point(i);
        if (i == 0 || i == 1e3-1) assert(added);
        points.propogate(i);
      }
    }
  }

  // check if center of mass is stable
  {
    std::string methods[] = { "rk4", "ef" };

    for (auto method : methods) {
      experiment_parameters p;
      p.Nt = 100;
      p.Np = 10;
      p.datapoints = 40;
      p.dt = 1e-3;
      p.a = 1.0;
      p.k = 1.0;
      p.l = 1.0;
      p.k_BT = 0.0;
      p.C = 0.0;
      p.space = "linear";
      p.method = method;
      p.datapoints = 40;

      MatrixXd Ri(3, p.Np);
      MatrixXd R(3*p.datapoints, p.Np);
      MatrixXd time(1, p.datapoints);
      polymer::generate(1.0, 0.3, 0.0, 0.3, Ri);
      one_polymer(p, Ri, R, time, false);
      MatrixXd CenterOfMass = MatrixXd::Zero(3, p.datapoints);
      // loop through datapoints
      for (int idp = 0; idp < p.datapoints; idp++) {
        // loop through dimensions
        for (int d = 0; d < 3; d++) {
          CenterOfMass(d, idp) = R.row(d+idp*3).mean();
        }
      }
      // calculate variance for all 3 dimensions of CenterOfMass
      // loop through dimensions
      for (int d = 0; d < 3; d++) {
        double mu = CenterOfMass.row(d).mean();
        double var = 0;
        for (int idp = 0; idp < p.datapoints; idp++)
          var += pow(CenterOfMass(d, idp) - mu, 2) / ((double)p.datapoints-1.0);
        assert(var < 1.0e-10);
      }
    }
  }
}
