#include "time_integration.h"
#include "stochastic.h"
#include <iostream>
#include <chrono>
#include <random>
#include <stdlib.h>

using namespace std;

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
  int n = space.size();
  double factora = log(a);
  double factorb = log(b);
  double delta = (factorb-factora)/((double)n-1.0);
  for (int i = 0; i < n; i++) {
    space[i] = exp(factora + delta*i);
  }
}

void time_integration::one_polymer(int Nt, 
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
    MatrixXd &time)
{
  one_polymer(Nt, dt, a, k, l, k_BT, C, Ri, Ndata, R, time, space, true);
}

void time_integration::one_polymer(int Nt, 
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
    bool gui_update) {

  int Np;
  std::chrono::time_point<std::chrono::system_clock> start, end; 
  MatrixXd Rold, Rnew, k1, k2, k3, k4, buf, xi_T, xi_A;
  double current_time;
  std::chrono::duration<double> elapsed;
  bool spring_forces, angle_forces, thermal_forces, active_forces;

  // Check if input is consistent
  assert(Ri.cols() > 0 && Ri.rows() == 3);
  assert(Nt > 0);
  assert(a >= 0);
  assert(k >= 0);
  assert(l > 0);
  assert(k_BT >= 0);
  assert(C >= 0);
  assert(Ndata > 1);
  assert(R.rows() == 3 * Ndata && R.cols() == Ri.cols());

  // Check if this simulation uses thermal forces
  if (k_BT == 0) thermal_forces = false;
  else thermal_forces = true; 
  // Check if this simulation uses active forces
  if (C == 0) active_forces = false;
  else active_forces = true; 
  // Check if this simulation uses spring forces
  if (k == 0) spring_forces = false;
  else spring_forces = true;
  // Check if this simulation uses angle forces
  if (a == 0) angle_forces = false;
  else angle_forces = true;
  // Number of particles
  Np = Ri.cols();
  // clock for reporting progress
  start = std::chrono::system_clock::now();
  // Do time integration
  // The scheme used here is RungeKutta4
  Rold = Ri; 
  Rnew = MatrixXd::Zero(3, Np); 
  k1 = MatrixXd::Zero(3, Np); 
  k2 = MatrixXd::Zero(3, Np); 
  k3 = MatrixXd::Zero(3, Np); 
  k4 = MatrixXd::Zero(3, Np); 
  buf = MatrixXd::Zero(3, Np);
  if (thermal_forces)
    xi_T = MatrixXd::Zero(3, Np);

  // the index for the datapoints
  time_integration::datapoints points(Nt, Ndata, dt, space);

  // I don't want to use all my mem for active forces :P
  vector<stochastic::normrnd_endless> ActiveForceGenerators;
  if (active_forces) {
    // todo: M could maybe be smaller
    int M = 100;
    stochastic::generate_active_init(C, dt, M, Np, ActiveForceGenerators);
    xi_A = MatrixXd::Zero(3, Np);
  }

  // loop through time
  for (int it = 0; it < Nt; it++) {
    if (spring_forces) spring::force(Rold, k, l, k1);
    if (angle_forces) angle::force(Rold, a, k1);
    k1 = dt * k1;
    buf = Rold + 0.5*k1;
    if (spring_forces) spring::force(buf, k, l, k2);
    if (angle_forces) angle::force(buf, a, k2);
    k2 = dt * k2;
    buf = Rold + 0.5*k2;
    if (spring_forces) spring::force(buf, k, l, k3);
    if (angle_forces) angle::force(buf, a, k3);
    k3 = dt * k3;
    buf = Rold + k3;
    if (spring_forces) spring::force(buf, k, l, k4);
    if (angle_forces) angle::force(buf, a, k4);
    k4 = dt * k4;
    Rnew = Rold + (k1 + 2*k2 + 2*k3 + k4)/6; 
    if (thermal_forces == true) {
      stochastic::generate_thermal(k_BT, xi_T);
      Rnew += sqrt(dt) * xi_T;
    }
    if (active_forces == true) {
      stochastic::generate_active(xi_A, ActiveForceGenerators);
      Rnew += sqrt(dt) * xi_A;
    }

    // check if we need to add another datapoint
    if (points.add_point(it)) {
      // calculate current time
      current_time = it * dt;
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
          printf("\r%.2f%%", round(100 * (double)it/(double)Nt)); 
          fflush(stdout);
          mtx.unlock();
        }
        // reset timer
        start = std::chrono::system_clock::now();
      }
      points.propogate(it);
    } 
    Rold = Rnew;
    k1.setZero(); 
    k2.setZero(); 
    k3.setZero(); 
    k4.setZero(); 
  }
  if (gui_update) printf("\n");
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
    int Nd = 40;
    int Np = 10;
    MatrixXd Ri(3, Np);
    MatrixXd R(3*Nd, Np);
    MatrixXd time(1, Nd);
    polymer::generate(1.0, 0.3, 0.0, 0.3, Ri);
    one_polymer(100, 0.001, 1.0, 1.0, 1.0, 0.0, 0.0, Ri, Nd, R, time, "linear", false);
    MatrixXd CenterOfMass = MatrixXd::Zero(3, Nd);
    // loop through datapoints
    for (int idp = 0; idp < Nd; idp++) {
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
      for (int idp = 0; idp < Nd; idp++)
        var += pow(CenterOfMass(d, idp) - mu, 2) / ((double)Nd- 1.0);
      assert(var < 1.0e-10);
    }
  }
}
