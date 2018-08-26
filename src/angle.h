#include "polymer.h"

#ifndef _ANGLE_
#define _ANGLE_
using namespace Eigen;
namespace angle {
  double f(MatrixXd &R, int n); 
  void Fn_n1(MatrixXd &R, int n, double a, MatrixXd &F); 
  void Fn_0(MatrixXd &R, int n, double a, MatrixXd &F); 
  void Fn_p1(MatrixXd &R, int n, double a, MatrixXd &F); 
  double potential(MatrixXd &R, double a); 
  void force(MatrixXd &R, double a, MatrixXd &F); 
  void force_numerical(MatrixXd &R, double a, MatrixXd &F); 
  void test(); 
}
#endif
