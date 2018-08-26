#include "polymer.h"

#ifndef _SPRING_
#define _SPRING_
using namespace Eigen;

namespace spring {
  double potential(MatrixXd &R, double k, double l); 
  void force_numerical(MatrixXd &R, double k, double l, MatrixXd &F);  
  void force(MatrixXd &R, double k, double l, MatrixXd &F); 
  void test(); 
}
#endif
