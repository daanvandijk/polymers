#include <Eigen/Dense>
#include <assert.h>

#ifndef _ROUSE_
#define _ROUSE_
namespace rouse {
  void transform(Eigen::MatrixXd &R, Eigen::MatrixXd &X);
  Eigen::MatrixXd transform(Eigen::MatrixXd &R);
  void invtransform(Eigen::MatrixXd &X, Eigen::MatrixXd &R);
  Eigen::MatrixXd invtransform(Eigen::MatrixXd &X);
  void test();
};
#endif
