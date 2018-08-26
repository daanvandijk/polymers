#include "angle.h"
#include <iostream>
#include <assert.h>

using namespace std;
using namespace polymer;

// define a function f such that:
// f(n) := 1 - cos theta_n
//       = 1 - {u_n \cdot u_{n+1} \over |u_n| \cdot |u_{n+1}|}
// with: 
// u_n := R_n - R_{n-1}
//
// Note that:
// f(n-1) := R_{n-2},R_{n-1},R_n \to [0,2]
// f(n) := R_{n-1},R_n,R_{n+1} \to [0,2]
// f(n+1) := R_n,R_{n+1},R_{n+2} \to [0,2]
double angle::f(MatrixXd &R, int n) {
  assert(n >= 1 && n <= R.cols() - 2);
  assert(R.rows() == 3);
  return 1.0 - (R.col(n)-R.col(n-1)).dot(R.col(n+1)-R.col(n)) / 
    sqrt((R.col(n)-R.col(n-1)).squaredNorm()) /
    sqrt((R.col(n+1)-R.col(n)).squaredNorm());
}

void angle::Fn_n1(MatrixXd &R, int n, double a, MatrixXd &F) {
  assert(R.cols() == F.cols());
  assert(R.rows() == F.rows());
  assert(n-2 >= 0 && n < R.cols());
  assert(R.rows() == 3);
  assert(a > 0);
  // Original matlab code:
  // force = (dot(Rnm1,Rnm1)-dot(Rnm2,Rnm1))*(Rnm1-Rn)/polymer.d(Rn,Rnm1)^3+...
  //    (Rnm2-Rnm1)/polymer.d(Rn,Rnm1) + ...
  //    dot(Rnm2-Rnm1,Rn) * (Rnm1-Rn) / polymer.d(Rn,Rnm1)^3;
  // force = force / polymer.d(Rnm1,Rnm2);
  //
  // This function works a little bit different than the original matlab function,
  // and the original definition in the pdf. Here, we incorporate the factor
  // from the function f immediately. This is convenient because of the data
  // structure now. 
  // loop through dimensions
  F.col(n) -= 
    a*f(R,n-1)*
    (
     (dot(R,n-1,n-1)-dot(R,n-2,n-1))*(R.col(n-1)-R.col(n))/pow(d(R,n,n-1),3)+
     (R.col(n-2)-R.col(n-1))/d(R,n,n-1)+
     (dot(R,n-2,n)-dot(R,n-1,n))*(R.col(n-1)-R.col(n))/pow(d(R,n,n-1),3)
    )/
    d(R,n-1,n-2); 
}

void angle::Fn_p1(MatrixXd &R, int n, double a, MatrixXd &F) {
  assert(R.cols() == F.cols());
  assert(R.rows() == F.rows());
  assert(n >= 0 && n+2 < R.cols());
  assert(R.rows() == 3);
  assert(a > 0);
  // Original matlab code:
  // force = polymer.Fn_n1(Rnp2, Rnp1, Rn);
  //
  // This function works a little bit different than the original matlab function,
  // and the original definition in the pdf. Here, we incorporate the factor
  // from the function f immediately. This is convenient because of the data
  // structure now. 
  //
  // Note that we must give an explicit definition of Fn_p1 now, in contrast
  // to the original Matlab code. But the following code is exactly the same as
  // Fn_m1 but most of the minus-signs are replaced with plus+signs. 
  // R_{n-2} becomes R_{n+2}
  // R_{n-1} becomes R_{n+1}
  // R_n stays R_n
  F.col(n) -= 
    a*f(R,n+1)*
    (
     (dot(R,n+1,n+1)-dot(R,n+2,n+1))*(R.col(n+1)-R.col(n))/pow(d(R,n,n+1),3)+
     (R.col(n+2)-R.col(n+1))/d(R,n,n+1)+
     (dot(R,n+2,n)-dot(R,n+1,n))*(R.col(n+1)-R.col(n))/pow(d(R,n,n+1),3)
    )/
    d(R,n+1,n+2); 
}

void angle::Fn_0(MatrixXd &R, int n, double a, MatrixXd &F) {
  assert(R.cols() == F.cols());
  assert(R.rows() == F.rows());
  assert(n-1 >= 0 && n+1 < R.cols());
  assert(R.rows() == 3);
  assert(a > 0);
  // Original matlab code:
  // force = 2*Rn - Rnm1 - Rnp1 ...
  //      + ((Rnp1 - Rn)/polymer.d(Rn,Rnp1)^2 ...
  //      + (Rnm1 - Rn)/polymer.d(Rn,Rnm1)^2) * ...
  //      (dot(Rnm1,Rnp1) + dot(Rn,Rn) - dot(Rn,Rnm1+Rnp1));
  // force = force / (polymer.d(Rn,Rnm1) * polymer.d(Rn,Rnp1));
  //
  // This function works a little bit different than the original matlab function,
  // and the original definition in the pdf. Here, we incorporate the factor
  // from the function f immediately. This is convenient because of the data
  // structure now. 
  F.col(n) -=
    a*f(R,n)*
    (
     2*R.col(n)-R.col(n-1)-R.col(n+1)+
     (
      (R.col(n+1)-R.col(n))/pow(d(R,n,n+1),2)+
      (R.col(n-1)-R.col(n))/pow(d(R,n,n-1),2)
     )*
     (dot(R,n-1,n+1)+dot(R,n,n)-dot(R,n,n-1)-dot(R,n,n+1))
    )/
    (d(R,n,n-1)*d(R,n,n+1));
}

void angle::force(MatrixXd &R, double a, MatrixXd &F) {
  // loop through particles
  for (int ip = 0; ip < R.cols(); ip++) {
    // check if theta_{ip-1} exists
    if (ip-2 >= 0 && ip < R.cols()) {
      // Fn_n1	
      Fn_n1(R, ip, a, F); 
    }
    // check if theta_{ip} exists
    if (ip-1 >= 0 && ip+1 < R.cols()) {
      // Fn_0
      Fn_0(R, ip, a, F); 
    }
    // check if theta_{ip+1} exists
    if (ip >= 0 && ip+2 < R.cols()) {
      // Fn_p1
      Fn_p1(R, ip, a, F); 
    }
  }
}

double angle::potential(MatrixXd &R, double a) {
  assert(a > 0);
  double flag = 0.0;
  // loop through angles theta_n 
  for (int n = 1; n < R.cols()-1; n++) {
    flag += 0.5 * a * pow(f(R, n),2);
  }
  return flag;
}

void angle::force_numerical(MatrixXd &R, double a, MatrixXd &F) {
  assert(a > 0);
  MatrixXd R2 = R;
  double delta = 1e-9;
  assert(compare(R, R2) < 1e-9);
  // loop through particles
  for (int ip = 0; ip < R.cols(); ip++) {
    // loop through dimensions
    for (int id = 0; id < R.rows(); id++) {
      R2(id, ip) += delta;
      F(id, ip) = -(potential(R2, a) - potential(R, a))/delta;
      R2(id, ip) = R(id, ip);
    }
  }
}

void angle::test() {
  cout << "Angle unit test\n";

  // test the f function
  // and test Fn_n1
  {
    MatrixXd s(3, 3);

    for (double theta = 0.0; theta < 3.0; theta += 0.01) {
      s(0, 0) = -1.0;
      s(1, 0) = 0.0;
      s(2, 0) = 0.0;

      s(0, 1) = 0.0;
      s(1, 1) = 0.0;
      s(2, 1) = 0.0;

      s(0, 2) = cos(theta);
      s(1, 2) = sin(theta);
      s(2, 2) = 0.0;

      assert(abs(f(s, 1) - (1 - cos(theta))) < 1e-9);
    }
  }

  // test force 
  // compare numerical differentation of the potential with the
  // analytical force
  {
    MatrixXd s(3, 10), F1 = MatrixXd::Zero(3, 10), F2 = MatrixXd::Zero(3, 10);
    generate(1.0, 0.3, 0.0, 0.3, s);
    force(s, 1.0, F1);
    force_numerical(s, 1.0, F2);
    assert(compare(F1, F2) < 1e-6);
  }
}
