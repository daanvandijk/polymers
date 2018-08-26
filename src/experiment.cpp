#include <iostream>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "time_integration.h"
#include "mex_miscellaneous.h"

using namespace std;

void _main();

const int n_inputs = 10;
const int n_outputs = 2;
string usage = "[R, t] = experiment(Nt, dt, a, k, l, k_BT, C, Ri, Nd, space)";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int dimensions, Np;
  int Nt, Nd;
  double dt, a, k, l, k_BT, C;
  std::string space;
  MatrixXd Ri, R, time;

  if (nrhs != n_inputs) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", usage.c_str());
  } 
  else if (nlhs != n_outputs) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", usage.c_str());
  }

  Nt = (int) *(double *) mxGetPr(prhs[0]);
  dt = *(double *) mxGetPr(prhs[1]);
  a =  *(double *) mxGetPr(prhs[2]);
  k = *(double *) mxGetPr(prhs[3]);
  l = *(double *) mxGetPr(prhs[4]);
  k_BT  = *(double *) mxGetPr(prhs[5]);
  C = *(double *) mxGetPr(prhs[6]);
  map_matrix_in(prhs, 7, Ri);
  dimensions = Ri.rows();
  Np = Ri.cols();
  Nd = (int) *(double *) mxGetPr(prhs[8]);
  space = (const char *) mxGetPr(prhs[9]);

  // Check if Ri is valid
  if (Np < 2) 
    mexErrMsgIdAndTxt("MATLAB:mexcpp:Np", "MEXCPP requires Np > 1.");
  else if (dimensions != 3) 
    mexErrMsgIdAndTxt("MATLAB:mexcpp:dimensions", "MEXCPP requires dim = 3.");

  cout << "Np: " << Np << endl;
  cout << "Nt: " << Nt << endl;
  cout << "dt: " << dt << endl;
  cout << "a: " << a << endl;
  cout << "k: " << k << endl;
  cout << "l: " << l << endl;
  cout << "k_B T: " << k_BT << endl;
  cout << "C: " << C << endl;
  cout << "Nd: " << Nd << endl;
  cout << "space: " << space << endl;
  cout << Ri << endl;

  // Init R matrix
  R = MatrixXd::Zero(dimensions*Nd, Np);
  // Init time matrix
  time = MatrixXd::Zero(1, Nd);

  // execute time integration
  time_integration::one_polymer(Nt, dt, a, k, l, k_BT, C, Ri, Nd, R, space, time);  

  // output
  map_matrix_out(plhs, 0, R);
  map_matrix_out(plhs, 1, time);

  return;
}

