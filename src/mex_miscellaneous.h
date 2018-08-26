#ifndef _MEX_MISCELLANEOUS_
#define _MEX_MISCELLANEOUS_
#include "mex.h"
#include "matrix.h"
#include <Eigen/Dense>

using namespace Eigen;

// Map an Eigen::MatrixXd to the n-th output of the mex function 
// Note that Matlab stores 2D matrices by row.
// So it first copies the first row, then the second, etc.
void map_matrix_out(mxArray *plhs[], int n, MatrixXd &m) {
  plhs[n] = mxCreateUninitNumericMatrix(m.rows(), m.cols(), mxDOUBLE_CLASS, mxREAL); 
  double *data = mxGetPr(plhs[n]);
  // loop through cols
  for (int j = 0; j < m.cols(); j++)
    // loop through rows
    for (int i = 0; i < m.rows(); i++)
      data[j * m.rows() + i] = m(i, j);
}

void map_matrix_out(mxArray *plhs[], int n, VectorXd &m) {
  plhs[n] = mxCreateUninitNumericMatrix(m.rows(), m.cols(), mxDOUBLE_CLASS, mxREAL); 
  double *data = mxGetPr(plhs[n]);
  // loop through cols
  for (int j = 0; j < m.cols(); j++)
    // loop through rows
    for (int i = 0; i < m.rows(); i++)
      data[j * m.rows() + i] = m(i, j);
}

// Map an Eigen::MatrixXd to the n-th cell 
void map_matrix_out(mxArray *cellp, int n, MatrixXd &m) {
  mxArray *p = mxCreateUninitNumericMatrix(m.rows(), m.cols(), mxDOUBLE_CLASS, mxREAL); 
  double *data = mxGetPr(p);
  // loop through cols
  for (int j = 0; j < m.cols(); j++)
    // loop through rows
    for (int i = 0; i < m.rows(); i++)
      data[j * m.rows() + i] = m(i, j);
  mxSetCell(cellp, n, p);
}

// Map a matlab matrix to an Eigen::MatrixXd
void map_matrix_in(const mxArray *prhs[], int n, MatrixXd &m) {
  size_t cols = mxGetN(prhs[n]);
  size_t rows = mxGetM(prhs[n]);
  double * _m = mxGetPr(prhs[n]);
  // Init m matrix
  m = MatrixXd(rows, cols);
  // loop throufh cols
  for (int j = 0; j < cols; j++) {
    // loop through rows
    for (int i = 0; i < rows; i++) {
      m(i, j) = _m[j * rows + i]; 
    }
  } 
}

void map_string(mxArray *plhs[], int n, const char *str) {
  plhs[n] = mxCreateString(str);
}
#endif
