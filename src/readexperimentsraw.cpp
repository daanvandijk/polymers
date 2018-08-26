#include <iostream>
#include <math.h>
#include <algorithm>
#include "mex_miscellaneous.h"
#include "mex.h" // these are from matlab
#include "matrix.h" // these are from matlab
#include "main.h"
#include "rousetransform.h"

using namespace std;

void _main();

const int n_inputs = 2;
const int n_outputs = 4;
string usage = "[parameters, time, R, X] = readexperimentsraw(name, max_trials)";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check for proper number of arguments 
  if (nrhs != n_inputs) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", usage.c_str());
  } 
  else if (nlhs != n_outputs) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", usage.c_str());
  } 

  // Check for proper input type 
  if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )  {
    mexErrMsgIdAndTxt( "MATLAB:mxmalloc:invalidInput", usage.c_str());
  }

  // Get number of characters in the input string.  Allocate enough
  // memory to hold the converted string. 
  size_t buflen = mxGetN(prhs[0]) + 1;
  char *filename = new char[buflen + 10];

  // Copy the string data into buf.  
  int status = mxGetString(prhs[0], filename, (mwSize)buflen);   
  int max_trials = (int) *(double *) mxGetPr(prhs[1]);

  cout << "Name: " << filename << endl;
  cout << "Max trials: " << max_trials << endl;
  sprintf(filename, "%s.dat", filename);
  ifstream is(filename);
  if (!is.is_open()) {
    mexErrMsgIdAndTxt("MATLAB:mexread:file", "Failed to open file.");
  }
  data_header header;
  cereal::PortableBinaryInputArchive archive(is);
  archive(header);

  stringstream out;
  {
    cereal::JSONOutputArchive ar( out );
    ar(cereal::make_nvp("parameters", header.p)); 
  }
  map_string(plhs, 0, out.str().c_str());

  streampos pos = is.tellg();
  for (int i = 0; i < min(header.size, max_trials); i++) {
    polymer_trial trial;
    archive(trial);

    if (i == 0) {
      // all of the experiment time matrices are the same
      map_matrix_out(plhs, 1, trial.time);

      mwSize ndim = 1;
      mwSize *dims = new mwSize[ndim];
      dims[0] = min(header.p.Ntrials, max_trials);
      plhs[2] = mxCreateCellArray(ndim, dims);
      plhs[3] = mxCreateCellArray(ndim, dims);
      delete[] dims;
    }

    // calculate rouse transform
    MatrixXd X = MatrixXd::Zero(trial.R.rows(), trial.R.cols());
    rouse::transform(trial.R, X);

    // map to output 
    map_matrix_out(plhs[2], i, trial.R);
    map_matrix_out(plhs[3], i, X);
  }
  is.close();

  delete[] filename;

  return;
}

