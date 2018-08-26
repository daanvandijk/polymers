#include <iostream>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "analysis.h"
#include "main.h"
#include "mex_miscellaneous.h"

using namespace std;

void _main();

const int n_inputs = 1;
const int n_outputs = 12;
string usage = "[parameters, time, cm_mean, cm_sigma, dR_mean, dR_sigma, theta_mean, theta_sigma, tau_mean, tau_sigma, mm_mean, mm_sigma] = readexperiments(name)";

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
  // memory to hold the converted string. */
  size_t buflen = mxGetN(prhs[0]) + 1;
  char *filename = new char[buflen + 10];

  // Copy the string data into buf.  
  int status = mxGetString(prhs[0], filename, (mwSize)buflen);   

  cout << "Name: " << filename << endl;
  sprintf(filename, "%s.dat", filename);
  ifstream is(filename);
  if (!is.is_open()) {
    mexErrMsgIdAndTxt("MATLAB:mexread:file", "Failed to open file.");
  }
  data_header header;
  MatrixXd time;
  cereal::PortableBinaryInputArchive archive(is);
  archive(header);
  analysis::ensemble ensemble(header.size);

  // We have to read the trials two times.
  // This is because we need to calculate variances,
  // and in order to do so, we need to first calculate mean values. 
  streampos pos = is.tellg();
  for (int i = 0; i < header.size; i++) {
    polymer_trial trial;
    archive(trial);

    if (i == 0) {
      time = trial.time;
    }

    // do analysis
    ensemble.add_trial(trial.R, trial.time);
  }
  // read trials again...
  is.seekg(pos);
  for (int i = 0; i < header.size; i++) {
    polymer_trial trial;
    archive(trial);

    // do analysis
    ensemble.add_trial(trial.R, trial.time);
  }
  is.close();

  stringstream out;
  {
    cereal::JSONOutputArchive ar( out );
    ar(cereal::make_nvp("parameters", header.p)); 
  }

  map_string(plhs, 0, out.str().c_str());
  map_matrix_out(plhs, 1, time);
  map_matrix_out(plhs, 2, ensemble.cm_mean);
  map_matrix_out(plhs, 3, ensemble.cm_sigma);
  map_matrix_out(plhs, 4, ensemble.dR_mean);
  map_matrix_out(plhs, 5, ensemble.dR_sigma);
  map_matrix_out(plhs, 6, ensemble.theta_mean);
  map_matrix_out(plhs, 7, ensemble.theta_sigma);
  map_matrix_out(plhs, 8, ensemble.tau_mean);
  map_matrix_out(plhs, 9, ensemble.tau_sigma);
  map_matrix_out(plhs, 10, ensemble.mm_mean);
  map_matrix_out(plhs, 11, ensemble.mm_sigma); 
  delete[] filename;

  return;
}

