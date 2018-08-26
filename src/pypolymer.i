/* pypolymer.i */
%module pypolymer
%{
/* Put header files here or function declarations like below */
#include "main.h"
%}
%include autodiff.i
%include analysis.i
%include polymer.i
%include rousetransform.i
%include "main.h"
/*extern void perform_experiments();*/
/*extern void list_experiments();*/
/*extern analysis::ensemble* read_experiment(char *title);*/
/*extern struct experiment_parameters;*/
