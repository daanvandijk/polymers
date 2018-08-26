%{
#include "rousetransform.h"
using namespace rouse;
%}
/*%include "rousetransform.h"*/
extern Eigen::MatrixXd transform(Eigen::MatrixXd &R);
extern Eigen::MatrixXd invtransform(Eigen::MatrixXd &X);
