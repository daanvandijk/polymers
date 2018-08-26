#include <iostream>
#include "polymer.h"
#include "spring.h"
#include "angle.h"
#include "time_integration.h"
#include "analysis.h"
#include "stochastic.h"
#include "rousetransform.h"

using namespace std;

int main(int argc, char **argv) {
  cout << "Unit tests\n";
  polymer::test();
  spring::test();
  angle::test();
  stochastic::test();
  time_integration::test();
  analysis::test();
  rouse::test();
  
  return 0;
}
