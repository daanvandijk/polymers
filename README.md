# Polymer simulation

## Goal
Perform long-time stable simulations of the Semiflexible chainmodel. 

## To-do list
- Publish compiled source on website
- Define Rouse timescale with respect to the experiments. 
- Configure comparison images.
- Add generated images to repository
- Add install and execute instructions
- Add results
- Don't use the auto keyword in C++, see [Eigen doc](https://eigen.tuxfamily.org/dox/TopicPitfalls.html)
- ... 

## Compiling
Compiling source code 
- The C++ code is compiled by g++, make sure it's installed.
- Make sure that Eigen and Cereal are installed.
  If you're using Ubuntu this can simply be done by the command `make library_install`.
  Otherwise these library can be found at [1] and [2].
- Make sure the correct `MATLABHOME` path is specified in the makefile

### Libraries
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Cereal](https://uscilab.github.io/cereal/)
