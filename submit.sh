#!/bin/bash
####################################
#
# Official submission execution code
# for consistent execution
#
# with one argument: the number of processes
#
####################################

# clean directory
rm -rf build/

# the execution according to the documentation
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make install
mpirun -n $1 ./numsim_parallel ../ini/lid_driven_cavity.txt 

# zip the solution
zip -r submission.zip src/ CMakeLists.txt