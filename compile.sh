#!/bin/bash

# load modules
module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

# clean directory
rm -rf build/

# the execution according to the documentation
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make install
