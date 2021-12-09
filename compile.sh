#!/bin/bash

# clean directory
rm -rf build/

# the execution according to the documentation
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make install
