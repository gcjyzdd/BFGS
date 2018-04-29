#!/bin/bash

if [ ! -d "build" ];then
    mkdir "build"
fi

cd build
cmake ..
make -j
cd ..
./build/bin/bfgs_test 

