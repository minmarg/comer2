#!/bin/bash

echo

MYHOMEDEF=${HOME}/local/comer2

read -ep "Enter COMER2 install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: $MYHOME
echo

if [ ! -d build ]; then mkdir build || exit 1; fi
cd build || exit 1


##cmake -DCMAKE_C_COMPILER=/home/mindaugas/bin/gcc -DCMAKE_CXX_COMPILER=/home/mindaugas/bin/c++ -DCMAKE_INSTALL_PREFIX=${HOME}/projects/share/cuda-comer ../src/
cmake -DCMAKE_INSTALL_PREFIX=${MYHOME} \
    -DCMAKE_CUDA_FLAGS="-gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_75,code=compute_75" \
    ../src/  ||  (cd ..; exit 1)

cmake --build . --config Release --target install  ||  (cd ..; exit 1)

cd ..

exit 0

