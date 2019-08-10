#!/bin/bash

echo

MYHOMEDEF=${HOME}/local/comer2

read -ep "Enter COMER2 install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: $MYHOME
echo

cd build || exit 1


cmake -DCMAKE_INSTALL_PREFIX=${MYHOME} \
    -DCMAKE_CUDA_FLAGS="-gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70" \
    ../src/  ||  (cd ..; exit 1)
#option: -DCMAKE_VERBOSE_MAKEFILE=ON

cmake --build . --config Release --target install  ||  (cd ..; exit 1)

cd ..

exit 0

