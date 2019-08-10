#!/bin/bash

echo

MYHOMEDEF=/mnt${HOME:1}/projects/share/cuda-comer

read -ep "Enter COMER2 install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: $MYHOME
echo

cd build || exit 1


##cmake -DCMAKE_C_COMPILER=/home/mindaugas/bin/gcc -DCMAKE_CXX_COMPILER=/home/mindaugas/bin/c++ -DCMAKE_INSTALL_PREFIX=${HOME}/projects/share/cuda-comer ../src/
#cmake -DCMAKE_INSTALL_PREFIX=/mnt${HOME:1}/projects/share/cuda-comer -DCMAKE_CUDA_FLAGS="-arch=sm_50" ../src/
#cmake -DCMAKE_INSTALL_PREFIX=/mnt${HOME:1}/projects/share/cuda-comer -DCMAKE_CUDA_FLAGS="-gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60" ../src/
cmake -DCMAKE_INSTALL_PREFIX=${MYHOME} \
    -DCMAKE_CUDA_FLAGS="-gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70" \
    ../src/  ||  (cd ..; exit 1)
#consider option: -DCMAKE_VERBOSE_MAKEFILE=ON

cmake --build . --config Release --target install  ||  (cd ..; exit 1)

cd ..

exit 0

