#!/bin/bash
g++ -std=c++0x -O3 -fopenmp PI.cc
export OMP_NUM_THREADS=32
./a.out
export OMP_NUM_THREADS=16
./a.out
export OMP_NUM_THREADS=8
./a.out
export OMP_NUM_THREADS=4
./a.out
export OMP_NUM_THREADS=2
./a.out
export OMP_NUM_THREADS=1
./a.out
