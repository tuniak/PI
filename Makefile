all:
	g++ -std=c++0x -fopenmp PI.cc
	export OMP_NUM_THREADS=2
	time ./a.out > 2s
