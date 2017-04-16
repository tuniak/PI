all:
	clear
	g++ -std=c++0x -fopenmp -O3 PI.cc
	export OMP_NUM_THREADS=16
	./a.out
