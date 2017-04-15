all:
	clear
	g++ -std=c++0x -O3  -fopenmp PI.cc
	export OMP_NUM_THREADS=12
	./a.out
