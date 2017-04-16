all:
	g++ -std=c++0x -O3 -fopenmp PI.cc
	./a.out
