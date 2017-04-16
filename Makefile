all:
	g++ -std=c++0x -O3 -fopenmp PI.cc
	{ time ./a.out ; } 2>> time.txt
