all:
	g++ -std=c++0x -O3 PI.cc
	{ time ./a.out ; } 2>> time.txt
