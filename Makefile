all:
	g++ -std=c++0x -Ofast PI.cc
	{ time ./a.out ; } 2>> time.txt
