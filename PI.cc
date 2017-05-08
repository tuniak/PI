#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

//mriezkou tohto CA je "cubic grid"
//aj ked kockova mriezka nema dostatocnu rotacnu symetriu, zavedenim hybnosti ako samostatnej veliciny, nezavislej na lattic velocity, sme ziskali dostatok stupnov volnosti
//predstavme si uzol ako kocku, kde bunky sidlia v rohoh kocky - mame tede osem buniek
#define A 1
#define B 2
#define C 4
#define D 8
#define E 16
#define F 32
#define G 64
#define H 128

#define allCells (A+B+C+D+E+F+G+H)
#define mirX (A+C+E+G)
#define dirX (B+D+F+H)
#define mirY (A+B+C+D)
#define dirY (E+F+G+H)
#define mirZ (C+D+G+H)
#define dirZ (A+B+E+F)

#define PI 3.141592 //M_PI

#define CIRC 128


//castica v bunke smeruje do uzlu, ktoreho sa tato bunka dotyka (teda do kocky s ktorou susedi rohom tejto bunky)
//siet teda netvori cela kockova mriezka, ale iba kocky, ktore sa dotykaju rohmi, teda "kazda druha kocka"

using namespace std;

//mass and momenta of 8 cells in one node (8 bits in char)
//int o = 0 ... obstacle is not present; int o != 0, obstacle in the node (int o = 1 is obstacle in x direction...)  

unsigned char cell[8] = { A, B, C, D, E, F, G, H };
typedef struct
{
	unsigned char m;
	unsigned char p[3];
//	int o;
} Node;
//na obstacle by stacil char, ale 64 bitove procesory preferuju pracu s 64 bitovymi strukturami, preto sme doplnili Node na 8 bytov (4x char + 1x int = 8 bytes)
//mohli sme stav uzlu reprezentovat aj ako char n[5], kde n[0] by bola hmotnost, n[1],n[2],n[3] by boli hybnosti a n[5] prekazka

// 3 directions, 4 pairs in each direction, 2 cells in each pair 
unsigned char Pair[3][4][2] =
{
	{
		{ A,E },{ B,F },{ C,G },{ D,H }
	},
	{
		{ A,C },{ B,D },{ E,G },{ F,H }
	},
	{
		{ A,B },{ C,D },{ E,F },{ G,H }
	}
};

//collision in the single node
void collision(Node &node)
{
	//mass

	//momentum

	//d,u ... index for downer and upper momentum
	int d, u;

	// l, r ... left/right cell in the pair
	// ml, mr ... particle in left/right cell is present (mass-left, mass-right)
	// lu, ld, ru, rd ... momenta of particles (left-upper, left-downer, right-upper, right-downer)
	unsigned char l, r, ml, mr, lu, ld, ru, rd, li, ri;

	for (int i = 0; i < 3; ++i)
	{
		d = (i + 1) % 3;
		u = (i + 2) % 3;
		for (int j = 0; j < 4; ++j)
		{
			l = Pair[i][j][0];
			r = Pair[i][j][1];

			ml = node.m&l;
			mr = node.m&r;
			ld = node.p[d] & l;
			lu = node.p[u] & l;
			li = node.p[i] & l;
			rd = node.p[d] & r;
			ru = node.p[u] & r;
			ri = node.p[i] & r;

			/* PAIR INTERACTIONS */
			if (!ml && !mr)
				continue;
			if (ml && mr)
			{
				// alternate momenta in direction of the pair-interaction
				if (!ld && !rd)
				{
					node.p[d] |= l;
					node.p[d] |= r;
				}
				else if (ld && rd)
				{
					node.p[d] ^= l;
					node.p[d] ^= r;
				}
				// alternate momenta in all other directions
				if (lu && !ru)
				{
					node.p[u] ^= l;
					node.p[u] |= r;
				}
				else if (!lu && ru)
				{
					node.p[u] |= l;
					node.p[u] ^= r;
				}
				if (li && !ri)
				{
					node.p[i] ^= l;
					node.p[i] |= r;
				}
				else if (!li && ri)
				{
					node.p[i] |= l;
					node.p[i] ^= r;
				}
			}
			else if (!ml && !rd)
			{
				node.m |= l;
				node.m ^= r;
				if (ru)
				{
					node.p[u] |= l;
					node.p[u] ^= r;
				}
				if (ri)
				{
					node.p[i] |= l;
					node.p[i] ^= r;
				}					
			}
			else if (!mr && !ld)
			{
				node.m |= r;
				node.m ^= l;
				if (lu)
				{
					node.p[u] |= r;
					node.p[u] ^= l;
				}
				if (li)
				{
					node.p[i] |= r;
					node.p[i] ^= l;
				}
			}
		}
	}
}

// we go through the whole grid and resolve collisions in nodes
void Collision(Node***array, int X, int Y, int Z, int start)
{
	int x, y, z;

#pragma omp parallel for private (x,y,z)
	for (x = start; x < X; x += 2)
	{
		for (y = start; y < Y; y += 2)
		{
			for (z = start; z < Z; z += 2)
			{
				collision(array[x][y][z]);
			}
		}
	}
}

int PeriodicBC(int n, int N)
{
	// could be also if (n==0)
	if (n<0) return	N - 1;
	else if (n<N) return n;
	// else n==N
	else return 0;
}

// Propagate particles from nodes in array to new nodes in new_array
void Propagation(Node***array, int X, int Y, int Z, int start)
{

#pragma omp parallel for
	//prejdeme kazdy uzol mriezky
	for (int x = start; x < X; x += 2)
	{
		for (int y = start; y < Y; y += 2)
		{
			for (int z = start; z < Z; z += 2)
			{
				// hmotnot uzla (m je char, cize je to 8 bitov - 1 bit pre kazdu bunku v uzli)
				// hybnost uzla, je to trojica charov (tri

				// we look at each cell in node
				for (int k = 1; k <= H; k <<= 1)
				{
					//if there is particle in the cell, we want to propagate it to corresponding node
					if (array[x][y][z].m & k)
					{
						// if particle from cell[c] propagates in positive direction of X
						int xN = k & dirX ? (x + 1) % X : (x - 1 + X) % X;
						int yN = k & dirY ? (y + 1) % Y : (y - 1 + Y) % Y;
						int zN = k & dirZ ? (z + 1) % Z : (z - 1 + Z) % Z;
						// else particle propagates in negative direction of x

#pragma omp atomic
						array[xN][yN][zN].m |= k & array[x][y][z].m;
#pragma omp atomic
						array[xN][yN][zN].p[0] |= k & array[x][y][z].p[0];
#pragma omp atomic
						array[xN][yN][zN].p[1] |= k & array[x][y][z].p[1];
#pragma omp atomic
						array[xN][yN][zN].p[2] |= k & array[x][y][z].p[2];
					}
				}
				//in the old array, we set the node to 0
				array[x][y][z].m = 0;
				array[x][y][z].p[0] = 0;
				array[x][y][z].p[1] = 0;
				array[x][y][z].p[2] = 0;
			}
		}
	}
}


void sphere_to_middle_flow(Node***a, int X, int Y, int Z, int R1, int R2in, int R2out, int start)
{
	//	random_device rd;
	//	mt19937 rng(rd());
	//	uniform_int_distribution<int> uni1(0,3);
	//	uniform_int_distribution<int> uni2(4,7);

	int x, y, z;
	int x1, y1, z1;
	int x2, y2, z2;
	int Rc;

	int is = cos(PI / 8) * R1;
	int isnot = sin(PI / 8) * R1;

	int c;
	char c1, c2;
#pragma omp parallel for private (x, y , z, x1, y1, z1, x2, y2, z2, Rc, c, c1, c2)
	for (x = start; x < X; x += 2)
	{
		x1 = (x - R1);
		x2 = x1*x1;
		for (y = start; y < Y; y += 2)
		{
			y1 = (y - R1);
			y2 = y1*y1;
			for (z = start; z < Z; z += 2)
			{
				z1 = (z - R1);
				z2 = z1*z1;
				Rc = x2 + y2 + z2;

				// we are on the sphere with thickness 1
				if (Rc < R2out && Rc > R2in)
				{
					//x_hp = (x > is);
					//x_mp = (x > )

					//case 1
					if (x1 > is)
					{
						for (c = 0; c < 8; c += 2)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[0] |= cell[c];
						}
					}
					//case 2
					else if (x1 < -is)
					{
						for (c = 1; c < 8; c += 2)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[0] |= cell[c];
						}
					}
					//case 3
					else if (y1 > is)
					{
						for (c = 0; c < 4; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[1] |= cell[c];
						}
					}
					//case 4
					else if (y1 < -is)
					{
						for (c = 4; c < 8; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[1] |= cell[c];
						}
					}
					//case 5
					else if (z1 > is)
					{
						for (c = 2; c < 4; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[2] |= cell[c];
						}
						for (c = 6; c < 8; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[2] |= cell[c];
						}
					}
					//case 6
					else if (z1 < -is)
					{
						for (c = 0; c < 2; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[2] |= cell[c];
						}
						for (c = 4; c < 6; ++c)
						{
							a[x][y][z].m |= cell[c];
							a[x][y][z].p[2] |= cell[c];
						}
					}
					else if (x1 < isnot && x1 > -isnot)
					{
						if (y1 > 0)
						{
							//case 7
							if (z1 > 0)
							{
								c1 = C;
								c2 = D;
							}
							//case 8
							else
							{
								c1 = A;
								c2 = B;
							}
						}
						else
						{
							//case 9
							if (z1 > 0)
							{
								c1 = G;
								c2 = H;
							}
							//case 10
							else
							{
								c1 = E;
								c2 = F;
							}
						}

						a[x][y][z].m |= c1;
						a[x][y][z].m |= c2;
						a[x][y][z].p[1] |= c1;
						a[x][y][z].p[1] |= c2;
						a[x][y][z].p[2] |= c1;
						a[x][y][z].p[2] |= c2;
					}
					else if (y1 < isnot && y1 > -isnot)
					{
						if (x1 > 0)
						{
							//case 11
							if (z1 > 0)
							{
								c1 = C;
								c2 = G;
							}
							//case 12
							else
							{
								c1 = A;
								c2 = E;
							}
						}
						else
						{
							//case 13
							if (z1 > 0)
							{
								c1 = D;
								c2 = H;
							}
							//case 14
							else
							{
								c1 = B;
								c2 = F;
							}
						}

						a[x][y][z].m |= c1;
						a[x][y][z].m |= c2;
						a[x][y][z].p[0] |= c1;
						a[x][y][z].p[0] |= c2;
						a[x][y][z].p[2] |= c1;
						a[x][y][z].p[2] |= c2;
					}
					else if (z1 < isnot && z1 > -isnot)
					{

						if (y1 > 0)
						{
							//15
							if (x1 > 0)
							{
								c1 = A;
								c2 = C;
							}
							//16
							else
							{
								c1 = B;
								c2 = D;
							}
						}
						else
						{
							//17
							if (x1 > 0)
							{
								c1 = E;
								c2 = G;
							}
							//18
							else
							{
								c1 = F;
								c2 = H;
							}
						}

						a[x][y][z].m |= c1;
						a[x][y][z].m |= c2;
						a[x][y][z].p[0] |= c1;
						a[x][y][z].p[0] |= c2;
						a[x][y][z].p[1] |= c1;
						a[x][y][z].p[1] |= c2;
					}
					if (x1 > 0)
					{
						if (y1 > 0)
						{
							if (z1 > 0)
							{
								a[x][y][z].m |= C;
								a[x][y][z].p[0] |= C;
								a[x][y][z].p[1] |= C;
								a[x][y][z].p[2] |= C;
							}
							else
							{
								a[x][y][z].m |= A;
								a[x][y][z].p[0] |= A;
								a[x][y][z].p[1] |= A;
								a[x][y][z].p[2] |= A;
							}
						}
						else
						{
							if (z1 > 0)
							{
								a[x][y][z].m |= G;
								a[x][y][z].p[0] |= G;
								a[x][y][z].p[1] |= G;
								a[x][y][z].p[2] |= G;
							}
							else
							{
								a[x][y][z].m |= E;
								a[x][y][z].p[0] |= E;
								a[x][y][z].p[1] |= E;
								a[x][y][z].p[2] |= E;
							}
						}
					}
					else if (x1 < 0)
					{
						if (y1 > 0)
						{
							if (z1 > 0)
							{
								a[x][y][z].m |= D;
								a[x][y][z].p[0] |= D;
								a[x][y][z].p[1] |= D;
								a[x][y][z].p[2] |= D;
							}
							else
							{
								a[x][y][z].m |= B;
								a[x][y][z].p[0] |= B;
								a[x][y][z].p[1] |= B;
								a[x][y][z].p[2] |= B;
							}
						}
						else
						{
							if (z1 > 0)
							{
								a[x][y][z].m |= H;
								a[x][y][z].p[0] |= H;
								a[x][y][z].p[1] |= H;
								a[x][y][z].p[2] |= H;
							}
							else
							{
								a[x][y][z].m |= F;
								a[x][y][z].p[0] |= F;
								a[x][y][z].p[1] |= F;
								a[x][y][z].p[2] |= F;
							}
						}
					}
				}
			}
		}
	}
}

// we allocate memory for the array, and set all elements to zero
Node*** allocate_grid_array(int X, int Y, int Z)
{
	Node***a = new Node**[X]();
	for (int x = 0; x < X; ++x)
	{
		a[x] = new Node*[Y]();
		for (int y = 0; y < Y; ++y)
		{
			a[x][y] = new Node[Z]();
		}
	}
	return a;
}

// set initial state of the grid
// in the middle of a plain z = const, we put one particle with p[0] = 1; p[1]=p[2]=0;
void set_initial(Node***a, int X, int Y, int Z)
{
	int Xb = X / 3;
	int Xe = 2 * Xb;
	int Yb = Y / 3;
	int Ye = 2 * Yb;
	int Zb = Z / 3;
	int Ze = 2 * Zb;
	for (int x = Xb; x < Xe; x += 2)
		for (int y = Yb; y < Ye; y += 2)
			for (int z = Zb; z < Ze; z += 2)
			{
				a[x][y][z].m = allCells;
				a[x][y][z].p[0] = 0;
				a[x][y][z].p[1] = 0;
				a[x][y][z].p[2] = 0;
			}
}

string write_sphere(int* S, int R, int X, int Y, int Z)
{
	ofstream out;
	string file_name = "sphere";
	out.open(file_name);


	int sX = S[0];
	int sY = S[1];
	int sZ = S[2];


	int x;
	int y1, y2, yS;
	int z1, z2, zS;
	double cos_theta;

	double phi, theta;

	double step = atan(30.0 / R);

	for (theta = 0; theta < PI; theta += step)
	{
		// 1) we take phi < PI, not 2*PI
		//    we use fact that cos(-phi) = cos(phi) and sin(-phi) = -sin(phi)
		//    so x(phi) = x(-phi), and y(phi) = -y(-phi) 
		zS = R*sin(theta);
		z1 = sZ - zS;
		z2 = sZ + zS;

		cos_theta = cos(theta);

		for (phi = 0; phi < PI; phi += step)
		{

			// cos(theta) is saved to variable, so we do not need compute it twice

			x = sX + R*cos(phi)*cos_theta;

			// yS is value of y from the "stredu" S
			yS = R*sin(phi)*cos_theta;
			y1 = sY - yS;
			y2 = sY + yS;


			out << x << "\t" << y1 << "\t" << z1 << endl;
			out << x << "\t" << y1 << "\t" << z2 << endl;
			out << x << "\t" << y2 << "\t" << z1 << endl;
			out << x << "\t" << y2 << "\t" << z2 << endl;

		}
	}
	out.close();
	return file_name;

}


//set nodes to be obstacles at the sphere (and inside) "so stredom s, polomerom R"
string set_and_write_sphere(Node***a, int* s, int R, int X, int Y, int Z)
{
	// we will save coordinates of the sphere, so we can plot it
	ofstream out;
	string file_name = "obstacle";
	out.open(file_name);


	int sX = s[0];
	int sY = s[1];
	int sZ = s[2];

	int x;
	int y1, y2, yS;
	int z1, z2, zS;
	double cos_theta;

	double phi, theta;

	double step = atan(1.0 / R);

	for (theta = 0; theta < PI; theta += step)
	{
		// 1) we take phi < PI, not 2*PI
		//    we use fact that cos(-phi) = cos(phi) and sin(-phi) = -sin(phi)
		//    so x(phi) = x(-phi), and y(phi) = -y(-phi) 
		// 2) we fill also inside of the sphere with obstacles
		//    we use for-cycle in y-direction

		zS = R*sin(theta);
		z1 = sZ - zS;
		z2 = sZ + zS;

		cos_theta = cos(theta);

		for (phi = 0; phi < PI; phi += step)
		{

			// cos(theta) is saved to variable, so we do not need compute it twice

			x = sX + R*cos(phi)*cos_theta;

			// yS is value of y from the "stredu" S
			yS = R*sin(phi)*cos_theta;
			y1 = sY - yS;
			y2 = sY + yS;


			out << x << "\t" << y1 << "\t" << z1 << endl;
			out << x << "\t" << y1 << "\t" << z2 << endl;
			out << x << "\t" << y2 << "\t" << z1 << endl;
			out << x << "\t" << y2 << "\t" << z2 << endl;

			for (int i = y1; i <= y2; ++i)
			{
				//		a[x][i][z1].o = 1;
				//		a[x][i][z2].o = 1;
			}

		}
	}
	out.close();
	return file_name;
}


/*
void set_plate(Node***a, int S, int R, int X, int Y, int Z)
{
int Sx = X / 2;
int Sy = Y / 2;
int z = S;

int R2 = R*R;

int x, y;
int x2, y2;

#pragma omp parallel for private (x,y,x2,y2)
for(x = 0; x < X; ++x)
{
x2 = (x - Sx) * (x - Sx);
for(y = 0; y < Y; ++y)
{
y2 = (y - Sy) * (y - Sy);
if (x2 + y2 < R2)
a[x][y][z].o = 1;
}
}
}
*/
string write_plate(int Sx, int Sy, int z, int X, int Y, int Z, int dx)
{
	ofstream out;
	string file_name = "plate";
	out.open(file_name);

	int R2 = pow(Sx / 2, 2);

	int step = dx / 2;
	int x, y;
	int x2, y2;

	for (x = 0; x < X; x += step)
	{
		x2 = pow(x - Sx, 2);
		for (y = 0; y < Y; y += step)
		{
			y2 = pow(y - Sy, 2);
			if (x2 + y2 < R2)
				out << x << "\t" << y << "\t" << z << endl;
		}
	}

	out.close();
	return file_name;
}
/*
void set_tunnel(Node***a, int X, int Y, int Z)
{
for (int z = 0; z < Z; ++z)
{
for (int y = 0; y < Y; ++y)
{
a[0][y][z].o = 1;
a[X-1][y][z].o = 1;
}
for (int x = 0; x < X; ++x)
{
a[x][0][z].o = 1;
a[x][Y-1][z].o = 1;
}
}
}
*/
/*
void set_sphere(Node***a, int S, int R, int X, int Y, int Z)
{
int Sx = X / 2;
int Sy = Y / 2;
int Sz = S;

int R2 = R * R;

int x, y, z;
int x1, y1, z1;
int x2, y2, z2;

#pragma omp parallel for private (x,y,z,x1,y1,z1,x2,y2,z2)
for(x = 0; x < X; ++x)
{
x1 = x - Sx;
x2 = x1 * x1;
for(y = 0; y < Y; ++y)
{
y1 = y - Sy;
y2 = y1 * y1;
for(z = 0; z < Sz + R + 1; ++z)
{
z1 = z - Sz;
z2 = z1 * z1;
if(x2 + y2 + z2 < R2)
a[x][y][z].o = 1;
}
}
}
}
*/

double **** allocate_velocity_array(int I, int J, int K)
{
	int i, j, k;
	double **** v = new double***[I];
	for (i = 0; i < I; ++i)
	{
		v[i] = new double**[J];
		for (j = 0; j < J; ++j)
		{
			v[i][j] = new double*[K];
			for (k = 0; k < K; ++k)
			{
				v[i][j][k] = new double[3]();
			}
		}
	}
	return v;
}

// we compute macroscopic limit of velocity, X,Y,Z is size of the grid; dx,dy,dz is size of area we use to mean velocity (vystredovat rychlost)

void compute_velocity(Node***array, double****v, double****mean, int dx, int dy, int dz, int I, int J, int K)
{
	ofstream out;
	out.open("Energy", ofstream::out | ofstream::app );
	double Energy = 0;
	double N = dx*dy*dz / 80;

	int i, j, k;
	int x, y, z;

	int c;
	int partnum;

#pragma omp parallel for private (i,j,k,x,y,z,c)
	for (i = 0; i < I; ++i)
	{
		for (j = 0; j < J; ++j)
		{
			for (k = 0; k < K; ++k)
			{
				partnum = 0;
				v[i][j][k][0] = 0;
				v[i][j][k][1] = 0;
				v[i][j][k][2] = 0;
				for (x = i*dx; x < (i + 1)*dx; x += 2)
				{
					for (y = j*dy; y < (j + 1)*dy; y += 2)
					{
						for (z = k*dz; z < (k + 1)*dz; z += 2)
						{
							for (c = 1; c <= H; c <<= 1)
							{
								if (c & array[x][y][z].m)
								{
									++partnum;
									if (c & array[x][y][z].p[0])
									{
										if (c & dirX)
										{
											++v[i][j][k][0];
										}
										else
										{
											--v[i][j][k][0];
										}
									}
									if (c & array[x][y][z].p[1])
									{
										if (c & dirY)
										{
											++v[i][j][k][1];
										}
										else
										{
											--v[i][j][k][1];
										}
									}
									if (c & array[x][y][z].p[2])
									{
										if (c & dirZ)
										{
											++v[i][j][k][2];
										}
										else
										{
											--v[i][j][k][2];
										}
									}
								}
							}
						}
					}
				}

				v[i][j][k][0] /= N;
				v[i][j][k][1] /= N;
				v[i][j][k][2] /= N;

				Energy += 0.0001 * partnum * (v[i][j][k][0]*v[i][j][k][0] + v[i][j][k][1]*v[i][j][k][1] + v[i][j][k][2]*v[i][j][k][2]);

				mean[i][j][k][0] += v[i][j][k][0];
				mean[i][j][k][1] += v[i][j][k][1];
				mean[i][j][k][2] += v[i][j][k][2];
			}
		}
	}
	cout << "Energy is " << Energy << endl;
	out << Energy << endl;
	out.close();
}


// correlation functions
// for constant R == X/4
// we taking discrete parts of angle theta and phi 


void null_covariance_tensor(double*****g, int I, int J, int K)
{
	int i, j, k, d, e;
#pragma omp parallel for private (i,j,k,d,e)
	for (i = 0; i<I; ++i)
		for (j = 0; j<J; ++j)
			for (k = 0; k<K; ++k)
				for (d = 0; d<3; ++d)
					for (e = 0; e<3; ++e)
						g[i][j][k][d][e] = 0;
}

void covariance_tensor(double****v, double*****g, int I, int J, int K)
{
	int i, j, k, d, e;
#pragma omp parallel for private (i,j,k,d,e)
	for (i = 0; i<I; ++i)
		for (j = 0; j<J; ++j)
			for (k = 0; k<K; ++k)
				for (d = 0; d<3; ++d)
					for (e = 0; e<3; ++e)
						g[i][j][k][d][e] += (v[i][j][k][d] * v[i][j][k][e]);
}

void finalize_covariance_tensor(double****mean, double*****g, int I, int J, int K, int p)
{
	int i, j, k, d, e;
#pragma omp parallel for private (i,j,k,d,e)
	for (i = 0; i<I; ++i)
		for (j = 0; j<J; ++j)
			for (k = 0; k<K; ++k)
				for (d = 0; d<3; ++d)
					for (e = 0; e<3; ++e)
						g[i][j][k][d][e] = (g[i][j][k][d][e] / p) - mean[i][j][k][d] * mean[i][j][k][e];
}

void Diagonalize(double**A_, double** Q, double** D_)
{
	// A must be a symmetric matrix.
	// returns Q and D such that 
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	double o[3], m[3];
	//quaternion
	double q[4] = { 0.0,0.0,0.0,1.0 };
	double jr[4];
	double sqw, sqx, sqy, sqz;
	double tmp1, tmp2, mq;
	double AQ[3][3];
	double thet, sgn, t, c;
	for (int i = 0; i < maxsteps; ++i)
	{
		// quat to matrix
		sqx = q[0] * q[0];
		sqy = q[1] * q[1];
		sqz = q[2] * q[2];
		sqw = q[3] * q[3];
		Q[0][0] = (sqx - sqy - sqz + sqw);
		Q[1][1] = (-sqx + sqy - sqz + sqw);
		Q[2][2] = (-sqx - sqy + sqz + sqw);
		tmp1 = q[0] * q[1];
		tmp2 = q[2] * q[3];
		Q[1][0] = 2.0 * (tmp1 + tmp2);
		Q[0][1] = 2.0 * (tmp1 - tmp2);
		tmp1 = q[0] * q[2];
		tmp2 = q[1] * q[3];
		Q[2][0] = 2.0 * (tmp1 - tmp2);
		Q[0][2] = 2.0 * (tmp1 + tmp2);
		tmp1 = q[1] * q[2];
		tmp2 = q[0] * q[3];
		Q[2][1] = 2.0 * (tmp1 + tmp2);
		Q[1][2] = 2.0 * (tmp1 - tmp2);

		// AQ = A * Q
		AQ[0][0] = Q[0][0] * A_[0][0] + Q[1][0] * A_[0][1] + Q[2][0] * A_[0][2];
		AQ[0][1] = Q[0][1] * A_[0][0] + Q[1][1] * A_[0][1] + Q[2][1] * A_[0][2];
		AQ[0][2] = Q[0][2] * A_[0][0] + Q[1][2] * A_[0][1] + Q[2][2] * A_[0][2];
		AQ[1][0] = Q[0][0] * A_[0][1] + Q[1][0] * A_[1][1] + Q[2][0] * A_[1][2];
		AQ[1][1] = Q[0][1] * A_[0][1] + Q[1][1] * A_[1][1] + Q[2][1] * A_[1][2];
		AQ[1][2] = Q[0][2] * A_[0][1] + Q[1][2] * A_[1][1] + Q[2][2] * A_[1][2];
		AQ[2][0] = Q[0][0] * A_[0][2] + Q[1][0] * A_[1][2] + Q[2][0] * A_[2][2];
		AQ[2][1] = Q[0][1] * A_[0][2] + Q[1][1] * A_[1][2] + Q[2][1] * A_[2][2];
		AQ[2][2] = Q[0][2] * A_[0][2] + Q[1][2] * A_[1][2] + Q[2][2] * A_[2][2];
		// D = Qt * AQ
		D_[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
		D_[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
		D_[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
		D_[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
		D_[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
		D_[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
		D_[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
		D_[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
		D_[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];
		o[0] = D_[1][2];
		o[1] = D_[0][2];
		o[2] = D_[0][1];
		m[0] = fabs(o[0]);
		m[1] = fabs(o[1]);
		m[2] = fabs(o[2]);

		k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2; // index of largest element of offdiag
		k1 = (k0 + 1) % 3;
		k2 = (k0 + 2) % 3;
		if (o[k0] == 0.0)
		{
			break;  // diagonal already
		}
		thet = (D_[k2][k2] - D_[k1][k1]) / (2.0*o[k0]);
		sgn = (thet > 0.0) ? 1.0 : -1.0;
		thet *= sgn; // make it positive
		t = sgn / (thet + ((thet < 1.E6) ? sqrt(thet*thet + 1.0) : thet)); // sign(T)/(|T|+sqrt(T^2+1))
		c = 1.0 / sqrt(t*t + 1.0); //  c= 1/(t^2+1) , t=s/c 
		if (c == 1.0)
		{
			break;  // no room for improvement - reached machine precision.
		}
		jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
		jr[k0] = sgn*sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
		jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
		jr[3] = sqrt(1.0f - jr[k0] * jr[k0]);
		if (jr[3] == 1.0)
		{
			break; // reached limits of floating point precision
		}
		q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
		q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
		q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
		q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
		mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		q[0] /= mq;
		q[1] /= mq;
		q[2] /= mq;
		q[3] /= mq;
	}
}




void diagonalize_and_write_covariance_tensor(double*****g, int t, int dx, int dy, int dz, int I, int J, int K)
{
	ofstream dout, tout;

	string tens = "tensor_" + to_string(t) + ".dat";
	dout.open("diagon.dat");
	tout.open(tens);

	int i, j, k;

	double**AA;

	double** Q = new double*[3];
	for (i = 0; i<3; ++i)
		Q[i] = new double[3]();

	double** DD = new double*[3];
	for (i = 0; i<3; ++i)
		DD[i] = new double[3]();

	for (i = 0; i<I; ++i)
		for (j = 0; j<J; ++j)
			for (k = 0; k<K; ++k)
			{
				AA = g[i][j][k];
				Diagonalize(AA, Q, DD);
				dout << dx / 2 + i*dx << "\t" << dy / 2 + j*dy << "\t" << dz / 2 + k*dz << "\t" << DD[0][0] << "\t" << DD[1][1] << "\t" << DD[2][2] << endl;
				tout << dx / 2 + i*dx << "\t" << dy / 2 + j*dy << "\t" << dz / 2 + k*dz << "\t" << AA[0][0] << "\t" << AA[0][1] << "\t" << AA[0][2] << "\t" << AA[1][0] << "\t" << AA[1][1] << "\t" << AA[1][2] << "\t" << AA[2][0] << "\t" << AA[2][1] << "\t" << AA[2][2] << endl;
			}

	dout.close();
	tout.close();
}

// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0

double ** XYZCorrelation(double****v, int I, int J, int K)
{

	double**Cor = new double*[I];
	for (int i = 0; i < I; ++i)
		Cor[i] = new double[3];

	int I2 = I / 2;
	int J2 = J / 2;
	int K2 = K / 2;

	double v_m = v[I2][J2][K2][0];

	int i;
	for (i = 0; i < I; i++)
	{
		Cor[i][0] = pow(v_m - v[i][J2][K2][0], 2);
		Cor[i][1] = pow(v_m - v[I2][i][K2][1], 2);
		Cor[i][2] = pow(v_m - v[I2][J2][i][2], 2);
	}
	return Cor;
}


double ** DiagonalCorrelation(double****v, int I, int J, int K, int N = 64)
{
	double theta, phi;

	int R = I / 2;

	theta = PI / 4;
	phi = PI / 4;

	double**Cor = new double*[R];

	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);

	int min = R - R*cos_theta;
	int max = R + R*cos_theta;


	double* v_s = v[I / 2][J / 2][K / 2];

	double* v_r = new double[3];

	double v_x, v_y, v_z;

	for (int i = min; i < max; i++)
	{


		v_x = pow(v_s[0] - v_r[0], 2) / 2;
		v_y = pow(v_s[1] - v_r[1], 2) / 2;
		v_z = 0;

	}
	return Cor;
}


void SRCorrelation(double****v, double**Cor, int I, int J, int K, int N = CIRC)
{
	int x, y;
	int i, j, k;

	int R = I / 4;

	double sin_theta, sin_phi, cos_theta, cos_phi;

	double theta = 0;
	double phi = 0;
	double step = PI / N;

	double* v_s = v[I / 2][J / 2][K / 2];
	double* v_r = new double[3]();

	//projections
	double v_x, v_y, v_z;


	//#pragma omp parallel for private (i,j,k, x,y, sin_theta,sin_phi, cos_theta,cos_phi, theta, phi, v_r, v_x, v_y, v_z)
	for (x = 0; x < 2 * N; ++x)
	{
		phi += step;
		sin_phi = sin(phi);
		cos_phi = cos(phi);

		for (y = 0; y < N; ++y)
		{
			theta += step;
			sin_theta = sin(theta);
			cos_theta = cos(theta);

			// coordinates of v_r
			i = R*cos_phi*cos_theta;
			j = R*sin_phi*cos_theta;
			k = R*sin_theta;

			v_r = v[i + I / 2][j + J / 2][k + K / 2];

			v_x = (v_s[0] - v_r[0])*cos_phi*cos_theta;
			v_y = (v_s[1] - v_r[1])*sin_phi*cos_theta;
			v_z = (v_s[2] - v_r[2])*sin_theta;

			Cor[x][y] += pow(v_x + v_y + v_z, 2);

		}
	}
}

string write_struct(double**str, int t, int N = CIRC)
{
	string file_name = "struct_" + to_string(t) + ".dat";

	ofstream out;
	out.open(file_name);

	double step = PI / N;
	int n, m;

	for (n = 0; n<2 * N; ++n)
		for (m = 0; m<N; ++m)
			out << step*n << "\t" << step*m << "\t" << str[n][m] << endl;

	out.close();
	return file_name;
}

void plot_struct(string file_name, double**str)
{
	const char* f = file_name.c_str();
	FILE * pipe = popen("gnuplot -persistent", "w");

	fprintf(pipe, "reset\n");
	fprintf(pipe, "set terminal pngcairo\n");
	fprintf(pipe, "set output '%s.png'\n", f);
	fprintf(pipe, "set view equal xyz\n");
	fprintf(pipe, "set xlabel \"phi\"\n");
	fprintf(pipe, "set ylabel \"theta\"\n");
	fprintf(pipe, "set zlabel \"S[(v_r-v_s)^2]\"\n");
	fprintf(pipe, "set xrange [-1:%i]\n", 2 * PI);
	fprintf(pipe, "set yrange [-1:%i]\n", PI);
	fprintf(pipe, "set zrange [-0.1:%i]\n", 0.1);

	fprintf(pipe, "splot \"%s\" u 1:2:3 with lines \n", f);
	pclose(pipe);
}

string write_velocity(double****v, int t, int I, int J, int K, int dx, int dy, int dz, string name = "velocity_")
{
	int i, j, k;
	ofstream out;
	string file_name = name + to_string(t) + ".dat";
	out.open(file_name);
	//#pragma omp parallel for private (i,j,k)
	for (i = 0; i < I; ++i)
		for (j = 0; j < J; ++j)
			for (k = 0; k < K; ++k)
				out << dx / 2 + i*dx << "\t" << dy / 2 + j*dy << "\t" << dz / 2 + k*dz << "\t" << v[i][j][k][0] << "\t" << v[i][j][k][1] << "\t" << v[i][j][k][2] << endl;
	out.close();
	return file_name;
}

void plot(string file_name, int I, int J, int K)
{
	const char* f = file_name.c_str();
	FILE * pipe = popen("gnuplot -persistent", "w");

	fprintf(pipe, "reset\n");
	fprintf(pipe, "set terminal pngcairo\n");
	fprintf(pipe, "set output '%s.png'\n", f);
	fprintf(pipe, "set view equal xyz\n");
	fprintf(pipe, "set xlabel \"X\"\n");
	fprintf(pipe, "set ylabel \"Y\"\n");
	fprintf(pipe, "set zlabel \"Z\"\n");
	fprintf(pipe, "set xrange [-1:%i]\n", I);
	fprintf(pipe, "set yrange [-1:%i]\n", J);
	fprintf(pipe, "set zrange [-1:%i]\n", K);
	fprintf(pipe, "set ticslevel 0\n");
	fprintf(pipe, "splot \"%s\" with vectors\n", f);
	pclose(pipe);
}

void plot(string file_name, int X, int Y, int Z, string obstacle)
{
	const char* f = file_name.c_str();
	const char* o = obstacle.c_str();
	const char* n = (file_name + obstacle).c_str();
	FILE * pipe = popen("gnuplot -persistent", "w");

	fprintf(pipe, "reset\n");
	fprintf(pipe, "set terminal pngcairo\n");
	fprintf(pipe, "set output '%s.png'\n", n);
	fprintf(pipe, "set view equal xyz\n");
	fprintf(pipe, "set xlabel \"X\"\n");
	fprintf(pipe, "set ylabel \"Y\"\n");
	fprintf(pipe, "set zlabel \"Z\"\n");
	fprintf(pipe, "set xrange [0:%i]\n", X);
	fprintf(pipe, "set yrange [0:%i]\n", Y);
	fprintf(pipe, "set zrange [0:%i]\n", Z);
	fprintf(pipe, "set ticslevel 0\n");
	fprintf(pipe, "splot \"%s\" with vectors, \"%s\" u 1:2:3 with points \n", f, o);
	pclose(pipe);
}
void null_struct(double**str, int N = CIRC)
{
	int n, m;
	for (n = 0; n<2 * N; ++n)
		for (m = 0; m<N; ++m)
			str[n][m] = 0;

}
void null_vel(double****v, int I, int J, int K)
{
	int i, j, k, l;

	for (i = 0; i < I; ++i)
		for (j = 0; j < J; ++j)
			for (k = 0; k < K; ++k)
				for (l = 0; l < 3; ++l)
					v[i][j][k][l] = 0;
}

double** alocate_SRC(int I, int J, int K, int N = CIRC)
{
	double**SRC = new double*[2 * N];
	for (int n = 0; n < 2 * N; ++n)
		SRC[n] = new double[N]();
	return SRC;
}

void finalize_mean(double****v, int I, int J, int K, int dx, int dy, int dz, int div)
{
	int i, j, k, l;
#pragma omp parallel for private (i,j,k,l)
	for (i = 0; i<I; ++i)
		for (j = 0; j<J; ++j)
			for (k = 0; k<K; ++k)
				for (l = 0; l<3; ++l)
					v[i][j][k][l] /= div;

}

void finalize_struct(double**str, double****mean, int I, int J, int K, int d, int N = CIRC)
{
	int x, y;
	int i, j, k;

	int R = I / 4;
	int I2 = I / 2;
	int J2 = J / 2;
	int K2 = K / 2;

	double sin_theta, sin_phi, cos_theta, cos_phi;

	double theta = 0;
	double phi = 0;
	double step = PI / N;

	double* v_s = mean[I2][J2][K2];
	double* v_r;

	//projections
	double v_x, v_y, v_z;


	//#pragma omp parallel for private (i,j,k,x,y,sin_theta,sin_phi,cos_theta,cos_phi,theta,phi,v_r,v_x,v_y,v_z)
	for (x = 0; x < 2 * N; ++x)
	{
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		phi += step;

		for (y = 0; y < N; ++y)
		{
			sin_theta = sin(theta);
			cos_theta = cos(theta);
			theta += step;

			// coordinates of v_r
			i = R*cos_phi*cos_theta;
			j = R*sin_phi*cos_theta;
			k = R*sin_theta;

			v_r = mean[i + I2][j + J2][k + K2];

			v_x = (v_s[0] - v_r[0])*cos_phi*cos_theta;
			v_y = (v_s[1] - v_r[1])*sin_phi*cos_theta;
			v_z = (v_s[2] - v_r[2])*sin_theta;

			str[x][y] = (str[x][y] / d) - pow(v_x + v_y + v_z, 2);
		}
	}
}

double***** allocate_Gamma(int I, int J, int K)
{
	int i, j, k, d;

	double***** g = new double****[I];
	for (i = 0; i < I; ++i)
	{
		g[i] = new double***[J];
		for (j = 0; j < J; ++j)
		{
			g[i][j] = new double**[K];
			for (k = 0; k < K; ++k)
			{
				g[i][j][k] = new double*[3];
				for (d = 0; d < 3; ++d)
				{
					g[i][j][k][d] = new double[3]();
				}
			}
		}
	}
	return g;
}

void set_speed(Node***a, int i, int j, int k, int dx, int dy, int dz, double vx, double vy, double vz)
{
	double n = 0;
	int mx = 0;
	int my = 0;
	int mz = 0;

	unsigned char sx = vx >= 0 ? dirX : mirX;
	unsigned char sy = vy >= 0 ? dirY : mirY;
	unsigned char sz = vz >= 0 ? dirZ : mirZ;

	double absx = vx >= 0 ? vx : -vx;
	double absy = vy >= 0 ? vy : -vy;
	double absz = vz >= 0 ? vz : -vz;

	for (int x = i; x < i + dx; x += 2)
	{
		for (int y = j; y < j + dy; y += 2)
		{
			for (int z = k; z < k + dz; z += 2)
			{
				n += 1;
				if (mx / n < absx)
				{
					++mx;
					a[x][y][z].m |= sx;
					a[x][y][z].p[0] |= sx;
				}
				if (my / n < absy)
				{
					++my;
					a[x][y][z].m |= sy;
					a[x][y][z].p[1] |= sy;
				}
				if (mz / n < absz)
				{
					++mz;
					a[x][y][z].m |= sz;
					a[x][y][z].p[2] |= sz;
				}
			}
		}
	}
}

void taylor_green_vortex(Node***a, int I, int J, int K, int X, int Y, int Z, int dx, int dy, int dz)
{
	double ampl = 1.0;
	double kx = 2 * PI / I;
	double ky = 2 * PI / J;
	double kz = 2 * PI / K;

	double ax = 1 * ampl;
	double ay = -1 * ampl;
	double az = 0;
	//double az = -2/ampl;

	double vx, vy, vz;

	for (int i = 0; i < I; ++i)
		for (int j = 0; j < J; ++j)
			for (int k = 0; k < K; ++k)
			{
				// continuum equation div v = 0 implies A*kx + B*ky + C*kz = 0
				// we chose A = 1, B = 1, C = -2 and kx=ky=kz
				vx = ax * sin(kx * i) * cos(ky * j) * cos(kz * k);
				vy = ay * cos(kx * i) * sin(ky * j) * cos(kz * k);
				vz = 0;
				//	vz = az * sin(kx * i) * sin(ky * j) * cos(kz * k);

				set_speed(a, i*dx, j*dy, k*dz, dx, dy, dz, vx, vy, vz);
			}
}

void total_speed(double****v, int I, int J, int K)
{
	double x = 0;
	double y = 0;
	double z = 0;
	for (int i = 0; i < I; ++i)
		for (int j = 0; j < J; ++j)
			for (int k = 0; k < K; ++k)
			{
				x += v[i][j][k][0];
				y += v[i][j][k][1];
				z += v[i][j][k][2];
			}
	cout << x << "    " << y << "    " << z << endl;
}

void flow_in_Z(Node***a, int X, int Y, int Z)
{
	int x, y;
	int z = 1;

#pragma omp parallel for private (x,y)
	for (x = 1; x < X - 1; ++x)
	{
		for (y = 1; y < Y - 1; ++y)
		{
			a[x][y][z].m |= dirZ;
			a[x][y][z].p[2] |= dirZ;
		}
	}
}

int main(int argc, char**argv)
{
	//start measure time
	//size of the grid
	int X = 320;
	int Y = X;
	int Z = X;

	int T = 200;

	//size of area we use to compute macroscopic velocity
	// PLEASE, use integer divisors of X,Y,Z
	int dx = 10;
	int dy = dx;
	int dz = dx;

	// number of areas along X,Y,Z axes
	int I = (X / dx);
	int J = (Y / dy);
	int K = (Z / dz);

	// radius and square of radius of the oblast
	//int R = X/2 - 1;
	//int R2out = R*R;
	//int R2in = (R-2)*(R-2);
	//++R;



	/* ALOCATE ARRAYS */
	Node *** array = allocate_grid_array(X, Y, Z);
	double****velocity = allocate_velocity_array(I, J, K);
	double****mean_vel = allocate_velocity_array(I, J, K);

	//double ** SRC = alocate_SRC(I,J,K);
	//double ***** Gamma = allocate_Gamma(I,J,K);

	/* NAMES OF FILES TO PLOT */
	string file_name, file_velocity, file_mean, file_struct, file_cov;

	// middle point and radius of the spherical obstacle
	//int S[3] = {X/2,Y/2,60};
	//int R = 30;


	/* SET TUNNEL BY OBSTACLES */
	//set_tunnel(array,X,Y,Z);
	//set_tunnel(even,X,Y,Z);

	/* RADIUS AND MIDDLE OF THE SPHERICAL OBSTACLE AND ROUND PLATE */
	//	int R = X/4;

	//	int S = 200;
	//	int Sp[3] = {X/2, Y/2, S};

	//	int R2out = R*R;

	/* SET OBSTACLES */
	//	set_sphere(array,S,R,X,Y,Z);
	//	set_plate(array,S,R,X,Y,Z);

	//	string obstacle;
	//	obstacle = write_sphere(Sp,R,X,Y,Z);
	//	obstacle = write_plate(2*R,2*R,S,X,Y,Z,dx);

//	set_initial(array,X,Y,Z);
	taylor_green_vortex(array, I, J, K, X, Y, Z, dx, dy, dz);

	int start;
	int div;

	time_t START = time(NULL);
	for (int t = 0; t <= T; ++t)
	{
		//cout << "som v kroku " << t << endl;

		start = t & 1;

		if (!(t % 2))
		{
			compute_velocity(array, velocity, mean_vel, dx, dy, dz, I, J, K);
			total_speed(velocity, I, J, K);
			cout << "krok " << t << " sekund " << time(NULL) - START << endl;
			//	SRCorrelation(velocity,SRC,I,J,K);		
			//	covariance_tensor(velocity,Gamma,I,J,K);
			file_name = write_velocity(velocity, t, I, J, K, dx, dy, dz);
			plot(file_name, X, Y, Z);
		}
		//if (!(t%100))
		//	cout << "krok " << t << " sekund " <<  time(NULL) - START << " sekund" << endl;
		/*
		if (t<=5000 && !(t%1000))
		{
		// prvnich 3 tisic kroku se tok ustaluje
		div = 100;
		//
		finalize_mean(mean_vel,I,J,K,dx,dy,dz,div);
		file_mean = write_velocity(mean_vel,t,I,J,K,dx,dy,dz,"mean");
		plot(file_mean,X,Y,Z);

		finalize_struct(SRC,mean_vel,I,J,K,div);
		file_struct = write_struct(SRC,t);
		plot_struct(file_struct,SRC);

		finalize_covariance_tensor(mean_vel,Gamma,I,J,K,div);
		diagonalize_and_write_covariance_tensor(Gamma,t,dx,dy,dz,I,J,K);

		null_vel(mean_vel,I,J,K);
		null_struct(SRC);
		null_covariance_tensor(Gamma,I,J,K);
		}
		if (t == 15000)
		{

		div = 1000;

		finalize_mean(mean_vel,I,J,K,dx,dy,dz,div);
		file_name = write_velocity(mean_vel, t, I,J,K, dx,dy,dz, "mean");
		plot(file_mean,X,Y,Z);

		finalize_struct(SRC,mean_vel,I,J,K,div);
		file_struct = write_struct(SRC,t);
		plot_struct(file_struct,SRC);

		finalize_covariance_tensor(mean_vel,Gamma,I,J,K,div);
		diagonalize_and_write_covariance_tensor(Gamma,t,dx,dy,dz,I,J,K);

		//null_vel(mean_vel,I,J,K);
		//null_struct(SRC);
		//null_covariance_tensor(Gamma,I,J,K);
		}
		*/
		//	sphere_to_middle_flow(array, X, Y, Z, R, R2in, R2out, start);
		//	flow_in_Z(array, X, Y, Z);
		Collision(array, X, Y, Z, start);
		Propagation(array, X, Y, Z, start);
	}
	total_speed(velocity, I, J, K);
	printf("trvalo to %lu sekund\n", time(NULL) - START);
	//	puts("este kreslim obrazky, chvilu strpenie...");
	return 0;
}
