#ifndef EWALD_H
#define EWALD_H
#define DIM 3
using namespace std;


void initialize();
void realspace();
void kspace();
void selfcontribution();
void writeout();
void free();

extern int nx, ny, nz, np, px, py, pz, repeat_x, repeat_y, repeat_z;
extern double Lx, Ly, Lz, xi, eta;

extern double *grid, *particle, *strength, *vel;

extern string outputfile;

#endif
