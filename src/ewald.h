#ifndef EWALD_H
#define EWALD_H
#define DIM 3
using namespace std;

double norm(double x,double y,double z);

double dp(double x,double y,double z,double xx,double yy,double zz);

void initialize();
void realspace();
void kspace();
void selfcontribution();
void writeout();
void free();

int nx, ny, nz, np;
double Lx, Ly, Lz, xi, eta;

double *grid, *particle, *strength, *vel;

char outputfile[256];

#endif
