#ifndef EWALD_H
#define EWALD_H
#define DIM 3
using namespace std;


void initialize();
void initialize_readinput(int N, int num_p, int P, double eta_in, 
  int rp, double L, double xi_in);
void realspace();
void realspaceOMP(int num_threads);
void kspace();
void kspaceParallel();
void selfcontribution();
void writeout();
void free();

extern int nx, ny, nz, np, px, py, pz, repeat_x, repeat_y, repeat_z;
extern double Lx, Ly, Lz, xi, eta;

extern double *grid, *particle, *strength, *vel;

extern string outputfile;

#endif
