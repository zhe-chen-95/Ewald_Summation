#include <string>
#include "ewald.h"

int nx, ny, nz, np, px, py, pz, repeat_x, repeat_y, repeat_z;
double Lx, Ly, Lz, xi, eta;

double *grid, *particle, *strength, *vel;

string outputfile;
