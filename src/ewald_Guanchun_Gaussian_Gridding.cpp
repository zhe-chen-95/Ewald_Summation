#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include "ewald.h"

using namespace std;

void initialize(){
  nx = 100;
  ny = 100;
  nz = 100;
  Lx = 1.0;
  Ly = 1.0;
  Lz = 1.0;
  eta = 1.0;
  xi = 5;
  np = 100;
  // grid = (double*)calloc(np*DIM,sizeof(double));
  particle = (double*) calloc(np*DIM,sizeof(double));
  strength = (double*) calloc(np*DIM,sizeof(double));
  vel = (double*) calloc(np*DIM,sizeof(double));
  // outputfile = "../results/vel.txt"
  // cout << "System initialized!" << '\n';
}

double* Gaussian_Gridding_type1(int Px, int Py, int Pz){
  double* H = (double*) calloc(DIM*nx*ny*nz, sizeof(double));
  double hx = Lx / nx, hy = Ly / ny, hz = Lz / nz;
  double hx_sq = hx * hx, hy_sq = hy * hy, hz_sq = hy * hy;
  int ig, jg, kg;
  double a = - 2 * xi * xi / eta;
  double xp, yp, zp, xp_o, yp_o, zp_o;
  int ip, jp, kp;
  double* E3_x = (double*) calloc(Px+1, sizeof(double));
  double* E3_y = (double*) calloc(Py+1, sizeof(double));
  double* E3_z = (double*) calloc(Pz+1, sizeof(double));
  double E1_x, E1_y, E1_z, E2_x, E2_y, E2_z;
  double* E2_xl = (double*) calloc(2*Px, sizeof(double));
  double* E2_yl = (double*) calloc(2*Py, sizeof(double));
  double* E2_zl = (double*) calloc(2*Pz, sizeof(double));
  double V0, Vx, Vy, Vz;
  for (long i = 0; i <= Px; i++) {
      E3_x[i] = exp(-a*i*i*hx_sq);
  }
  for (long i = 0; i <= Py; i++) {
      E3_y[i] = exp(-a*i*i*hy_sq);
  }
  for (long i = 0; i <= Pz; i++) {
      E3_z[i] = exp(-a*i*i*hz_sq);
  }
  for (long n = 0; n < np; n++){
    xp = particle[DIM*n+0];
    yp = particle[DIM*n+1];
    zp = particle[DIM*n+2];
    ip = xp/hx; jp = yp/hy; kp = zp/hz;
    xp_o = xp - ip*hx; yp_o = yp - jp*hy; zp_o = zp - kp*hz;
    
    E1_x = exp(-a*xp_o*xp_o);
    E1_y = exp(-a*yp_o*yp_o);
    E1_z = exp(-a*zp_o*zp_o);
    E2_x = exp(2*a*xp_o*hx);
    E2_y = exp(2*a*yp_o*hy);
    E2_z = exp(2*a*zp_o*hz);
    for (long i = - Px+1; i <= Px; i++) {
        E2_xl[i+Px-1] = pow(E2_x, i);
    }
    for (long j = -Py+1; j <= Py; j++) {
        E2_yl[j+Py-1] = pow(E2_y, j);
    }
    for (long k = -Pz+1; k <= Pz; k++) {
        E2_zl[k+Pz-1] = pow(E2_z, k);
    }
    V0 = E1_x * E1_y * E1_z;
    for (long i = - Px+1; i <= Px; i++){
      Vx = V0 * E2_xl[i+Px-1] * E3_x[abs(i)];
      for (long j = -Py+1; j <= Py; j++){
        Vy = Vx * E2_yl[j+Py-1] * E3_y[abs(j)];
        for (long k = -Pz+1; k <= Pz; k++){
            Vz = Vy * E2_zl[k+Pz-1] * E3_z[abs(k)];
            ig = (ip+i+nx) % nx; jg = (jp+j+ny) % ny; kg = (kp+k+nz) % nz;
            for (long m = 0; m < DIM; m++){
              H[kg + nz*(jg + ny*(ig + m*nx))] += Vz * strength[DIM*n+m];
            }
          }
        }
      }
    }
  return H;
}

double* Gaussian_Gridding_type2(int Px, int Py, int Pz, double* H){
  double* vel_F = (double*) calloc(np*DIM, sizeof(double));
  double hx = Lx / nx, hy = Ly / ny, hz = Lz / nz;
  double hx_sq = hx * hx, hy_sq = hy * hy, hz_sq = hy * hy;
  int ig, jg, kg;
  double a = - 2 * xi * xi / eta;
  double xp, yp, zp, xp_o, yp_o, zp_o;
  int ip, jp, kp;
  double* E3_x = (double*) calloc(Px+1, sizeof(double));
  double* E3_y = (double*) calloc(Py+1, sizeof(double));
  double* E3_z = (double*) calloc(Pz+1, sizeof(double));
  double E1_x, E1_y, E1_z, E2_x, E2_y, E2_z;
  double* E2_xl = (double*) calloc(2*Px, sizeof(double));
  double* E2_yl = (double*) calloc(2*Py, sizeof(double));
  double* E2_zl = (double*) calloc(2*Pz, sizeof(double));
  double V0, Vx, Vy, Vz;
  for (long i = 0; i <= Px; i++) {
      E3_x[i] = exp(-a*i*i*hx_sq);
  }
  for (long i = 0; i <= Py; i++) {
      E3_y[i] = exp(-a*i*i*hy_sq);
  }
  for (long i = 0; i <= Pz; i++) {
      E3_z[i] = exp(-a*i*i*hz_sq);
  }
  for (long n = 0; n < np; n++){
    xp = particle[DIM*n+0];
    yp = particle[DIM*n+1];
    zp = particle[DIM*n+2];
    ip = xp/hx; jp = yp/hy; kp = zp/hz;
    xp_o = xp - ip*hx; yp_o = yp - jp*hy; zp_o = zp - kp*hz;
    
    E1_x = exp(-a*xp_o*xp_o);
    E1_y = exp(-a*yp_o*yp_o);
    E1_z = exp(-a*zp_o*zp_o);
    E2_x = exp(2*a*xp_o*hx);
    E2_y = exp(2*a*yp_o*hy);
    E2_z = exp(2*a*zp_o*hz);
    for (long i = - Px+1; i <= Px; i++) {
        E2_xl[i+Px-1] = pow(E2_x, i);
    }
    for (long j = -Py+1; j <= Py; j++) {
        E2_yl[j+Py-1] = pow(E2_y, j);
    }
    for (long k = -Pz+1; k <= Pz; k++) {
        E2_zl[k+Pz-1] = pow(E2_z, k);
    }
    V0 = E1_x * E1_y * E1_z;
    for (long i = - Px+1; i <= Px; i++){
      Vx = V0 * E2_xl[i+Px-1] * E3_x[abs(i)];
      for (long j = -Py+1; j <= Py; j++){
        Vy = Vx * E2_yl[j+Py-1] * E3_y[abs(j)];
        for (long k = -Pz+1; k <= Pz; k++){
            Vz = Vy * E2_zl[k+Pz-1] * E3_z[abs(k)];
            ig = (ip+i+nx) % nx; jg = (jp+j+ny) % ny; kg = (kp+k+nz) % nz;
            for (long m = 0; m < DIM; m++){
              vel_F[DIM*n+m] += Vz * H[kg + ny*(jg + nz*(ig + m*nx))];
            }      
        }
      }
    }
  }
  return vel_F;
}

int main(int argc, char* argv[]){
  initialize();
  for (int i = 0; i < np; i++) {
    particle[i*DIM+0] = drand48();
    particle[i*DIM+1] = drand48();
    particle[i*DIM+2] = drand48();
    strength[i*DIM+0] = drand48();
    strength[i*DIM+1] = drand48();
    strength[i*DIM+2] = drand48();
  }
  printf("%f %f \n", particle[0], particle[1]);
  long tt = clock();
  double* H = Gaussian_Gridding_type1(24, 24, 24);
  printf("Type1 Time:    %f     Value: %f %f \n", (clock()-tt)*1.0/CLOCKS_PER_SEC, H[0], H[nx*ny*nz/2]);
  tt = clock();
  double* vel_F = Gaussian_Gridding_type2(24, 24, 24, H);
  printf("Type2 Time:    %f     Value: %f %f \n", (clock()-tt)*1.0/CLOCKS_PER_SEC, vel_F[0], vel_F[3*10]);
  return 0;
}
