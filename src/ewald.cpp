#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <complex>
#include "ewald.h"

using namespace std;

void initialize(){
  nx = 10;
  ny = 10;
  nz = 10;
  Lx = 10.0;
  Ly = 10.0;
  Lz = 10.0;
  eta = 1.0;
  xi = 1,0;
  np = 100;
  px = 7;
  py = 7;
  pz = 7;
  grid = (double*)calloc(np*DIM,sizeof(double));
  particle = (double*)calloc(np*DIM,sizeof(double));
  strength = (double*)calloc(np*DIM,sizeof(double));
  vel = (double*)calloc(np*DIM,sizeof(double));
  outputfile = "../results/vel.txt";
  cout << "System initialized! # of particles: " << np << '\n';
}


double* realfunc(double x, double y, double z, double xi, double *st1, double *st2){
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double e1 = exp(xi*xi*r2);
  double v[3][2], tmp1[3], tmp2[3];
  double coef[2] = {2*(xi*e1/sqrt(M_PI)/r2+erfc(xi*r)/2/r2/r,4*xi/sqrt(M_PI)*e1};
  double *st = st1;
  tmp1 = {
    (r2+x*x)*st[0]+(x*y)*st[1]+(x*z)*st[2]),
    (y*x)*st[0]+(r2+y*y)*st[1]+(y*z)*st[2],
    (z*x)*st[0]+(z*y)*st[1]+(r2+z*z)*st[2]
  };
  st = st2;
  tmp2 = {
    (r2+x*x)*st[0]+(x*y)*st[1]+(x*z)*st[2]),
    (y*x)*st[0]+(r2+y*y)*st[1]+(y*z)*st[2],
    (z*x)*st[0]+(z*y)*st[1]+(r2+z*z)*st[2]
  };
  v =  {
    coef[0]*tmp1[0]-coef[1]*st1[0],
    coef[0]*tmp1[1]-coef[1]*st1[1],
    coef[0]*tmp1[2]-coef[1]*st1[2],
    coef[0]*tmp2[0]-coef[1]*st2[0],
    coef[0]*tmp2[1]-coef[1]*st2[1],
    coef[0]*tmp2[2]-coef[1]*st2[2]
  };
  return v;



}

void realspace(){
  double rx, ry, rz;
  double *v;
  long tt = clock();
  for (int i = 0; i < np; i++){
    for (int j = i+1; j < np; j++){
      for (int px = -repeat_x; px < repeat_x; px++){
        for (int py = -repeat_y; py < repeat_y; py++){
          for (int pz = -repeat_z; pz < repeat_z; pz++){
            if (px == 0 && py == 0 && pz == 0){
              if (i != j) {
                rx = particle[DIM*j+0]-particle[DIM*i+0];
                ry = particle[DIM*j+1]-particle[DIM*i+1];
                rz = particle[DIM*j+2]-particle[DIM*i+2];
                v = realfunc(rx, ry, rz, xi, &(strength[DIM*i]), &(strength[DIM*j]));
                vel[DIM*i+0] += v[0];
                vel[DIM*i+1] += v[1];
                vel[DIM*i+2] += v[2];

                vel[DIM*j+0] += v[3];
                vel[DIM*j+1] += v[4];
                vel[DIM*j+2] += v[5];
              }
            }
            else{
              rx = particle[DIM*j+0]+Lx*px-particle[DIM*i+0];
              ry = particle[DIM*j+1]+Ly*py-particle[DIM*i+1];
              rz = particle[DIM*j+2]+Lz*pz-particle[DIM*i+2];
              v = realfunc(rx, ry, rz, xi, &(strength[DIM*i]), &(strength[DIM*j]));
              vel[DIM*i+0] += v[0];
              vel[DIM*i+1] += v[1];
              vel[DIM*i+2] += v[2];

              vel[DIM*j+0] += v[3];
              vel[DIM*j+1] += v[4];
              vel[DIM*j+2] += v[5];
            }
          }
        }
      }
    }
  }
  printf("Real space part finished with %ds\n",(clock()-tt)*1.0/CLOCK_PER_SEC);
}

double* Gaussian_Gridding_type1(){

}

double* Gaussian_Gridding_type2(double *Hx){

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

void kspace(){
  double Hx[3*nx*ny*nz];
  complex<double> Hx_hat[3*(nx/2+1)*(ny/2+1)*(nz/2+1)];
  Hx = Gaussian_Gridding_type1();
  Hx_hat =  FFT3D(Hx);

  for (int idim = 0; idim<3; idim++){
    for (int i = 0; i<nx/2+1; i++){
      for (int j = 0; j<ny/2+1; j++){
        for (int k = 0; k<nz/2+1; k++){
          Hx_hat[idim*(nx/2+1)*(ny/2+1)*(nz/2+1)+i*(ny/2+1)*(nz/2+1)+j*(nx/2+1)+k] =
        }
      }
    }
  }

}
void selfcontribution(){
  long tt = clock();
  double tmp = (4*xi)/sqrt(M_PI);
  for (long i = 0; i < np; i++) {
    vel[DIM*i+0] -= tmp*strength[DIM*i+0];
    vel[DIM*i+1] -= tmp*strength[DIM*i+1];
    vel[DIM*i+2] -= tmp*strength[DIM*i+2];
  }
  printf("Self Contribution part finished with %ds\n",(clock()-tt)*1.0/CLOCK_PER_SEC);

}

void writeout(){
  ofstream output(outputfile);
  output << "#Velocity obtained from Ewald summation" << endl;
  output << np << endl;
  for(long i = 0; i < np; i += 1){
      output << i << " " << vel[DIM*i+0] << " " <<
      vel[DIM*i+1] << " " <<
      vel[DIM*i+2]<< "\n" ;
    }
  outfile.close();
  outfile.clear();
  printf("Write files into %s",outputfile);
}

void freeall(){
  free(grid);
  free(particle);
  free(strength);
  free(vel);
  cout << "Dynamical variable destoyed!" << '\n';
}
