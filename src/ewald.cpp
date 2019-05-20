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
  eta = 0.5;
  xi = 1.0;
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
  double v[6], tmp1[3], tmp2[3];
  double coef[2] = {2*(xi*e1/sqrt(M_PI)/r2+erfc(xi*r)/2/r2/r,4*xi/sqrt(M_PI)*e1};
  double *st = st1;
  tmp1[0] = (r2+x*x)*st[0]+(x*y)*st[1]+(x*z)*st[2];
  tmp1[1] = (y*x)*st[0]+(r2+y*y)*st[1]+(y*z)*st[2];
  tmp1[2] = (z*x)*st[0]+(z*y)*st[1]+(r2+z*z)*st[2];
  st = st2;
  tmp2[0] = (r2+x*x)*st[0]+(x*y)*st[1]+(x*z)*st[2];
  tmp2[1] = (y*x)*st[0]+(r2+y*y)*st[1]+(y*z)*st[2];
  tmp2[2] = (z*x)*st[0]+(z*y)*st[1]+(r2+z*z)*st[2];
  v[0] = coef[0]*tmp1[0]-coef[1]*st1[0];
  v[1] = coef[0]*tmp1[1]-coef[1]*st1[1];
  v[2] = coef[0]*tmp1[2]-coef[1]*st1[2];
  v[3] = coef[0]*tmp2[0]-coef[1]*st2[0];
  v[4] = coef[0]*tmp2[1]-coef[1]*st2[1];
  v[5] = coef[0]*tmp2[2]-coef[1]*st2[2];
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
                vel[DIM*i+0] += v[3];
                vel[DIM*i+1] += v[4];
                vel[DIM*i+2] += v[5];

                vel[DIM*j+0] += v[0];
                vel[DIM*j+1] += v[1];
                vel[DIM*j+2] += v[2];
              }
            }
            else{
              rx = particle[DIM*j+0]+Lx*px-particle[DIM*i+0];
              ry = particle[DIM*j+1]+Ly*py-particle[DIM*i+1];
              rz = particle[DIM*j+2]+Lz*pz-particle[DIM*i+2];
              v = realfunc(rx, ry, rz, xi, &(strength[DIM*i]), &(strength[DIM*j]));
              vel[DIM*i+0] += v[3];
              vel[DIM*i+1] += v[4];
              vel[DIM*i+2] += v[5];

              vel[DIM*j+0] += v[0];
              vel[DIM*j+1] += v[1];
              vel[DIM*j+2] += v[2];
            }
          }
        }
      }
    }
  }
  printf("Real space part finished with %ds\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
}



double* Gaussian_Gridding_type1(){
  long tt = clock();
  double* H = (double*) calloc(DIM*nx*ny*nz, sizeof(double));
  double hx = Lx / nx, hy = Ly / ny, hz = Lz / nz;
  double hx_sq = hx * hx, hy_sq = hy * hy, hz_sq = hy * hy;
  int ig, jg, kg;
  double a = - 2 * xi * xi / eta;
  double xp, yp, zp, xp_o, yp_o, zp_o;
  int ip, jp, kp;
  double* E3_x = (double*) calloc(px+1, sizeof(double));
  double* E3_y = (double*) calloc(py+1, sizeof(double));
  double* E3_z = (double*) calloc(pz+1, sizeof(double));
  double E1_x, E1_y, E1_z, E2_x, E2_y, E2_z;
  double* E2_xl = (double*) calloc(2*px, sizeof(double));
  double* E2_yl = (double*) calloc(2*py, sizeof(double));
  double* E2_zl = (double*) calloc(2*pz, sizeof(double));
  double V0, Vx, Vy, Vz;
  for (long i = 0; i <= px; i++) {
    E3_x[i] = exp(-a*i*i*hx_sq);
  }
  for (long i = 0; i <= py; i++) {
    E3_y[i] = exp(-a*i*i*hy_sq);
  }
  for (long i = 0; i <= pz; i++) {
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
    for (long i = - px+1; i <= px; i++) {
      E2_xl[i+px-1] = pow(E2_x, i);
    }
    for (long j = -py+1; j <= py; j++) {
      E2_yl[j+py-1] = pow(E2_y, j);
    }
    for (long k = -pz+1; k <= pz; k++) {
      E2_zl[k+pz-1] = pow(E2_z, k);
    }
    V0 = E1_x * E1_y * E1_z;
    for (long i = - px+1; i <= px; i++){
      Vx = V0 * E2_xl[i+px-1] * E3_x[abs(i)];
      for (long j = -py+1; j <= py; j++){
        Vy = Vx * E2_yl[j+py-1] * E3_y[abs(j)];
        for (long k = -pz+1; k <= pz; k++){
          Vz = Vy * E2_zl[k+pz-1] * E3_z[abs(k)];
          ig = (ip+i+nx) % nx; jg = (jp+j+ny) % ny; kg = (kp+k+nz) % nz;
          for (long m = 0; m < DIM; m++){
            H[kg + nz*(jg + ny*(ig + m*nx))] += Vz * strength[DIM*n+m];
          }
        }
      }
    }
  }
  printf("Gaussian_Gridding_type1 finished with %ds\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
  return H;
}

void Gaussian_Gridding_type2(double* H){
  long tt = clock();
  double scale_factor = hx*hy*hz * pow(2*xi*xi/(M_PI*eta), 1.5);
  double* vel_F = (double*) calloc(np*DIM, sizeof(double));
  double hx = Lx / nx, hy = Ly / ny, hz = Lz / nz;
  double hx_sq = hx * hx, hy_sq = hy * hy, hz_sq = hy * hy;
  int ig, jg, kg;
  double a = - 2 * xi * xi / eta;
  double xp, yp, zp, xp_o, yp_o, zp_o;
  int ip, jp, kp;
  double* E3_x = (double*) calloc(px+1, sizeof(double));
  double* E3_y = (double*) calloc(py+1, sizeof(double));
  double* E3_z = (double*) calloc(pz+1, sizeof(double));
  double E1_x, E1_y, E1_z, E2_x, E2_y, E2_z;
  double* E2_xl = (double*) calloc(2*px, sizeof(double));
  double* E2_yl = (double*) calloc(2*py, sizeof(double));
  double* E2_zl = (double*) calloc(2*pz, sizeof(double));
  double V0, Vx, Vy, Vz;
  for (long i = 0; i <= px; i++) {
    E3_x[i] = exp(-a*i*i*hx_sq);
  }
  for (long i = 0; i <= py; i++) {
    E3_y[i] = exp(-a*i*i*hy_sq);
  }
  for (long i = 0; i <= pz; i++) {
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
    for (long i = - px+1; i <= px; i++) {
      E2_xl[i+px-1] = pow(E2_x, i);
    }
    for (long j = -py+1; j <= py; j++) {
      E2_yl[j+py-1] = pow(E2_y, j);
    }
    for (long k = -pz+1; k <= pz; k++) {
      E2_zl[k+pz-1] = pow(E2_z, k);
    }
    V0 = E1_x * E1_y * E1_z;
    for (long i = - px+1; i <= px; i++){
      Vx = V0 * E2_xl[i+px-1] * E3_x[abs(i)];
      for (long j = -py+1; j <= py; j++){
        Vy = Vx * E2_yl[j+py-1] * E3_y[abs(j)];
        for (long k = -pz+1; k <= pz; k++){
          Vz = Vy * E2_zl[k+pz-1] * E3_z[abs(k)];
          ig = (ip+i+nx) % nx; jg = (jp+j+ny) % ny; kg = (kp+k+nz) % nz;
          for (long m = 0; m < DIM; m++){
            vel[DIM*n+m] += scale_factor * Vz * H[kg + ny*(jg + nz*(ig + m*nx))];
          }
        }
      }
    }
  }
  printf("Gaussian_Gridding_type2 finished with %ds\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
}

complex<double> FFT3D(double *H){
  cufftHandle plan;
  cufftDoubleReal *data;
  cufftDoubleComplex *data1;
  int n[DIM] = {nx, ny, nz};
  cudaMalloc((void**)&data, sizeof(cufftDoubleReal)*nx*ny*nz*3);
  cudaMalloc((void**)&data1, sizeof(cufftDoubleComplex)*nx*ny*(nz/2+1)*3);
  cudaMemcpy(data,H,nx*ny*nz*3*sizeof(cufftDoubleReal));
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return;
  }
  /* Create a 3D FFT plan. */
  if (cufftPlanMany(&plan, DIM, n,
    NULL, 1, nx*ny*nz, // *inembed, istride, idist
    NULL, 1, nx*ny*nz, // *onembed, ostride, odist
    CUFFT_D2Z, 3) != CUFFT_SUCCESS){
      fprintf(stderr, "CUFFT error: Plan creation failed");
      return;
    }
    /* Use the CUFFT plan to transform the signal in place. */
    if (cufftExecD2Z(plan, data, data1) != CUFFT_SUCCESS){
      fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
      return;
    }
    if (cudaDeviceSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize\n");
      return;
    }
    cufftDestroy(plan);
    cudaFree(data);
    cudaFree(data1);
}

complex<double> IFFT3D(complex<double> *H){

}

void kspace(){
  long tt = clock();
  double Hx[3*nx*ny*nz];
  complex<double> Hx_hat[3*(nx)*(ny)*(nz/2+1)];
  Hx = Gaussian_Gridding_type1();
  Hx_hat =  FFT3D(Hx);
  complex<double> Hx_tilde[3*(nx)*(ny)*(nz/2+1)];
  double kx, ky, kz, k2, e1;
  double kx0=2*M_PI/Lx, ky0=2*M_PI/Ly, kz0=2*M_PI/Lz;
  for (int i = 0; i<nx; i++){
    for (int j = 0; j<ny; j++){
      for (int k = 0; k<nz/2+1; k++){
        if (i<nx/2){
          kx = i*kx0;
        }
        else{
          kx = (i-nx)*kx0;
        }
        if (j<ny/2){
          ky = j*ky0;
        }
        else{
          ky = (j-ny)*ky0;
        }
        if (k<nz/2){
          kz = k*kz0;
        }
        else{
          kz = (k-nz)*kz0;
        }
        k2 = kx*kx + ky*ky + kz*kz;
        e1 = exp(-(1-eta)*k2/4/xi/xi);
        Hx_tilde[0*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k] = e1*
        8*M_PI*(1+k2/4/xi/xi)/k2/k2*
        ((k2-kx*kx)*Hx_hat[0*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (-kx*ky)*Hx_hat[1*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (-kx*kz)*Hx_hat[2*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]);

        Hx_tilde[1*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k] = e1*
        8*M_PI*(1+k2/4/xi/xi)/k2/k2*
        ((-ky*kx)*Hx_hat[0*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (k2-ky*ky)*Hx_hat[1*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (-ky*kz)*Hx_hat[2*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]);

        Hx_tilde[2*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k] = e1*
        8*M_PI*(1+k2/4/xi/xi)/k2/k2*
        ((-kz*kx)*Hx_hat[0*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (-kz*ky)*Hx_hat[1*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]+
        (k2-kz*kz)*Hx_hat[2*(nx)*(ny)*(nz/2+1)+i*(ny)*(nz/2+1)+j*(nz/2+1)+k]);
      }
    }
  }
  Hx = IFFT3D(Hx_tilde);
  Gaussian_Gridding_type2(Hx);
  free(Hx);
  printf("k-pace part finished with %ds\n",(clock()-tt)*1.0/CLOCK_PER_SEC);
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
