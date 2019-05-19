#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
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


void kspace(){

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
