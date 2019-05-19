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
  outputfile = "../results/vel.txt"
  cout << "System initialized!" << '\n';
}

void realspace(){
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
