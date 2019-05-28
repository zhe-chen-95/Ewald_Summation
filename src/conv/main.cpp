#include <cmath>
#include <vector>
#include <ctime>
#include <cstdio>
#include <string>
#include <iostream>
#include "ewald.h"

using namespace std;

// nvcc -std=c++11 -O0 -g  main.cpp ewald.cu global.cpp -lcufft -o ewald
// nvcc -std=c++11 -O3   main.cpp ewald.cu global.cpp -lcufft -o ewald
int main(int argc, char* argv[]){
  int num_threads = 1;
  // cout << "Please input num of threads(default=1): " << endl;
	// cin >> num_threads;
  sscanf(argv[1], "%d", &num_threads);
  int N, num_p, P, rp;
  // double L = 1.0, xi_in = 0.2, eta_in = 0.025;
  float L, xi_in, eta_in;
  sscanf(argv[2], "%d", &N);
  sscanf(argv[3], "%d", &num_p);
  // sscanf(argv[3], "%f", &xi_in);
  // sscanf(argv[4], "%f", &eta_in);
  sscanf(argv[4], "%d", &P);
  sscanf(argv[5], "%d", &rp);
  sscanf(argv[6], "%f", &L);
  sscanf(argv[7], "%f", &xi_in);
  sscanf(argv[8], "%f", &eta_in);

  // sscanf(argv[7], "%f", &L);
  printf("Size N = %d with %d particles \n", N, num_p);
  printf("Gaussian P = %d, Real layer rp = %d \n", P, rp);
  printf("Length L = %f, xi = %f, eta = %f", L, xi_in, eta_in);

  long tt = clock();
  cout << "\n---------------initialize---------------" << endl;
  // initialize();
  initialize_readinput(N, num_p, P, eta_in, rp, L, xi_in);
  // cout << "\n---------------Real Space---------------" << endl;
  // realspace();
  cout << "\n---------------Real Space OMP---------------" << endl;
  realspaceOMP(num_threads);
  // cout << "\n---------------K-Space---------------" << endl;
  // kspace();
  cout << "\n---------------K-Space Parallel---------------" << endl;
  kspaceParallel();
  cout << "\n---------------Self Contribution---------------" << endl;
  selfcontribution();
  writeout();
  cout << "\n---------------In Total---------------" << endl;
  printf("Total time is %f(s)\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
  return 0;
}
