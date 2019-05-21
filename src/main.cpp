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
  cout << "Please input num of threads(default=1): " << endl;
	cin >> num_threads;
  long tt = clock();
  cout << "\n---------------initialize---------------" << endl;
  initialize();
  cout << "\n---------------Real Space---------------" << endl;
  realspace();
  cout << "\n---------------Real Space OMP---------------" << endl;
  realspaceOMP(num_threads);
  cout << "\n---------------K-Space---------------" << endl;
  kspace();
  cout << "\n---------------K-Space Parallel---------------" << endl;
  kspaceParallel();
  cout << "\n---------------Self Contribution---------------" << endl;
  selfcontribution();
  writeout();
  cout << "\n---------------In Total---------------" << endl;
  printf("Total time is %f(s)\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
  return 0;
}
