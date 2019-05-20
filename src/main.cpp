#include <cmath>
#include <vector>
#include <ctime>
#include <cstdio>
#include <string>
#include "ewald.h"

using namespace std;

// nvcc -std=c++11 -O0 -g  main.cpp ewald.cu global.cpp -lcufft -o ewald
// nvcc -std=c++11 -O3   main.cpp ewald.cu global.cpp -lcufft -o ewald
int main(int argc, char* argv[])
{
  long tt = clock();
  initialize();
  realspace();
  kspace();
  selfcontribution();
  writeout();
  printf("Total time is %f(s)\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
  return 0;
}
