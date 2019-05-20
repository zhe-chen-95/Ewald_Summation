#include <cmath>
#include <vector>
#include <ctime>
#include <cstdio>
#include <string>
#include "ewald.h"

using namespace std;


int main(int argc, char* argv[])
{
  long tt = clock();
  initialize();
  realspace();
  kspace();
  selfcontribution();
  writeout();
  printf("Total time is %ds\n",(clock()-tt)*1.0/CLOCKS_PER_SEC);
  return 0;
}
