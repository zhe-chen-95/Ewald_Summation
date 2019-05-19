#include <cmath>
#include <vector>
#include <ctime>
#include "ewald.h"

using namespace std;

const double Pi = 3.1415927;

// Main implementation
int main(int argc, char* argv[])
{
  long tt = clock();
  initialize();
  realspace();
  kspace();
  selfcontribution();
  writeout();
  printf("Total time is %ds\n",(clock()-tt)*1.0/CLOCK_PER_SEC);
  return 0;
}
