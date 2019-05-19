#include <cmath>
#include <vector>
#include <ctime>
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
  printf("Total time is %ds\n",(clock()-tt)*1.0/CLOCK_PER_SEC);
  return 0;
}
