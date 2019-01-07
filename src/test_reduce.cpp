#include "classgroupSquarer.h"
#include "relic_int.h"
#include <iostream>
int main(int argc, char* argv[]) {
  relic_int a, b, c, d;

  // 3rd
  a = 76;
  b = 217;
  c = 155;
  d = -31;

  Classgroup C(a, b, c, d);
  ClassgroupSquarer D1(C);
  Classgroup R1 = C;

  D1.reduced();

  int count = 0;
  std::cout << "Start = " << R1 << "\n";
  std::cout << "R1 = " << D1.result() << "\n";

  ClassgroupSquarer D2(C);
  D2.fast_reduce();
  std::cout << "R2 = " << D2.result() << "\n";

  return 0;
}
