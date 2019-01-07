#include "classgroup.h"
#include <iostream>
int main(int argc, char* argv[]) {
  bn_t aa, bb, cc;

  int64_t a_int = -1;
  int64_t b_int = 1;

  bn_new(aa);
  bn_new(bb);
  bn_new(cc);

  bn_set_dig(aa, a_int);
  bn_set_dig(bb, b_int);

  bn_div(cc, aa, bb);
  bn_print(aa);
  bn_print(bb);
  bn_print(cc);

  relic_int a(-1);
  relic_int b(1);

  std::cout << "a,b = " << a << "," << b << "\n";

  auto div = a / b;
  std::cout << "for -1/1 div = " << div << "\n";

  std::cout << "compile check classgroup\n";
  Classgroup c(2, 1, 3);
  auto d = c.discriminant();
  std::cout << "Disc = " << d << "\n";  // Should be -23 add check later

  Classgroup c2(11, 49, 55);
  Classgroup n = normalized(c2);
  std::cout << "Normalized = " << n << "\n";  // 11,5,1
  Classgroup r2 = c2.reduced();
  std::cout << "Reduced = " << r2 << "\n";  // 1,1,5

  Classgroup t2_1_6(2, 1, 6);
  Classgroup t4_1_3(4, 1, 3);
  auto p = t2_1_6 * t4_1_3;
  std::cout << "mult = " << p << "\n";

  return 0;
}
