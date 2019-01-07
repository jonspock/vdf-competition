#include "classgroup.h"
#include "classgroupSquarer.h"
#include "timer.h"
#include <iostream>

int main(int argc, char* argv[]) {
  relic_int D;

  D = "-117486081361718979337668786659568125756677310916264821465703349385729177030045693203384244824767784173284"
      "73277751253005897383252417809294541920707334384671";

  std::cout << "D = " << D << " " << D.bitSize() << "\n";

  relic_int one(1);
  relic_int two(2);
  Classgroup g = from_ab_discriminant(two, one, D);
  while ((g.a_bits() < g.c_bits()) || (g.b_bits() < g.c_bits())) { g = pow(g, 2); }
  std::cout << "g = " << g << "\n";
  Classgroup g2 = pow(g, 2);

  bn_t a;
  bn_t b;
  bn_t w;
  bn_t s;
  bn_t abs_a;
  bn_t xa, yb, r, q;

  bn_new(a);
  bn_new(w);
  bn_new(s);
  bn_new(abs_a);
  bn_new(xa);
  bn_new(yb);
  bn_new(r);
  bn_new(q);

  g2.get_a(a);
  g2.get_b(b);
  bn_gcd(w, a, b);  // w is +

  std::cout << "Input w,a,b = ";
  bn_print(w);
  bn_print(a);
  bn_print(b);

  my_bn_div(s, a, w);

  std::cout << "Original s = ";
  bn_print(s);

  bn_abs(abs_a, a);
  int original_a_sign = bn_sign(a);

  std::cout << "abs (a) = ";
  bn_print(abs_a);

  bn_vdf_div(s, r, q,  // 4 temp vals
             original_a_sign, abs_a, w);

  std::cout << "New s = ";
  bn_print(s);

  return 0;
}
