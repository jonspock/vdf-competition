#include "relic_int.h"
#include <iostream>
int main(int argc, char* argv[]) {
  const int loops = 8;

  relic_int a, b;

  // 3rd loop breaks
  //  a.setHex("6C4E77EE1AABB301D6055808277A0B6D2DF4D1B357F7C9C6AB80BA3AB6D5C1A8");
  //  b.setHex("-2CC92A915DE86CA6D6FC3CA6EC1E3042D728806E11F784A6DE96A7C0FD4C7F11");

  // 3rd
  a.setHex("6C4E77EE1AABB301");
  b.setHex("-2CC92A915DE86CA");

  for (int i = 0; i < loops; i++) {
    // b = a + b; --> causes issues
    a = a + b;

    relic_int r, s, t;

    std::cout << "------------a,b (hex....) = " << a << " ------ " << b << "-----------------\n";
    std::cout << "------------a,b (decimal) = " << a.getDec() << " ------ " << b.getDec() << "-----------------\n";
    std::tie(r, s, t) = extended_gcd(a, b);  // will ignore t
    std::cout << "r,s = " << r << "," << s << "\n";

    bn_t aa, bb;
    bn_new(aa);
    bn_new(bb);

    a.get_val(aa);
    b.get_val(bb);

    bn_t rr, ss;
    bn_new(rr);
    bn_new(ss);
    /*
    std::cout << "a,b = "; bn_print(aa);
    std::cout << " ";     bn_print(bb);
    std::cout << "\n";
    */
    // std::cout << "Top a,b Sizes = " << aa->used << "," << bb->used << "\n";

    bn_gcdext(aa, bb, rr, ss);

    std::cout << "r,s = ";
    bn_print(rr);
    std::cout << " ";
    bn_print(ss);
    std::cout << "\n";
  }
  return 0;
}
