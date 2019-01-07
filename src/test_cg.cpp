#include "classgroup.h"
#include <iostream>
int main(int argc, char* argv[]) {
  std::cout << "compile check classgroup\n";
  Classgroup c(2, 1, 3);
  auto d = c.discriminant();
  std::cout << "Disc = " << d << "\n";  // Should be -23 add check later

  // Inverse
  relic_int D = -7 - 8;
  relic_int one(1);
  relic_int two(2);
  Classgroup e_id = identity_for_discriminant(D);
  Classgroup e_gen = from_ab_discriminant(two, one, D);
  std::cout << "e_gen = " << e_gen << "\n";
  Classgroup e_gen_inv = inverse(e_gen);
  std::cout << "e_gen_inv = " << e_gen_inv << "\n";
  Classgroup e_prod = e_gen * e_gen_inv;
  std::cout << "e_prod = " << e_prod << "\n";  // Should be 1,1,4

  // Seems OK
  for (int i = 0; i < 10; i++) {
    relic_int D = -7 - 8 * i;
    Classgroup e0 = from_ab_discriminant(two, one, D);
    Classgroup e1 = inverse(e0);

    int power = 9;
    Classgroup p0 = pow(e0, power);
    Classgroup p1 = pow(e1, power);

    std::cout << "e0 = " << e0 << " e0^" << power << " = " << p0 << "\n";
    std::cout << "e1 = " << e1 << " e1^" << power << " = " << p1 << "\n";
  }

  Classgroup c2(11, 49, 55);
  Classgroup n = normalized(c2);
  std::cout << "Normalized = " << n << "\n";  // 11,5,1
  Classgroup r2 = c2.reduced();
  std::cout << "Reduced = " << r2 << "\n";  // 1,1,5

  Classgroup t2_1_6(2, 1, 6);
  Classgroup t4_1_3(4, 1, 3);
  auto p = t2_1_6 * t4_1_3;
  std::cout << "mult = " << p << "\n";  // should be 3,1,4

  Classgroup t12_11_3(12, 11, 3);
  Classgroup t93_109_32(93, 109, 32);

  std::cout << "disc = " << t12_11_3.discriminant() << " vs  -23\n";
  std::cout << "disc = " << t93_109_32.discriminant() << " vs  -23\n";

  Classgroup t = t12_11_3 * t93_109_32;
  std::cout << "t = " << t << "\n";  // (1, -15, 62))
  std::cout << "t disc = " << t.discriminant() << " vs  -23\n";

  Classgroup t1 = t93_109_32 * t12_11_3;
  std::cout << "t1 = " << t1 << "\n";  // (1, -15, 62))

  // the normalize and reduce example from the paper
  Classgroup f(195751, 1212121, 1876411);
  std::cout << "f.normalized() = " << normalized(f) << "\n";  // (195751, 37615, 1807))
  std::cout << "f.reduced() = " << f.reduced() << "\n";       // (1, 1, 1))

  return 0;
}
