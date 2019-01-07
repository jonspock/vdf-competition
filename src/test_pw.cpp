#include "proof_wesolowski.h"
#include "timer.h"
#include <iostream>

int main(int argc, char* argv[]) {
  const int num_iters = 10000;
  timer time_calc;

  std::cout << "benchmark small 512 bit prime\n";
  // relic_int oldD = -77447;
  relic_int D(
      "-117486081361718979337668786659568125756677310916264821465703349385729177030045693203384244824767784173284732777"
      "51253005897383252417809294541920707334384671");

  std::cout << "D = " << D << " " << D.bitSize() << "\n";

  auto [L, k, tmp] = approximate_parameters(num_iters);

  std::cout << "L,k,tmp= " << L << "," << k << "," << tmp << "\n";

  int powers = ceil(num_iters / double(k * L)) + 1 + num_iters;

  relic_int one(1);
  relic_int two(2);
  Classgroup g = from_ab_discriminant(two, one, D);
  std::cout << "Starting g = " << g << " {" << g.a_bits() << "," << g.b_bits() << "," << g.c_bits() << "}\n";

  std::vector<Classgroup> C;
  Classgroup x = g;
  for (int i = 0; i < powers; i++) {
    x = pow(x, 2);  //?
    C.push_back(x);
  }
  Classgroup y = C[powers - 1];
  bool send_conn = false;

  auto identity = identity_for_discriminant(D);

  relic_int Bhash = hash_prime(g.Serialize(), y.Serialize());

  time_calc.start();
  auto proof = eval_optimized(identity, g, Bhash, num_iters, k, L, &C[0]);
  time_calc.stop();
  std::cout << "Time Elapsed for eval_optimized: " << time_calc.duration_in_milliseconds() << " ms total\n";

  time_calc.start();
  bool ok = verify_proof(Bhash, g, y, proof, num_iters);
  time_calc.stop();

  std::cout << "Time Elapsed for verify: " << time_calc.duration_in_milliseconds() << " ms total\n";
  std::cout << " verify = " << ok << "\n";

  return 0;
}
