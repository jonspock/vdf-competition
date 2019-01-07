#include "classgroup.h"
#include "timer.h"
#include <iostream>

// bn_mul_basic is too slow
// Speed Karat 4 is optimum

int main(int argc, char* argv[]) {
  const int num_iters = 1000000;
  timer time_calc;
  relic_int D;
  D = "-117486081361718979337668786659568125756677310916264821465703349385729177030045693203384244824767784173284732777"
      "51253005897383252417809294541920707334384671";

  relic_int G = D;

  const int EL = 16;

  bn_t a[EL];
  bn_t b[EL];
  bn_t m[EL];

  bn_t sum;
  bn_new(sum);

  for (int i = 0; i < EL; i++) {
    bn_new(a[i]);
    bn_new(b[i]);
    bn_new(m[i]);
    G.get_val(b[i]);
    bn_set_dig(a[i],i);
  }

  bn_zero(sum);
  time_calc.start();
  for (int j = 0; j < num_iters; j++) {
    for (int i = 0; i < EL; i++) bn_mul_comba(m[i], a[i], b[i]);
    for (int i = 0; i < EL; i++) bn_add(sum, sum, m[i]);
  }
  time_calc.stop();
  std::cout << "Time Elapsed for Comba: " << time_calc.duration_in_milliseconds(num_iters) << " ms per iteration\t"
            << time_calc.duration_in_milliseconds() << " ms total \n";
  std::cout << "sum = ";
  bn_print(sum);


  bn_zero(sum);
  time_calc.start();
  for (int j = 0; j < num_iters; j++) {
    for (int i = 0; i < EL; i++) bn_mul_simp(m[i], a[i], b[i]);
    for (int i = 0; i < EL; i++) bn_add(sum, sum, m[i]);
  }
  time_calc.stop();
  std::cout << "Time Elapsed for Simp: " << time_calc.duration_in_milliseconds(num_iters) << " ms per iteration\t"
            << time_calc.duration_in_milliseconds() << " ms total \n";
  std::cout << "sum = ";
  bn_print(sum);

  /*
time_calc.start();
for (int i = 0; i < num_iters; i++) bn_mul_karat_imp(m[i], g[i], d[i], 4);
time_calc.stop();
std::cout << "Time Elapsed for Karat 4: " << time_calc.duration_in_milliseconds(num_iters) << " ms per iteration\t" <<
time_calc.duration_in_milliseconds() << " ms total \n";
  */
  return 0;
}
