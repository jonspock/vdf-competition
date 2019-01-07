#include "classgroup.h"
#include "timer.h"
#include <iostream>
#include <omp.h>

#define EL 2

int main(int argc, char* argv[]) {
  const int num_iters = 100000;
  timer time_calc;
  relic_int D;
  D = "-117486081361718979337668786659568125756677310916264821465703349385729177030045693203384244824767784173284732777"
      "51253005897383252417809294541920707334384671";

  relic_int G = D;

  G = G*D;
  G = G*G;
  
  std::cout << "G size = " <<  G.bitSize() << "\n";

  bn_t g[EL];
  bn_t d[EL];
  bn_t m[EL];

  bn_t sum;
    
  for (int i = 0; i < EL; i++) {
    bn_new(g[i]);
    bn_new(d[i]);
    bn_new(m[i]);
    G.get_val(g[i]);
    G.get_val(d[i]);
  }

    
  bn_new(sum);
  bn_zero(sum);

  time_calc.start();
  for (int j = 0; j < num_iters; j++) {
    for (int i = 0; i < EL; i++) bn_mul_comba(m[i],g[i],d[i]);
    for (int i = 0; i < EL; i++) bn_add(sum, sum, m[i]);
  }
  time_calc.stop();
  std::cout << "Time Elapsed for Comba: " << time_calc.duration_in_milliseconds(num_iters) << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";
  std::cout << "sum = ";
  bn_print(sum);

  bn_zero(sum);
  time_calc.start();
  for (int j = 0; j < num_iters; j++) {
#pragma omp for
    for (int i = 0; i < EL; i++) bn_mul_comba(m[i],g[i],d[i]);
    for (int i = 0; i < EL; i++) bn_add(sum, sum, m[i]);
  }
  time_calc.stop();
  std::cout << "Time Elapsed for Comba: " << time_calc.duration_in_milliseconds(num_iters) << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";
  std::cout << "sum = ";
  bn_print(sum);
    
  return 0;
}
