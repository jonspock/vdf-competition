#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "relic_int.h"
#include "type_info.h"

using namespace std;

// Quick way to enable/disable asserts
//#define ASSERT(t) assert(t)
#define ASSERT(t)

int main(int argc, char* argv[]) {
  int status;
  int num = 10;
  double rinp;
  relic_int x;
  relic_int y;
  relic_int diff;
  relic_int ls;
  relic_int rs;

  relic_int a, b;
  relic_int c;

  relic_int inc = 1;
  cout << "incremental value = " << inc << "\n";

  // relic_int neg = -1;
  // cout << "negative value = " << neg << "\n";

  x = y = diff = ls = rs = 0;

  for (int j = 0; j < num; j++) {
    // Create uniform random nunber
    rinp = rand() & 0xfff;  // for now
    y = x;
    x = rinp;

    diff = x - y;
    ls = diff << 2;
    rs = ls;
    rs >>= 2;
  }

  /// integer math
  relic_int t_one = 1;
  relic_int one = 1;
  std::cout << "************* VALUES ****************\n";

  //
  a = 10;
  b = 1;
  c = 2;

  std::cout << "a = " << a << " b = " << b << " c = " << c << "\n";
  std::cout << "************* Real values ****************\n";

  MIXED_VALUES(a, b, 1);
  THIS_MIXED_VALUES(a, b, 1);

  DMIXED_VALUES(a, c, 1);
  DMIXED_THIS_VALUES(a, c, 1);

  FUNC_VALUES(a, 1, 8);

  relic_int sat_a;

  c = 64;
  sat_a = 1;
  b = sat_a;  // save
  VALUE_NAME_SHIFT(sat_a, 1, >>=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 8, >>=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 1, <<=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 8, <<=);

  sat_a = 1;
  b = sat_a;  // save
  VALUE_NAME_SHIFT(sat_a, 1, >>=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 8, >>=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 1, <<=);
  sat_a = b;  // reset
  VALUE_NAME_SHIFT(sat_a, 8, <<=);

  VALUE_NAME_OP(a, b, /);
  VALUE_NAME_OP(a, c, /);
  // VALUE_NAME_OP(a, -c, /);

  // OTHER_FUNCS(a, b)

  VALUE_NAME_FUNC2(a, b, pow);
  // VALUE_NAME_FUNC2(a, b, inverse);
  // Seg faults VALUE_NAME_FUNC(a, isPrime);

  return (0);
}
