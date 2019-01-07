#pragma once
#include "classgroup.h"

class ClassgroupSquarer {

private:
  bn_t a;
  bn_t b;
  bn_t c;
  bn_t d;
  
  bn_t e;
  bn_t g;
  bn_t q;
  bn_t r;
  bn_t s;
  bn_t m;

  bn_t C[2]; // for reduced

  bn_t tu;
  bn_t sc;
  bn_t cb;
  bn_t ab;
  bn_t ra;
  bn_t ma;

  bn_t a2;
  bn_t ra2;
  bn_t bar;
  bn_t arrbr;

  // for solve_linear_congruence
  bn_t gg;
  bn_t dd;
  bn_t rr;

  // normalized/reduced
  bn_t aaa;
  bn_t bbb;
  bn_t ccc;
  bn_t csmb;
  bn_t cssmbs;

  // alt_reduce
  bn_t bsqr;
  bn_t bbmd;
  bn_t div_tmp;

  bn_t fba;
  bn_t fbb;
  bn_t fbc;
  bn_t fca;
  bn_t fcb;
  bn_t fcc;
  
  void init();
  
public:  
  ClassgroupSquarer(const Classgroup& g);
  ClassgroupSquarer(const bn_t& aa, const bn_t& bb, const bn_t& cc);
  void square();
  void normalized();
  void normalize();
  void reduced();
  bool test_reduction();
  void fast_reduce();
  Classgroup result();
};
