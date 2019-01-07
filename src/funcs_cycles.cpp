#include "funcs.h"

#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

static inline uint64_t get_cycles()
{
  uint64_t t;
           // editor's note: "=A" is unsafe for this in x86-64
  __asm volatile ("rdtsc" : "=A"(t));
  return t;
}

void func_rr(int64_t u, int64_t w, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& rr) {
  uint64_t c1,c2;
  int64_t aa,ab,ac;
  c1 = get_cycles();
  aa = u * u; ab = u * w; ac = w * w;
  bn_mul_dig(aaa, a, aa); // since sqr
  bn_mul_sdig(bbb, b, ab);
  bn_mul_dig(ccc, c, ac); // since sqr
  bn_add(rr, aaa, bbb);
  c2 = get_cycles();
  std::cout << "Cycles used in func_rr = " << (int)(c2-c1) << "\n";
}

void func_gg(int64_t u, int64_t v, int64_t x, int64_t w, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& gg) {
  uint64_t c1,c2;
  int64_t ba,bb,bc;
  c1 = get_cycles();
  ba = u * v << 1; bb = u * x + v * w; bc = w * x << 1;
  bn_mul_sdig(aaa, a, ba);
  bn_mul_sdig(bbb, b, bb);
  bn_mul_sdig(ccc, c, bc);
  bn_add(gg, aaa, bbb);
  c2 = get_cycles();
  std::cout << "Cycles used in func_gg = " << (int)(c2-c1) << "\n";
}

void func_dd(int64_t v, int64_t x, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& dd) {
  uint64_t c1,c2;
  int64_t ca,cb,cc;
  c1 = get_cycles();
  ca = v * v; cb = v * x; cc = x * x;
  bn_mul_dig(aaa, a, ca); // since sqr
  bn_mul_sdig(bbb, b, cb);
  bn_mul_dig(ccc, c, cc); // since sqr
  bn_add(dd, aaa, bbb);
  c2 = get_cycles();
  std::cout << "Cycles used in func_dd = " << (int)(c2-c1) << "\n";
}      
  
