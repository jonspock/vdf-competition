#include "funcs.h"

void func_rr(int64_t u, int64_t w, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& rr) {
  int64_t aa,ab,ac;
  aa = u * u; ab = u * w; ac = w * w;
  bn_mul_dig(aaa, a, aa); // since sqr
  bn_mul_sdig(bbb, b, ab);
  bn_mul_dig(ccc, c, ac); // since sqr
  bn_add(rr, aaa, bbb);
}

void func_gg(int64_t u, int64_t v, int64_t x, int64_t w, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& gg) {
  int64_t ba,bb,bc;
  ba = u * v << 1; bb = u * x + v * w; bc = w * x << 1;
  bn_mul_sdig(aaa, a, ba);
  bn_mul_sdig(bbb, b, bb);
  bn_mul_sdig(ccc, c, bc);
  bn_add(gg, aaa, bbb);
}

void func_dd(int64_t v, int64_t x, const bn_t& a,  const bn_t& b, const bn_t& c,
             bn_t& aaa, bn_t& bbb, bn_t& ccc,
             bn_t& dd) {
  int64_t ca,cb,cc;
  ca = v * v; cb = v * x; cc = x * x;
  bn_mul_dig(aaa, a, ca); // since sqr
  bn_mul_sdig(bbb, b, cb);
  bn_mul_dig(ccc, c, cc); // since sqr
  bn_add(dd, aaa, bbb);
}      
  
