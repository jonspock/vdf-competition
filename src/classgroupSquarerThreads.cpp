#include "classgroupSquarer.h"
#include <thread>

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


void ClassgroupSquarer::init() {
  bn_new(a);
  bn_new(b);
  bn_new(c);
  bn_new(d);
  bn_new(e);
  bn_new(g);

  bn_new(q);
  bn_new(r);
  bn_new(s);
  bn_new(m);

  bn_new(tu);
  bn_new(sc);

  bn_new(gg);
  bn_new(dd);
  bn_new(rr);

  // in normalized
  bn_new(ab);
  bn_new(a2);
  bn_new(ra);
  bn_new(ra2);
  bn_new(bar);
  bn_new(arrbr);
  bn_new(ma);

  bn_new(cb);
  bn_new(aaa);
  bn_new(bbb);
  bn_new(ccc);
  bn_new(csmb);
  bn_new(cssmbs);

  bn_new(bbmd);
  bn_new(bsqr);
  bn_new(div_tmp);

  bn_new(C[0]);
  bn_new(C[1]);

  bn_new(fba);
  bn_new(fbb);
  bn_new(fbc);

  bn_new(fca);
  bn_new(fcb);
  bn_new(fcc);

}

ClassgroupSquarer::ClassgroupSquarer(const Classgroup& x) {
  init();
  x.get_a(a);
  x.get_b(b);
  x.get_c(c);
  x.get_d(d);
}

ClassgroupSquarer::ClassgroupSquarer(const bn_t& aa, const bn_t& bb, const bn_t& cc) {
  init();
  // Initial values
  bn_copy(a, aa);
  bn_copy(b, bb);
  bn_copy(c, cc);
}

void ClassgroupSquarer::square() {
  //  bn_zero(aaa);
  //  bn_zero(bbb);
  bn_copy(aaa, a);
  bn_copy(bbb, b);
  bn_gcdext(bbb, aaa, g, s);

  bn_div(q, c, g);
  bn_mul(tu, q, s);
  bn_mod(tu, tu, a);

  bn_mul(m, b, tu);
  bn_sub(m, m, c);
  bn_div(m, m, a);

  bn_mul(sc, tu, a);  // a*tu
  bn_dbl(a2, sc);     // 2*a*tu

  // a -> a*a
  bn_copy(aaa, a);
  bn_sqr(a, aaa);  // a = a^2
  //  B = B - 2*a*tu
  bn_sub(b, b, a2);
  // c
  bn_sqr(c, tu);
  bn_sub(c, c, m);
#ifdef OLD_METHOD
  reduced();
#else
  fast_reduce();
#endif
}

void ClassgroupSquarer::reduced() {
  normalized();

  bn_copy(C[0], c);
  bn_copy(C[1], a);
  int IC = 0;
  int IA = 1;

  int CMP = bn_cmp(C[IA], C[IC]);

  // Roughly  reduce  ~74% 
  // bn_div_fast 26%
  // bn_mul      18%
  // bn_add      19%
  // bn_sub       6%
  // bn_dbl       4%
  
  // while (a > c or (a == c and b < 0)) {
  while ((CMP == CMP_GT) || (CMP == CMP_EQ && bn_sign(b) == BN_NEG)) {
    bn_add(cb, C[IC], b);
    bn_dbl(ccc, C[IC]);  // Seems faster than lsh with 1 or bn_add
    //  s = (c + b) / (c + c);  //
    // Sizes of cb and ccc are quite close
    bn_div_fast(s, rr, cb, ccc);  // cb, ccc modified but ok since not used after

    // Combining mul + add/sub may not help since both outputs are needed
    // Usually s is small and C[IC] large
    bn_mul(sc, s, C[IC]);
    bn_sub(csmb, sc, b);
    // Usually s is small and csmb large
    bn_mul(cssmbs, s, csmb);
    // New b = 2*s*c-b; == ((c*s-b) + sc)
    bn_add(b, csmb, sc);
    // save c (as IA because will swap)
    // c = c * s * s - b * s + a;  = (c*s - b)*s + a;
    bn_add(C[IA], cssmbs, C[IA]);
    // did have swap first, but moved to lessen stall dependency (swap IC/IA instead)
    CMP = bn_cmp(C[IC], C[IA]);
    std::swap(IA, IC); // not needed for a few cycles 
  }
  // copy final result to a,c
  bn_copy(c, C[IC]);
  bn_copy(a, C[IA]);
  normalized();
}

void ClassgroupSquarer::normalized() {
  bn_neg(ma, a);

  if ((bn_cmp(ma, b) == CMP_LT) && ((bn_cmp(b, a) == CMP_LT) || (bn_cmp(b, a) == CMP_EQ))) return;
  bn_sub(ab, a, b);
  bn_dbl(a2, a);
  bn_div(r, ab, a2);
  // my_bn_div2(r, rr, ab, a2); // ab,a2 can be modified since not used again

  bn_mul(ra, r, a);
  bn_dbl(ra2, ra);

  bn_add(bar, b, ra);
  bn_mul(arrbr, bar, r);
  bn_add(c, arrbr, c);

  bn_add(b, b, ra2);
}

void ClassgroupSquarer::normalize() {

  bn_add(ma, b, c);
  bn_dbl(a2, c);
   
  bn_div(m, ma, a2); //???

  bn_copy(aaa, c);
  bn_copy(ccc, a);
  
  bn_dbl(a2, m);
  bn_neg(bbb, b);
  bn_mul(gg, c, a2);
  bn_add(bbb, bbb, gg); // ok
  
  bn_mul(gg, b, m);
  bn_sub(ccc, ccc, gg);
  
  bn_sqr(ma, m); //denom*denom
  
  bn_mul(gg, c, ma);
  bn_add(ccc, ccc, gg);
  
  bn_copy(a, aaa);
  bn_copy(b, bbb);
  bn_copy(c, ccc);;
}

Classgroup ClassgroupSquarer::result() {
  relic_int aa(a);
  relic_int bb(b);
  relic_int cc(c);
  relic_int dd(d);
  Classgroup cg(aa, bb, cc, dd);
  return cg;
}

/**
Copyright 2018 Chia Network Inc
Modifications copyright (C) 2019 Akashnil Dutta

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
**/


#define LOG2(X) (63 - __builtin_clzll((X)))

const int64_t THRESH = 1UL<<31;
const int64_t EXP_THRESH = 31;
//const int64_t DOUBLE_THRESH = 1UL<<31;
//const int64_t MAXV = ((1UL<<63) - 1);

inline uint64_t signed_shift(uint64_t op, int shift) {
    if (shift > 0) return op << shift;
    if (shift <= -64) return 0;
    return op >> (-shift);
}

// Return an approximation x of the large mpz_t op by an int64_t and the exponent e adjustment.
// We must have (x * 2^e) / op = constant approximately.
inline int64_t mpz_get_si_2exp (signed long int *exp, const bn_t& op) {
  uint64_t size = op->used;
  uint64_t last = op->dp[size - 1];
  uint64_t ret;
  int lg2 = LOG2(last) + 1;
  *exp = lg2; ret = signed_shift(last, 63 - *exp);
  if (size > 1) {
    *exp += (size-1) * 64;
    uint64_t prev = op->dp[size - 2];
    ret += signed_shift(prev, -1 - lg2);
  }
  if (op->sign) return - ((int64_t)ret);
  return ret;
}

// Test if f is reduced. If it almost is but a, c are swapped, 
// then just swap them to make it reduced.
bool ClassgroupSquarer::test_reduction() {
  
  int a_b = bn_cmp_abs(a, b); // check
  int c_b = bn_cmp_abs(c, b);

  if (a_b == CMP_LT || c_b == CMP_LT) return false;

  int a_c = bn_cmp(a, c);

  if (a_c == CMP_GT) {
    // do better swap
    bn_copy(bbb, c);
    bn_copy(c, a);
    bn_copy(a, bbb);
    bn_neg(b,b);
  }

  if (a_c == CMP_EQ && b->sign) {
    bn_neg(b,b);
  }
    
  return true;
}

// This is a replacement of the original reduce procedure.
void ClassgroupSquarer::fast_reduce() {

    int64_t u, v, w, x, u_, v_, w_, x_;
    int64_t delta, gamma;
    int64_t a_, b_, c_;
    int64_t iaa, ibb, icc;
    int64_t aa, ab, ac, ba, bb, bc, ca, cb, cc;
    signed long int a_exp, b_exp, c_exp, max_exp, min_exp;
    
//    while ((CMP == CMP_GT) || (CMP == CMP_EQ && bn_sign(b) == BN_NEG)) {
    while (!test_reduction()) {
      iaa = mpz_get_si_2exp(&a_exp, a);
      ibb = mpz_get_si_2exp(&b_exp, b);
      icc = mpz_get_si_2exp(&c_exp, c);

      max_exp = a_exp;
      min_exp = a_exp;
      
      if (max_exp < b_exp) max_exp = b_exp;
      if (min_exp > b_exp) min_exp = b_exp;
      
      if (max_exp < c_exp) max_exp = c_exp;
      if (min_exp > c_exp) min_exp = c_exp;
      
      if (max_exp - min_exp > EXP_THRESH) {
        normalize();
      } else {
        max_exp++; // for safety vs overflow

        // Ensure a, b, c are shifted so that a : b : c ratios are same as f.a : f.b : f.c
        // a, b, c will be used as approximations to f.a, f.b, f.c
        iaa >>= (max_exp - a_exp);
        ibb >>= (max_exp - b_exp);
        icc >>= (max_exp - c_exp);
        
        u_ = 1; v_ = 0; w_ = 0; x_ = 1;
        
        // We must be very careful about overflow in the following steps
        do {
          u = u_; v = v_; w = w_; x = x_;
          // Ensure that delta = floor ((b+c) / 2c)
          delta = ibb >= 0 ? (ibb+icc) / (icc<<1) : - (-ibb+icc) / (icc<<1);
          a_ = icc;
          c_ = icc * delta;
          b_ = -ibb + (c_ << 1);
          gamma = ibb - c_;
          c_ = iaa - delta * gamma;
          
          iaa = a_; ibb = b_; icc = c_;
          
          u_ = v;
          v_ = -u + delta * v;
          w_ = x;
          x_ = -w + delta * x;
          // The condition (abs(v_) | abs(x_)) <= THRESH protects against overflow
        } while ((abs(v_) | abs(x_)) <= THRESH && iaa > icc && icc > 0);
        
        if ((abs(v_) | abs(x_)) <= THRESH) {
          u = u_; v = v_; w = w_; x = x_;
        }

        // The following operations take about 40% of the overall runtime.
        std::thread th1(func_rr, u,w, std::ref(a),std::ref(b),std::ref(c), std::ref(aaa), std::ref(bbb), std::ref(ccc), std::ref(rr));
        std::thread th2(func_gg, u,v,x,w, std::ref(a), std::ref(b), std::ref(c), std::ref(fba), std::ref(fbb), std::ref(fbc), std::ref(gg));
        std::thread th3(func_dd, v,x, std::ref(a), std::ref(b), std::ref(c), std::ref(fca), std::ref(fcb), std::ref(fcc), std::ref(dd));

        th1.join();
        th2.join();
        th3.join();
        
        bn_add(a, rr, ccc);
        bn_add(b, gg, fbc);
        bn_add(c, dd, fcc);
      }
    }
}

