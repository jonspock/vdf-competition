#include "classgroupSquarer.h"

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
  bn_new(u);
  bn_new(k);
  bn_new(m);
  bn_new(l);
  bn_new(w);

  bn_new(k_temp);
  bn_new(constant_factor);
  bn_new(tu);
  bn_new(sc);
  bn_new(tcf);
  bn_new(tm);

  bn_new(zero);

  bn_new(qs1);
  bn_new(qt1);

  bn_new(qd);

  bn_new(tuk);
  bn_new(nom);
  bn_new(tmc);

  bn_new(ju);
  bn_new(kt);
  bn_new(ls);
  bn_new(kpl);
  bn_new(kl);
  bn_new(jm);

  bn_new(ktls);

  bn_new(gg);
  bn_new(dd);
  bn_new(ee);

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
  bn_new(sc2);
  bn_new(sb);
  bn_new(cs);
  bn_new(csmb);
  bn_new(cssmbs);

  bn_new(xa);
  bn_new(yb);

  bn_new(xx);
  bn_new(yy);
  bn_new(qq);
  bn_new(rr);

  // Initial values
  bn_zero(zero);

  bn_new(R[0]);
  bn_new(R[1]);
  bn_new(S[0]);
  bn_new(S[1]);
  bn_new(T[0]);
  bn_new(T[1]);
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

void ClassgroupSquarer::old_square() {
  bn_gcd(w, a, b);
  my_bn_div(s, a, w);  // t = s
  my_bn_div(u, b, w);

  // Solve mod 1 inputs
  bn_mul(tu, s, u);
  bn_mul(sc, s, c);
  bn_sqr(a, s);  // New 'a'

  // Solve_mod 1 - common factor here s*u & s*s ?
  ext_gcd(tu, a, g, d, e);

  my_bn_div(q, sc, g);
  my_bn_div(constant_factor, a, g);

  bn_mul(qd, q, d);
  bn_mod(k, qd, a);

  // Solve mod 2 - depends on above - inputs
  bn_mul(tcf, s, constant_factor);

  // Solve mod 2 - shortcut for extended_gcd
  //  m = (tu * k - hu - sc) / (st);
  bn_mul(tuk, tu, k);
  bn_sub(nom, tuk, sc);
  my_bn_div(m, nom, a);

  //  l = (tm + c) / u;
  bn_mul(tm, s, m);
  bn_add(tmc, tm, c);
  my_bn_div(l, tmc, u);

  // For b
  bn_mul(ju, w, u);
  bn_add(kpl, k, l);
  bn_mul(ktls, kpl, s);
  bn_sub(b, ju, ktls);

  // For c
  bn_mul(kl, k, l);
  bn_mul(jm, w, m);
  bn_sub(c, kl, jm);

  reduced();
}

void ClassgroupSquarer::solve_linear_congruence(bn_t& mu) {
  ext_gcd(b, a, gg, dd, ee);  // ee ignored (or a,b?)
  my_bn_div(q, c, gg);
  bn_mul(mu, q, dd);
  bn_mod(mu, mu, a);
}

void ClassgroupSquarer::square() {
  //  bn_zero(aaa);
  //  bn_zero(bbb);
  bn_copy(aaa, a);
  bn_copy(bbb, b);
  bn_gcdext(bbb, aaa, g, d);

  my_bn_div(q, c, g);
  bn_mul(tu, q, d);
  bn_mod(tu, tu, a);

  bn_mul(m, b, tu);
  bn_sub(m, m, c);
  my_bn_div(m, m, a);

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
  reduced();
}

void ClassgroupSquarer::reduced() {
  normalized();

  bn_copy(C[0], c);
  bn_copy(C[1], a);
  int IC = 0;
  int IA = 1;

  int CMP = bn_cmp(C[IA], C[IC]);

  while ((CMP == CMP_GT) || (CMP == CMP_EQ && bn_sign(b) == BN_NEG)) {
    // while (a > c or (a == c and b < 0)) {
    bn_add(cb, C[IC], b);
    bn_dbl(ccc, C[IC]);  // Seems faster than lsh with 1 or bn_add
    //  s = (c + b) / (c + c);  //
    my_bn_div2(s, rr, cb, ccc);  // cb, ccc modified but ok since not used after
    bn_mul(sc, s, C[IC]);
    bn_sub(csmb, sc, b);
    bn_mul(cssmbs, csmb, s);
    // New b = 2*s*c-b; == ((cs-b) + sc)
    bn_add(b, csmb, sc);
    // save to C's other element
    bn_add(C[IA], cssmbs, C[IA]);
    // c = c * s * s - b * s + a;  = (c*s - b)*s + a;
    // bn_copy(a,C[IC]); // Swap old C to a
    std::swap(IA, IC);
    CMP = bn_cmp(C[IA], C[IC]);
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
  my_bn_div(r, ab, a2);
  // my_bn_div2(r, rr, ab, a2); // ab,a2 can be modified since not used again

  bn_mul(ra, r, a);
  bn_dbl(ra2, ra);

  bn_add(bar, b, ra);
  bn_mul(arrbr, bar, r);
  bn_add(c, arrbr, c);

  bn_add(b, b, ra2);
}

void ClassgroupSquarer::ext_gcd(const bn_t& aa, const bn_t& bb, bn_t& r0, bn_t& s0, bn_t& t0) {
  // Return r, s, t such that gcd(a, b) = r = a * s + b * t
  bn_copy(R[0], aa);
  bn_copy(R[1], bb);
  bn_set_dig(S[0], 1);
  bn_zero(S[1]);
  bn_zero(T[0]);
  bn_set_dig(T[1], 1);

  if (bn_cmp(R[0], R[1]) == CMP_GT) {
    bn_copy(R[0], bb);
    bn_copy(R[1], aa);
    bn_zero(S[0]);
    bn_zero(T[1]);
    bn_set_dig(S[1], 1);
    bn_set_dig(T[0], 1);
  }

  // typically about 140?
  int IA = 0;
  int IB = 1;
  while (bn_cmp(R[IB], zero) == CMP_GT) {
    // Get Quotient and Remainder in one shot
    my_bn_div_rem_temp_a(q, r, R[IA], R[IB]);  // Allow R[IA] to be overriden/used since it will be replaced

    bn_mul(qt1, q, T[IB]);
    bn_mul(qs1, q, S[IB]);

    bn_sub(T[IA], T[IA], qt1);
    bn_sub(S[IA], S[IA], qs1);

    bn_copy(R[IA], r);

    std::swap(IA, IB);
  }

  if (IA) {
    bn_copy(r0, R[1]);
    bn_copy(s0, S[1]);
    bn_copy(t0, T[1]);
  } else {
    bn_copy(r0, R[0]);
    bn_copy(s0, S[0]);
    bn_copy(t0, T[0]);
  }
}

Classgroup ClassgroupSquarer::result() {
  relic_int aa(a);
  relic_int bb(b);
  relic_int cc(c);
  relic_int dd(d);
  Classgroup cg(aa, bb, cc, dd);
  return cg;
}
