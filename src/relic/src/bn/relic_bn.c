/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (C) 2007-2017 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * RELIC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RELIC. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * Implementation of the multiple precision addition and subtraction functions.
 *
 * @ingroup bn
 */

#include "relic_bn_low.h"
#include "relic_core.h"
#include <inttypes.h>
#include <errno.h>

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Adds two multiple precision integers, where a >= b.
 *
 * @param[out] c        - the result.
 * @param[in] a         - the first multiple precision integer to add.
 * @param[in] b         - the second multiple precision integer to add.
 */
static void bn_add_imp(bn_t c, const bn_t a, const bn_t b) {
  int max, min;
  dig_t carry;

  max = a->used;
  min = b->used;

  /* Grow the result. */
  bn_grow(c, max);

  if (a->used == b->used) {
    carry = bn_addn_low(c->dp, a->dp, b->dp, max);
  } else {
    carry = bn_addn_low(c->dp, a->dp, b->dp, min);
    carry = bn_add1_low(c->dp + min, a->dp + min, carry, max - min);
  }
  if (carry) {
    bn_grow(c, max + 1);
    c->dp[max] = carry;
  }
  c->used = max + carry;
  bn_trim(c);
}

/**
 * Subtracts two multiple precision integers, where a >= b.
 *
 * @param[out] c        - the result.
 * @param[in] a         - the first multiple precision integer to subtract.
 * @param[in] b         - the second multiple precision integer to subtract.
 */
static void bn_sub_imp(bn_t c, const bn_t a, const bn_t b) {
  int max, min;
  dig_t carry;

  max = a->used;
  min = b->used;

  /* Grow the destination to accomodate the result. */
  bn_grow(c, max);

  if (a->used == b->used) {
    carry = bn_subn_low(c->dp, a->dp, b->dp, min);
  } else {
    carry = bn_subn_low(c->dp, a->dp, b->dp, min);
    carry = bn_sub1_low(c->dp + min, a->dp + min, carry, max - min);
  }
  c->used = max;
  bn_trim(c);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_add(bn_t c, const bn_t a, const bn_t b) {
  int sa, sb;

  sa = a->sign;
  sb = b->sign;

  if (sa == sb) {
    /* If the signs are equal, copy the sign and add. */
    c->sign = sa;
    if (bn_cmp_abs(a, b) == CMP_LT) {
      bn_add_imp(c, b, a);
    } else {
      bn_add_imp(c, a, b);
    }
  } else {
    /* If the signs are different, subtract. */
    if (bn_cmp_abs(a, b) == CMP_LT) {
      c->sign = sb;
      bn_sub_imp(c, b, a);
    } else {
      c->sign = sa;
      bn_sub_imp(c, a, b);
    }
  }
}

void bn_add_dig(bn_t c, const bn_t a, dig_t b) {
  dig_t carry;

  bn_grow(c, a->used);

  if (a->sign == BN_POS) {
    carry = bn_add1_low(c->dp, a->dp, b, a->used);
    if (carry) {
      bn_grow(c, a->used + 1);
      c->dp[a->used] = carry;
    }
    c->used = a->used + carry;
    c->sign = BN_POS;
  } else {
    /* If a < 0 && |a| >= b, compute c = -(|a| - b). */
    if (a->used > 1 || a->dp[0] >= b) {
      carry = bn_sub1_low(c->dp, a->dp, b, a->used);
      c->used = a->used;
      c->sign = BN_NEG;
    } else {
      /* If a < 0 && |a| < b. */
      if (a->used == 1) {
        c->dp[0] = b - a->dp[0];
      } else {
        c->dp[0] = b;
      }
      c->used = 1;
      c->sign = BN_POS;
    }
  }
  bn_trim(c);
}

void bn_sub(bn_t c, const bn_t a, const bn_t b) {
  int sa, sb;

  sa = a->sign;
  sb = b->sign;

  if (sa != sb) {
    /* If the signs are different, copy the sign of the first number and
     * add. */
    c->sign = sa;
    if (bn_cmp_abs(a, b) == CMP_LT) {
      bn_add_imp(c, b, a);
    } else {
      bn_add_imp(c, a, b);
    }
  } else {
    /* If the signs are equal, adjust the sign and subtract. */
    if (bn_cmp_abs(a, b) != CMP_LT) {
      c->sign = sa;
      bn_sub_imp(c, a, b);
    } else {
      c->sign = (sa == BN_POS) ? BN_NEG : BN_POS;
      bn_sub_imp(c, b, a);
    }
  }
}

void bn_sub_dig(bn_t c, const bn_t a, dig_t b) {
  dig_t carry;

  bn_grow(c, a->used);

  /* If a < 0, compute c = -(|a| + b). */
  if (a->sign == BN_NEG) {
    carry = bn_add1_low(c->dp, a->dp, b, a->used);
    if (carry) {
      bn_grow(c, a->used + 1);
      c->dp[a->used] = carry;
    }
    c->used = a->used + carry;
    c->sign = BN_NEG;
  } else {
    /* If a > 0 && |a| >= b, compute c = (|a| - b). */
    if (a->used > 1 || a->dp[0] >= b) {
      carry = bn_sub1_low(c->dp, a->dp, b, a->used);
      c->used = a->used;
      c->sign = BN_POS;
    } else {
      /* If a > 0 && a < b. */
      if (a->used == 1) {
        c->dp[0] = b - a->dp[0];
      } else {
        c->dp[0] = b;
      }
      c->used = 1;
      c->sign = BN_NEG;
    }
  }
  bn_trim(c);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int bn_cmp_abs(const bn_t a, const bn_t b) {
  if (a->used > b->used) { return CMP_GT; }

  if (a->used < b->used) { return CMP_LT; }

  return bn_cmpn_low(a->dp, b->dp, a->used);
}

int bn_cmp_dig(const bn_t a, dig_t b) {
  if (a->sign == BN_NEG) { return CMP_LT; }

  if (a->used > 1) { return CMP_GT; }

  return bn_cmp1_low(a->dp[0], b);
}

int bn_cmp(const bn_t a, const bn_t b) {
  if (a->sign == BN_POS && b->sign == BN_NEG) { return CMP_GT; }
  if (a->sign == BN_NEG && b->sign == BN_POS) { return CMP_LT; }

  if (a->sign == BN_NEG) { return bn_cmp_abs(b, a); }

  return bn_cmp_abs(a, b);
}

void bn_div_dig(bn_t c, const bn_t a, dig_t b) {
  bn_t q;
  dig_t r;

  if (b == 0) { THROW(ERR_NO_VALID); }

  if (b == 1 || bn_is_zero(a) == 1) {
    if (c != NULL) { bn_copy(c, a); }
    return;
  }

  TRY {
    bn_new(q);
    int size = a->used;
    const dig_t *ap = a->dp;

    bn_div1_low(q->dp, &r, ap, size, b);

    if (c != NULL) {
      q->used = a->used;
      q->sign = a->sign;
      bn_trim(q);
      bn_copy(c, q);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_div_rem_dig(bn_t c, dig_t *d, const bn_t a, dig_t b) {
  bn_t q;
  dig_t r;


  if (b == 0) { THROW(ERR_NO_VALID); }

  if (b == 1 || bn_is_zero(a) == 1) {
    if (d != NULL) { *d = 0; }
    if (c != NULL) { bn_copy(c, a); }
    return;
  }

  TRY {
    bn_new(q);
    int size = a->used;
    const dig_t *ap = a->dp;

    bn_div1_low(q->dp, &r, ap, size, b);

    if (c != NULL) {
      q->used = a->used;
      q->sign = a->sign;
      bn_trim(q);
      bn_copy(c, q);
    }

    if (d != NULL) { *d = r; }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
void bn_div_exact(bn_t c, bn_t rr, bn_t a, bn_t b) {
  int sign;

  TRY {
    bn_new_size(c, a->used + 1);
    bn_zero(c);
    bn_zero(rr);
    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);
    // make a/b + for tdiv_qr
    a->sign = BN_POS;
    b->sign = BN_POS;
    bn_divn_low(c->dp, rr->dp, a->dp, a->used, b->dp, b->used);
    // mpn_divexact(c->dp, a->dp, a->used, b->dp, b->used);
    /* We have the quotient in q -> ignore r */
    c->used = a->used - b->used + 1;
    c->sign = sign;
    bn_trim(c);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_impf(bn_t c, bn_t rr, bn_t a, bn_t b) {
  int sign;

  /* If |a| < |b|, we're done. */
  if (bn_cmp_abs(a, b) == CMP_LT) {
    if (bn_sign(a) == bn_sign(b)) {
      bn_zero(c);
    } else {
      bn_set_neg_dig(c, 1);
    }
    return;
  }

  TRY {
    bn_new_size(c, a->used + 1);
    bn_zero(c);
    bn_zero(rr);

    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);
    // make a/b + for tdiv_qr
    a->sign = BN_POS;
    b->sign = BN_POS;
    bn_divn_low(c->dp, rr->dp, a->dp, a->used, b->dp, b->used);
    // mpn_tdiv_qr(c->dp, rr->dp, 0, a->dp, a->used, b->dp, b->used);

    /* We have the quotient in q and the remainder in r. */
    c->used = a->used - b->used + 1;
    c->sign = sign;
    bn_trim(c);
    if ((sign == BN_POS) || bn_is_zero(rr)) {
      return;
    } else {
      bn_sub_dig(c, c, 1);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[out] d		- the remainder.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_impf_with_temps(bn_t c, bn_t d, bn_t xx, bn_t yy, bn_t rr, bn_t qq, const bn_t a, const bn_t b) {
  int sign;

  /* If |a| < |b|, we're done. */
  if (bn_cmp_abs(a, b) == CMP_LT) {
    if (bn_sign(a) == bn_sign(b)) {
      if (c != NULL) { bn_zero(c); }
      if (d != NULL) { bn_copy(d, a); }
    } else {
      if (c != NULL) {
        bn_set_dig(c, 1);
        bn_neg(c, c);
      }
      if (d != NULL) { bn_add(d, a, b); }
    }
    return;
  }

  TRY {
    bn_new_size(qq, a->used + 1);
    bn_zero(qq);
    bn_zero(rr);
    bn_abs(xx, a);
    bn_abs(yy, b);

    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);

    bn_divn_low(qq->dp, rr->dp, xx->dp, a->used, yy->dp, b->used);
    // mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);

    /* We have the quotient in q and the remainder in r. */
    if (c != NULL) {
      qq->used = a->used - b->used + 1;
      qq->sign = sign;
      bn_trim(qq);
      if ((sign == BN_POS) || bn_is_zero(rr)) {
        bn_copy(c, qq);
      } else {
        bn_sub_dig(c, qq, 1);
      }
    }

    if (d != NULL) {
      rr->used = b->used;
      rr->sign = b->sign;
      bn_trim(rr);
      if ((sign == BN_POS) || bn_is_zero(rr)) {
        bn_copy(d, rr);
      } else {
        bn_sub(d, b, rr);
      }
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}
/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_impf1_const(bn_t c, bn_t xx, bn_t yy, bn_t rr, bn_t qq, const bn_t a, const bn_t b) {
  int sign;

  /* If |a| < |b|, we're done. */
  if (bn_cmp_abs(a, b) == CMP_LT) {
    if (bn_sign(a) == bn_sign(b)) {
      bn_zero(c);
    } else {
      bn_set_neg_dig(c, 1);
    }
    return;
  }

  TRY {
    bn_new_size(qq, a->used + 1);
    bn_zero(qq);
    bn_zero(rr);
    bn_abs(xx, a);
    bn_abs(yy, b);

    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);
    bn_divn_low(qq->dp, rr->dp, xx->dp, a->used, yy->dp, b->used);
    // mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);

    /* We have the quotient in q and the remainder in r. */
    qq->used = a->used - b->used + 1;
    qq->sign = sign;
    bn_trim(qq);
    if ((sign == BN_POS) || bn_is_zero(rr)) {
      bn_copy(c, qq);
    } else {
      bn_sub_dig(c, qq, 1);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_div_with_temps(bn_t c, bn_t xx, bn_t yy, bn_t rr, bn_t qq, const bn_t a, const bn_t b) {
  if (bn_is_zero(b)) { THROW(ERR_NO_VALID); }
  bn_div_impf1_const(c, xx, yy, qq, rr, a, b);
}

void bn_div_fast(bn_t c,
                 // temporary
                 bn_t rr,
                 // avoid copy since not const and can be changed
                 bn_t a, bn_t b) {
  if (bn_is_zero(b)) { THROW(ERR_NO_VALID); }
  bn_div_impf(c, rr, a, b);
}

void bn_div_rem_with_temps(bn_t c, bn_t d,
                           // temporaries
                           bn_t xx, bn_t yy, bn_t rr, bn_t qq, const bn_t a, const bn_t b) {
  if (bn_is_zero(b)) { THROW(ERR_NO_VALID); }
  bn_div_impf_with_temps(c, d, xx, yy, rr, qq, a, b);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_gcd_basic(bn_t c, const bn_t a, const bn_t b) {
  bn_t u, v;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    return;
  }


  TRY {
    bn_new(u);
    bn_new(v);

    bn_abs(u, a);
    bn_abs(v, b);
    while (!bn_is_zero(v)) {
      bn_copy(c, v);
      bn_mod(v, u, v);
      bn_copy(u, c);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

// a,b inputs - destroyed, c,d,e outputs
void bn_gcdext(bn_t u, bn_t v, bn_t g, bn_t s) {
  mp_size_t sn;
  //  bn_zero(s);
  //  bn_zero(g);
  //  bn_trim(u);
  //  bn_trim(v);
  mp_size_t size = mpn_gcdext(g->dp, s->dp, &sn, u->dp, u->used, v->dp, v->used);
  g->used = size;
  g->sign = BN_POS;

  if (sn < 0) {
    s->used = -sn;
    s->sign = (u->sign) ? BN_POS : BN_NEG;
  } else {
    s->used = sn;
    s->sign = u->sign;
  }
  bn_trim(s);
  bn_trim(g);
}

void bn_gcd_ext_basic(bn_t c, bn_t d, bn_t e, const bn_t a, const bn_t b) {
  bn_t u, v, x_1, y_1, q, r;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    bn_zero(d);
    if (e != NULL) { bn_set_dig(e, 1); }
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    bn_set_dig(d, 1);
    if (e != NULL) { bn_zero(e); }
    return;
  }

  TRY {
    bn_new(u);
    bn_new(v);
    bn_new(x_1);
    bn_new(y_1);
    bn_new(q);
    bn_new(r);

    bn_abs(u, a);
    bn_abs(v, b);

    bn_zero(x_1);
    bn_set_dig(y_1, 1);

    if (e != NULL) {
      bn_set_dig(d, 1);
      bn_zero(e);

      while (!bn_is_zero(v)) {
        bn_div_rem(q, r, u, v);

        bn_copy(u, v);
        bn_copy(v, r);

        bn_mul(c, q, x_1);
        bn_sub(r, d, c);
        bn_copy(d, x_1);
        bn_copy(x_1, r);

        bn_mul(c, q, y_1);
        bn_sub(r, e, c);
        bn_copy(e, y_1);
        bn_copy(y_1, r);
      }
    } else {
      bn_set_dig(d, 1);

      while (!bn_is_zero(v)) {
        bn_div_rem(q, r, u, v);

        bn_copy(u, v);
        bn_copy(v, r);

        bn_mul(c, q, x_1);
        bn_sub(r, d, c);
        bn_copy(d, x_1);
        bn_copy(x_1, r);
      }
    }
    bn_copy(c, u);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_lehme(bn_t c, const bn_t a, const bn_t b) {
  bn_t x, y, u, v, t0, t1, t2, t3;
  dig_t _x, _y, q, _q, t, _t;
  dis_t _a, _b, _c, _d;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    return;
  }

  /*
   * Taken from Handbook of Hyperelliptic and Elliptic Cryptography.
   */
  TRY {
    bn_new(x);
    bn_new(y);
    bn_new(u);
    bn_new(v);
    bn_new(t0);
    bn_new(t1);
    bn_new(t2);
    bn_new(t3);

    if (bn_cmp(a, b) == CMP_GT) {
      bn_abs(x, a);
      bn_abs(y, b);
    } else {
      bn_abs(x, b);
      bn_abs(y, a);
    }
    while (y->used > 1) {
      bn_rsh(u, x, bn_bits(x) - BN_DIGIT);
      _x = u->dp[0];
      bn_rsh(v, y, bn_bits(x) - BN_DIGIT);
      _y = v->dp[0];
      _a = _d = 1;
      _b = _c = 0;
      t = 0;
      if (_y != 0) {
        q = _x / _y;
        t = _x % _y;
      }
      if (t >= ((dig_t)1 << (BN_DIGIT / 2))) {
        while (1) {
          _q = _y / t;
          _t = _y % t;
          if (_t < ((dig_t)1 << (BN_DIGIT / 2))) { break; }
          _x = _y;
          _y = t;
          t = _a - q * _c;
          _a = _c;
          _c = t;
          t = _b - q * _d;
          _b = _d;
          _d = t;
          t = _t;
          q = _q;
        }
      }
      if (_b == 0) {
        bn_mod(t0, x, y);
        bn_copy(x, y);
        bn_copy(y, t0);
      } else {
        bn_rsh(u, x, bn_bits(x) - 2 * BN_DIGIT);
        bn_rsh(v, y, bn_bits(x) - 2 * BN_DIGIT);
        if (_a < 0) {
          bn_mul_dig(t0, u, -_a);
          bn_neg(t0, t0);
        } else {
          bn_mul_dig(t0, u, _a);
        }
        if (_b < 0) {
          bn_mul_dig(t1, v, -_b);
          bn_neg(t1, t1);
        } else {
          bn_mul_dig(t1, v, _b);
        }
        if (_c < 0) {
          bn_mul_dig(t2, u, -_c);
          bn_neg(t2, t2);
        } else {
          bn_mul_dig(t2, u, _c);
        }
        if (_d < 0) {
          bn_mul_dig(t3, v, -_d);
          bn_neg(t3, t3);
        } else {
          bn_mul_dig(t3, v, _d);
        }
        bn_add(u, t0, t1);
        bn_add(v, t2, t3);
        bn_rsh(t0, u, bn_bits(u) - BN_DIGIT);
        _x = t0->dp[0];
        bn_rsh(t1, v, bn_bits(u) - BN_DIGIT);
        _y = t1->dp[0];
        t = 0;
        if (_y != 0) {
          q = _x / _y;
          t = _x % _y;
        }
        if (t >= ((dig_t)1 << BN_DIGIT / 2)) {
          while (1) {
            _q = _y / t;
            _t = _y % t;
            if (_t < ((dig_t)1 << BN_DIGIT / 2)) { break; }
            _x = _y;
            _y = t;
            t = _a - q * _c;
            _a = _c;
            _c = t;
            t = _b - q * _d;
            _b = _d;
            _d = t;
            t = _t;
            q = _q;
          }
        }
        if (_a < 0) {
          bn_mul_dig(t0, x, -_a);
          bn_neg(t0, t0);
        } else {
          bn_mul_dig(t0, x, _a);
        }
        if (_b < 0) {
          bn_mul_dig(t1, y, -_b);
          bn_neg(t1, t1);
        } else {
          bn_mul_dig(t1, y, _b);
        }
        if (_c < 0) {
          bn_mul_dig(t2, x, -_c);
          bn_neg(t2, t2);
        } else {
          bn_mul_dig(t2, x, _c);
        }
        if (_d < 0) {
          bn_mul_dig(t3, y, -_d);
          bn_neg(t3, t3);
        } else {
          bn_mul_dig(t3, y, _d);
        }
        bn_add(x, t0, t1);
        bn_add(y, t2, t3);
      }
    }
    bn_gcd_ext_dig(c, u, v, x, y->dp[0]);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_ext_lehme(bn_t c, bn_t d, bn_t e, const bn_t a, const bn_t b) {
  bn_t x, y, u, v, t0, t1, t2, t3, t4;
  dig_t _x, _y, q, _q, t, _t;
  dis_t _a, _b, _c, _d;
  int swap;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    bn_zero(d);
    if (e != NULL) { bn_set_dig(e, 1); }
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    bn_set_dig(d, 1);
    if (e != NULL) { bn_zero(e); }
    return;
  }

  /*
   * Taken from Handbook of Hyperelliptic and Elliptic Cryptography.
   */
  TRY {
    bn_new(x);
    bn_new(y);
    bn_new(u);
    bn_new(v);
    bn_new(t0);
    bn_new(t1);
    bn_new(t2);
    bn_new(t3);
    bn_new(t4);

    if (bn_cmp(a, b) != CMP_LT) {
      bn_abs(x, a);
      bn_abs(y, b);
      swap = 0;
    } else {
      bn_abs(x, b);
      bn_abs(y, a);
      swap = 1;
    }

    bn_zero(t4);
    bn_set_dig(d, 1);

    while (y->used > 1) {
      bn_rsh(u, x, bn_bits(x) - BN_DIGIT);
      _x = u->dp[0];
      bn_rsh(v, y, bn_bits(x) - BN_DIGIT);
      _y = v->dp[0];
      _a = _d = 1;
      _b = _c = 0;
      t = 0;
      if (_y != 0) {
        q = _x / _y;
        t = _x % _y;
      }
      if (t >= ((dig_t)1 << (BN_DIGIT / 2))) {
        while (1) {
          _q = _y / t;
          _t = _y % t;
          if (_t < ((dig_t)1 << (BN_DIGIT / 2))) { break; }
          _x = _y;
          _y = t;
          t = _a - q * _c;
          _a = _c;
          _c = t;
          t = _b - q * _d;
          _b = _d;
          _d = t;
          t = _t;
          q = _q;
        }
      }
      if (_b == 0) {
        bn_div_rem(t1, t0, x, y);
        bn_copy(x, y);
        bn_copy(y, t0);
        bn_mul(t1, t1, d);
        bn_sub(t1, t4, t1);
        bn_copy(t4, d);
        bn_copy(d, t1);
      } else {
        bn_rsh(u, x, bn_bits(x) - 2 * BN_DIGIT);
        bn_rsh(v, y, bn_bits(x) - 2 * BN_DIGIT);
        if (_a < 0) {
          bn_mul_dig(t0, u, -_a);
          bn_neg(t0, t0);
        } else {
          bn_mul_dig(t0, u, _a);
        }
        if (_b < 0) {
          bn_mul_dig(t1, v, -_b);
          bn_neg(t1, t1);
        } else {
          bn_mul_dig(t1, v, _b);
        }
        if (_c < 0) {
          bn_mul_dig(t2, u, -_c);
          bn_neg(t2, t2);
        } else {
          bn_mul_dig(t2, u, _c);
        }
        if (_d < 0) {
          bn_mul_dig(t3, v, -_d);
          bn_neg(t3, t3);
        } else {
          bn_mul_dig(t3, v, _d);
        }
        bn_add(u, t0, t1);
        bn_add(v, t2, t3);
        bn_rsh(t0, u, bn_bits(u) - BN_DIGIT);
        _x = t0->dp[0];
        bn_rsh(t1, v, bn_bits(u) - BN_DIGIT);
        _y = t1->dp[0];
        t = 0;
        if (_y != 0) {
          q = _x / _y;
          t = _x % _y;
        }
        if (t >= ((dig_t)1 << BN_DIGIT / 2)) {
          while (1) {
            _q = _y / t;
            _t = _y % t;
            if (_t < ((dig_t)1 << BN_DIGIT / 2)) { break; }
            _x = _y;
            _y = t;
            t = _a - q * _c;
            _a = _c;
            _c = t;
            t = _b - q * _d;
            _b = _d;
            _d = t;
            t = _t;
            q = _q;
          }
        }
        if (_a < 0) {
          bn_mul_dig(t0, x, -_a);
          bn_neg(t0, t0);
        } else {
          bn_mul_dig(t0, x, _a);
        }
        if (_b < 0) {
          bn_mul_dig(t1, y, -_b);
          bn_neg(t1, t1);
        } else {
          bn_mul_dig(t1, y, _b);
        }
        if (_c < 0) {
          bn_mul_dig(t2, x, -_c);
          bn_neg(t2, t2);
        } else {
          bn_mul_dig(t2, x, _c);
        }
        if (_d < 0) {
          bn_mul_dig(t3, y, -_d);
          bn_neg(t3, t3);
        } else {
          bn_mul_dig(t3, y, _d);
        }
        bn_add(x, t0, t1);
        bn_add(y, t2, t3);

        if (_a < 0) {
          bn_mul_dig(t0, t4, -_a);
          bn_neg(t0, t0);
        } else {
          bn_mul_dig(t0, t4, _a);
        }
        if (_b < 0) {
          bn_mul_dig(t1, d, -_b);
          bn_neg(t1, t1);
        } else {
          bn_mul_dig(t1, d, _b);
        }
        if (_c < 0) {
          bn_mul_dig(t2, t4, -_c);
          bn_neg(t2, t2);
        } else {
          bn_mul_dig(t2, t4, _c);
        }
        if (_d < 0) {
          bn_mul_dig(t3, d, -_d);
          bn_neg(t3, t3);
        } else {
          bn_mul_dig(t3, d, _d);
        }
        bn_add(t4, t0, t1);
        bn_add(d, t2, t3);
      }
    }
    bn_gcd_ext_dig(c, u, v, x, y->dp[0]);
    if (!swap) {
      bn_mul(t0, t4, u);
      bn_mul(t1, d, v);
      bn_add(t4, t0, t1);
      bn_mul(x, b, t4);
      bn_sub(x, c, x);
      bn_div(d, x, a);
      if (bn_sign(x) != bn_sign(a)) { bn_add_dig(d, d, 1); }
    } else {
      bn_mul(t0, t4, u);
      bn_mul(t1, d, v);
      bn_add(d, t0, t1);
      bn_mul(x, a, d);
      bn_sub(x, c, x);
      bn_div(t4, x, b);
      if (bn_sign(x) != bn_sign(b)) { bn_add_dig(t4, t4, 1); }
    }
    if (e != NULL) { bn_copy(e, t4); }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_stein(bn_t c, const bn_t a, const bn_t b) {
  bn_t u, v, t;
  int shift;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    return;
  }


  TRY {
    bn_new(u);
    bn_new(v);
    bn_new(t);

    bn_abs(u, a);
    bn_abs(v, b);

    shift = 0;
    while (bn_is_even(u) && bn_is_even(v)) {
      bn_hlv(u, u);
      bn_hlv(v, v);
      shift++;
    }
    while (!bn_is_zero(u)) {
      while (bn_is_even(u)) { bn_hlv(u, u); }
      while (bn_is_even(v)) { bn_hlv(v, v); }
      bn_sub(t, u, v);
      bn_abs(t, t);
      bn_hlv(t, t);
      if (bn_cmp(u, v) != CMP_LT) {
        bn_copy(u, t);
      } else {
        bn_copy(v, t);
      }
    }
    bn_lsh(c, v, shift);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_ext_stein(bn_t c, bn_t d, bn_t e, const bn_t a, const bn_t b) {
  bn_t x, y, u, v, _a, _b, _e;
  int shift, found;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    bn_zero(d);
    if (e != NULL) { bn_set_dig(e, 1); }
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    bn_set_dig(d, 1);
    if (e != NULL) { bn_zero(e); }
    return;
  }

  TRY {
    bn_new(x);
    bn_new(y);
    bn_new(u);
    bn_new(v);
    bn_new(_a);
    bn_new(_b);
    bn_new(_e);

    bn_abs(x, a);
    bn_abs(y, b);

    /* g = 1. */
    shift = 0;
    /* While x and y are both even, x = x/2 and y = y/2, g = 2g. */
    while (bn_is_even(x) && bn_is_even(y)) {
      bn_hlv(x, x);
      bn_hlv(y, y);
      shift++;
    }

    bn_copy(u, x);
    bn_copy(v, y);

    /* u = x, y = v, A = 1, B = 0, C = 0, D = 1. */
    bn_set_dig(_a, 1);
    bn_zero(_b);
    bn_zero(d);
    bn_set_dig(_e, 1);

    found = 0;
    while (!found) {
      /* While u is even, u = u/2. */
      while ((u->dp[0] & 0x01) == 0) {
        bn_hlv(u, u);
        /* If A = B = 0 (mod 2) then A = A/2, B = B/2. */
        if ((_a->dp[0] & 0x01) == 0 && (_b->dp[0] & 0x01) == 0) {
          bn_hlv(_a, _a);
          bn_hlv(_b, _b);
        } else {
          /* Otherwise A = (A + y)/2, B = (B - x)/2. */
          bn_add(_a, _a, y);
          bn_hlv(_a, _a);
          bn_sub(_b, _b, x);
          bn_hlv(_b, _b);
        }
      }
      /* While v is even, v = v/2. */
      while ((v->dp[0] & 0x01) == 0) {
        bn_hlv(v, v);
        /* If C = D = 0 (mod 2) then C = C/2, D = D/2. */
        if ((d->dp[0] & 0x01) == 0 && (_e->dp[0] & 0x01) == 0) {
          bn_hlv(d, d);
          bn_hlv(_e, _e);
        } else {
          /* Otherwise C = (C + y)/2, D = (D - x)/2. */
          bn_add(d, d, y);
          bn_hlv(d, d);
          bn_sub(_e, _e, x);
          bn_hlv(_e, _e);
        }
      }
      /* If u >= v then u = u - v, A = A - C, B = B - D. */
      if (bn_cmp(u, v) != CMP_LT) {
        bn_sub(u, u, v);
        bn_sub(_a, _a, d);
        bn_sub(_b, _b, _e);
      } else {
        /* Otherwise, v = v - u, C = C - a, D = D - B. */
        bn_sub(v, v, u);
        bn_sub(d, d, _a);
        bn_sub(_e, _e, _b);
      }
      /* If u = 0 then d = C, e = D and return (d, e, g * v). */
      if (bn_is_zero(u)) {
        bn_lsh(c, v, shift);
        found = 1;
      }
    }
    if (e != NULL) { bn_copy(e, _e); }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_ext_mid(bn_t c, bn_t d, bn_t e, bn_t f, const bn_t a, const bn_t b) {
  bn_t p, q, r, s, t, u, v, x, w, y, z;

  if (bn_is_zero(a)) {
    bn_abs(c, b);
    bn_zero(d);
    bn_zero(e);
    return;
  }

  if (bn_is_zero(b)) {
    bn_abs(c, a);
    bn_set_dig(d, 1);
    bn_set_dig(e, 1);
    return;
  }

  TRY {
    bn_new(p);
    bn_new(q);
    bn_new(r);
    bn_new(s);
    bn_new(t);
    bn_new(u);
    bn_new(v);
    bn_new(x);
    bn_new(w);
    bn_new(y);
    bn_new(z);

    if (bn_cmp_abs(a, b) == CMP_GT) {
      bn_abs(u, a);
      bn_abs(v, b);
    } else {
      bn_abs(u, b);
      bn_abs(v, a);
    }

    bn_srt(p, u);

    bn_set_dig(x, 1);
    bn_zero(t);

    int wait = 0;
    while (!bn_is_zero(v)) {
      bn_div_rem(q, r, u, v);

      bn_copy(u, v);
      bn_copy(v, r);

      bn_mul(s, q, x);
      bn_sub(s, t, s);
      bn_copy(t, x);
      bn_copy(x, s);

      if (wait) {
        bn_copy(e, r);
        bn_neg(f, x);
        wait = 0;
      }
      if (bn_cmp(u, p) == CMP_GT) {
        bn_copy(c, r);
        bn_neg(d, x);
        bn_copy(w, u);
        bn_neg(y, t);
        wait = 1;
      }
    }
    /* Compute r as the norm of vector (w, y). */
    bn_sqr(s, w);
    bn_sqr(t, y);
    bn_add(t, t, s);

    /* Compute q as the norm of vector (e, f). */
    bn_sqr(r, e);
    bn_sqr(q, f);
    bn_add(q, q, r);

    /* Output (e, f) as the vector of smaller norm. */
    if (bn_cmp(t, q) == CMP_LT) {
      bn_copy(e, w);
      bn_copy(f, y);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_gcd_dig(bn_t c, const bn_t a, dig_t b) {
  dig_t _u, _v, _t = 0;

  if (bn_is_zero(a)) {
    bn_set_dig(c, b);
    return;
  }

  if (b == 0) {
    bn_abs(c, a);
    return;
  }

  bn_mod_dig(&(c->dp[0]), a, b);
  _v = c->dp[0];
  _u = b;
  while (_v != 0) {
    _t = _v;
    _v = _u % _v;
    _u = _t;
  }
  bn_set_dig(c, _u);
}

void bn_gcd_ext_dig(bn_t c, bn_t d, bn_t e, const bn_t a, const dig_t b) {
  bn_t u, v, x1, y1, q, r;
  dig_t _v, _q, _t, _u;

  if (d == NULL && e == NULL) {
    bn_gcd_dig(c, a, b);
    return;
  }

  if (bn_is_zero(a)) {
    bn_set_dig(c, b);
    bn_zero(d);
    if (e != NULL) { bn_set_dig(e, 1); }
    return;
  }

  if (b == 0) {
    bn_abs(c, a);
    bn_set_dig(d, 1);
    if (e != NULL) { bn_zero(e); }
    return;
  }

  TRY {
    bn_new(u);
    bn_new(v);
    bn_new(x1);
    bn_new(y1);
    bn_new(q);
    bn_new(r);

    bn_abs(u, a);
    bn_set_dig(v, b);

    bn_zero(x1);
    bn_set_dig(y1, 1);
    bn_set_dig(d, 1);

    if (e != NULL) { bn_zero(e); }

    bn_div_rem(q, r, u, v);

    bn_copy(u, v);
    bn_copy(v, r);

    bn_mul(c, q, x1);
    bn_sub(r, d, c);
    bn_copy(d, x1);
    bn_copy(x1, r);

    if (e != NULL) {
      bn_mul(c, q, y1);
      bn_sub(r, e, c);
      bn_copy(e, y1);
      bn_copy(y1, r);
    }

    _v = v->dp[0];
    _u = u->dp[0];
    while (_v != 0) {
      _q = _u / _v;
      _t = _u % _v;

      _u = _v;
      _v = _t;

      bn_mul_dig(c, x1, _q);
      bn_sub(r, d, c);
      bn_copy(d, x1);
      bn_copy(x1, r);

      if (e != NULL) {
        bn_mul_dig(c, y1, _q);
        bn_sub(r, e, c);
        bn_copy(e, y1);
        bn_copy(y1, r);
      }
    }
    bn_set_dig(c, _u);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_init(bn_t a, int digits) {
  /* Verify if the number of digits is sane. */
  if (digits > BN_SIZE) {
    THROW(ERR_NO_PRECI);
  } else {
    digits = BN_SIZE;
  }
  if (a != NULL) {
    a->used = 0;
    a->alloc = digits;
    a->sign = BN_POS;
  }
}

void bn_clean(bn_t a) {
  if (a != NULL) {
    a->used = 0;
    a->sign = BN_POS;
  }
}

#ifdef USE_GROW
void bn_grow(bn_t a, int digits) {
  if (digits > BN_SIZE) { THROW(ERR_NO_PRECI) }
  (void)a;
}
#endif

void bn_trim(bn_t a) {
  while (a->used > 0 && a->dp[a->used - 1] == 0) { --(a->used); }
  /* Zero can't be negative. */
  if (a->used <= 0) {
    a->used = 1;
    a->dp[0] = 0;
    a->sign = BN_POS;
  }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_mod_2b(bn_t c, const bn_t a, int b) {
  int i, first, d;

  if (b <= 0) {
    bn_zero(c);
    return;
  }

  if (b >= (int)(a->used * BN_DIGIT)) {
    bn_copy(c, a);
    return;
  }

  bn_copy(c, a);

  SPLIT(b, d, b, BN_DIG_LOG);

  first = (d) + (b == 0 ? 0 : 1);
  for (i = first; i < c->used; i++) c->dp[i] = 0;

  c->dp[d] &= MASK(b);

  bn_trim(c);
}

void bn_mod_dig(dig_t *c, const bn_t a, dig_t b) { bn_div_rem_dig(NULL, c, a, b); }

void bn_mod_basic(bn_t c, const bn_t a, const bn_t m) { bn_div_rem(NULL, c, a, m); }

void bn_mod_pre_barrt(bn_t u, const bn_t m) {
  bn_set_2b(u, m->used * 2 * BN_DIGIT);
  bn_div(u, u, m);
}

void bn_mod_barrt(bn_t c, const bn_t a, const bn_t m, const bn_t u) {
  unsigned long mu;
  bn_t q, t;


  if (bn_cmp(a, m) == CMP_LT) {
    bn_copy(c, a);
    return;
  }
  TRY {
    bn_new(q);
    bn_new(t);
    bn_zero(t);

    mu = m->used;

    bn_rsh(q, a, (mu - 1) * BN_DIGIT);

    if (mu > ((dig_t)1) << (BN_DIGIT - 1)) {
      bn_mul(t, q, u);
    } else {
      if (q->used > u->used) {
        bn_muld_low(t->dp, q->dp, q->used, u->dp, u->used, mu, q->used + u->used);
      } else {
        bn_muld_low(t->dp, u->dp, u->used, q->dp, q->used, mu - (u->used - q->used), q->used + u->used);
      }
      t->used = q->used + u->used;
      bn_trim(t);
    }

    bn_rsh(q, t, (mu + 1) * BN_DIGIT);

    if (q->used > m->used) {
      bn_muld_low(t->dp, q->dp, q->used, m->dp, m->used, 0, q->used + 1);
    } else {
      bn_muld_low(t->dp, m->dp, m->used, q->dp, q->used, 0, mu + 1);
    }
    t->used = mu + 1;
    bn_trim(t);

    bn_mod_2b(q, t, BN_DIGIT * (mu + 1));
    bn_mod_2b(t, a, BN_DIGIT * (mu + 1));
    bn_sub(t, t, q);

    if (bn_sign(t) == BN_NEG) {
      bn_set_dig(q, (dig_t)1);
      bn_lsh(q, q, (mu + 1) * BN_DIGIT);
      bn_add(t, t, q);
    }

    while (bn_cmp(t, m) != CMP_LT) { bn_sub(t, t, m); }

    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mod_pre_monty(bn_t u, const bn_t m) {
  dig_t x, b;
  b = m->dp[0];

  if ((b & 0x01) == 0) { THROW(ERR_NO_VALID); }

  x = (((b + 2) & 4) << 1) + b; /* here x*a==1 mod 2**4 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**8 */
#if WORD != 8
  x *= 2 - b * x; /* here x*a==1 mod 2**16 */
#endif
#if WORD == 64 || WORD != 8 || WORD == 16
  x *= 2 - b * x; /* here x*a==1 mod 2**32 */
#endif
#if WORD == 64
  x *= 2 - b * x; /* here x*a==1 mod 2**64 */
#endif
  /* u = -1/m0 (mod 2^BN_DIGIT) */
  bn_set_dig(u, -x);
}

void bn_mod_monty_conv(bn_t c, const bn_t a, const bn_t m) {
  if (bn_sign(a) == BN_NEG) {
    bn_add(c, m, a);
  } else {
    bn_copy(c, a);
  }
  bn_lsh(c, c, m->used * BN_DIGIT);
  bn_mod(c, c, m);
}

void bn_mod_monty_back(bn_t c, const bn_t a, const bn_t m) {
  bn_t u;

  TRY {
    bn_new(u);

    bn_mod_pre_monty(u, m);
    bn_mod_monty(c, a, m, u);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mod_monty_basic(bn_t c, const bn_t a, const bn_t m, const bn_t u) {
  int digits, i;
  dig_t r, u0, *tmp;
  bn_t t;

  digits = 2 * m->used;

  TRY {
    bn_new_size(t, digits);
    bn_zero(t);
    bn_copy(t, a);

    u0 = u->dp[0];
    tmp = t->dp;

    for (i = 0; i < m->used; i++, tmp++) {
      r = (dig_t)(*tmp * u0);
      *tmp = bn_mula_low(tmp, m->dp, r, m->used);
    }
    if (bn_addn_low(t->dp, t->dp, tmp, m->used)) { bn_subn_low(t->dp, t->dp, m->dp, m->used); }
    t->used = m->used;
    bn_trim(t);

    if (bn_cmp_abs(t, m) != CMP_LT) { bn_sub(t, t, m); }

    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mod_monty_comba(bn_t c, const bn_t a, const bn_t m, const bn_t u) {
  int digits;
  bn_t t;

  digits = 2 * m->used;

  TRY {
    bn_new_size(t, digits);
    bn_zero(t);

    bn_modn_low(t->dp, a->dp, a->used, m->dp, m->used, u->dp[0]);
    t->used = m->used;

    bn_trim(t);
    if (bn_cmp_abs(t, m) != CMP_LT) { bn_sub(t, t, m); }
    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mod_pre_pmers(bn_t u, const bn_t m) {
  int bits;

  bits = bn_bits(m);

  bn_set_2b(u, bits);
  bn_sub(u, u, m);
}

void bn_mod_pmers(bn_t c, const bn_t a, const bn_t m, const bn_t u) {
  bn_t q, t, r;
  int bits;

  TRY {
    bn_new(q);
    bn_new(t);
    bn_new(r);

    bn_copy(t, a);

    bits = bn_bits(m);

    bn_rsh(q, t, bits);
    bn_mod_2b(r, t, bits);

    while (!bn_is_zero(q)) {
      if (u->used == 1) {
        bn_mul_dig(t, q, u->dp[0]);
      } else {
        bn_mul(t, q, u);
      }
      bn_rsh(q, t, bits);
      bn_mod_2b(t, t, bits);

      bn_add(r, r, t);
    }
    while (bn_cmp_abs(r, m) != CMP_LT) { bn_sub(r, r, m); }

    bn_copy(c, r);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Multiplies two multiple precision integers using recursive Karatsuba
 * multiplication.
 *
 * @param[out] c			- the result.
 * @param[in] a				- the first multiple precision integer.
 * @param[in] b				- the second multiple precision integer.
 * @param[in] level			- the number of Karatsuba steps to apply.
 */
void bn_mul_karat_imp(bn_t c, const bn_t a, const bn_t b, int level) {
  int h;
  bn_t a0, a1, b0, b1, a0b0, a1b1;
  bn_t t;
  const dig_t *tmpa, *tmpb;
  dig_t *t0;

  /* Compute half the digits of a or b. */

  // Should truncate to digits boundary
  h = MIN(a->used, b->used) >> 1;

  TRY {
    /* Allocate the temp variables. */
    bn_new(a0);
    bn_new(a1);
    bn_new(b0);
    bn_new(b1);
    bn_new(a0b0);
    bn_new(a1b1);
    bn_new(t);

    a0->used = b0->used = h;
    a1->used = a->used - h;
    b1->used = b->used - h;

    tmpa = a->dp;
    tmpb = b->dp;

    // No need for full copying here

    /* a = a1 || a0 */
    t0 = a0->dp;
    for (int i = 0; i < h; i++, t0++, tmpa++) *t0 = *tmpa;
    t0 = a1->dp;
    for (int i = 0; i < a1->used; i++, t0++, tmpa++) *t0 = *tmpa;

    /* b = b1 || b0 */
    t0 = b0->dp;
    for (int i = 0; i < h; i++, t0++, tmpb++) *t0 = *tmpb;
    t0 = b1->dp;
    for (int i = 0; i < b1->used; i++, t0++, tmpb++) *t0 = *tmpb;

    bn_trim(a0);
    bn_trim(b0);
    bn_trim(a1);
    bn_trim(b1);

    // Look at Parallel/Threaded calcs for this part

    /* a0b0 = a0 * b0 and a1b1 = a1 * b1 */
    if (level <= 1) {
      bn_mul_comba(a0b0, a0, b0);
      bn_mul_comba(a1b1, a1, b1);
    } else {
      bn_mul_karat_imp(a0b0, a0, b0, level - 1);
      bn_mul_karat_imp(a1b1, a1, b1, level - 1);
    }

    /* t = (a1 + a0) */
    bn_add(a1, a1, a0);
    /* t2 = (b1 + b0) */
    bn_add(b1, b1, b0);

    /* t = (a1 + a0)*(b1 + b0) */
    if (level <= 1) {
      bn_mul_comba(t, a1, b1);
    } else {
      bn_mul_karat_imp(t, a1, b1, level - 1);
    }
    /* t2 = (a0*b0 + a1*b1) */
    bn_sub(t, t, a0b0);

    /* t = (a1 + a0)*(b1 + b0) - (a0*b0 + a1*b1) */
    bn_sub(t, t, a1b1);

    /* t = (a1 + a0)*(b1 + b0) - (a0*b0 + a1*b1) << h digits */
    bn_lsh(t, t, h * BN_DIGIT);

    /* t2 = a1 * b1 << 2*h digits */
    bn_lsh(a1b1, a1b1, 2 * h * BN_DIGIT);

    /* t = t + a0*b0 */
    bn_add(t, t, a0b0);

    /* c = t + a1*b1 */
    bn_add(t, t, a1b1);

    t->sign = a->sign ^ b->sign;
    bn_copy(c, t);
    bn_trim(c);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_mul_dig(bn_t c, const bn_t a, dig_t b) {
  bn_grow(c, a->used + 1);
  c->sign = a->sign;
  c->dp[a->used] = bn_mul1_low(c->dp, a->dp, b, a->used);
  c->used = a->used + 1;
  bn_trim(c);
}

void bn_mul_sdig(bn_t c, const bn_t a, int64_t b) {
  dig_t ab = (b < 0) ? -b : b;
  bn_grow(c, a->used + 1);
  c->sign = (b < 0) ? ((a->sign) ? 0 : 1) : a->sign;
  c->dp[a->used] = bn_mul1_low(c->dp, a->dp, ab, a->used);
  c->used = a->used + 1;
  bn_trim(c);
}

void bn_mul_simp(bn_t c, const bn_t a, const bn_t b) {
  int i;
  dig_t carry;
  //  bn_clear(c);
  c->used = a->used + b->used;
  for (i = 0; i < a->used; i++) {
    // carry = mpn_addmul_1(c->dp+i, b->dp, b->used, *(a->dp+i));
    carry = mpn_mul_1(c->dp + i, b->dp, b->used, *(a->dp + i));
    // carry = bn_mula_low(c->dp + i, b->dp, *(a->dp + i), b->used);
    *(c->dp + i + b->used) = carry;
  }
  c->sign = a->sign ^ b->sign;
  bn_trim(c);
}

void bn_mul_simp_new(bn_t c, const bn_t a, const bn_t b) {
  int i;
  // bn_clear(c);
  c->used = a->used + b->used;
  dig_t *bp = (dig_t *)&b->dp[0];
  dig_t *ap = (dig_t *)&a->dp[0];
  for (i = 0; i < a->used; i++) {
    dig_t carry = 0;
    // This should be fixed but can be OPENMPed....
    for (int j = 0; j < b->used; j++, bp++, ap++) {
      dbl_t r = (dbl_t)(*bp) * (dbl_t)(*ap + j) + (dbl_t)(carry);
      *c->dp = (dig_t)r;
      carry = (dig_t)(r >> (dbl_t)BN_DIGIT);
    }
    // carry = bn_mula_low(c->dp + i, b->dp, *(a->dp + i), b->used);
    *(c->dp + i + b->used) = carry;
  }
  c->sign = a->sign ^ b->sign;
  bn_trim(c);
}

void bn_mul_basic(bn_t c, const bn_t a, const bn_t b) {
  int i;
  bn_t t;
  dig_t carry;

  TRY {
    /* We need a temporary variable so that c can be a or b. */
    bn_new_size(t, a->used + b->used);
    bn_zero(t);
    t->used = a->used + b->used;

    for (i = 0; i < a->used; i++) {
      carry = bn_mula_low(t->dp + i, b->dp, *(a->dp + i), b->used);
      *(t->dp + i + b->used) = carry;
    }
    t->sign = a->sign ^ b->sign;
    bn_trim(t);

    /* Swap c and t. */
    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mul_comba(bn_t t, const bn_t a, const bn_t b) {
  int digits;
  // bn_t t;

  TRY {
    digits = a->used + b->used;

    /* We need a temporary variable so that c can be a or b. */
    bn_new_size(t, digits);
    t->used = digits;

    if (a->used == b->used) {
      bn_muln_low(t->dp, a->dp, b->dp, a->used);
    } else {
      if (a->used > b->used) {
        bn_muld_low(t->dp, a->dp, a->used, b->dp, b->used, 0, a->used + b->used);
      } else {
        bn_muld_low(t->dp, b->dp, b->used, a->dp, a->used, 0, a->used + b->used);
      }
    }

    t->sign = a->sign ^ b->sign;
    bn_trim(t);
    // bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_mul_karat(bn_t c, const bn_t a, const bn_t b) { bn_mul_karat_imp(c, a, b, BN_KARAT); }

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Size of precomputation table.
 */
#define RELIC_TABLE_SIZE 64

static void swap_cond(dig_t *c, dig_t *a, int digits, dig_t cond) {
  dig_t mask, t;

  mask = -cond;
  for (int i = 0; i < digits; i++) {
    t = (a[i] ^ c[i]) & mask;
    a[i] ^= t;
    c[i] ^= t;
  }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

#if BN_MXP == BASIC || !defined(STRIP)

void bn_mxp_basic(bn_t c, const bn_t a, const bn_t b, const bn_t m) {
  int i, l;
  bn_t t, u, r;

  if (bn_is_zero(b)) {
    bn_set_dig(c, 1);
    return;
  }

  TRY {
    bn_new(t);
    bn_new(u);
    bn_new(r);

    bn_mod_pre(u, m);

    l = bn_bits(b);

#if BN_MOD == MONTY
    bn_mod_monty_conv(t, a, m);
#else
    bn_copy(t, a);
#endif

    bn_copy(r, t);

    for (i = l - 2; i >= 0; i--) {
      bn_sqr(r, r);
      bn_mod(r, r, m, u);
      if (bn_get_bit(b, i)) {
        bn_mul(r, r, t);
        bn_mod(r, r, m, u);
      }
    }

#if BN_MOD == MONTY
    bn_mod_monty_back(c, r, m);
#else
    bn_copy(c, r);
#endif
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

#endif

#if BN_MXP == SLIDE || !defined(STRIP)

void bn_mxp_slide(bn_t c, const bn_t a, const bn_t b, const bn_t m) {
  bn_t tab[RELIC_TABLE_SIZE], t, u, r;
  int i, j, l, w = 1;
  uint8_t win[RELIC_BN_BITS];

  TRY {
    /* Find window size. */
    i = bn_bits(b);
    if (i <= 21) {
      w = 2;
    } else if (i <= 32) {
      w = 3;
    } else if (i <= 128) {
      w = 4;
    } else if (i <= 256) {
      w = 5;
    } else {
      w = 6;
    }

    for (i = 1; i < (1 << w); i += 2) { bn_new(tab[i]); }

    bn_new(t);
    bn_new(u);
    bn_new(r);
    bn_mod_pre(u, m);

#if BN_MOD == MONTY
    bn_set_dig(r, 1);
    bn_mod_monty_conv(r, r, m);
    bn_mod_monty_conv(t, a, m);
#else /* BN_MOD == BARRT || BN_MOD == RADIX */
    bn_set_dig(r, 1);
    bn_copy(t, a);
#endif

    bn_copy(tab[1], t);
    bn_sqr(t, tab[1]);
    bn_mod(t, t, m, u);
    /* Create table. */
    for (i = 1; i < 1 << (w - 1); i++) {
      bn_mul(tab[2 * i + 1], tab[2 * i - 1], t);
      bn_mod(tab[2 * i + 1], tab[2 * i + 1], m, u);
    }

    l = RELIC_BN_BITS + 1;
    bn_rec_slw(win, &l, b, w);
    for (i = 0; i < l; i++) {
      if (win[i] == 0) {
        bn_sqr(r, r);
        bn_mod(r, r, m, u);
      } else {
        for (j = 0; j < util_bits_dig(win[i]); j++) {
          bn_sqr(r, r);
          bn_mod(r, r, m, u);
        }
        bn_mul(r, r, tab[win[i]]);
        bn_mod(r, r, m, u);
      }
    }
    bn_trim(r);
#if BN_MOD == MONTY
    bn_mod_monty_back(c, r, m);
#else
    bn_copy(c, r);
#endif
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

#endif

#if BN_MXP == MONTY || !defined(STRIP)

void bn_mxp_monty(bn_t c, const bn_t a, const bn_t b, const bn_t m) {
  bn_t tab[2], u;
  dig_t mask;
  int t;

  TRY {
    bn_new(u);
    bn_mod_pre(u, m);

    bn_new(tab[0]);
    bn_new(tab[1]);

#if BN_MOD == MONTY
    bn_set_dig(tab[0], 1);
    bn_mod_monty_conv(tab[0], tab[0], m);
    bn_mod_monty_conv(tab[1], a, m);
#else
    bn_set_dig(tab[0], 1);
    bn_copy(tab[1], a);
#endif

    for (int i = bn_bits(b) - 1; i >= 0; i--) {
      int j = bn_get_bit(b, i);
      swap_cond(tab[0]->dp, tab[1]->dp, BN_DIGS, j ^ 1);
      mask = -(j ^ 1);
      t = (tab[0]->used ^ tab[1]->used) & mask;
      tab[0]->used ^= t;
      tab[1]->used ^= t;
      bn_mul(tab[0], tab[0], tab[1]);
      bn_mod(tab[0], tab[0], m, u);
      bn_sqr(tab[1], tab[1]);
      bn_mod(tab[1], tab[1], m, u);
      swap_cond(tab[0]->dp, tab[1]->dp, BN_DIGS, j ^ 1);
      mask = -(j ^ 1);
      t = (tab[0]->used ^ tab[1]->used) & mask;
      tab[0]->used ^= t;
      tab[1]->used ^= t;
    }

#if BN_MOD == MONTY
    bn_mod_monty_back(c, tab[0], m);
#else
    bn_copy(c, tab[0]);
#endif
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

#endif

void bn_mxp_dig(bn_t c, const bn_t a, dig_t b, const bn_t m) {
  int i, l;
  bn_t t, u, r;

  if (b == 0) {
    bn_set_dig(c, 1);
    return;
  }

  TRY {
    bn_new(t);
    bn_new(u);
    bn_new(r);

    bn_mod_pre(u, m);

    l = util_bits_dig(b);

#if BN_MOD == MONTY
    bn_mod_monty_conv(t, a, m);
#else
    bn_copy(t, a);
#endif

    bn_copy(r, t);

    for (i = l - 2; i >= 0; i--) {
      bn_sqr(r, r);
      bn_mod(r, r, m, u);
      if (b & ((dig_t)1 << i)) {
        bn_mul(r, r, t);
        bn_mod(r, r, m, u);
      }
    }

#if BN_MOD == MONTY
    bn_mod_monty_back(c, r, m);
#else
    bn_copy(c, r);
#endif
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}
/**
 * Implementation of the prime number generation and testing functions.
 *
 * Strong prime generation is based on Gordon's Algorithm, taken from Handbook
 * of Applied Cryptography.
 *
 * @ingroup bn
 */

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Number of trial division tests.
 */
#define BASIC_TESTS ((int)(sizeof(primes) / sizeof(dig_t)))

/**
 * Small prime numbers table.
 */
static const dig_t primes[] = {
    0x0002, 0x0003, 0x0005, 0x0007, 0x000B, 0x000D, 0x0011, 0x0013, 0x0017, 0x001D, 0x001F, 0x0025, 0x0029,
    0x002B, 0x002F, 0x0035, 0x003B, 0x003D, 0x0043, 0x0047, 0x0049, 0x004F, 0x0053, 0x0059, 0x0061, 0x0065,
    0x0067, 0x006B, 0x006D, 0x0071, 0x007F, 0x0083, 0x0089, 0x008B, 0x0095, 0x0097, 0x009D, 0x00A3, 0x00A7,
    0x00AD, 0x00B3, 0x00B5, 0x00BF, 0x00C1, 0x00C5, 0x00C7, 0x00D3, 0x00DF,
#if WORD > 8
    0x00E3, 0x00E5, 0x00E9, 0x00EF, 0x00F1, 0x00FB, 0x0101, 0x0107, 0x010D, 0x010F, 0x0115, 0x0119, 0x011B,
    0x0125, 0x0133, 0x0137,

    0x0139, 0x013D, 0x014B, 0x0151, 0x015B, 0x015D, 0x0161, 0x0167, 0x016F, 0x0175, 0x017B, 0x017F, 0x0185,
    0x018D, 0x0191, 0x0199, 0x01A3, 0x01A5, 0x01AF, 0x01B1, 0x01B7, 0x01BB, 0x01C1, 0x01C9, 0x01CD, 0x01CF,
    0x01D3, 0x01DF, 0x01E7, 0x01EB, 0x01F3, 0x01F7, 0x01FD, 0x0209, 0x020B, 0x021D, 0x0223, 0x022D, 0x0233,
    0x0239, 0x023B, 0x0241, 0x024B, 0x0251, 0x0257, 0x0259, 0x025F, 0x0265, 0x0269, 0x026B, 0x0277, 0x0281,
    0x0283, 0x0287, 0x028D, 0x0293, 0x0295, 0x02A1, 0x02A5, 0x02AB, 0x02B3, 0x02BD, 0x02C5, 0x02CF,

    0x02D7, 0x02DD, 0x02E3, 0x02E7, 0x02EF, 0x02F5, 0x02F9, 0x0301, 0x0305, 0x0313, 0x031D, 0x0329, 0x032B,
    0x0335, 0x0337, 0x033B, 0x033D, 0x0347, 0x0355, 0x0359, 0x035B, 0x035F, 0x036D, 0x0371, 0x0373, 0x0377,
    0x038B, 0x038F, 0x0397, 0x03A1, 0x03A9, 0x03AD, 0x03B3, 0x03B9, 0x03C7, 0x03CB, 0x03D1, 0x03D7, 0x03DF,
    0x03E5, 0x03F1, 0x03F5, 0x03FB, 0x03FD, 0x0407, 0x0409, 0x040F, 0x0419, 0x041B, 0x0425, 0x0427, 0x042D,
    0x043F, 0x0443, 0x0445, 0x0449, 0x044F, 0x0455, 0x045D, 0x0463, 0x0469, 0x047F, 0x0481, 0x048B,

    0x0493, 0x049D, 0x04A3, 0x04A9, 0x04B1, 0x04BD, 0x04C1, 0x04C7, 0x04CD, 0x04CF, 0x04D5, 0x04E1, 0x04EB,
    0x04FD, 0x04FF, 0x0503, 0x0509, 0x050B, 0x0511, 0x0515, 0x0517, 0x051B, 0x0527, 0x0529, 0x052F, 0x0551,
    0x0557, 0x055D, 0x0565, 0x0577, 0x0581, 0x058F, 0x0593, 0x0595, 0x0599, 0x059F, 0x05A7, 0x05AB, 0x05AD,
    0x05B3, 0x05BF, 0x05C9, 0x05CB, 0x05CF, 0x05D1, 0x05D5, 0x05DB, 0x05E7, 0x05F3, 0x05FB, 0x0607, 0x060D,
    0x0611, 0x0617, 0x061F, 0x0623, 0x062B, 0x062F, 0x063D, 0x0641, 0x0647, 0x0649, 0x064D, 0x0653,

    0x0655, 0x065B, 0x0665, 0x0679, 0x067F, 0x0683, 0x0685, 0x069D, 0x06A1, 0x06A3, 0x06AD, 0x06B9, 0x06BB,
    0x06C5, 0x06CD, 0x06D3, 0x06D9, 0x06DF, 0x06F1, 0x06F7, 0x06FB, 0x06FD, 0x0709, 0x0713, 0x071F, 0x0727,
    0x0737, 0x0745, 0x074B, 0x074F, 0x0751, 0x0755, 0x0757, 0x0761, 0x076D, 0x0773, 0x0779, 0x078B, 0x078D,
    0x079D, 0x079F, 0x07B5, 0x07BB, 0x07C3, 0x07C9, 0x07CD, 0x07CF, 0x07D3, 0x07DB, 0x07E1, 0x07EB, 0x07ED,
    0x07F7, 0x0805, 0x080F, 0x0815, 0x0821, 0x0823, 0x0827, 0x0829, 0x0833, 0x083F, 0x0841, 0x0851,

    0x0853, 0x0859, 0x085D, 0x085F, 0x0869, 0x0871, 0x0883, 0x089B, 0x089F, 0x08A5, 0x08AD, 0x08BD, 0x08BF,
    0x08C3, 0x08CB, 0x08DB, 0x08DD, 0x08E1, 0x08E9, 0x08EF, 0x08F5, 0x08F9, 0x0905, 0x0907, 0x091D, 0x0923,
    0x0925, 0x092B, 0x092F, 0x0935, 0x0943, 0x0949, 0x094D, 0x094F, 0x0955, 0x0959, 0x095F, 0x096B, 0x0971,
    0x0977, 0x0985, 0x0989, 0x098F, 0x099B, 0x09A3, 0x09A9, 0x09AD, 0x09C7, 0x09D9, 0x09E3, 0x09EB, 0x09EF,
    0x09F5, 0x09F7, 0x09FD, 0x0A13, 0x0A1F, 0x0A21, 0x0A31, 0x0A39, 0x0A3D, 0x0A49, 0x0A57, 0x0A61,

    0x0A63, 0x0A67, 0x0A6F, 0x0A75, 0x0A7B, 0x0A7F, 0x0A81, 0x0A85, 0x0A8B, 0x0A93, 0x0A97, 0x0A99, 0x0A9F,
    0x0AA9, 0x0AAB, 0x0AB5, 0x0ABD, 0x0AC1, 0x0ACF, 0x0AD9, 0x0AE5, 0x0AE7, 0x0AED, 0x0AF1, 0x0AF3, 0x0B03,
    0x0B11, 0x0B15, 0x0B1B, 0x0B23, 0x0B29, 0x0B2D, 0x0B3F, 0x0B47, 0x0B51, 0x0B57, 0x0B5D, 0x0B65, 0x0B6F,
    0x0B7B, 0x0B89, 0x0B8D, 0x0B93, 0x0B99, 0x0B9B, 0x0BB7, 0x0BB9, 0x0BC3, 0x0BCB, 0x0BCF, 0x0BDD, 0x0BE1,
    0x0BE9, 0x0BF5, 0x0BFB, 0x0C07, 0x0C0B, 0x0C11, 0x0C25, 0x0C2F, 0x0C31, 0x0C41, 0x0C5B, 0x0C5F,

    0x0C61, 0x0C6D, 0x0C73, 0x0C77, 0x0C83, 0x0C89, 0x0C91, 0x0C95, 0x0C9D, 0x0CB3, 0x0CB5, 0x0CB9, 0x0CBB,
    0x0CC7, 0x0CE3, 0x0CE5, 0x0CEB, 0x0CF1, 0x0CF7, 0x0CFB, 0x0D01, 0x0D03, 0x0D0F, 0x0D13, 0x0D1F, 0x0D21,
    0x0D2B, 0x0D2D, 0x0D3D, 0x0D3F, 0x0D4F, 0x0D55, 0x0D69, 0x0D79, 0x0D81, 0x0D85, 0x0D87, 0x0D8B, 0x0D8D,
    0x0DA3, 0x0DAB, 0x0DB7, 0x0DBD, 0x0DC7, 0x0DC9, 0x0DCD, 0x0DD3, 0x0DD5, 0x0DDB, 0x0DE5, 0x0DE7, 0x0DF3,
    0x0DFD, 0x0DFF, 0x0E09, 0x0E17, 0x0E1D, 0x0E21, 0x0E27, 0x0E2F, 0x0E35, 0x0E3B, 0x0E4B, 0x0E57,
#endif
};

#if BN_MOD == PMERS

/**
 * Computes c = a ^ b mod m.
 *
 * @param c				- the result.
 * @param a				- the basis.
 * @param b				- the exponent.
 * @param m				- the modulus.
 */
static void bn_exp(bn_t c, const bn_t a, const bn_t b, const bn_t m) {
  int i, l;
  bn_t t;

  TRY {
    bn_new(t);

    l = bn_bits(b);

    bn_copy(t, a);

    for (i = l - 2; i >= 0; i--) {
      bn_sqr(t, t);
      bn_mod(t, t, m);
      if (bn_get_bit(b, i)) {
        bn_mul(t, t, a);
        bn_mod(t, t, m);
      }
    }

    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

#endif

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

dig_t bn_get_prime(int pos) {
  if (pos >= BASIC_TESTS) { return 0; }
  return primes[pos];
}

int bn_is_prime(const bn_t a) {
  int result;

  result = 0;
  if (!bn_is_prime_basic(a)) { goto end; }

  if (!bn_is_prime_rabin(a)) { goto end; }

  result = 1;
end:
  return result;
}

int bn_is_prime_basic(const bn_t a) {
  dig_t t;
  int i, result;

  result = 1;

  if (bn_cmp_dig(a, 1) == CMP_EQ) { return 0; }

  /* Trial division. */
  for (i = 0; i < BASIC_TESTS; i++) {
    bn_mod_dig(&t, a, primes[i]);
    if (t == 0 && bn_cmp_dig(a, primes[i]) != CMP_EQ) {
      result = 0;
      break;
    }
  }
  return result;
}

int bn_is_prime_rabin(const bn_t a) {
  bn_t t, n1, y, r;
  int i, s, j, result, b, tests = 0;

  tests = 0;
  result = 1;

  if (bn_cmp_dig(a, 1) == CMP_EQ) { return 0; }

  TRY {
    /*
     * These values are taken from Table 4.4 inside Handbook of Applied
     * Cryptography.
     */
    b = bn_bits(a);
    if (b >= 1300) {
      tests = 2;
    } else if (b >= 850) {
      tests = 3;
    } else if (b >= 650) {
      tests = 4;
    } else if (b >= 550) {
      tests = 5;
    } else if (b >= 450) {
      tests = 6;
    } else if (b >= 400) {
      tests = 7;
    } else if (b >= 350) {
      tests = 8;
    } else if (b >= 300) {
      tests = 9;
    } else if (b >= 250) {
      tests = 12;
    } else if (b >= 200) {
      tests = 15;
    } else if (b >= 150) {
      tests = 18;
    } else {
      tests = 27;
    }

    bn_new(t);
    bn_new(n1);
    bn_new(y);
    bn_new(r);

    /* r = (n - 1)/2^s. */
    bn_sub_dig(n1, a, 1);
    s = 0;
    while (bn_is_even(n1)) {
      s++;
      bn_rsh(n1, n1, 1);
    }
    bn_lsh(r, n1, s);

    for (i = 0; i < tests; i++) {
      /* Fix the basis as the first few primes. */
      bn_set_dig(t, primes[i]);

      /* y = b^r mod a. */
#if BN_MOD != PMERS
      bn_mxp(y, t, r, a);
#else
      bn_exp(y, t, r, a);
#endif

      if (bn_cmp_dig(y, 1) != CMP_EQ && bn_cmp(y, n1) != CMP_EQ) {
        j = 1;
        while ((j <= (s - 1)) && bn_cmp(y, n1) != CMP_EQ) {
          bn_sqr(y, y);
          bn_mod(y, y, a);

          /* If y == 1 then composite. */
          if (bn_cmp_dig(y, 1) == CMP_EQ) {
            result = 0;
            break;
          }
          ++j;
        }

        /* If y != n1 then composite. */
        if (bn_cmp(y, n1) != CMP_EQ) {
          result = 0;
          break;
        }
      }
    }
  }
  CATCH_ANY {
    result = 0;
    THROW(ERR_CAUGHT);
  }
  return result;
}

/**
 * @file
 *
 * Implementation of the multiple precision integer recoding functions.
 *
 * @ingroup bn
 */

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Returns a maximum of eight contiguous bits from a multiple precision integer.
 *
 * @param[in] a				- the multiple precision integer.
 * @param[in] from			- the first bit position.
 * @param[in] to			- the last bit position, inclusive.
 * @return the bits in the chosen positions.
 */
static char get_bits(const bn_t a, int from, int to) {
  int f, t;
  dig_t mf, mt;

  SPLIT(from, f, from, BN_DIG_LOG);
  SPLIT(to, t, to, BN_DIG_LOG);

  if (f == t) {
    /* Same digit. */

    mf = MASK(from);
    mt = MASK(to + 1);

    if (to + 1 == BN_DIGIT) { mt = DMASK; }

    mf = mf ^ mt;

    return ((a->dp[f] & (mf)) >> from);
  } else {
    mf = MASK(BN_DIGIT - from) << from;
    mt = MASK(to + 1);

    return ((a->dp[f] & mf) >> from) | ((a->dp[t] & mt) << (BN_DIGIT - from));
  }
}

/**
 * Constant C for the partial reduction modulo (t^m - 1)/(t - 1).
 */
#define MOD_C 8

/**
 * Constant 2^C.
 */
#define MOD_2TC (1 << MOD_C)

/**
 * Mask to calculate reduction modulo 2^C.
 */
#define MOD_CMASK (MOD_2TC - 1)

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_rec_win(uint8_t *win, int *len, const bn_t k, int w) {
  int i, j, l;

  l = bn_bits(k);

  if (*len < CEIL(l, w)) { THROW(ERR_NO_BUFFER); }

  j = 0;
  for (i = 0; i < l - w; i += w) { win[j++] = get_bits(k, i, i + w - 1); }
  win[j++] = get_bits(k, i, bn_bits(k) - 1);
  *len = j;
}

void bn_rec_slw(uint8_t *win, int *len, const bn_t k, int w) {
  int i, j, l, s;

  l = bn_bits(k);

  if (*len < CEIL(l, w)) { THROW(ERR_NO_BUFFER); }

  i = l - 1;
  j = 0;
  while (i >= 0) {
    if (!bn_get_bit(k, i)) {
      i--;
      win[j++] = 0;
    } else {
      s = MAX(i - w + 1, 0);
      while (!bn_get_bit(k, s)) { s++; }
      win[j++] = get_bits(k, s, i);
      i = s - 1;
    }
  }
  *len = j;
}

void bn_rec_naf(int8_t *naf, int *len, const bn_t k, int w) {
  int i, l;
  bn_t t;
  dig_t t0, mask;
  int8_t u_i;

  if (*len < (bn_bits(k) + 1)) { THROW(ERR_NO_BUFFER); }

  TRY {
    bn_new(t);
    bn_abs(t, k);

    mask = MASK(w);
    l = (1 << w);

    i = 0;
    if (w == 2) {
      while (!bn_is_zero(t)) {
        if (!bn_is_even(t)) {
          bn_get_dig(&t0, t);
          u_i = 2 - (t0 & mask);
          if (u_i < 0) {
            bn_add_dig(t, t, -u_i);
          } else {
            bn_sub_dig(t, t, u_i);
          }
          *naf = u_i;
        } else {
          *naf = 0;
        }
        bn_hlv(t, t);
        i++;
        naf++;
      }
    } else {
      while (!bn_is_zero(t)) {
        if (!bn_is_even(t)) {
          bn_get_dig(&t0, t);
          u_i = t0 & mask;
          if (u_i > l / 2) { u_i = (int8_t)(u_i - l); }
          if (u_i < 0) {
            bn_add_dig(t, t, -u_i);
          } else {
            bn_sub_dig(t, t, u_i);
          }
          *naf = u_i;
        } else {
          *naf = 0;
        }
        bn_hlv(t, t);
        i++;
        naf++;
      }
    }
    *len = i;
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_tnaf_get(uint8_t *t, int8_t *beta, int8_t *gama, int8_t u, int w) {
  if (u == -1) {
    switch (w) {
      case 2:
      case 3:
        *t = 2;
        break;
      case 4:
        *t = 10;
        break;
      case 5:
      case 6:
        *t = 26;
        break;
      case 7:
      case 8:
        *t = 90;
        break;
    }
  } else {
    switch (w) {
      case 2:
        *t = 2;
        break;
      case 3:
      case 4:
      case 5:
        *t = 6;
        break;
      case 6:
      case 7:
        *t = 38;
        break;
      case 8:
        *t = 166;
        break;
    }
  }

  beta[0] = 1;
  gama[0] = 0;

  if (w >= 3) {
    beta[1] = 1;
    gama[1] = (int8_t)-u;
  }

  if (w >= 4) {
    beta[1] = -3;
    beta[2] = -1;
    beta[3] = 1;
    gama[1] = gama[2] = gama[3] = (int8_t)u;
  }

  if (w >= 5) {
    beta[4] = -3;
    beta[5] = -1;
    beta[6] = beta[7] = 1;
    gama[4] = gama[5] = gama[6] = (int8_t)(2 * u);
    gama[7] = (int8_t)(-3 * u);
  }

  if (w >= 6) {
    beta[1] = beta[8] = beta[14] = 3;
    beta[2] = beta[9] = beta[15] = 5;
    beta[3] = -5;
    beta[4] = beta[10] = beta[11] = -3;
    beta[5] = beta[12] = -1;
    beta[6] = beta[7] = beta[13] = 1;
    gama[1] = gama[2] = 0;
    gama[3] = gama[4] = gama[5] = gama[6] = (int8_t)(2 * u);
    gama[7] = gama[8] = gama[9] = (int8_t)(-3 * u);
    gama[10] = (int8_t)(4 * u);
    gama[11] = gama[12] = gama[13] = (int8_t)(-u);
    gama[14] = gama[15] = (int8_t)(-u);
  }

  if (w >= 7) {
    beta[3] = beta[22] = beta[29] = 7;
    beta[4] = beta[16] = beta[23] = -5;
    beta[5] = beta[10] = beta[17] = beta[24] = -3;
    beta[6] = beta[11] = beta[18] = beta[25] = beta[30] = -1;
    beta[7] = beta[12] = beta[14] = beta[19] = beta[26] = beta[31] = 1;
    beta[8] = beta[13] = beta[20] = beta[27] = 3;
    beta[9] = beta[21] = beta[28] = 5;
    beta[15] = -7;
    gama[3] = 0;
    gama[4] = gama[5] = gama[6] = (int8_t)(-3 * u);
    gama[11] = gama[12] = gama[13] = (int8_t)(4 * u);
    gama[14] = (int8_t)(-6 * u);
    gama[15] = gama[16] = gama[17] = gama[18] = (int8_t)u;
    gama[19] = gama[20] = gama[21] = gama[22] = (int8_t)u;
    gama[23] = gama[24] = gama[25] = gama[26] = (int8_t)(-2 * u);
    gama[27] = gama[28] = gama[29] = (int8_t)(-2 * u);
    gama[30] = gama[31] = (int8_t)(5 * u);
  }

  if (w == 8) {
    beta[10] = beta[17] = beta[48] = beta[55] = beta[62] = 7;
    beta[11] = beta[18] = beta[49] = beta[56] = beta[63] = 9;
    beta[12] = beta[22] = beta[29] = -3;
    beta[36] = beta[43] = beta[50] = -3;
    beta[13] = beta[23] = beta[30] = beta[37] = -1;
    beta[44] = beta[51] = beta[58] = -1;
    beta[14] = beta[24] = beta[31] = beta[38] = 1;
    beta[45] = beta[52] = beta[59] = 1;
    beta[15] = beta[32] = beta[39] = beta[46] = beta[53] = beta[60] = 3;
    beta[16] = beta[40] = beta[47] = beta[54] = beta[61] = 5;
    beta[19] = beta[57] = 11;
    beta[20] = beta[27] = beta[34] = beta[41] = -7;
    beta[21] = beta[28] = beta[35] = beta[42] = -5;
    beta[25] = -11;
    beta[26] = beta[33] = -9;
    gama[10] = gama[11] = (int8_t)(-3 * u);
    gama[12] = gama[13] = gama[14] = gama[15] = (int8_t)(-6 * u);
    gama[16] = gama[17] = gama[18] = gama[19] = (int8_t)(-6 * u);
    gama[20] = gama[21] = gama[22] = (int8_t)(8 * u);
    gama[23] = gama[24] = (int8_t)(8 * u);
    gama[25] = gama[26] = gama[27] = gama[28] = (int8_t)(5 * u);
    gama[29] = gama[30] = gama[31] = gama[32] = (int8_t)(5 * u);
    gama[33] = gama[34] = gama[35] = gama[36] = (int8_t)(2 * u);
    gama[37] = gama[38] = gama[39] = gama[40] = (int8_t)(2 * u);
    gama[41] = gama[42] = gama[43] = gama[44] = (int8_t)(-1 * u);
    gama[45] = gama[46] = gama[47] = gama[48] = (int8_t)(-1 * u);
    gama[49] = (int8_t)(-1 * u);
    gama[50] = gama[51] = gama[52] = gama[53] = (int8_t)(-4 * u);
    gama[54] = gama[55] = gama[56] = gama[57] = (int8_t)(-4 * u);
    gama[58] = gama[59] = gama[60] = (int8_t)(-7 * u);
    gama[61] = gama[62] = gama[63] = (int8_t)(-7 * u);
  }
}

void bn_rec_tnaf_mod(bn_t r0, bn_t r1, const bn_t k, int u, int m) {
  bn_t t, t0, t1, t2, t3;

  TRY {
    bn_new(t);
    bn_new(t0);
    bn_new(t1);
    bn_new(t2);
    bn_new(t3);

    /* (a0, a1) = (1, 0). */
    bn_set_dig(t0, 1);
    bn_zero(t1);
    /* (b0, b1) = (0, 0). */
    bn_zero(t2);
    bn_zero(t3);
    /* (r0, r1) = (k, 0). */
    bn_abs(r0, k);
    bn_zero(r1);

    for (int i = 0; i < m; i++) {
      if (!bn_is_even(r0)) {
        /* r0 = r0 - 1. */
        bn_sub_dig(r0, r0, 1);
        /* (b0, b1) = (b0 + a0, b1 + a1). */
        bn_add(t2, t2, t0);
        bn_add(t3, t3, t1);
      }

      bn_hlv(t, r0);
      /* r0 = r1 + mu * r0 / 2. */
      if (u == -1) {
        bn_sub(r0, r1, t);
      } else {
        bn_add(r0, r1, t);
      }
      /* r1 = - r0 / 2. */
      bn_neg(r1, t);

      bn_dbl(t, t1);
      /* a1 = a0 + mu * a1. */
      if (u == -1) {
        bn_sub(t1, t0, t1);
      } else {
        bn_add(t1, t0, t1);
      }
      /* a0 = - 2 * a1. */
      bn_neg(t0, t);
    }

    /*r 0 = r0 + b0, r1 = r1 + b1. */
    bn_add(r0, r0, t2);
    bn_add(r1, r1, t3);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_tnaf(int8_t *tnaf, int *len, const bn_t k, int8_t u, int m, int w) {
  int i, l;
  bn_t tmp, r0, r1;
  int8_t beta[1 << (w - 2)], gama[1 << (w - 2)];
  uint8_t t_w;
  dig_t t0, t1, mask;
  int s, t, u_i;

  if (*len < (bn_bits(k) + 1)) { THROW(ERR_NO_BUFFER); }

  TRY {
    bn_new(r0);
    bn_new(r1);
    bn_new(tmp);

    bn_rec_tnaf_get(&t_w, beta, gama, u, w);
    bn_abs(tmp, k);
    bn_rec_tnaf_mod(r0, r1, tmp, u, m);

    mask = MASK(w);
    l = 1 << w;

    i = 0;
    while (!bn_is_zero(r0) || !bn_is_zero(r1)) {
      while ((r0->dp[0] & 1) == 0) {
        tnaf[i++] = 0;
        /* tmp = r0. */
        bn_hlv(tmp, r0);
        /* r0 = r1 + mu * r0 / 2. */
        if (u == -1) {
          bn_sub(r0, r1, tmp);
        } else {
          bn_add(r0, r1, tmp);
        }
        /* r1 = - r0 / 2. */
        bn_copy(r1, tmp);
        r1->sign = tmp->sign ^ 1;
      }
      /* If r0 is odd. */
      if (w == 2) {
        t0 = r0->dp[0];
        if (bn_sign(r0) == BN_NEG) { t0 = l - t0; }
        t1 = r1->dp[0];
        if (bn_sign(r1) == BN_NEG) { t1 = l - t1; }
        u_i = 2 - ((t0 - 2 * t1) & mask);
        tnaf[i++] = u_i;
        if (u_i < 0) {
          bn_add_dig(r0, r0, -u_i);
        } else {
          bn_sub_dig(r0, r0, u_i);
        }
      } else {
        /* t0 = r0 mod_s 2^w. */
        t0 = r0->dp[0];
        if (bn_sign(r0) == BN_NEG) { t0 = l - t0; }
        /* t1 = r1 mod_s 2^w. */
        t1 = r1->dp[0];
        if (bn_sign(r1) == BN_NEG) { t1 = l - t1; }
        /* u = r0 + r1 * (t_w) mod_s 2^w. */
        u_i = (t0 + t_w * t1) & mask;

        if (u_i >= (l / 2)) {
          /* If u < 0, s = -1 and u = -u. */
          u_i = (int8_t)(u_i - l);
          tnaf[i++] = u_i;
          u_i = (int8_t)(-u_i >> 1);
          t = -beta[u_i];
          s = -gama[u_i];
        } else {
          /* If u > 0, s = 1. */
          tnaf[i++] = u_i;
          u_i = (int8_t)(u_i >> 1);
          t = beta[u_i];
          s = gama[u_i];
        }
        /* r0 = r0 - s * beta_u. */
        if (t > 0) {
          bn_sub_dig(r0, r0, t);
        } else {
          bn_add_dig(r0, r0, -t);
        }
        /* r1 = r1 - s * gama_u. */
        if (s > 0) {
          bn_sub_dig(r1, r1, s);
        } else {
          bn_add_dig(r1, r1, -s);
        }
      }
      /* tmp = r0. */
      bn_hlv(tmp, r0);
      /* r0 = r1 + mu * r0 / 2. */
      if (u == -1) {
        bn_sub(r0, r1, tmp);
      } else {
        bn_add(r0, r1, tmp);
      }
      /* r1 = - r0 / 2. */
      bn_copy(r1, tmp);
      r1->sign = tmp->sign ^ 1;
    }
    *len = i;
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_rtnaf(int8_t *tnaf, int *len, const bn_t k, int8_t u, int m, int w) {
  int i, l;
  bn_t tmp, r0, r1;
  int8_t beta[1 << (w - 2)], gama[1 << (w - 2)];
  uint8_t t_w;
  dig_t t0, t1, mask;
  int s, t, u_i;

  if (*len < (bn_bits(k) + 1)) { THROW(ERR_NO_BUFFER); }

  TRY {
    bn_new(r0);
    bn_new(r1);
    bn_new(tmp);

    bn_rec_tnaf_get(&t_w, beta, gama, u, w);
    bn_abs(tmp, k);
    bn_rec_tnaf_mod(r0, r1, tmp, u, m);
    mask = MASK(w);
    l = CEIL(m + 2, (w - 1));

    i = 0;
    while (i < l) {
      /* If r0 is odd. */
      if (w == 2) {
        t0 = r0->dp[0];
        if (bn_sign(r0) == BN_NEG) { t0 = (1 << w) - t0; }
        t1 = r1->dp[0];
        if (bn_sign(r1) == BN_NEG) { t1 = (1 << w) - t1; }
        u_i = ((t0 - 2 * t1) & mask) - 2;
        tnaf[i++] = u_i;
        if (u_i < 0) {
          bn_add_dig(r0, r0, -u_i);
        } else {
          bn_sub_dig(r0, r0, u_i);
        }
      } else {
        /* t0 = r0 mod_s 2^w. */
        t0 = r0->dp[0];
        if (bn_sign(r0) == BN_NEG) { t0 = (1 << w) - t0; }
        /* t1 = r1 mod_s 2^w. */
        t1 = r1->dp[0];
        if (bn_sign(r1) == BN_NEG) { t1 = (1 << w) - t1; }
        /* u = r0 + r1 * (t_w) mod_s 2^w. */
        u_i = ((t0 + t_w * t1) & mask) - (1 << (w - 1));
        if (u_i < 0) {
          /* If u < 0, s = -1 and u = -u. */
          tnaf[i++] = u_i;
          u_i = (int8_t)(-u_i >> 1);
          t = -beta[u_i];
          s = -gama[u_i];
        } else {
          /* If u > 0, s = 1. */
          tnaf[i++] = u_i;
          u_i = (int8_t)(u_i >> 1);
          t = beta[u_i];
          s = gama[u_i];
        }
        /* r0 = r0 - s * beta_u. */
        if (t > 0) {
          bn_sub_dig(r0, r0, t);
        } else {
          bn_add_dig(r0, r0, -t);
        }
        /* r1 = r1 - s * gama_u. */
        if (s > 0) {
          bn_sub_dig(r1, r1, s);
        } else {
          bn_add_dig(r1, r1, -s);
        }
      }
      for (int j = 0; j < (w - 1); j++) {
        /* tmp = r0. */
        bn_hlv(tmp, r0);
        /* r0 = r1 + mu * r0 / 2. */
        if (u == -1) {
          bn_sub(r0, r1, tmp);
        } else {
          bn_add(r0, r1, tmp);
        }
        /* r1 = - r0 / 2. */
        bn_copy(r1, tmp);
        r1->sign = tmp->sign ^ 1;
      }
    }
    s = r0->dp[0];
    t = r1->dp[0];
    if (bn_sign(r0) == BN_NEG) { s = -s; }
    if (bn_sign(r1) == BN_NEG) { t = -t; }
    if (s != 0 && t != 0) {
      for (int j = 0; j < (1 << (w - 2)); j++) {
        if (beta[j] == s && gama[j] == t) {
          tnaf[i++] = 2 * j + 1;
          break;
        }
      }
      for (int j = 0; j < (1 << (w - 2)); j++) {
        if (beta[j] == -s && gama[j] == -t) {
          tnaf[i++] = -(2 * j + 1);
          break;
        }
      }
    } else {
      if (t != 0) {
        tnaf[i++] = t;
      } else {
        tnaf[i++] = s;
      }
    }
    *len = i;
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_reg(int8_t *naf, int *len, const bn_t k, int n, int w) {
  int i, l;
  bn_t t;
  dig_t t0, mask;
  int8_t u_i;

  mask = MASK(w);
  l = CEIL(n, (w - 1));

  if (*len < l) { THROW(ERR_NO_BUFFER); }

  TRY {
    bn_new(t);
    bn_abs(t, k);

    i = 0;
    if (w == 2) {
      for (i = 0; i < l; i++, naf++) {
        u_i = (t->dp[0] & mask) - 2;
        t->dp[0] -= u_i;
        *naf = u_i;
        bn_hlv(t, t);
      }
      bn_get_dig(&t0, t);
      *naf = t0;
    } else {
      for (i = 0; i < l; i++, naf++) {
        u_i = (t->dp[0] & mask) - (1 << (w - 1));
        t->dp[0] -= u_i;
        *naf = u_i;
        bn_rsh(t, t, w - 1);
      }
      bn_get_dig(&t0, t);
      *naf = t0;
    }
    *len = l + 1;
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_jsf(int8_t *jsf, int *len, const bn_t k, const bn_t l) {
  bn_t n0, n1;
  dig_t l0, l1;
  int8_t u0, u1, d0, d1;
  int i, j, offset;

  if (*len < (2 * bn_bits(k) + 1)) { THROW(ERR_NO_BUFFER); }

  TRY {
    bn_new(n0);
    bn_new(n1);

    bn_abs(n0, k);
    bn_abs(n1, l);

    i = bn_bits(k);
    j = bn_bits(l);
    offset = MAX(i, j) + 1;

    i = 0;
    d0 = d1 = 0;
    while (!(bn_is_zero(n0) && d0 == 0) || !(bn_is_zero(n1) && d1 == 0)) {
      bn_get_dig(&l0, n0);
      bn_get_dig(&l1, n1);
      /* For reduction modulo 8. */
      l0 = (l0 + d0) & MASK(3);
      l1 = (l1 + d1) & MASK(3);

      if (l0 % 2 == 0) {
        u0 = 0;
      } else {
        u0 = 2 - (l0 & MASK(2));
        if ((l0 == 3 || l0 == 5) && ((l1 & MASK(2)) == 2)) { u0 = (int8_t)-u0; }
      }
      jsf[i] = u0;
      if (l1 % 2 == 0) {
        u1 = 0;
      } else {
        u1 = 2 - (l1 & MASK(2));
        if ((l1 == 3 || l1 == 5) && ((l0 & MASK(2)) == 2)) { u1 = (int8_t)-u1; }
      }
      jsf[i + offset] = u1;

      if (d0 + d0 == 1 + u0) { d0 = (int8_t)(1 - d0); }
      if (d1 + d1 == 1 + u1) { d1 = (int8_t)(1 - d1); }

      i++;
      bn_hlv(n0, n0);
      bn_hlv(n1, n1);
    }
    *len = i;
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_rec_glv(bn_t k0, bn_t k1, const bn_t k, const bn_t n, const bn_t *v1, const bn_t *v2) {
  bn_t t, b1, b2;
  int r1, r2, bits;

  TRY {
    bn_new(b1);
    bn_new(b2);
    bn_new(t);

    bn_abs(t, k);
    bits = bn_bits(n);

    bn_mul(b1, t, v1[0]);
    r1 = bn_get_bit(b1, bits);
    bn_rsh(b1, b1, bits + 1);
    bn_add_dig(b1, b1, r1);

    bn_mul(b2, t, v2[0]);
    r2 = bn_get_bit(b2, bits);
    bn_rsh(b2, b2, bits + 1);
    bn_add_dig(b2, b2, r2);

    bn_mul(k0, b1, v1[1]);
    bn_mul(k1, b2, v2[1]);
    bn_add(k0, k0, k1);
    bn_sub(k0, t, k0);

    bn_mul(k1, b1, v1[2]);
    bn_mul(t, b2, v2[2]);
    bn_add(k1, k1, t);

    bn_neg(k1, k1);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}
/**
 * @file
 *
 * Implementation of the multiple precision arithmetic shift functions.
 *
 * @ingroup bn
 */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_dbl(bn_t c, const bn_t a) {
  dig_t carry;

  bn_grow(c, a->used + 1);

  c->used = a->used;
  carry = bn_lsh1_low(c->dp, a->dp, c->used);

  /* If there is an additional carry. */
  if (carry != 0) {
    c->dp[c->used] = carry;
    (c->used)++;
  }

  c->sign = a->sign;
}

void bn_hlv(bn_t c, const bn_t a) {
  bn_grow(c, a->used);

  c->used = a->used;
  bn_rsh1_low(c->dp, a->dp, c->used);

  c->sign = a->sign;
  bn_trim(c);
}

void bn_lsh(bn_t c, const bn_t a, int bits) {
  int digits;
  dig_t carry;

  bn_copy(c, a);
  if (bits <= 0) { return; }

  SPLIT(bits, digits, bits, BN_DIG_LOG);

  if (bits > 0) {
    if (bn_bits(c) + bits > c->used * (int)BN_DIGIT) { bn_grow(c, c->used + digits + 1); }
  } else {
    bn_grow(c, c->used + digits);
  }

  if (digits > 0) { bn_lshd_low(c->dp, a->dp, a->used, digits); }
  c->used = a->used + digits;
  c->sign = a->sign;

  if (bits > 0) {
    if (c != a) {
      carry = bn_lshb_low(c->dp + digits, a->dp, a->used, bits);
    } else {
      carry = bn_lshb_low(c->dp + digits, c->dp + digits, c->used - digits, bits);
    }
    if (carry != 0) {
      c->dp[c->used] = carry;
      (c->used)++;
    }
  }
  bn_trim(c);
}

void bn_rsh(bn_t c, const bn_t a, int bits) {
  int digits = 0;

  if (bits <= 0) {
    bn_copy(c, a);
    return;
  }

  SPLIT(bits, digits, bits, BN_DIG_LOG);

  if (digits > 0) { bn_rshd_low(c->dp, a->dp, a->used, digits); }
  c->used = a->used - digits;
  c->sign = a->sign;

  if (c->used > 0 && bits > 0) {
    if (digits == 0 && c != a) {
      bn_rshb_low(c->dp, a->dp + digits, a->used - digits, bits);
    } else {
      bn_rshb_low(c->dp, c->dp, c->used, bits);
    }
  }
  bn_trim(c);
}

/**
 * @file
 *
 * Implementation of the multiple precision arithmetic squaring
 * functions.
 *
 * @ingroup bn
 */

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Computes the square of a multiple precision integer using recursive Karatsuba
 * squaring.
 *
 * @param[out] c			- the result.
 * @param[in] a				- the multiple precision integer to square.
 * @param[in] level			- the number of Karatsuba steps to apply.
 */
static void bn_sqr_karat_imp(bn_t c, const bn_t a, int level) {
  int h;
  bn_t a0, a1, a0a0, a1a1, t;
  const dig_t *tmpa;
  dig_t *t0;

  /* Compute half the digits of a or b. */
  h = a->used >> 1;

  TRY {
    /* Allocate the temp variables. */
    bn_new(a0);
    bn_new(a1);
    bn_new(a0a0);
    bn_new(a1a1);
    bn_new(t);

    a0->used = h;
    a1->used = a->used - h;

    tmpa = a->dp;

    /* a = a1 || a0 */
    t0 = a0->dp;
    for (int i = 0; i < h; i++, t0++, tmpa++) *t0 = *tmpa;
    t0 = a1->dp;
    for (int i = 0; i < a1->used; i++, t0++, tmpa++) *t0 = *tmpa;

    bn_trim(a0);

    if (level <= 1) {
      /* a0a0 = a0 * a0 and a1a1 = a1 * a1 */
#if BN_SQR == BASIC
      bn_sqr_basic(a0a0, a0);
      bn_sqr_basic(a1a1, a1);
#elif BN_SQR == COMBA
      bn_sqr_comba(a0a0, a0);
      bn_sqr_comba(a1a1, a1);
#elif BN_SQR == MULTP
      bn_mul_comba(a0a0, a0, a0);
      bn_mul_comba(a1a1, a1, a1);
#endif
    } else {
      bn_sqr_karat_imp(a0a0, a0, level - 1);
      bn_sqr_karat_imp(a1a1, a1, level - 1);
    }

    /* t = (a1 + a0) */
    bn_add(t, a1, a0);

    if (level <= 1) {
      /* t = (a1 + a0)*(a1 + a0) */
#if BN_SQR == BASIC
      bn_sqr_basic(t, t);
#elif BN_SQR == COMBA
      bn_sqr_comba(t, t);
#elif BN_SQR == MULTP
      bn_mul_comba(t, t, t);
#endif
    } else {
      bn_sqr_karat_imp(t, t, level - 1);
    }

    /* t2 = (a0*a0 + a1*a1) */
    bn_add(a0, a0a0, a1a1);
    /* t = (a1 + a0)*(b1 + b0) - (a0*a0 + a1*a1) */
    bn_sub(t, t, a0);

    /* t = (a1 + a0)*(a1 + a0) - (a0*a0 + a1*a1) << h digits */
    bn_lsh(t, t, h * BN_DIGIT);

    /* t2 = a1 * b1 << 2*h digits */
    bn_lsh(a1a1, a1a1, 2 * h * BN_DIGIT);

    /* t = t + a0*a0 */
    bn_add(t, t, a0a0);
    /* c = t + a1*a1 */
    bn_add(t, t, a1a1);

    t->sign = BN_POS;
    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_sqr_basic(bn_t c, const bn_t a) {
  int i, digits;
  bn_t t;

  digits = 2 * a->used;

  TRY {
    bn_new_size(t, digits);
    bn_zero(t);
    t->used = digits;

    for (i = 0; i < a->used; i++) { bn_sqra_low(t->dp + (2 * i), a->dp + i, a->used - i); }

    t->sign = BN_POS;
    bn_trim(t);
    bn_copy(c, t);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_sqr_comba(bn_t c, const bn_t a) {
  int digits;
  digits = 2 * a->used;

  TRY {
    bn_new_size(c, digits);
    c->used = digits;

    // bn_sqrn_low(c->dp, a->dp, a->used);
    // mpn_mul_n(c->dp, a->dp, a->dp, a->used);
    mpn_sqr(c->dp, a->dp, a->used);

    c->sign = BN_POS;
    bn_trim(c);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_sqr_karat(bn_t c, const bn_t a) { bn_sqr_karat_imp(c, a, BN_KARAT); }

/**
 * @file
 *
 * Implementation of the multiple precision integer square root extraction.
 *
 * @ingroup bn
 */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_srt(bn_t c, bn_t a) {
  bn_t h, l, m, t;
  int bits, cmp;

  if (bn_sign(a) == BN_NEG) { THROW(ERR_NO_VALID); }

  bits = bn_bits(a);
  bits += (bits % 2);

  TRY {
    bn_new(h);
    bn_new(l);
    bn_new(m);
    bn_new(t);

    bn_set_2b(h, bits >> 1);
    bn_set_2b(l, (bits >> 1) - 1);

    /* Trivial binary search approach. */
    do {
      bn_add(m, h, l);
      bn_hlv(m, m);
      bn_sqr(t, m);
      cmp = bn_cmp(t, a);

      if (cmp == CMP_GT) {
        bn_copy(h, m);
      } else if (cmp == CMP_LT) {
        bn_copy(l, m);
      }
      bn_sub(t, h, l);
    } while (bn_cmp_dig(t, 1) == CMP_GT && cmp != CMP_EQ);

    bn_copy(c, m);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_copy(bn_t c, const bn_t a) {
  int i;

  if (c->dp == a->dp) return;

  bn_grow(c, a->used);

  for (i = 0; i < a->used; i++) { c->dp[i] = a->dp[i]; }

  c->used = a->used;
  c->sign = a->sign;
}

void bn_abs(bn_t c, const bn_t a) {
  if (c->dp != a->dp) { bn_copy(c, a); }
  c->sign = BN_POS;
}

void bn_neg(bn_t c, const bn_t a) {
  if (c->dp != a->dp) { bn_copy(c, a); }
  if (!bn_is_zero(c)) { c->sign = a->sign ^ 1; }
}

int bn_sign(const bn_t a) { return a->sign; }

void bn_clear(bn_t a) {
  a->sign = BN_POS;
  a->used = 1;
  for (int i = 0; i < a->alloc; i++) a->dp[i] = 0;
}

void bn_zero(bn_t a) {
  a->sign = BN_POS;
  a->used = 1;
#define ZERO_ALL
#ifdef ZERO_ALL
  // dv_zero(a->dp, a->alloc);
  for (int i = 0; i < a->alloc; i++) a->dp[i] = 0;
#else
  // Don't need to zero everything, do we?
  a->dp[0] = 0;
#endif
}

int bn_is_zero(const bn_t a) {
  if (a->used == 0) { return 1; }
  if ((a->used == 1) && (a->dp[0] == 0)) { return 1; }
  return 0;
}

int bn_is_even(const bn_t a) {
  if (bn_is_zero(a)) { return 1; }
  if ((a->dp[0] & 0x01) == 0) { return 1; }
  return 0;
}

int bn_bits(const bn_t a) {
  int bits;

  if (a->used == 0) { return 0; }

  /* Bits in lower digits. */
  bits = (a->used - 1) * BN_DIGIT;

  return bits + util_bits_dig(a->dp[a->used - 1]);
}

int bn_get_bit(const bn_t a, int bit) {
  int d;

  SPLIT(bit, d, bit, BN_DIG_LOG);

  if (d >= a->used) {
    return 0;
  } else {
    return (a->dp[d] >> bit) & (dig_t)1;
  }
}

void bn_set_bit(bn_t a, int bit, int value) {
  int d;

  SPLIT(bit, d, bit, BN_DIG_LOG);

  if (value == 1) {
    a->dp[d] |= ((dig_t)1 << bit);
    if ((d + 1) > a->used) { a->used = d + 1; }
  } else {
    a->dp[d] &= ~((dig_t)1 << bit);
    bn_trim(a);
  }
}

int bn_ham(const bn_t a) {
  int c = 0;

  for (int i = 0; i < bn_bits(a); i++) { c += bn_get_bit(a, i); }

  return c;
}

void bn_get_dig(dig_t *c, const bn_t a) { *c = a->dp[0]; }

void bn_set_dig(bn_t a, dig_t digit) {
  bn_zero(a);
  a->dp[0] = digit;
  a->used = 1;
  a->sign = BN_POS;
}
void bn_set_neg_dig(bn_t a, dig_t digit) {
  bn_zero(a);
  a->dp[0] = digit;
  a->used = 1;
  a->sign = BN_NEG;
}

void bn_set_2b(bn_t a, int b) {
  int i, d;

  SPLIT(b, d, b, BN_DIG_LOG);

  bn_grow(a, d + 1);
  for (i = 0; i < d; i++) a->dp[i] = 0;
  a->used = d + 1;
  a->dp[d] = ((dig_t)1 << b);
  a->sign = BN_POS;
}

void bn_print(const bn_t a) {
  int i;

  if (a->sign == BN_NEG) { util_print("-"); }
  if (a->used == 0) {
    util_print("0\n");
  } else {
#if WORD == 64
    util_print_dig(a->dp[a->used - 1], 0);
    for (i = a->used - 2; i >= 0; i--) { util_print_dig(a->dp[i], 1); }
#else
    util_print_dig(a->dp[a->used - 1], 0);
    for (i = a->used - 2; i >= 0; i--) { util_print_dig(a->dp[i], 1); }
#endif
    util_print("\n");
  }
}

int bn_size_str(const bn_t a, int radix) {
  int digits = 0;
  bn_t t;

  /* Binary case requires the bits, a sign and the null terminator. */
  if (radix == 2) { return bn_bits(a) + (a->sign == BN_NEG ? 1 : 0) + 1; }

  /* Check the radix. */
  if (radix < 2 || radix > 64) { THROW(ERR_NO_VALID); }

  if (bn_is_zero(a)) { return 2; }

  if (a->sign == BN_NEG) { digits++; }

  TRY {
    bn_new(t);
    bn_copy(t, a);

    t->sign = BN_POS;

    while (!bn_is_zero(t)) {
      bn_div_dig(t, t, (dig_t)radix);
      digits++;
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }

  return digits + 1;
}

void bn_read_str(bn_t a, const char *str, int len, int radix) {
  int sign, i, j;
  char c;

  bn_zero(a);

  if (radix < 2 || radix > 64) { THROW(ERR_NO_VALID) }

  j = 0;
  if (str[0] == '-') {
    j++;
    sign = BN_NEG;
  } else {
    sign = BN_POS;
  }

  while (str[j] && j < len) {
    c = (char)((radix < 36) ? TOUPPER(str[j]) : str[j]);
    for (i = 0; i < 64; i++) {
      if (c == util_conv_char(i)) { break; }
    }

    if (i < radix) {
      bn_mul_dig(a, a, (dig_t)radix);
      bn_add_dig(a, a, (dig_t)i);
    } else {
      break;
    }
    j++;
  }

  a->sign = sign;
}

void bn_write_str(char *str, int len, const bn_t a, int radix) {
  bn_t t;
  dig_t d;
  int digits, l, i, j;
  char c;

  l = bn_size_str(a, radix);
  if (len < l) { THROW(ERR_NO_BUFFER); }

  if (radix < 2 || radix > 64) { THROW(ERR_NO_VALID) }

  if (bn_is_zero(a) == 1) {
    *str++ = '0';
    *str = '\0';
    return;
  }

  TRY {
    bn_new(t);
    bn_copy(t, a);

    j = 0;
    if (t->sign == BN_NEG) {
      str[j] = '-';
      j++;
      t->sign = BN_POS;
    }

    digits = 0;
    while (!bn_is_zero(t)) {
      bn_div_rem_dig(t, &d, t, (dig_t)radix);
      str[j] = util_conv_char(d);
      digits++;
      j++;
    }

    /* Reverse the digits of the string. */
    i = 0;
    if (str[0] == '-') { i = 1; }

    j = l - 2;
    while (i < j) {
      c = str[i];
      str[i] = str[j];
      str[j] = c;
      ++i;
      --j;
    }

    str[l - 1] = '\0';
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

int bn_size_bin(const bn_t a) {
  dig_t d;
  int digits;

  digits = (a->used - 1) * (BN_DIGIT / 8);
  d = a->dp[a->used - 1];

  while (d != 0) {
    d = d >> 8;
    digits++;
  }
  return digits;
}

void bn_read_bin(bn_t a, const uint8_t *bin, int len) {
  int i, j;
  dig_t d = (BN_DIGIT / 8);
  int digs = (len % d == 0 ? len / d : len / d + 1);

  bn_grow(a, digs);
  bn_zero(a);
  a->used = digs;

  for (i = 0; i < digs - 1; i++) {
    d = 0;
    for (j = (BN_DIGIT / 8) - 1; j >= 0; j--) {
      d = d << 8;
      d |= bin[len - 1 - (i * (BN_DIGIT / 8) + j)];
    }
    a->dp[i] = d;
  }
  d = 0;
  for (j = (BN_DIGIT / 8) - 1; j >= 0; j--) {
    if ((int)(i * (BN_DIGIT / 8) + j) < len) {
      d = d << 8;
      d |= bin[len - 1 - (i * (BN_DIGIT / 8) + j)];
    }
  }
  a->dp[i] = d;

  a->sign = BN_POS;
  bn_trim(a);
}

void bn_write_bin(uint8_t *bin, int len, const bn_t a) {
  int size, k;
  dig_t d;

  size = bn_size_bin(a);

  if (len < size) { THROW(ERR_NO_BUFFER); }

  k = 0;
  for (int i = 0; i < a->used - 1; i++) {
    d = a->dp[i];
    for (int j = 0; j < (int)(BN_DIGIT / 8); j++) {
      bin[len - 1 - k++] = d & 0xFF;
      d = d >> 8;
    }
  }

  d = a->dp[a->used - 1];
  while (d != 0) {
    bin[len - 1 - k++] = d & 0xFF;
    d = d >> 8;
  }

  while (k < len) { bin[len - 1 - k++] = 0; }
}

int bn_size_raw(const bn_t a) { return a->used; }

void bn_read_raw(bn_t a, const dig_t *raw, int len) {
  TRY {
    bn_grow(a, len);
    a->used = len;
    a->sign = BN_POS;
    for (int i = 0; i < len; i++) { a->dp[i] = *raw++; }
    bn_trim(a);
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

void bn_write_raw(dig_t *raw, int len, const bn_t a) {
  int i, size;

  size = a->used;

  if (len < size) { THROW(ERR_NO_BUFFER); }

  for (i = 0; i < size; i++) { raw[i] = a->dp[i]; }
  for (; i < len; i++) { raw[i] = 0; }
}



/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/
static bn_t x;
static bn_t y;
static bn_t r;
static int init = 0;
static bn_t q;

/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[out] d		- the remainder.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_imp(bn_t c, bn_t d, const bn_t a, const bn_t b) {
  int sign;
  // bn_t q;

  if (!init) {
    bn_new(x);
    bn_new(y);
    bn_new(r);
    bn_new(q);
    init = 1;
  }

  /* If |a| < |b|, we're done. */
  if (bn_cmp_abs(a, b) == CMP_LT) {
    if (bn_sign(a) == bn_sign(b)) {
      if (c != NULL) { bn_zero(c); }
      if (d != NULL) { bn_copy(d, a); }
    } else {
      if (c != NULL) {
        bn_set_dig(c, 1);
        bn_neg(c, c);
      }
      if (d != NULL) { bn_add(d, a, b); }
    }
    return;
  }

  TRY {
    bn_new_size(q, a->used + 1);
    bn_zero(q);
    bn_zero(r);
    bn_abs(x, a);
    bn_abs(y, b);

    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);

    bn_divn_low(q->dp, r->dp, x->dp, a->used, y->dp, b->used);
    // mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);

    /* We have the quotient in q and the remainder in r. */
    if (c != NULL) {
      q->used = a->used - b->used + 1;
      q->sign = sign;
      bn_trim(q);
      if ((sign == BN_POS) || bn_is_zero(r)) {
        bn_copy(c, q);
      } else {
        bn_sub_dig(c, q, 1);
      }
    }

    if (d != NULL) {
      r->used = b->used;
      r->sign = b->sign;
      bn_trim(r);
      if ((sign == BN_POS) || bn_is_zero(r)) {
        bn_copy(d, r);
      } else {
        bn_sub(d, b, r);
      }
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}
/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_imp1(bn_t c, const bn_t a, const bn_t b) {
  int sign;

  if (!init) {
    bn_new(x);
    bn_new(y);
    bn_new(r);
    bn_new(q);
    init = 1;
  }

  /* If |a| < |b|, we're done. */
  if (bn_cmp_abs(a, b) == CMP_LT) {
    if (bn_sign(a) == bn_sign(b)) {
      bn_zero(c);
    } else {
      bn_set_neg_dig(c, 1);
    }
    return;
  }

  TRY {
    bn_new_size(q, a->used + 1);
    bn_zero(q);
    bn_zero(r);
    bn_abs(x, a);
    bn_abs(y, b);

    /* Find the sign. */
    sign = (a->sign == b->sign ? BN_POS : BN_NEG);
    bn_divn_low(q->dp, r->dp, x->dp, a->used, y->dp, b->used);
    // mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);

    /* We have the quotient in q and the remainder in r. */
    q->used = a->used - b->used + 1;
    q->sign = sign;
    bn_trim(q);
    if ((sign == BN_POS) || bn_is_zero(r)) {
      bn_copy(c, q);
    } else {
      bn_sub_dig(c, q, 1);
    }
  }
  CATCH_ANY { THROW(ERR_CAUGHT); }
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_div(bn_t c, const bn_t a, const bn_t b) {
  if (bn_is_zero(b)) { THROW(ERR_NO_VALID); }
  bn_div_imp1(c, a, b);
}

void bn_div_rem(bn_t c, bn_t d, const bn_t a, const bn_t b) {
  if (bn_is_zero(b)) { THROW(ERR_NO_VALID); }
  bn_div_imp(c, d, a, b);
}
