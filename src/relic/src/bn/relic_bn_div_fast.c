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
 * Implementation of the multiple precision division functions.
 *
 * @ingroup bn
 */

#include "relic_core.h"
#include "relic_bn_low.h"

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
void bn_div_exact(bn_t c,
                     bn_t rr, 
                     bn_t a, bn_t b) {
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
    //mpn_divexact(c->dp, a->dp, a->used, b->dp, b->used);
    /* We have the quotient in q -> ignore r */
    c->used = a->used - b->used + 1;
    c->sign = sign;
    bn_trim(c);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
}

/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_impf(bn_t c,
                       bn_t rr, 
                       bn_t a, bn_t b) {
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
    //mpn_tdiv_qr(c->dp, rr->dp, 0, a->dp, a->used, b->dp, b->used);
        
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
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
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
static void bn_div_impf_with_temps(bn_t c, bn_t d,
                                  bn_t xx, bn_t yy, bn_t rr, bn_t qq,
                                  const bn_t a, const bn_t b) {
	int sign;

	/* If |a| < |b|, we're done. */
	if (bn_cmp_abs(a, b) == CMP_LT) {
		if (bn_sign(a) == bn_sign(b)) {
			if (c != NULL) {
				bn_zero(c);
			}
			if (d != NULL) {
				bn_copy(d, a);
			}
		} else {
			if (c != NULL) {
				bn_set_dig(c, 1);
				bn_neg(c, c);
			}
			if (d != NULL) {
				bn_add(d, a, b);
			}
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
    //mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);
        
        
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
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
}
/**
 * Divides two multiple precision integers, computing the quotient and the
 * remainder.
 *
 * @param[out] c		- the quotient.
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_impf1_const(bn_t c,
                              bn_t xx, bn_t yy, bn_t rr, bn_t qq,
                              const bn_t a, const bn_t b) {
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
    //mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);
        
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
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_div_with_temps(bn_t c,
            bn_t xx, bn_t yy, bn_t rr, bn_t qq,
            const bn_t a, const bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_impf1_const(c,
                    xx,yy,qq,rr,
                    a, b);
}

void bn_div_fast(bn_t c,
                 // temporary
                 bn_t rr,
                 // avoid copy since not const and can be changed
                 bn_t a, bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_impf(c, rr, a, b);
}

void bn_div_rem_with_temps(bn_t c, bn_t d,
                           // temporaries
                           bn_t xx, bn_t yy, bn_t rr, bn_t qq,
                           const bn_t a, const bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_impf_with_temps(c, d,
                        xx,yy,rr,qq,
                        a, b);
}


