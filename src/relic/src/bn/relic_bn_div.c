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
 * @param[in] a			- the dividend.
 * @param[in] b			- the the divisor.
 */
static void bn_div_imp2(bn_t c,
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
static void bn_div_imp(bn_t c, bn_t d, const bn_t a, const bn_t b) {
	int sign;
    //bn_t q;

  if (!init) {
    bn_null(x);
    bn_null(y);
    bn_null(r);
    bn_null(q);
    bn_new(x);
    bn_new(y);
    bn_new(r);
    bn_new(q);
    init = 1;
  }

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
		bn_new_size(q, a->used + 1);
		bn_zero(q);
		bn_zero(r);
		bn_abs(x, a);
		bn_abs(y, b);

		/* Find the sign. */
		sign = (a->sign == b->sign ? BN_POS : BN_NEG);

    bn_divn_low(q->dp, r->dp, x->dp, a->used, y->dp, b->used);
    //mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);
        
        
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
static void bn_div_imp1(bn_t c, const bn_t a, const bn_t b) {
	int sign;

  if (!init) {
    bn_null(x);
    bn_null(y);
    bn_null(r);
    bn_null(q);
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
    //mpn_tdiv_qr(q->dp, r->dp, 0, x->dp, a->used, y->dp, b->used);
        
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
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void bn_div(bn_t c, const bn_t a, const bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_imp1(c, a, b);
}

void bn_div2(bn_t c,
                // temporary
                bn_t rr,
                // avoid copy since not const and can be changed
                bn_t a, bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_imp2(c, rr, a, b);
}

void bn_div_rem(bn_t c, bn_t d, const bn_t a, const bn_t b) {
	if (bn_is_zero(b)) {
		THROW(ERR_NO_VALID);
	}
	bn_div_imp(c, d, a, b);
}


void bn_div_dig(bn_t c, const bn_t a, dig_t b) {
	bn_t q;
	dig_t r;

	bn_null(q);

	if (b == 0) {
		THROW(ERR_NO_VALID);
	}

	if (b == 1 || bn_is_zero(a) == 1) {
		if (c != NULL) {
			bn_copy(c, a);
		}
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
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(q);
	}
}

void bn_div_rem_dig(bn_t c, dig_t *d, const bn_t a, dig_t b) {
	bn_t q;
	dig_t r;

	bn_null(q);

	if (b == 0) {
		THROW(ERR_NO_VALID);
	}

	if (b == 1 || bn_is_zero(a) == 1) {
		if (d != NULL) {
			*d = 0;
		}
		if (c != NULL) {
			bn_copy(c, a);
		}
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

		if (d != NULL) {
			*d = r;
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(q);
	}
}

