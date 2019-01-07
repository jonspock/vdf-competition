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
 * Implementation of the multiple precision integer memory management routines.
 *
 * @ingroup bn
 */

#include <errno.h>
#include "relic_core.h"

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
	if (digits > BN_SIZE) {
		THROW(ERR_NO_PRECI)
	}
	(void)a;
}
#endif

void bn_trim(bn_t a) {
	while (a->used > 0 && a->dp[a->used - 1] == 0) {
		--(a->used);
	}
	/* Zero can't be negative. */
	if (a->used <= 0) {
		a->used = 1;
		a->dp[0] = 0;
		a->sign = BN_POS;
	}
}
