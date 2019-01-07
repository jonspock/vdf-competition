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
 * Implementation of the library basic functions.
 *
 * @ingroup relic
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "relic_core.h"
#include "relic_types.h"
#include "relic_err.h"
#include "relic_arch.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

/**
 * If multi-threading is enabled, assigns each thread a local copy of the data.
 */
#if MULTI == PTHREAD
#define thread 	__thread
#else
#define thread /* */
#endif

/**
 * Default library context.
 */
thread ctx_t first_ctx;

/**
 * Active library context.
 */
thread ctx_t *core_ctx = NULL;

#if MULTI == OPENMP
#pragma omp threadprivate(first_ctx, core_ctx)
#endif

int core_init(void) {
	if (core_ctx == NULL) {
		core_ctx = &(first_ctx);
	}

#if defined(CHECK) && defined(TRACE)
	core_ctx->trace = 0;
#endif

#ifdef CHECK
	core_ctx->reason[ERR_NO_MEMORY] = MSG_NO_MEMORY;
	core_ctx->reason[ERR_NO_PRECI] = MSG_NO_PRECI;
	core_ctx->reason[ERR_NO_FILE] = MSG_NO_FILE;
	core_ctx->reason[ERR_NO_READ] = MSG_NO_READ;
	core_ctx->reason[ERR_NO_VALID] = MSG_NO_VALID;
	core_ctx->reason[ERR_NO_BUFFER] = MSG_NO_BUFFER;
	core_ctx->reason[ERR_NO_FIELD] = MSG_NO_FIELD;
	core_ctx->reason[ERR_NO_CURVE] = MSG_NO_CURVE;
	core_ctx->reason[ERR_NO_CONFIG] = MSG_NO_CONFIG;
	core_ctx->last = NULL;
#endif /* CHECK */

#ifdef OVERH
	core_ctx->over = 0;
#endif

	core_ctx->code = STS_OK;
	return STS_OK;
}

int core_clean(void) {
	core_ctx = NULL;
	return STS_OK;
}

ctx_t *core_get(void) {
	return core_ctx;
}

void core_set(ctx_t *ctx) {
	core_ctx = ctx;
}
