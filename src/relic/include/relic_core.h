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
 * @defgroup relic Core functions
 */

/**
 * @file
 *
 * Interface of the library core functions.
 *
 * @ingroup relic
 */

#ifndef RELIC_CORE_H
#define RELIC_CORE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "relic_err.h"
#include "relic_bn.h"
#include "relic_conf.h"

#if MULTI != RELIC_NONE
#include <math.h>

#if MULTI == OPENMP
#include <omp.h>
#elif MULTI == PTHREAD
#include <pthread.h>
#endif /* OPENMP */

#endif /* MULTI != RELIC_NONE */

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/**
 * Indicates that the function executed correctly.
 */
#define STS_OK			0

/**
 * Indicates that an error occurred during the function execution.
 */
#define STS_ERR			1

/**
 * Indicates that a comparison returned that the first argument was lesser than
 * the second argument.
 */
#define CMP_LT			-1

/**
 * Indicates that a comparison returned that the first argument was equal to
 * the second argument.
 */
#define CMP_EQ			0

/**
 * Indicates that a comparison returned that the first argument was greater than
 * the second argument.
 */
#define CMP_GT			1

/**
 * Indicates that two incomparable elements are not equal.
 */
#define CMP_NE			2

/**
 * Optimization identifer for the case where a coefficient is 0.
 */
#define OPT_ZERO		0

/**
 * Optimization identifer for the case where a coefficient is 1.
 */
#define OPT_ONE			1

/**
 * Optimization identifer for the case where a coefficient is 1.
 */
#define OPT_TWO			2

/**
 * Optimization identifer for the case where a coefficient is small.
 */
#define OPT_DIGIT		3

/**
 * Optimization identifier for the case where a coefficient is -3.
 */
#define OPT_MINUS3		4

/**
 * Optimization identifier for the case where the coefficient is random
 */
#define RELIC_OPT_NONE		5

/**
 * Maximum number of terms to describe a sparse object.
 */
#define MAX_TERMS		16

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/**
 * Library context.
 */
typedef struct _ctx_t {
	/** The value returned by the last call, can be STS_OK or STS_ERR. */
	int code;

#ifdef CHECK
	/** The state of the last error caught. */
	sts_t *last;
	/** Error state to be used outside try-catch blocks. */
	sts_t error;
	/** Error number to be used outside try-catch blocks. */
	err_t number;
	/** The error message respective to the last error. */
	char *reason[ERR_MAX];
	/** A flag to indicate if the last error was already caught. */
	int caught;
#endif /* CHECK */

#if BENCH > 0
	/** Stores the time measured before the execution of the benchmark. */
	bench_t before;
	/** Stores the time measured after the execution of the benchmark. */
	bench_t after;
	/** Stores the sum of timings for the current benchmark. */
	long long total;
#ifdef OVERH
	/** Benchmarking overhead to be measured and subtracted from benchmarks. */
	long long over;
#endif
#endif

	/** Flag to indicate if PRNG is seed. */
	int seeded;
	/** Counter to keep track of number of calls since last seeding. */
	int counter;
} ctx_t;

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initializes the library.
 *
 * @return STS_OK if no error occurs, STS_ERR otherwise.
 */
int core_init(void);

/**
 * Finalizes the library.
 *
 * @return STS_OK if no error occurs, STS_ERR otherwise.
 */
int core_clean(void);

/**
 * Returns a pointer to the current library context.
 *
 * @return a pointer to the library context.
 */
ctx_t *core_get(void);

/**
 * Switched the library context to a new context.
 *
 * @param[in] ctx					- the new library context.
 */
void core_set(ctx_t *ctx);

#endif /* !RELIC_CORE_H */
