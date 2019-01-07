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
 * @defgroup arch Architecture-dependent utilities
 */

/**
 * @file
 *
 * Interface of architecture-dependent functions.
 *
 * @ingroup arch
 */

#ifndef RELIC_ARCH_H
#define RELIC_ARCH_H

#include "relic_types.h"

#if ARCH == AVR
#include <avr/pgmspace.h>
#endif

/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/**
 * Chooses a proper way to store a string in the target architecture.
 *
 * @param[in] STR		- the string to store.
 */
#if ARCH == AVR
#define STRING(STR)			PSTR(STR)
#else
#define STRING(STR)			STR
#endif

/**
 * Fetches a constant string to be used by the library.
 *
 * @param[out] STR		- the resulting prepared parameter.
 * @param[in] ID		- the parameter represented as a string.
 * @param[in] L			- the length of the string.
 */
#if ARCH == AVR
#define FETCH(STR, ID, L)	arch_copy_rom(STR, STRING(ID), L);
#else
#define FETCH(STR, ID, L)	memcpy(STR, ID, L);
#endif

#endif /* !RELIC_ARCH_H */
