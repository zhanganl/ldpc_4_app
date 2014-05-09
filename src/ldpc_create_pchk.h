#ifndef LDPC_CREATE_PCHK__
#define LDPC_CREATE_PCHK__

#include <stdbool.h>
#include "ldpc_matrix_sparse.h"

typedef enum make_method_enum
{
	Evencol, 	/* Uniform number of bits per column, with number specified */
	Evenboth 	/* Uniform (as possible) over both columns and rows */
} make_method; 


typedef enum SessionType_enum
{
	TypeLDGM,
	TypeSTAIRS,
	TypeTRIANGLE
} SessionType;

static unsigned long seed;

/**
 * Initialize the PRNG with a seed between 1 and 0x7FFFFFFE
 * (2^^31-2) inclusive.
 */
static inline void ldpc_srand (unsigned long s)
{
	if ((s >= 1) && (s <= 0x7FFFFFFE))
		seed = s;
	else 
		exit(-1);
}

/**
 * Returns a random integer between 0 and maxv-1 inclusive.
 *	16807		multiplier constant (7^^5)
 *	0x7FFFFFFF	modulo constant (2^^31-1)
 * The inner PRNG produces a value between 1 and 0x7FFFFFFE
 * (2^^31-2) inclusive.
 * This value is then scaled between 0 and maxv-1 inclusive.
 */
static inline unsigned long
ldpc_rand (unsigned long	maxv)
{
	unsigned long	hi, lo;
	lo = 16807 * (seed & 0xFFFF);
	hi = 16807 * (seed >> 16);
	lo += (hi & 0x7FFF) << 16;
	lo += hi >> 15;
	if (lo > 0x7FFFFFFF)
		lo -= 0x7FFFFFFF;
	seed = (long) lo;
	return ((unsigned long)
			((double)seed * (double)maxv / (double)0x7FFFFFFF));
}

mod2sparse* CreatePchkMatrix (int nbRows, int nbCols, make_method makeMethod, int leftDegree, int seed, bool no4cycle, SessionType type);

#endif

