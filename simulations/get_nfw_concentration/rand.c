#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

/* Minimal random number generator of Park and Miller. Returns a uniform random
 * deviate between 0.0 and 1.0. Set or reset seed to any integer value (except 
 * the unlikely value MASK) to initialize the sequence; Adapted from Numerical 
 * Recipies section 7.1 */

float ran0(long *seed)
{
	long k;
	float ans;

	*seed ^= MASK;
	k=(*seed)/IQ;
	*seed = IA*(*seed-k*IQ)-IR*k;
	if (*seed < 0) *seed += IM;
	ans = AM*(*seed);
	*seed ^= MASK;
	return ans;
}

