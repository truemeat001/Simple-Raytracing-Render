#pragma once
#ifndef MATHF_H
#define MATHF_H
#include <math.h>

#define __m 0x100000000LL
#define __c 0xB16
#define __a 0x5DEECE66DLL

#define M_PI 3.14159265358979323846

static unsigned long long seed = 1;

double drand48(void)
{
	seed = (__a * seed + __c) & 0xFFFFFFFFFFFFLL;
	unsigned int x = seed >> 16;
	return ((double)x / (double)__m);
}

void srand48(unsigned int i)
{
	seed = (((long long int)i << 16) | rand());
}
#endif // !MATHF_H
