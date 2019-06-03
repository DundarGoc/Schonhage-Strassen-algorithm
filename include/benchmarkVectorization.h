#pragma once

#include <assert.h>
#include <flint/nmod_poly.h>
#include <sys/time.h>
#include <time.h>

double GetTime(void);
void PrintPolynomial(ulong *poly, ulong length);
void PrintCharNTimes(char character, ulong numberOfTimes);
void PrintResults(double data[12][13]);
void BenchmarkButterfly(unsigned int modBits, ulong maximumPossibleLength, flint_rand_t state);
