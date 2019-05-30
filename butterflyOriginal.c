#include "include/butterfly.h"


ulong nmodAdd(ulong a, ulong b, nmod_t mod)
{
	ulong neg = mod.n - a;
	return (neg > b) ? a + b : b - neg;
}

ulong nmodSub(ulong a, ulong b, nmod_t mod)
{
	ulong diff = a - b;
	return (a < b) ? mod.n + diff : diff;
}


void ButterflyInPlace(ulong *op1, ulong *op2, ulong n, nmod_t modFLINT)
{
	for(; n; --n, ++op1, ++op2)
	{
		ulong x = *op1;
		ulong y = *op2;
		*op1 = nmodAdd(y, x, modFLINT);
		*op2 = nmodSub(y, x, modFLINT);
	}
}
