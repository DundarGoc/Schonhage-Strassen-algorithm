#pragma once

#include <assert.h>
#include <flint/nmod_poly.h>
#include <sys/time.h>
#include <time.h>


mp_limb_t nmodAdd(mp_limb_t a, mp_limb_t b, nmod_t mod);
mp_limb_t nmodSub(mp_limb_t a, mp_limb_t b, nmod_t mod);
void ButterflyInPlace(ulong *op1, ulong *op2, ulong n, nmod_t modFLINT);
