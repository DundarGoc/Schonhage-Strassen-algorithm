#include "include/SSA.h"

/*
   A rough estimate of the L1 cache size in bytes. This is to take advantage of the cache by switching from the recursive Truncated Fourier Transform to an iterative version once the input size
      becomes small enough. David Harvey, the person writing zn_poly which my program is based on, states that if this is a bit on the small side then it's probably not a big deal. If it's on the big
      side, that might start to seriously degrade performance. I haven't tested this yet but David's a smart fellow so it's probably true.
 */
#define CACHE_SIZE_ESTIMATE 32768


/*
   My recreation of the zn_poly library but with all irrelevant parts removed. Many functions have also been replaced with similar (often identical) functions from FLINT.
 */


/*
   For integers a >= 1 and b >= 1, returns ceil(a / b).
 */
#define CEIL_DIV(a, b) ((((a) - 1) / (b)) + 1)

void PMFVectorFFTIterative(ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	assert(lgK <= lgM + 1);
	assert((t << lgK) < (2UL << lgM));

	ulong pmfLength = 1UL << lgM;
	ulong numberOfCells = 1UL << lgK;

	if(lgK == 0)
	{
		return;
	}

	ptrdiff_t half = skip << (lgK - 1);
	ulong *end = pmfVector + skip * numberOfCells;

	for(ulong r = pmfLength >> (lgK - 1); r <= pmfLength; r <<= 1, half >>= 1, t <<= 1)
	{
		for(ulong *start = pmfVector, s = t; s < pmfLength; s += r, start += skip)
		{
			for(ulong *p = start; p < end; p += 2 * half)
			{
				PMFButterfly(p, p + half, pmfLength, modFLINT);
				p[half] += pmfLength + s;
			}
		}
	}
}
void PMFVectorTFTDivideAndConquer(ulong numberOfOutputCoefficients, ulong z, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	assert(lgK <= lgM + 1);
	assert((t << lgK) < (2UL << lgM));
	assert(numberOfOutputCoefficients >= 1 && numberOfOutputCoefficients <= (1UL << lgK));
	assert(z >= 1 && z <= (1UL << lgK));

	ulong numberOfCells = 1UL << lgK;
	ulong pmfLength = 1UL << lgM;

	if(numberOfCells == 1)
	{
		return;
	}

	if(numberOfOutputCoefficients == numberOfCells && z == numberOfCells)
	{
		PMFVectorFFTIterative(t, lgK, lgM, skip, pmfVector, modFLINT);
		return;
	}

	--lgK;
	numberOfCells = 1UL << lgK;
	ulong numberOfColumns = numberOfCells;
	ulong *p = pmfVector;
	ptrdiff_t half = skip * numberOfCells;
	ulong z2 = FLINT_MIN(z, numberOfColumns);
	long i = 0;

	if(numberOfOutputCoefficients <= numberOfColumns)
	{
		for(; i < (long)(z - numberOfColumns); ++i, p += skip)
		{
			PMFAdd(p, p + half, pmfLength, modFLINT);
		}

		PMFVectorTFTDivideAndConquer(numberOfOutputCoefficients, z2, t << 1, lgK, lgM, skip, pmfVector, modFLINT);
	}
	else
	{
		ulong s = t;
		ulong r = pmfLength >> lgK;

		for(; i < (long)(z - numberOfColumns); ++i, p += skip, s += r)
		{
			PMFButterfly(p, p + half, pmfLength, modFLINT);
			p[half] += pmfLength + s;
		}

		for(; i < (long)z2; ++i, p += skip, s += r)
		{
			_nmod_vec_set(p + half, p, pmfLength + 1);
			p[half] += s;
		}

		PMFVectorTFTDivideAndConquer(numberOfColumns, z2, t << 1, lgK, lgM, skip, pmfVector, modFLINT);

		pmfVector += half;
		PMFVectorTFTDivideAndConquer(numberOfOutputCoefficients - numberOfColumns, z2, t << 1, lgK, lgM, skip, pmfVector, modFLINT);
	}
}
void PMFVectorTFTHuge(ulong numberOfOutputCoefficients, ulong numberOfInputCoefficients, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	ulong lgT = lgK / 2;
	ulong pmfLength = 1UL << lgM;

	assert(lgK <= lgM + 1);
	assert(t << lgK < 2 * pmfLength);
	assert(lgT > 0 && lgT < lgK);
	assert(numberOfOutputCoefficients >= 1 && numberOfOutputCoefficients <= (1UL << lgK));
	assert(numberOfInputCoefficients >= 1 && numberOfInputCoefficients <= (1UL << lgK));

	ulong numberOfRows = 1UL << lgT;
	ulong lgU = lgK - lgT;
	ulong numberOfColumns = 1UL << lgU;
	ptrdiff_t skip_U = skip << lgU;
	ulong *pmfVectorOriginal = pmfVector;
	ulong numberOfOutputCoefficientsOnLastRow = numberOfOutputCoefficients & (numberOfColumns - 1);
	ulong numberOfCompleteOutputRows = numberOfOutputCoefficients >> lgU;
	ulong numberOfOutputRows = numberOfCompleteOutputRows + (numberOfOutputCoefficientsOnLastRow > 0);
	ulong numberOfCompleteInputRows = numberOfInputCoefficients >> lgU;
	ulong numberOfInputCoefficientsOnLastRow = numberOfInputCoefficients & (numberOfColumns - 1);
	ulong numberOfInputColumns = numberOfCompleteInputRows ? numberOfColumns : numberOfInputCoefficientsOnLastRow;
	ulong r = pmfLength  >> (lgK - 1);
	ulong s = t;

	// Transform the first numberOfInputCoefficientsOnLastRow columns.
	for(ulong iColumn = 0; iColumn < numberOfInputCoefficientsOnLastRow; ++iColumn, pmfVector += skip, s += r)
	{
		PMFVectorTFT(numberOfOutputRows, numberOfCompleteInputRows + 1, s, lgT, lgM, skip_U, pmfVector, modFLINT);
	}

	// Transform the rest of the columns.
	for(ulong iColumn = numberOfInputCoefficientsOnLastRow; iColumn < numberOfInputColumns; ++iColumn, pmfVector += skip, s += r)
	{
		PMFVectorTFT(numberOfOutputRows, numberOfCompleteInputRows, s, lgT, lgM, skip_U, pmfVector, modFLINT);
	}

	pmfVector = pmfVectorOriginal;
	t *= numberOfRows;


	// Transform first numberOfCompleteOutputRows rows.
	for(ulong iRow = 0; iRow < numberOfCompleteOutputRows; ++iRow, pmfVector += skip_U)
	{
		PMFVectorTFT(numberOfColumns, numberOfInputColumns, t, lgU, lgM, skip, pmfVector, modFLINT);
	}

	// Transform last partial row if there is one.
	if(numberOfOutputCoefficientsOnLastRow)
	{
		PMFVectorTFT(numberOfOutputCoefficientsOnLastRow, numberOfInputColumns, t, lgU, lgM, skip, pmfVector, modFLINT);
	}
}
void PMFVectorTFT(ulong numberOfOutputCoefficients, ulong numberOfInputCoefficients, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	ulong pmfVectorLength = 1UL << lgK;
	ulong pmfLength = 1UL << lgM;

	assert(lgK <= lgM + 1);
	assert(t * pmfVectorLength < 2 * pmfLength);
	assert(numberOfOutputCoefficients >= 1 && numberOfOutputCoefficients <= pmfVectorLength);
	assert(numberOfInputCoefficients >= 1 && numberOfInputCoefficients <= pmfVectorLength);

	if(pmfVectorLength <= 2 || 2 * pmfVectorLength * pmfLength * sizeof(ulong) <= CACHE_SIZE_ESTIMATE)
	{
		PMFVectorTFTDivideAndConquer(numberOfOutputCoefficients, numberOfInputCoefficients, t, lgK, lgM, skip, pmfVector, modFLINT);
	}
	else
	{
		PMFVectorTFTHuge(numberOfOutputCoefficients, numberOfInputCoefficients, t, lgK, lgM, skip, pmfVector, modFLINT);
	}
}

void PMFVectorIFFTIterative(ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	assert(lgK <= lgM + 1);
	assert(t * (1UL << lgK) < (2UL << lgM));

	ulong pmfLength = 1UL << lgM;

	if(lgK == 0)
	{
		return;
	}

	ulong r = pmfLength;
	ulong r_last = pmfLength >> (lgK - 1);
	t <<= (lgK - 1);
	ptrdiff_t half = skip;
	ulong *end = pmfVector + (skip << lgK);

	for(; r >= r_last; r >>= 1, half <<= 1, t >>= 1)
	{
		for(ulong *start = pmfVector, s = t; s < pmfLength; s += r, start += skip)
		{
			for(ulong *p = start; p < end; p += 2 * half)
			{
				p[half] += pmfLength - s;
				PMFButterfly(p + half, p, pmfLength, modFLINT);
			}
		}
	}
}
void PMFVectorITFTDivideAndConquer(ulong n, int fwd, ulong z, ulong t, ulong lgM, ulong lgK, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	assert(lgK <= lgM + 1);
	assert((t << lgK) < (2UL << lgM));
	assert(z >= 1 && z <= (1UL << lgK));
	assert(n + fwd >= 1 && n + fwd <= (1UL << lgK));
	assert(n <= z);

	ulong pmfVectorLength = 1UL << lgK;
	ulong pmfLength = 1UL << lgM;

	if(pmfVectorLength == 1)
	{
		return;
	}

	if(n == pmfVectorLength)
	{
		PMFVectorIFFTIterative(t, lgK, lgM, skip, pmfVector, modFLINT);
		return;
	}

	pmfVectorLength >>= 1;
	--lgK;

	ulong U = pmfVectorLength;
	ptrdiff_t half = skip << lgK;

	if(n + fwd <= U)
	{
		long zU2 = FLINT_MIN(z, U);
		long last_zero_fwd_bfly = FLINT_MAX(z - zU2, n);
		long last_zero_cross_bfly = FLINT_MIN(z - zU2, n);
		long i = zU2 - 1;
		ulong *p = pmfVector + skip * i;

		for(; i >= last_zero_fwd_bfly; --i, p -= skip)
		{
			_nmod_vec_scalar_mul_nmod(p + 1, p + 1, pmfLength, nmod_inv(2, modFLINT), modFLINT);
		}

		for(; i >= (long)n; --i, p -= skip)
		{
			PMFAdd(p, p + half, pmfLength, modFLINT);
			_nmod_vec_scalar_mul_nmod(p + 1, p + 1, pmfLength, nmod_inv(2, modFLINT), modFLINT);
		}

		PMFVectorITFTDivideAndConquer(n, fwd, zU2, t << 1, lgM, lgK, skip, pmfVector, modFLINT);

		for(; i >= last_zero_cross_bfly; --i, p -= skip)
		{
			PMFAdd(p, p, pmfLength, modFLINT);
		}

		for(; i >= 0; --i, p -= skip)
		{
			PMFAdd(p, p, pmfLength, modFLINT);
			PMFSub(p, p + half, pmfLength, modFLINT);
		}
	}
	else
	{
		PMFVectorIFFTIterative(t << 1, lgK, lgM, skip, pmfVector, modFLINT);

		long i = U - 1;
		ulong r = pmfLength >> lgK;
		ulong s = t + r * i;
		ulong *p = pmfVector + skip * i;
		long last_zero_cross_bfly = z - U;
		long last_cross_bfly = n - U;

		for(; i >= last_zero_cross_bfly; --i, s -= r, p -= skip)
		{
			_nmod_vec_set(p + half, p, pmfLength + 1);
			p[half] += s;
			PMFAdd(p, p, pmfLength, modFLINT);
		}

		for(; i >= last_cross_bfly; --i, s -= r, p -= skip)
		{
			PMFSub(p + half, p, pmfLength, modFLINT);
			PMFSub(p, p + half, pmfLength, modFLINT);
			p[half] += pmfLength + s;
		}

		pmfVector += half;
		PMFVectorITFTDivideAndConquer(n - U, fwd, U, t << 1, lgM, lgK, skip, pmfVector, modFLINT);

		for(; i >= 0; --i, s -= r, p -= skip)
		{
			p[half] += pmfLength - s;
			PMFButterfly(p + half, p, pmfLength, modFLINT);
		}
	}
}
void PMFVectorITFTHuge(ulong lgT, ulong n, int fwd, ulong z, ulong t, nmod_t modFLINT, ulong lgK, ptrdiff_t skip, ulong lgM, ulong *pmfVector)
{
	ulong pmfLength = 1UL << lgM;

	assert(lgK <= lgM + 1);
	assert(t * (1UL << lgK) < 2 * pmfLength);
	assert(z >= 1 && z <= (1UL << lgK));
	assert(n + fwd >= 1 && n + fwd <= (1UL << lgK));
	assert(n <= z);
	assert(lgT > 0 && lgT < lgK);

	ulong lgU = lgK - lgT;
	ulong U = 1UL << lgU;
	ptrdiff_t skip_U = skip << lgU;
	ulong *pmfVectorOriginal = pmfVector;
	ulong nU = n & (U - 1);
	ulong nT = n >> lgU;
	ulong zU = z & (U - 1);
	ulong zT = z >> lgU;
	ulong zU2 = zT ? U : zU;
	ulong mU1 = FLINT_MIN(zU, nU);
	ulong mU2 = FLINT_MAX(zU, nU);
	int fwd2 = nU || fwd;
	ulong r = pmfLength >> (lgK - 1);
	ulong s, i;
	ulong tT = t << lgT;

	for(i = 0; i < nT; ++i, pmfVector += skip_U)
	{
		PMFVectorITFT(U, 0, U, tT, lgM, lgU, skip, pmfVector, modFLINT);
	}

	for(i = nU, pmfVector = pmfVectorOriginal + (skip * nU), s = t + (r * nU); i < mU2; ++i, pmfVector += skip, s += r)
	{
		PMFVectorITFT(nT, fwd2, zT + 1, s, lgM, lgT, skip_U, pmfVector, modFLINT);
	}

	for(; i < zU2; ++i, pmfVector += skip, s += r)
	{
		PMFVectorITFT(nT, fwd2, zT, s, lgM, lgT, skip_U, pmfVector, modFLINT);
	}

	if(fwd2)
	{
		pmfVector = pmfVectorOriginal + nT * skip_U;
		PMFVectorITFT(nU, fwd, zU2, tT, lgM, lgU, skip, pmfVector, modFLINT);

		for(i = 0, pmfVector = pmfVectorOriginal, s = t; i < mU1; ++i, pmfVector += skip, s += r)
		{
			PMFVectorITFT(nT + 1, 0, zT + 1, s, lgM, lgT, skip_U, pmfVector, modFLINT);
		}

		for(; i < nU; ++i, pmfVector += skip, s += r)
		{
			PMFVectorITFT(nT + 1, 0, zT, s, lgM, lgT, skip_U, pmfVector, modFLINT);
		}
	}
}
void PMFVectorITFT(ulong n, int fwd, ulong z, ulong t, ulong lgM, ulong lgK, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT)
{
	ulong pmfLength = 1UL << lgM;
	ulong pmfVectorLength = 1UL << lgK;

	assert(lgK <= lgM + 1);
	assert(t * pmfVectorLength < 2 * pmfLength);
	assert(z <= pmfVectorLength);
	assert(n <= z);
	assert(n + fwd <= pmfVectorLength);

	if(pmfVectorLength <= 2 || 2 * pmfVectorLength * pmfLength * sizeof(ulong) <= CACHE_SIZE_ESTIMATE)
	{
		PMFVectorITFTDivideAndConquer(n, fwd, z, t, lgM, lgK, skip, pmfVector, modFLINT);
	}
	else
	{
		PMFVectorITFTHuge(lgK / 2, n, fwd, z, t, modFLINT, lgK, skip, lgM, pmfVector);
	}
}

void ButterflyInPlace(ulong *op1, ulong *op2, ulong n, nmod_t modFLINT)
{
	if(modFLINT.n == FLINT_BITS)
	{
		for(; n; --n, ++op1, ++op2)
		{
			ulong x = *op1;
			ulong y = *op2;
			*op1 = nmod_add(y, x, modFLINT);
			*op2 = nmod_sub(y, x, modFLINT);
		}
	}
	else
	{
		for(; n; --n, ++op1, ++op2)
		{
			ulong x = *op1;
			ulong y = *op2;
			*op1 = _nmod_add(y, x, modFLINT);
			*op2 = _nmod_sub(y, x, modFLINT);
		}

	}
}
void PMFButterfly(ulong *op1, ulong *op2, ulong pmfLength, nmod_t modFLINT)
{
	ulong relativeBias = op2[0] - op1[0];

	if(relativeBias & pmfLength)
	{
		relativeBias &= (pmfLength - 1);
		ButterflyInPlace(op1 + 1, op2 + 1 + pmfLength - relativeBias, relativeBias, modFLINT);
		ButterflyInPlace(op2 + 1, op1 + 1 + relativeBias, pmfLength - relativeBias, modFLINT);
	}
	else
	{
		relativeBias &= (pmfLength - 1);
		ButterflyInPlace(op2 + 1 + pmfLength - relativeBias, op1 + 1, relativeBias, modFLINT);
		ButterflyInPlace(op1 + 1 + relativeBias, op2 + 1, pmfLength - relativeBias, modFLINT);
	}
}
void PMFAdd(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT)
{
	ulong relativeBias = op2[0] - op1[0];

	if(relativeBias & pmfLength)
	{
		relativeBias &= (pmfLength - 1);
		_nmod_vec_add(op1 + 1, op1 + 1, op2 + 1 + pmfLength - relativeBias, relativeBias, modFLINT);
		_nmod_vec_sub(op1 + 1 + relativeBias, op1 + 1 + relativeBias, op2 + 1, pmfLength - relativeBias, modFLINT);
	}
	else
	{
		relativeBias &= (pmfLength - 1);
		_nmod_vec_sub(op1 + 1, op1 + 1, op2 + 1 + pmfLength - relativeBias, relativeBias, modFLINT);
		_nmod_vec_add(op1 + 1 + relativeBias, op1 + 1 + relativeBias, op2 + 1, pmfLength - relativeBias, modFLINT);
	}
}
void PMFSub(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT)
{
	ulong relativeBias = op2[0] - op1[0];

	if(relativeBias & pmfLength)
	{
		relativeBias &= (pmfLength - 1);
		_nmod_vec_sub(op1 + 1, op1 + 1, op2 + 1 + pmfLength - relativeBias, relativeBias, modFLINT);
		_nmod_vec_add(op1 + 1 + relativeBias, op1 + 1 + relativeBias, op2 + 1, pmfLength - relativeBias, modFLINT);
	}
	else
	{
		relativeBias &= (pmfLength - 1);
		_nmod_vec_add(op1 + 1, op1 + 1, op2 + 1 + pmfLength - relativeBias, relativeBias, modFLINT);
		_nmod_vec_sub(op1 + 1 + relativeBias, op1 + 1 + relativeBias, op2 + 1, pmfLength - relativeBias, modFLINT);
	}
}
ulong GetFudgeFactorFromPointwiseMultiplication(ulong lgM, nmod_t modFLINT, ulong wordMod)
{
	ulong pmfLength = 1UL << lgM;
	return GetFudgeFactorFromMultiplication(pmfLength, pmfLength, modFLINT, wordMod);
}
void PMFVectorPointwiseMultiplication(ulong n, nmod_t modFLINT, ulong lgM, ulong *data1, ulong *data2, ptrdiff_t skip, ulong wordMod)
{
	assert(modFLINT.n & 1);

	ulong pmfLength = 1UL << lgM;
	const ulong *p1 = data1;
	const ulong *p2 = data2;
	ulong *p3 = data1;

	ulong i = 0;
	ulong *temp = (ulong *)flint_malloc(sizeof(ulong) * 2 * pmfLength);
	ulong fudge2 = GetFudgeFactorFromMultiplication(pmfLength, pmfLength, modFLINT, wordMod);
	ulong fudge1 = GetFudgeFactorFromMultiplication(pmfLength / 2, pmfLength / 2, modFLINT, wordMod);
	ulong fudge = nmod_div(fudge1, fudge2, modFLINT);

	for(; i < 2 && i < n; ++i, p3 += skip, p1 += skip, p2 += skip)
	{
		p3[0] = p1[0] + p2[0];
		MultiplicationInternal(temp, p1 + 1, pmfLength / 2, p2 + 1, pmfLength / 2, 1, modFLINT, wordMod);
		_nmod_vec_scalar_mul_nmod(p3 + 1, temp, pmfLength - 1, fudge, modFLINT);
		p3[pmfLength] = 0;
	}


	temp[2 * pmfLength - 1] = 0;

	for(; i < n; ++i, p3 += skip, p1 += skip, p2 += skip)
	{
		p3[0] = p1[0] + p2[0];

		MultiplicationInternal(temp, p1 + 1, pmfLength, p2 + 1, pmfLength, 1, modFLINT, wordMod);
		_nmod_vec_sub(p3 + 1, temp, temp + pmfLength, pmfLength, modFLINT);
	}


	flint_free(temp);
}

void TransformPolynomialToPMFVector(const ulong *poly, size_t polyLength, ulong x, ulong pmfLength, nmod_t modFLINT, ulong *pmfVector)
{
	for(; polyLength >= pmfLength / 2; polyLength -= pmfLength / 2, poly += pmfLength / 2, pmfVector += pmfLength)
	{
		*pmfVector++ = 0;
		_nmod_vec_scalar_mul_nmod(pmfVector, poly, pmfLength / 2, x, modFLINT);
		_nmod_vec_zero(pmfVector + pmfLength / 2, pmfLength / 2);
	}

	if(polyLength)
	{
		*pmfVector++ = 0;
		_nmod_vec_scalar_mul_nmod(pmfVector, poly, polyLength, x, modFLINT);
		_nmod_vec_zero(pmfVector + polyLength, pmfLength - polyLength);
	}
}
void NegateOrCopy(ulong *res, const ulong *op, ulong n, int negate, nmod_t modFLINT)
{
	if(negate)
	{
		_nmod_vec_neg(res, op, n, modFLINT);
	}
	else
	{
		assert(op != NULL);
		_nmod_vec_set(res, op, n);
	}
}
void SchonhageStrassenCombineChunk(ulong *res, size_t n, ulong *op1, ulong *op2, ulong pmfLength, nmod_t modFLINT)
{
	n = FLINT_MIN(n, pmfLength / 2);

	if(op1 == NULL && op2 == NULL)
	{
		_nmod_vec_zero(res, n);
		return;
	}

	ulong s1 = ULONG_MAX;
	int neg1 = -1;

	if(op1)
	{
		s1 = (pmfLength / 2 - op1[0]) & (2 * pmfLength - 1);
		neg1 = (s1 >= pmfLength);

		if(neg1)
		{
			s1 -= pmfLength;
		}
	}

	ulong s2 = ULONG_MAX;
	int neg2 = -1;

	if(op2)
	{
		s2 = (-op2[0]) & (2 * pmfLength - 1);
		neg2 = (s2 >= pmfLength);

		if(neg2)
		{
			s2 -= pmfLength;
		}
	}

	if(s1 > s2)
	{
		MP_PTR_SWAP(op1, op2);
		ulong s_temp = s1;
		s1 = s2;
		s2 = s_temp;
		int neg_temp = neg1;
		neg1 = neg2;
		neg2 = neg_temp;
	}

	++op1;
	++op2;

	if(s2 == ULONG_MAX)
	{
		if(n <= pmfLength - s1)
		{
			NegateOrCopy(res, op1 + s1, n, neg1, modFLINT);
		}
		else
		{
			NegateOrCopy(res, op1 + s1, pmfLength - s1, neg1, modFLINT);
			NegateOrCopy(res + pmfLength - s1, op1, n - pmfLength + s1, !neg1, modFLINT);
		}

		return;
	}

	if(n <= pmfLength - s2)
	{
		SkipSignedAdd(res, 1, n, op2 + s2, neg2, op1 + s1, neg1, modFLINT);
		return;
	}

	res = SkipSignedAdd(res, 1, pmfLength - s2, op2 + s2, neg2, op1 + s1, neg1, modFLINT);
	n -= (pmfLength - s2);

	if(n <= s2 - s1)
	{
		SkipSignedAdd(res, 1, n, op2, !neg2, op1 + s1 + pmfLength - s2, neg1, modFLINT);
		return;
	}

	res = SkipSignedAdd(res, 1, s2 - s1, op2, !neg2, op1 + s1 + pmfLength - s2, neg1, modFLINT);
	n -= (s2 - s1);

	SkipSignedAdd(res, 1, (n >= s1) ? s1 : n, op2 + s2 - s1, !neg2, op1, !neg1, modFLINT);
}
void TransformPMFVectorToPolynomial(ulong *poly, ulong polyLength, ulong m3, ulong pmfLength, nmod_t modFLINT, ptrdiff_t skip, ulong *pmfVector)
{
	if(m3 == 0)
	{
		_nmod_vec_zero(poly, polyLength);
		return;
	}

	ulong pmfVectorLength = FLINT_MIN(polyLength, pmfLength / 2);
	SchonhageStrassenCombineChunk(poly, pmfVectorLength, NULL, pmfVector, pmfLength, modFLINT);
	poly += pmfVectorLength;
	polyLength -= pmfVectorLength;

	ulong *ptr1 = pmfVector;
	ulong *ptr2 = pmfVector + skip;
	ulong i;

	for(i = 1; i < m3 && polyLength >= pmfLength / 2; ++i, polyLength -= pmfLength / 2, poly += pmfLength / 2, ptr1 += skip, ptr2 += skip)
	{
		SchonhageStrassenCombineChunk(poly, polyLength, ptr1, ptr2, pmfLength, modFLINT);
	}

	if(i < m3)
	{
		SchonhageStrassenCombineChunk(poly, polyLength, ptr1, ptr2, pmfLength, modFLINT);
		return;
	}

	SchonhageStrassenCombineChunk(poly, polyLength, ptr1, NULL, pmfLength, modFLINT);

	if(polyLength > pmfLength / 2)
	{
		_nmod_vec_zero(poly + pmfLength / 2, polyLength - pmfLength / 2);
	}
}
ulong GetFudgeFactorFromTFT(size_t poly1Length, size_t poly2Length, nmod_t modFLINT, ulong wordMod)
{
	ulong lgM, m1, m2, m3, pmfLength;

	for(lgM = 1;; ++lgM)
	{
		pmfLength = 1UL << lgM;
		m1 = CEIL_DIV(poly1Length, pmfLength / 2);
		m2 = CEIL_DIV(poly2Length, pmfLength / 2);
		m3 = m1 + m2 - 1;

		if(m3 <= 2 * pmfLength)
		{
			break;
		}
	}

	ulong lgK = (m3 > pmfLength) ? (lgM + 1) : lgM;
	ulong pmfVectorLength = 1UL << lgK;
	ulong fudgeFromTFT = nmod_inv(pmfVectorLength, modFLINT);
	ulong fudgeFromPointwiseMultiplication = GetFudgeFactorFromPointwiseMultiplication(lgM, modFLINT, wordMod);

	return nmod_mul(fudgeFromTFT, fudgeFromPointwiseMultiplication, modFLINT);
}
void MultiplicationSchonhageStrassen(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, ulong x, nmod_t modFLINT, ulong wordMod)
{
	assert(modFLINT.n & 1);
	assert(poly2Length >= 1);
	assert(poly1Length >= poly2Length);

	ulong lgM, m1, m2, m3, pmfLength;
	ulong polyResultLength = poly1Length + poly2Length - 1;

	for(lgM = 1;; ++lgM)
	{
		pmfLength = 1UL << lgM;
		m1 = CEIL_DIV(poly1Length, pmfLength / 2);
		m2 = CEIL_DIV(poly2Length, pmfLength / 2);
		m3 = m1 + m2 - 1;

		if(m3 <= 2 * pmfLength)
		{
			break;
		}
	}

	ulong lgK = (m3 > pmfLength) ? (lgM + 1) : lgM;
	ulong pmfVectorLength = 1UL << lgK;
	ptrdiff_t pmfLengthWithBias = pmfLength + 1;
	ulong *pmfVector1 = (ulong *)flint_malloc(sizeof(ulong) * pmfLengthWithBias * pmfVectorLength);
	ulong *pmfVector2 = (ulong *)flint_malloc(sizeof(ulong) * pmfLengthWithBias * pmfVectorLength);

	TransformPolynomialToPMFVector(poly1, poly1Length, 1, pmfLength, modFLINT, pmfVector1);
	TransformPolynomialToPMFVector(poly2, poly2Length, x, pmfLength, modFLINT, pmfVector2);

	PMFVectorTFT(m3, m1, 0, lgK, lgM, pmfLengthWithBias, pmfVector1, modFLINT);
	PMFVectorTFT(m3, m2, 0, lgK, lgM, pmfLengthWithBias, pmfVector2, modFLINT);

	PMFVectorPointwiseMultiplication(m3, modFLINT, lgM, pmfVector1, pmfVector2, pmfLengthWithBias, wordMod);

	PMFVectorITFT(m3, 0, m3, 0, lgM, lgK, pmfLengthWithBias, pmfVector1, modFLINT);

	TransformPMFVectorToPolynomial(polyResult, polyResultLength, m3, pmfLength, modFLINT, pmfLengthWithBias, pmfVector1);

	flint_free(pmfVector1);
	flint_free(pmfVector2);
}

ulong *SkipSignedAdd(ulong *vectorResult, ptrdiff_t skipSize, size_t n, const ulong *vector1, int negateVector1, const ulong *vector2, int negateVector2, nmod_t modFLINT)
{
	if(negateVector1)
	{
		if(negateVector2)
		{
			for(; n > 0; --n, vectorResult += skipSize, ++vector1, ++vector2)
			{
				ulong x = nmod_add(*vector1, *vector2, modFLINT);
				*vectorResult = nmod_neg(x, modFLINT);
			}
		}
		else
		{
			for(; n > 0; --n, vectorResult += skipSize, ++vector1, ++vector2)
			{
				*vectorResult = nmod_sub(*vector2, *vector1, modFLINT);
			}
		}
	}
	else
	{
		if(negateVector2)
		{
			for(; n > 0; --n, vectorResult += skipSize, ++vector1, ++vector2)
			{
				*vectorResult = nmod_sub(*vector1, *vector2, modFLINT);
			}
		}
		else
		{
			for(; n > 0; --n, vectorResult += skipSize, ++vector1, ++vector2)
			{
				*vectorResult = nmod_add(*vector1, *vector2, modFLINT);
			}
		}
	}

	return vectorResult;
}

ulong GetFudgeFactorFromMultiplication(ulong poly1Length, size_t poly2Length, nmod_t modFLINT, ulong wordMod)
{
	assert(poly2Length >= 1);
	assert(poly1Length >= poly2Length);

	if(!(modFLINT.n & 1))
	{
		return 1;
	}


	ulong modBits = FLINT_BIT_COUNT(modFLINT.n);
	int ks1Used = poly2Length < tuningModulo[modBits].thresholdKS2;
	int ks2Used = poly2Length < tuningModulo[modBits].thresholdKS4;
	int ks4Used = poly2Length < tuningModulo[modBits].thresholdFFT;
	int ksUsed = ks1Used || ks2Used || ks4Used;
	if(ksUsed)
	{
		return nmod_neg(wordMod, modFLINT);
	}



	return GetFudgeFactorFromTFT(poly1Length, poly2Length, modFLINT, wordMod);
}
void MultiplicationInternal(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, int useREDCIfOddModulo, nmod_t modFLINT, ulong wordMod)
{
	assert(poly2Length >= 1);
	assert(poly1Length >= poly2Length);

	int moduloIsOdd = (modFLINT.n & 1);
	int useREDC = useREDCIfOddModulo && moduloIsOdd;

	ulong modBits = FLINT_BIT_COUNT(modFLINT.n);

	if(poly2Length < tuningModulo[modBits].thresholdKS2)
	{
		_nmod_poly_mul_KS(polyResult, poly1, poly1Length, poly2, poly2Length, 0, modFLINT);
	}
	else if(poly2Length < tuningModulo[modBits].thresholdKS4)
	{
		_nmod_poly_mul_KS2(polyResult, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(!moduloIsOdd || poly2Length < tuningModulo[modBits].thresholdFFT)
	{
		_nmod_poly_mul_KS4(polyResult, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else
	{
		ulong x = useREDCIfOddModulo ? 1 : GetFudgeFactorFromTFT(poly1Length, poly2Length, modFLINT, wordMod);
		MultiplicationSchonhageStrassen(polyResult, poly1, poly1Length, poly2, poly2Length, x, modFLINT, wordMod);
	}

	if(useREDC)
	{
		ulong polyResultLength = poly1Length + poly2Length - 1;
		ulong wordModInverse = nmod_inv(wordMod, modFLINT);
		_nmod_vec_scalar_mul_nmod(polyResult, polyResult, polyResultLength, wordModInverse, modFLINT);
		_nmod_vec_neg(polyResult, polyResult, polyResultLength, modFLINT);
	}
}
void SSA(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, nmod_t modFLINT)
{
	assert(modFLINT.n > 1);

	// Not sure if it's the correct term but word refers to word size i.e. 2^FLINT_BITS, where FLINT_BITS=64 or 32 bits depending system.
	ulong wordMod = nmod_pow_ui(2, FLINT_BITS, modFLINT);

	MultiplicationInternal(polyResult, poly1, poly1Length, poly2, poly2Length, 0, modFLINT, wordMod);
}
