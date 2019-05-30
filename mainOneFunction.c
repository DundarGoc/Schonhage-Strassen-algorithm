#include <flint/nmod_poly.h>
#include <sys/time.h>
#include <time.h>
#include <zn_poly/zn_poly.h>

/*
   Colored output.
 */
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define NUMBER_OF_ROWS_RESULT 12
#define NUMBER_OF_COLUMNS_RESULT 13
#define LOWER_LIMIT 0.8
#define UPPER_LIMIT 1.2

/*
   Choose multiplication function to use in Makefile as the first argument.
 #0 - MultiplicationSSA
 #1 - MultiplicationFLINT
 #2 - zn_poly
 */

int chooseFunction = 0;
double timeElapsed;
ulong tuneMod[65][4];
double measurementTimeInSeconds = 0.01;

/*
   Simple custom timer.
 */
double GetTime(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);

	return t.tv_sec + t.tv_usec * 1e-6;
}


//#############################################################################################################################
//Declaration cemetary
void PMFButterfly(ulong *op1, ulong *op2, ulong pmfLength, nmod_t modFLINT);
void MultiplicationInternal(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, int useREDCIfOddModulo, nmod_t modFLINT, ulong wordMod);
void PMFAdd(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT);
void PMFVectorTFT(ulong numberOfOutputCoefficients, ulong numberOfInputCoefficients, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFSub(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT);
void PMFVectorITFT(ulong n, int fwd, ulong z, ulong t, ulong lgM, ulong lgK, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
ulong *SkipSignedAdd(ulong *vectorResult, ptrdiff_t skipSize, size_t n, const ulong *vector1, int negateVector1, const ulong *vector2, int negateVector2, nmod_t modFLINT);
ulong GetFudgeFactorFromMultiplication(ulong poly1Length, size_t poly2Length, nmod_t modFLINT, ulong wordMod);



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
	assert(lgK <= lgM + 1);
	assert(t * (1UL << lgK) < 2 * (1UL << lgM));
	assert(numberOfOutputCoefficients >= 1 && numberOfOutputCoefficients <= (1UL << lgK));
	assert(numberOfInputCoefficients >= 1 && numberOfInputCoefficients <= (1UL << lgK));


	PMFVectorTFTDivideAndConquer(numberOfOutputCoefficients, numberOfInputCoefficients, t, lgK, lgM, skip, pmfVector, modFLINT);

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
		PMFVectorITFTDivideAndConquer(1UL << lgK, 0, 1UL << lgK, t << 1, lgM, lgK, skip, pmfVector, modFLINT);

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
		pmfVector -= half;

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
	assert(lgK <= lgM + 1);
	assert(t * (1UL << lgK) < 2 * (1UL << lgM));
	assert(z <= (1UL << lgK));
	assert(n <= z);
	assert(n + fwd <= (1UL << lgK));


	PMFVectorITFTDivideAndConquer(n, fwd, z, t, lgM, lgK, skip, pmfVector, modFLINT);
}


void NussbaumerSplit(const ulong *op, nmod_t modFLINT, ulong lgK, ulong lgM, ulong *pmfVector, ptrdiff_t skip)
{
	assert(lgK >= 2);
	assert(lgM + 1 >= lgK);

	ulong pmfLength = 1UL << lgM;
	ulong pmfVectorLength = 1UL << lgK;
	ulong *dest = pmfVector + 1;
	ptrdiff_t half = skip << (lgK - 2);
	ulong r = pmfLength >> (lgK - 1);

	for(ulong j = 0, s = 0; j < pmfVectorLength / 4; ++j, dest += skip, s += r)
	{
		const ulong *src = op + j;

		dest[-1] = 0;
		dest[-1 + half] = 2 * s;
		dest[-1 + 2 * half] = s;
		dest[-1 + 3 * half] = 3 * s;

		for(ulong i = 0; i < pmfLength / 2; ++i, src += pmfVectorLength / 2)
		{
			ulong x0 = src[0];
			ulong x1 = src[pmfVectorLength / 4];
			ulong x2 = src[pmfLength * pmfVectorLength / 4];
			ulong x3 = src[pmfLength * pmfVectorLength / 4 + pmfVectorLength / 4];

			dest[i] = nmod_add(x0, x1, modFLINT);
			dest[i + half] = nmod_sub(x0, x1, modFLINT);
			dest[i + 2 * half] = nmod_sub(x0, x3, modFLINT);
			dest[i + 3 * half] = nmod_add(x0, x3, modFLINT);
			dest[i + pmfLength / 2] = nmod_add(x2, x3, modFLINT);
			dest[i + half + pmfLength / 2] = nmod_sub(x2, x3, modFLINT);
			dest[i + 2 * half + pmfLength / 2] = nmod_add(x2, x1, modFLINT);
			dest[i + 3 * half + pmfLength / 2] = nmod_sub(x2, x1, modFLINT);
		}
	}
}
void NussbaumerFFT(nmod_t modFLINT, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector)
{
	assert(lgK >= 2);
	assert(lgM + 1 >= lgK);

	if(lgK == 2)
	{
		return;
	}

	ulong pmfLength = 1UL << lgM;
	ulong s, r = pmfLength >> (lgK - 3);
	ptrdiff_t half = skip << (lgK - 3);
	ulong *end = pmfVector + (skip << lgK);
	ulong *start;
	ulong *p;

	for(; r <= pmfLength; r <<= 1, half >>= 1)
	{
		for(s = 0, start = pmfVector; s < pmfLength; s += r, start += skip)
		{
			for(p = start; p < end; p += 2 * half)
			{
				PMFButterfly(p, p + half, pmfLength, modFLINT);
				p[half] += pmfLength + s;
			}
		}
	}
}
void NussbaumerIFFT(ulong pmfLength, nmod_t modFLINT, ulong lgK, ulong *pmfVector, ptrdiff_t skip)
{
	ulong s;
	ulong r = pmfLength;
	ulong r_last = pmfLength >> (lgK - 1);
	ptrdiff_t half = skip;
	ulong *end = pmfVector + (skip << lgK);
	ulong *p;
	ulong *start;

	for(; r >= r_last; r >>= 1, half <<= 1)
	{
		for(s = 0, start = pmfVector; s < pmfLength; s += r, start += skip)
		{
			for(p = start; p < end; p += 2 * half)
			{
				p[half] += pmfLength - s;
				PMFButterfly(p + half, p, pmfLength, modFLINT);
			}
		}
	}
}
void NussbaumerCombine(ulong *res, ulong pmfLength, nmod_t modFLINT, ulong *pmfVector, ptrdiff_t skip, ulong pmfVectorLength)
{
	ulong *src1 = pmfVector + 1;
	ulong *src2 = pmfVector + skip * pmfVectorLength / 2 + 1;

	for(ulong i = 0; i < pmfVectorLength / 2; ++i, src1 += skip, src2 += skip)
	{
		ulong *dest = res + i;
		ulong s1 = (-src1[-1]) & (2 * pmfLength - 1);
		int neg1 = (s1 >= pmfLength);

		if(neg1)
		{
			s1 -= pmfLength;
		}

		ulong s2 = (-1 - src2[-1]) & (2 * pmfLength - 1);
		int neg2 = (s2 >= pmfLength);

		if(neg2)
		{
			s2 -= pmfLength;
		}

		ulong *x1;
		ulong *x2;

		if(s1 < s2)
		{
			x1 = src1;
			x2 = src2;
		}
		else
		{
			x1 = src2;
			x2 = src1;
			ulong s_temp = s1;
			s1 = s2;
			s2 = s_temp;
			int neg_temp = neg1;
			neg1 = neg2;
			neg2 = neg_temp;
		}

		dest = SkipSignedAdd(dest, pmfVectorLength / 2, pmfLength - s2, x2 + s2, neg2, x1 + s1, neg1, modFLINT);
		dest = SkipSignedAdd(dest, pmfVectorLength / 2, s2 - s1, x2, !neg2, x1 + s1 + pmfLength - s2, neg1, modFLINT);
		SkipSignedAdd(dest, pmfVectorLength / 2, s1, x2 + s2 - s1, !neg2, x1, !neg1, modFLINT);
	}
}
void NussbaumerPointwiseMultiplication(nmod_t modFLINT, ulong pmfLength, ulong pmfVectorLength, ulong *pmfVector1, ulong *pmfVector2, ptrdiff_t skip, ulong wordMod)
{
	ulong i;
	ulong *dest = pmfVector1;
	const ulong *src1 = pmfVector1;
	const ulong *src2 = pmfVector2;
	ulong *temp = (ulong *)flint_malloc(sizeof(ulong) * 2 * pmfLength);

	temp[2 * pmfLength - 1] = 0;

	for(i = 0; i < pmfVectorLength; ++i, dest += skip, src1 += skip, src2 += skip)
	{
		dest[0] = src1[0] + src2[0];

		MultiplicationInternal(temp, src1 + 1, pmfLength, src2 + 1, pmfLength, 1, modFLINT, wordMod);
		_nmod_vec_sub(dest + 1, temp, temp + pmfLength, pmfLength, modFLINT);
	}

	flint_free(temp);
}
ulong GetFudgeFactorFromNussbaumer(ulong lgL, nmod_t modFLINT, ulong wordMod)
{
	ulong lgK = (lgL / 2) + 1;
	ulong lgM = lgL - lgK + 1;
	ulong pmfVectorLength = 1UL << lgK;
	ulong pmfLength = 1UL << lgM;
	ulong fudge = GetFudgeFactorFromMultiplication(pmfLength, pmfLength, modFLINT, wordMod);

	return nmod_div(fudge, pmfVectorLength, modFLINT);
}
void NussbaumerMultiplication(ulong *pmfVectorResult, const ulong *pmfVector1, const ulong *pmfVector2, ulong lgM, ulong lgK, nmod_t modFLINT, ulong *nussbaumerVector1, ptrdiff_t skip,
                              ulong *nussbaumerVector2, ulong wordMod)
{
	assert(lgM + 1 >= lgK);

	ulong pmfLength = 1UL << lgM;
	ulong pmfVectorLength = 1UL << lgK;
	NussbaumerSplit(pmfVector1, modFLINT, lgK, lgM, nussbaumerVector1, skip);
	NussbaumerFFT(modFLINT, lgK, lgM, skip, nussbaumerVector1);
	NussbaumerSplit(pmfVector2, modFLINT, lgK, lgM, nussbaumerVector2, skip);
	NussbaumerFFT(modFLINT, lgK, lgM, skip, nussbaumerVector2);
	NussbaumerPointwiseMultiplication(modFLINT, pmfLength, pmfVectorLength, nussbaumerVector1, nussbaumerVector2, skip, wordMod);
	NussbaumerIFFT(pmfLength, modFLINT, lgK, nussbaumerVector1, skip);
	NussbaumerCombine(pmfVectorResult, pmfLength, modFLINT, nussbaumerVector1, skip, pmfVectorLength);
}


void ButterflyInPlace(ulong *op1, ulong *op2, ulong n, nmod_t modFLINT)
{
	for(; n; --n, ++op1, ++op2)
	{
		ulong x = *op1;
		ulong y = *op2;
		*op1 = nmod_add(y, x, modFLINT);
		*op2 = nmod_sub(y, x, modFLINT);
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

void TransformPolynomialToPMFVector(const ulong *poly, size_t polyLength, ulong x, ulong pmfLength, nmod_t modFLINT, ptrdiff_t skip, ulong *pmfVector)
{
	for(; polyLength >= pmfLength / 2; polyLength -= pmfLength / 2, poly += pmfLength / 2, pmfVector += skip)
	{
		pmfVector[0] = 0;
		_nmod_vec_scalar_mul_nmod(pmfVector + 1, poly, pmfLength / 2, x, modFLINT);
		_nmod_vec_zero(pmfVector + 1 + pmfLength / 2, pmfLength / 2);
	}

	if(polyLength)
	{
		pmfVector[0] = 0;
		_nmod_vec_scalar_mul_nmod(pmfVector + 1, poly, polyLength, x, modFLINT);
		_nmod_vec_zero(pmfVector + 1 + polyLength, pmfLength - polyLength);
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

	TransformPolynomialToPMFVector(poly1, poly1Length, 1, pmfLength, modFLINT, pmfLengthWithBias, pmfVector1);
	TransformPolynomialToPMFVector(poly2, poly2Length, x, pmfLength, modFLINT, pmfLengthWithBias, pmfVector2);

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


	ulong modBits = FLINT_CLOG2(modFLINT.n);
	int ks1Used = poly2Length < tuneMod[modBits][0];
	int ks2Used = poly2Length < tuneMod[modBits][1];
	int ks4Used = poly2Length < tuneMod[modBits][2];
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

	ulong modBits = FLINT_CLOG2(modFLINT.n);

	if(poly2Length < tuneMod[modBits][0])
	{
		_nmod_poly_mul_KS(polyResult, poly1, poly1Length, poly2, poly2Length, 0, modFLINT);
	}
	else if(poly2Length < tuneMod[modBits][1])
	{
		_nmod_poly_mul_KS2(polyResult, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(!moduloIsOdd || poly2Length < tuneMod[modBits][2])
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
void Multiplication(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, nmod_t modFLINT)
{
	assert(modFLINT.n > 1);

	// Not sure if it's the correct term but word refers to word size i.e. 2^FLINT_BITS, where FLINT_BITS=64 or 32 bits depending system.
	ulong wordMod = nmod_pow_ui(2, FLINT_BITS, modFLINT);

	MultiplicationInternal(polyResult, poly1, poly1Length, poly2, poly2Length, 0, modFLINT, wordMod);
}



void MultiplicationSSA(ulong *poly1, ulong *poly2, ulong poly1Length, ulong poly2Length, ulong *polyResult, nmod_t modFLINT)
{
	if(poly1Length >= poly2Length)
	{
		Multiplication(polyResult, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else
	{
		Multiplication(polyResult, poly2, poly2Length, poly1, poly1Length, modFLINT);
	}
}

//#############################################################################################################################
void MultiplicationFLINT(ulong *poly1, ulong *poly2, ulong poly1Length, ulong poly2Length, ulong *polyResult, nmod_t modFLINT)
{
	if(poly1Length >= poly2Length)
	{
		_nmod_poly_mul(polyResult, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else
	{
		_nmod_poly_mul(polyResult, poly2, poly2Length, poly1, poly1Length, modFLINT);
	}
}

//######################################################################################################################
void MultiplicationZN_POLY(ulong *poly1, ulong *poly2, ulong poly1Length, ulong poly2Length, ulong *polyResult, zn_mod_t modZN)
{
	if(poly1Length >= poly2Length)
	{
		zn_array_mul(polyResult, poly1, poly1Length, poly2, poly2Length, modZN);
	}
	else
	{
		zn_array_mul(polyResult, poly2, poly2Length, poly1, poly1Length, modZN);
	}
}
//######################################################################################################################

void PrintCharNTimes(char character, ulong numberOfTimes)
{
	for(ulong i = 0; i < numberOfTimes; ++i)
	{
		putchar(character);
	}
}

void WriteResultsToFile(double **data)
{
	ulong numberOfRows = NUMBER_OF_ROWS_RESULT;
	ulong numberOfColumns = NUMBER_OF_COLUMNS_RESULT;

	FILE *f = fopen("data.txt", "w");
	fprintf(f, "\\begin{table}[ht]\n");
	fprintf(f, "\\begin{center}\n");
	fprintf(f, "\\begin{tabular}{c|*{13}{R}}\n");
	fprintf(f, "\\multicolumn{1}{c|}{$k$} &\n");
	fprintf(f, "\\multicolumn{13}{c}{$n / 1024$} \\\\\n");
	fprintf(f, "\\multicolumn{1}{c|}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{1} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{2} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{4} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{8} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{16} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{32} &\n");
	fprintf(f, "\\multicolumn{1}{c}{} &\n");
	fprintf(f, "\\multicolumn{1}{c}{64} \n");
	fprintf(f, "\\\\");
	fprintf(f, "\\hline\n");

	for(ulong i = 0; i < numberOfRows; ++i)
	{

		fprintf(f, "%li ", (i + 1) * 5);

		for(ulong j = 0; j < numberOfColumns; ++j)
		{
			fprintf(f, "& %.2lf ", data[i][j]);

			if(j == numberOfColumns - 1)
			{
				fprintf(f, " \\\\");
			}
		}


		if(i != numberOfRows - 1)
		{
			fprintf(f, "\n");
		}
	}

	fprintf(f, "\n\\end{tabular}\n");
	fprintf(f, "\\end{center}\n");
	fprintf(f, "\\end{table}");

	fclose(f);
}

void PrintResults(double **data)
{
	double totalAverage = 0;
	double modAverage[NUMBER_OF_ROWS_RESULT];
	double degreeAverage[NUMBER_OF_COLUMNS_RESULT];

	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		modAverage[i] = 0;

		for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
		{
			modAverage[i] += data[i][j] / NUMBER_OF_COLUMNS_RESULT;
		}
	}

	for(ulong i = 0; i < NUMBER_OF_COLUMNS_RESULT; ++i)
	{
		degreeAverage[i] = 0;

		for(ulong j = 0; j < NUMBER_OF_ROWS_RESULT; ++j)
		{
			degreeAverage[i] += data[j][i] / NUMBER_OF_ROWS_RESULT;
		}
	}

	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
		{
			totalAverage += data[i][j] / (double)(NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT);
		}
	}

	putchar('\n');
	PrintCharNTimes(' ', 7);

	for(ulong i = 0; i < NUMBER_OF_COLUMNS_RESULT; ++i)
	{
		if(degreeAverage[i] < LOWER_LIMIT)
		{
			printf(ANSI_COLOR_RED "%.2lf  " ANSI_COLOR_RESET, degreeAverage[i]);
		}
		else if(degreeAverage[i] > UPPER_LIMIT)
		{
			printf(ANSI_COLOR_GREEN "%.2lf  " ANSI_COLOR_RESET, degreeAverage[i]);
		}
		else
		{
			printf(ANSI_COLOR_BLUE "%.2lf  " ANSI_COLOR_RESET, degreeAverage[i]);
		}
	}

	putchar('\n');
	PrintCharNTimes(' ', 7);
	PrintCharNTimes('_', 6 * NUMBER_OF_COLUMNS_RESULT - 2);
	putchar('\n');

	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		if(modAverage[i] < LOWER_LIMIT)
		{
			printf(ANSI_COLOR_RED "%.2lf  |" ANSI_COLOR_RESET, modAverage[i]);
		}
		else if(modAverage[i] > UPPER_LIMIT)
		{
			printf(ANSI_COLOR_GREEN "%.2lf  |" ANSI_COLOR_RESET, modAverage[i]);
		}
		else
		{
			printf(ANSI_COLOR_BLUE "%.2lf  |" ANSI_COLOR_RESET, modAverage[i]);
		}


		for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
		{
			if(data[i][j] < LOWER_LIMIT)
			{
				printf(ANSI_COLOR_RED "%.2lf  " ANSI_COLOR_RESET, data[i][j]);
			}
			else if(data[i][j] > UPPER_LIMIT)
			{
				printf(ANSI_COLOR_GREEN "%.2lf  " ANSI_COLOR_RESET, data[i][j]);
			}
			else
			{
				printf(ANSI_COLOR_BLUE "%.2lf  " ANSI_COLOR_RESET, data[i][j]);
			}
		}

		putchar('\n');
	}

	printf("\nThe numbers show how much faster the second implementation is compared to the first one.\n");
	printf("Average performance increase: ");

	if(totalAverage < LOWER_LIMIT)
	{
		printf(ANSI_COLOR_RED "%.2lf\n\n" ANSI_COLOR_RESET, totalAverage);
	}
	else if(totalAverage > UPPER_LIMIT)
	{
		printf(ANSI_COLOR_GREEN "%.2lf\n\n" ANSI_COLOR_RESET, totalAverage);
	}
	else
	{
		printf(ANSI_COLOR_BLUE "%.2lf\n\n" ANSI_COLOR_RESET, totalAverage);
	}
}

void BenchmarkMultiplication(ulong modBits, ulong maxPolyBits, flint_rand_t state)
{
	long maximumPossibleLength = 1024L << maxPolyBits / 2;

	if(maxPolyBits & 1)
	{
		maximumPossibleLength += maximumPossibleLength / 2;
	}


	ulong poly1Length = n_randint(state, maximumPossibleLength) + 1;
	ulong poly2Length = n_randint(state, maximumPossibleLength) + 1;
	ulong mod = n_randbits(state, modBits);

	ulong polyResultLength = poly1Length + poly2Length - 1;
	ulong *poly1 = (ulong *)flint_calloc(poly1Length, sizeof(ulong));
	ulong *poly2 = (ulong *)flint_calloc(poly2Length, sizeof(ulong));
	ulong *polyResult = (ulong *)flint_calloc(polyResultLength, sizeof(ulong));
	nmod_t modFLINT;
	nmod_init(&modFLINT, mod);
	zn_mod_t modZN;
	zn_mod_init(modZN, mod);

	for(ulong i = 0; i < poly1Length; ++i)
	{
		poly1[i] = n_randlimb(state);
	}

	for(ulong i = 0; i < poly2Length; ++i)
	{
		poly2[i] = n_randlimb(state);
	}

	_nmod_vec_reduce(poly1, poly1, poly1Length, modFLINT);
	_nmod_vec_reduce(poly2, poly2, poly2Length, modFLINT);


	double t = GetTime();

	if(chooseFunction == 0)
	{
		MultiplicationSSA(poly1, poly2, poly1Length, poly2Length, polyResult, modFLINT);
	}
	else if(chooseFunction == 1)
	{
		MultiplicationFLINT(poly1, poly2, poly1Length, poly2Length, polyResult, modFLINT);
	}
	else if(chooseFunction == 2)
	{
		MultiplicationZN_POLY(poly1, poly2, poly1Length, poly2Length, polyResult, modZN);
	}

	timeElapsed += GetTime() - t;

	flint_free(poly1);
	flint_free(poly2);
	flint_free(polyResult);
}

int main(int argc, char *argv[])
{
	if(argc != 3)
	{
		printf("Wrong number of arguments. Setting measurement time to 0.01 and function to SSA.\n");
	}
	else
	{
		measurementTimeInSeconds = atof(argv[1]);
		chooseFunction = atoi(argv[2]);
	}

	FILE *f;

	double data[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
	double dataRunOne[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
	ulong numberOfIterationsTable[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
	setbuf(stdout, NULL);

	flint_rand_t state;
	flint_randinit(state);
	flint_randseed(state, time(0), time(0));

	if((f = fopen("tuneMod.txt", "r")) != NULL)
	{
		for(ulong i = 0; i < 65; ++i)
		{
			for(ulong j = 0; j < 4; ++j)
			{
				if(!fscanf(f, "%lu", &tuneMod[i][j]))
				{
					printf("Failed to read integer.\n");
				}
			}
		}

		fclose(f);
	}
	else
	{
		printf("Failed to open tuneMod.txt.");
	}

	int currentlyOnRunOne = (f = fopen("data.txt", "r")) == NULL;
	int currentlyOnRunTwo = !currentlyOnRunOne;
	if(currentlyOnRunTwo)
	{
		if(!fread(dataRunOne, sizeof(double), NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT, f))
		{
			printf("Failed to read dataRunOne.\n");
		}

		fclose(f);

		if((f = fopen("numberOfIterations.txt", "r")) != NULL)
		{
			if(!fread(numberOfIterationsTable, sizeof(ulong), NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT, f))
			{
				printf("Failed to read dataRunOne.\n");
			}

			fclose(f);

		}
		else
		{
			printf("\nFailed to open numberOfIterations.txt\n");
		}
	}

	for(ulong modBits = 5, row = 0; modBits <= 60; modBits += 5, ++row)
	{
		for(long idx = 0, column = 0; idx < 13; ++idx, ++column)
		{
			long maxPolyBits = 1024 * (1L << idx / 2);
			if(idx & 1)
			{
				maxPolyBits += maxPolyBits / 2;
			}

			timeElapsed = 0;
			if(currentlyOnRunOne)
			{
				double startTime = GetTime();
				ulong numberOfIterations = 0;
				while(GetTime() - startTime < measurementTimeInSeconds)
				{
					BenchmarkMultiplication(modBits, maxPolyBits, state);
					++numberOfIterations;
				}
				numberOfIterationsTable[modBits / 5 - 1][idx] = numberOfIterations;
				data[modBits / 5 - 1][idx] = timeElapsed;
			}
			else
			{
				for(ulong i = 0; i < numberOfIterationsTable[modBits / 5 - 1][idx]; ++i)
				{
					BenchmarkMultiplication(modBits, maxPolyBits, state);
				}
				data[modBits / 5 - 1][idx] = timeElapsed;
			}
		}
	}


	if(currentlyOnRunOne)
	{
		f = fopen("data.txt", "w");
		fwrite(data, sizeof(double), NUMBER_OF_COLUMNS_RESULT * NUMBER_OF_ROWS_RESULT, f);
		fclose(f);

		f = fopen("numberOfIterations.txt", "w");
		fwrite(numberOfIterationsTable, sizeof(ulong), NUMBER_OF_COLUMNS_RESULT * NUMBER_OF_ROWS_RESULT, f);
		fclose(f);
	}
	else
	{
		double *performanceResultsEntries = (double *) malloc(sizeof(double) * NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT);
		double **performanceResults = (double **) malloc(sizeof(double *) * NUMBER_OF_ROWS_RESULT);
		for(ulong ix = 0, jx = 0; ix < NUMBER_OF_ROWS_RESULT; ++ix, jx += NUMBER_OF_COLUMNS_RESULT)
		{
			performanceResults[ix] = performanceResultsEntries + jx;
		}

		for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
		{
			for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
			{
				performanceResults[i][j] = dataRunOne[i][j] / data[i][j];
			}
		}

		PrintResults(performanceResults);

		WriteResultsToFile(performanceResults);
		flint_free(performanceResults);
		flint_free(performanceResultsEntries);
	}
}
