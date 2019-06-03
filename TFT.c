#include "include/TFT.h"

#define CACHE_SIZE_ESTIMATE 32768

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
