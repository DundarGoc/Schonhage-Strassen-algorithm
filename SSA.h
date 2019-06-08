#pragma once

#include <flint/nmod_poly.h>
#include <assert.h>

typedef struct
{
	ulong thresholdKS2;
	ulong thresholdKS4;
	ulong thresholdFFT;
}
tuningModulo_t;

extern tuningModulo_t tuningModulo[];

void PMFVectorFFTIterative(ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorTFTDivideAndConquer(ulong numberOfOutputCoefficients, ulong z, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorTFTHuge(ulong numberOfOutputCoefficients, ulong numberOfInputCoefficients, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorTFT(ulong numberOfOutputCoefficients, ulong numberOfInputCoefficients, ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorIFFTIterative(ulong t, ulong lgK, ulong lgM, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorITFTDivideAndConquer(ulong n, int fwd, ulong z, ulong t, ulong lgM, ulong lgK, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);
void PMFVectorITFTHuge(ulong lgT, ulong n, int fwd, ulong z, ulong t, nmod_t modFLINT, ulong lgK, ptrdiff_t skip, ulong lgM, ulong *pmfVector);
void PMFVectorITFT(ulong n, int fwd, ulong z, ulong t, ulong lgM, ulong lgK, ptrdiff_t skip, ulong *pmfVector, nmod_t modFLINT);

void ButterflyInPlace(ulong *op1, ulong *op2, ulong n, nmod_t modFLINT);
void PMFButterfly(ulong *op1, ulong *op2, ulong pmfLength, nmod_t modFLINT);
void PMFAdd(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT);
void PMFSub(ulong *op1, const ulong *op2, ulong pmfLength, nmod_t modFLINT);
ulong GetFudgeFactorFromPointwiseMultiplication(ulong lgM, nmod_t modFLINT, ulong wordMod);
void PMFVectorPointwiseMultiplication(ulong n, nmod_t modFLINT, ulong lgM, ulong *data1, ulong *data2, ptrdiff_t skip, ulong wordMod);

void TransformPolynomialToPMFVector(const ulong *poly, size_t polyLength, ulong x, ulong pmfLength, nmod_t modFLINT, ulong *pmfVector);
void NegateOrCopy(ulong *res, const ulong *op, ulong n, int negate, nmod_t modFLINT);
void SchonhageStrassenCombineChunk(ulong *res, size_t n, ulong *op1, ulong *op2, ulong pmfLength, nmod_t modFLINT);
void TransformPMFVectorToPolynomial(ulong *poly, ulong polyLength, ulong m3, ulong pmfLength, nmod_t modFLINT, ptrdiff_t skip, ulong *pmfVector);
ulong GetFudgeFactorFromTFT(size_t poly1Length, size_t poly2Length, nmod_t modFLINT, ulong wordMod);
void MultiplicationSchonhageStrassen(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, ulong x, nmod_t modFLINT, ulong wordMod);

ulong *SkipSignedAdd(ulong *vectorResult, ptrdiff_t skipSize, size_t n, const ulong *vector1, int negateVector1, const ulong *vector2, int negateVector2, nmod_t modFLINT);

ulong GetFudgeFactorFromMultiplication(ulong poly1Length, size_t poly2Length, nmod_t modFLINT, ulong wordMod);
void MultiplicationInternal(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, int useREDCIfOddModulo, nmod_t modFLINT, ulong wordMod);
void SSA(ulong *polyResult, const ulong *poly1, size_t poly1Length, const ulong *poly2, size_t poly2Length, nmod_t modFLINT);
