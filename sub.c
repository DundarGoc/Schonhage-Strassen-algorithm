#include "immintrin.h"

void
test(
    int c[],
    int const a[],
    int const b[]
    )
{
  __m256i ar = _mm256_loadu_si256((__m256i const*)a);
  __m256i br = _mm256_loadu_si256((__m256i const*)b);

  __m256i cr = _mm256_sub_epi32(ar, br);

  _mm256_storeu_si256((__m256i *)c, cr);
}
