#include <immintrin.h>  // portable to all x86 compilers
#include <stdio.h>

void print128_num(__m128 var);

int main()
{
  // high element first, opposite of C array order.  Use _mm_setr_ps if you want "little endian" element order in the source.
  __m128 vector1 = _mm_set_ps(4, 3, 2, 1);
  __m128 vector2 = _mm_set_ps(7, 8, 9, 0);

  // result = vector1 + vector 2
  __m128 sum = _mm_add_ps(vector1, vector2);

  // vector1 is now (1, 2, 3, 4) (above shuffle reversed it)
  vector1 = _mm_shuffle_ps(vector1, vector1, _MM_SHUFFLE(0,1,2,3));
 
  print128_num(sum);
}

void print128_num(__m128 var) 
{
  int64_t *v64val = (int64_t*) &var;
  printf("%.16lx %.16lx\n", v64val[1], v64val[0]);
}
