#include "stdio.h"


void
test(
    int c[],
    int const a[],
    int const b[]
    );

int
main()
{
  int const a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int const b[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int c[8];
  
  test(c,a,b);

  printf("%i %i %i %i %i %i %i %i\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);

  return 0;
}
