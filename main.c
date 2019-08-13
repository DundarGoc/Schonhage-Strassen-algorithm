#include <stdio.h>
#include <stdlib.h>

int main(){

  int *a = malloc(sizeof(int) * 8);
  int *b = malloc(sizeof(int) * 8);
  int *c = malloc(sizeof(int) * 8);

  for(int i=0; i<8; ++i){
    a[i]=8-i;
    b[i]=i;
  }

  for(int i=0; i<8; ++i){
    c[i]=a[i]-b[i];
  }

  free(a);
  free(b);
  free(c);

}
