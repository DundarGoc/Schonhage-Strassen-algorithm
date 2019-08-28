void
sub_novec(
    int c[],
    int const a[],
    int const b[]
    )
{
  for(int i=0; i<8; ++i){
    c[i]=a[i]-b[i];
  }
}
