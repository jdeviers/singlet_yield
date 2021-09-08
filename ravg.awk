#!/usr/bin/awk -f

BEGIN{
  rsum = 0.;
  rcount = 0.
}
{
  rsum = $1;
  rcount += 99.;           # N(N-1) offdiag elements
  print rcount,rsum/rcount
}
