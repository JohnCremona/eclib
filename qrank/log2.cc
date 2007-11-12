#include <math.h>
#include "log2.h"


// Function which returns 1 and sets e such that 2**e=n if n is a power of 2.
// If the "roundup" flag is set and n is not a power of 2 it increases n to
// the next power of 2 (and returns 0)

int intlog2(long& n, long& e, int roundup)
{ 
  e = 0;
  if (n<1) {if(roundup) n=1; return 0;}
  long m=n;
  while (m) { m >>= 1; e++; }
  e--;
  m=1<<e;
//                     at this point m=2^e <= n < 2^(e+1)
  if(m==n) return 1;
  if(roundup) {n=m<<1; e++;}
  return 0;
}

