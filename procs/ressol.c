\// The following code is borrowed from LiDIA file ressol.c
// simplified to work under libg++
//
// LiDIA - a library for computational number theory
//   Copyright (c) 1994, 1995 by the LiDIA Group
//
// File        : ressol.c 
// Author      : Thomas Denny (TD), Andreas M"uller (AM)
// Last change : AM, Feb 7 1995, initial version
//		 Oliver Morsch, Mar 1996, bigint-version
//		 MM, Sep 24 1996, cleaning up
//


#include<LiDIA/bigint.h>

//
// Input:
//    prime number p
//    0 <= a < p
//
// Output:
//    square root r of a mod p
//    0 <= r < p
// 
// Condition: p = 2^s * (2k+1) + 1
//	      where s fits in a long
//
// Algorithm: Shanks-Ressol
//

long jacobi(const bigint& a, long bigint& p)
{
  return legendre(a,p);
}

power_mod(bigint& ans, const bigint& base, const bigint& expo, const bigint& m)
{
  bigint ans(1), b(base), e(expo);
  while(e>0)
  {
    if(odd(e)) {ans*=b; ans%=m;}
    b*=b; b%=m;
    e>>=1;
  }
}

void ressol(bigint & x, const bigint & a, const bigint & p)
{
  bigint v, one(1), two(2);
   
  if ( p == two ) { x=a;  return;   }
 
  // (1) p = 3 mod 4 
  if ( (p%4 == 3 ) {
    if (jacobi(a, p) == 1) {
      v = (p + one)/4;
      power_mod(x, a, v, p);
      return;
    }
    else
      {
	cerr<<"In ressol: a = " << a << " is not a quadratic residue";
	cerr<< " mod p = "<<p<<endl;
	exit(1);
      }
  }

  bigint r, n, c, z, k;
  long   s, t;

  // (2) initialisation:
  // compute k, s : p = 2^s (2k+1) + 1 
  // (MM: shift_right as part of the loop)
  k = p - one;
  s = divide_out(k,two);
  
  k = (k-one)/two;
  
  // (3) initial values
  power_mod(r, a, k, p);
  
  square(n, r);
  n %= p; // n = (r * r) % p;
  n *= a;
  n %= p; // n = (n * a) % p;
  r *= a;
  r %= p; // r = (r * a) % p;
  
  if ( n==one ) {
    x=r;
    return;
  }
  
  // (4) non-quadratic residue
  z=two;
  while ( jacobi(z,p) == 1 ) z+=one;
  
  v=(k*two)+one; // v = (k << 1) + 1;
  power_mod(c, z, v, p);
  
  // (5) iteration
  while (n > one)
    {
      k  = n;
      t  = s;
      s  = 0;
      
      while ( k!=one )
	{
	  square(k, k);
	  k%= p;         // k = (k * k) % p; 
	  s++; 
	}
      
      t -= s; 
      if (t == 0)
	{
	  cerr<<"In ressol: a = " << a << " is not a quadratic residue";
	  cerr<< " mod p = "<<p<<endl;
	  exit(1);
	}
       
      v = one<<(t-1);     // v = 1 << (t-1);
      power_mod(c, c, v, p);
      r *= c;
      r %= p;             // r = (r * c) % p;
      square(c, c);
      c %= p;             // c = (c * c) % p;
      n *= c;
      n %= p;             // n = (n * c) % p;
    } 
  
  x=r;
  return;
} 
