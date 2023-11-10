// xmod.h: declarations of basic modular arithmetic functions 
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
// All are inline functions, there's no .cc file

#ifndef _ECLIB_XMOD_H
#define _ECLIB_XMOD_H      1

#include <eclib/interface.h>

// undefine this to use int/long/longlong arithmetic only

//#define USE_DMOD

// We'll define BIGPRIME to be one of the following, the largest prime
// possible for the given strategy without overlow.

// We allow all residues a with |a|<p, for speed

// This condition is ASSUMED ON INPUT to the xmodmul functions

// In fact we only use PRIME27 and PRIME30 which are ints.

#define PRIME27   134217689  // = largest p such that p < 2^27
#define PRIME63  6074000981  // = largest p such that (p/2)^2 < 2^63. 
#define PRIME30  1073741789  // = largest p such that p < 2^30. 
#define PRIME31a 2147483647  // = largest p such that p < 2^31. 
#define PRIME31b      92681  // = largest p such that (p/2)^2 < 2^31. 

//////////////////////////////////////////////////////////////////////
#ifdef USE_DMOD

// Strategy: whether base type is 32- or 64-bit, modular operations
// are done via doubles which have 52 bits for the mantissa, so any
// modulus < 2^27 is ok, and the default modulus PRIME27 is the
// largest such prime.

#include "math.h" // for floor()
#define XMOD_METHOD "doubles"

const int BIGPRIME =  PRIME27;
const double BIGPRIME_D =  PRIME27;
const double BIGPRIME_D_INV =  1/BIGPRIME_D;
const int BIGPRIME_I =  PRIME27;

inline long xmod(long a, double m) {return (long)(a-(m*floor((a/m)+0.5)));}
inline int  xmod(int  a, double m) {return (int)(a-(m*floor((a/m)+0.5)));}
inline long xmod(long a, long m)   {return (long)(a-(m*floor((a/m)+0.5)));}
inline int  xmod(int  a, int m)    {return (int)(a-(m*floor((a/m)+0.5)));}
inline int  xmod(long a, int m)    {return (int)(a-(m*floor((a/m)+0.5)));}

inline long xmod0(long a) {return (long)(a-(BIGPRIME_D*floor((a*BIGPRIME_D_INV))));}
inline int  xmod0(int a)  {return (int)(a-(BIGPRIME_D*floor((a*BIGPRIME_D_INV))));}

inline long mod0(long a) {return (long)(a-(BIGPRIME_D*floor((a*BIGPRIME_D_INV)+0.5)));}
inline int  mod0(int a)  {return (int)(a-(BIGPRIME_D*floor((a*BIGPRIME_D_INV)+0.5)));}

inline long xmodmul(long a, long b, double m)
{ double c = (double)a * (double)b;
  return (long)(c-m*floor((c/m)+0.5));
}
inline int xmodmul(int a, int b, double m)
{ double c = (double)a * (double)b;
  return (int)(c-m*floor((c/m)+0.5));
}
inline long xmodmul0(long a, long b)
{ double c = (double)a * (double)b;
  return (long)(c-BIGPRIME_D*floor((c*BIGPRIME_D_INV)+0.5));
}
inline int xmodmul0(int a, int b)
{ double c = (double)a * (double)b;
  return (int)(c-BIGPRIME_D*floor((c*BIGPRIME_D_INV)+0.5));
}
inline long xmodmul(long a, long b, long m)
{ double c = (double)a * (double)b;
  return (long)(c-m*floor((c/m)+0.5));
}
inline int xmodmul(int a, int b, int m)
{ double c = (double)a * (double)b;
  return (int)(c-m*floor((c/m)+0.5));
}


//////////////////////////////////////////////////////////////////////
#else // use some int/long/longlong combination

// Strategy: modular multiplication for ints and longs is done via
// int64_ts which have 63 bits
// modular addition is simply xmod(a+b,m) since 2*m<2^31

#define XMOD_METHOD "ints and longs"
const int BIGPRIME = PRIME30;
const int HALF_BIGPRIME = BIGPRIME>>1; // = 536870894;
const int TWO_BIGPRIME = 2147483578; // 2*BIGPRIME
const int64_t INV_BIGPRIME = 4294967436LL; // = 2^32+140 = [2^62/p]

inline int xmod(int a, int m) {return a%m;}
inline long xmod(long a, long m) {return a%m;}

inline int xmod(long a, int m) { return (int)(a%(long)m);}
inline long xmod(int a, long m) { return (long)((long)(a)%m);}

inline int xmod0(int a) {return a%BIGPRIME;}
inline int mod0(int a)
{a%=BIGPRIME;
 if(a>0)
   while(a>HALF_BIGPRIME) a-=BIGPRIME; 
 else
   while(-a>HALF_BIGPRIME) a+=BIGPRIME; 
 return a;
}

inline long xmod0(long a) {return (a%(long)BIGPRIME);}
inline long mod0(long a) {return mod0((int)(a%BIGPRIME));}

inline int xmodmul0(int a, int b)
{
  return ((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(BIGPRIME) ))%BIGPRIME;
}

inline long xmodmul0(long a, long b)
{
  return (long)(((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(BIGPRIME) ))%BIGPRIME);
}

// This special version only works modulo BIGPRIME, not a general modulus:
// It should work faster (no divisions)!  Thanks to David Harvey.

inline int xmm0(int a, int b)
{
  if (a==1) return b;
  if (a==-1) return -b;
  if (b==1) return a;
  if (b==-1) return -a;
  // check:
  //  int r2 = (a*(int64_t)b) % BIGPRIME;
  if(a<0) a+=BIGPRIME;
  if(b<0) b+=BIGPRIME;
  int64_t ab = a*(int64_t)b;
  int64_t r = ab-((INV_BIGPRIME*(ab>>30))>>32)*BIGPRIME;
  r -= ( ((r>=TWO_BIGPRIME)?BIGPRIME:0) + ((r>=BIGPRIME)?BIGPRIME:0) );
  if (r>HALF_BIGPRIME) r-=BIGPRIME;
  // check:
  // if (r!=r2)
  //   {
  //     cout << "Problem with "<<a<<"*"<<b<<" (mod "<<BIGPRIME
  //          <<"): computed "<<r<<", not "<<r2<<endl;
  //     return r2;
  //   }
  return (int)r;
}

inline long xmm0(long a, long b)
{
  return (a*(int64_t)b) % BIGPRIME;
}

inline int xmodmul(int a, int b, int m)
{
  if (m==BIGPRIME) return xmm0(a,b);
  return ((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(m) ))%m;
}

inline int xmodmul(int a, int b, long m)
{
  return (int)(((long)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(m) ))%m);
}

inline int xmodmul(long a, long b, int m)
{
  return ((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(m) ))%m;
}

inline long xmodmul(long a, long b, long m)
{
  return ((long)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(m) ))%m;
}

#endif // ifdef USE_DMOD

#if(1)
const int DEFAULT_MODULUS = BIGPRIME;
// table of inversers of residues<20 modulo BIGPRIME:
static int table_invs[20] = {0,1, 536870895, 357913930, 805306342, 214748358, 178956965, 920350105, 402653171, 477218573, 107374179, 97612890, 626349377, 330382089, 997045947, 71582786, 738197480, 442128972, 775480181, 226050903};
#else
const int DEFAULT_MODULUS = 1073741783; //BIGPRIME-6;
// table of inversers of residues<20 modulo 1073741783:
static int table_invs[20] = {0, 1, 536870892, 357913928, 268435446,
 644245070, 178956964, 460175050, 134217723, 835132498, 322122535,
 780903115, 89478482, 743359696, 230087525, 930576212, 603979753,
 884257939, 417566249, 395589078};
#endif

inline long invmod0(long aa)
{
  long a=aa;

  // if |a| is small, use look-up table:
  if ((a>0)&&(a<20)) return table_invs[a];
  long ma=-a;
  if ((ma>0)&&(ma<20)) return -table_invs[ma];
  // if a = BIGPRIME-ma with ma small, use look-up table:
  ma+=BIGPRIME; // = BIGPRIME-a
  if ((ma>0)&&(ma<20)) return -table_invs[ma];
  // if a = -BIGPRIME+ma with ma small, use look-up table:
  ma=a-BIGPRIME;
  if ((ma>0)&&(ma<20)) return table_invs[ma];

 // General code, use Euclidean Algorithm:
 long x=0,oldx=1,newx,b=BIGPRIME,c,q;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a==1)  {return oldx;}
 if (a==-1) {return -oldx;}
 cout << "invmod0 called with " << a << " -- not invertible!\n";
 return 0;
}

inline int invmod0(int aa)
{
  int a=aa;

  // if |a| is small, use look-up table:
  if ((a>0)&&(a<20)) return table_invs[a];
  int ma=-a;
  if ((ma>0)&&(ma<20)) return -table_invs[ma];
  // if a = BIGPRIME-ma with ma small, use look-up table:
  ma+=BIGPRIME; // = BIGPRIME-a
  if ((ma>0)&&(ma<20)) return -table_invs[ma];
  // if a = -BIGPRIME+ma with ma small, use look-up table:
  ma=a-BIGPRIME;
  if ((ma>0)&&(ma<20)) return table_invs[ma];

 // General code, use Euclidean Algorithm:
 int x=0,oldx=1,newx,b=BIGPRIME,c,q;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a==1)  {return oldx;}
 if (a==-1) {return -oldx;}
 cout << "invmod0 called with " << a << " -- not invertible!\n";
 return 0;
}


#endif // ifndef _XMOD_H
