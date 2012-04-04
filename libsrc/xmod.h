// xmod.h: declarations of basic modular arithmetic functions 
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 

// All are inline functions, there's no .cc file

#ifndef _XMOD_H
#define _XMOD_H      1

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
const int HALF_BIGPRIME = BIGPRIME>>1;

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

inline int xmodmul(int a, int b, int m)
{
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

inline int xmodmul0(int a, int b)
{
  return ((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(BIGPRIME) ))%BIGPRIME;
}

inline long xmodmul0(long a, long b) 
{
  return (long)(((int)( ( (int64_t)(a)*(int64_t)(b) ) % (int64_t)(BIGPRIME) ))%BIGPRIME);
}

#endif // ifdef USE_DMOD

inline long invmod0(long aa)
{long x=0,oldx=1,newx,a=aa,b=BIGPRIME,c,q;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a==1)  {return oldx;}
 if (a==-1) {return -oldx;}
 cout << "invmod0 called with " << a << " -- not invertible!\n";
 abort();
 return 0;
}

inline int invmod0(int aa)
{int x=0,oldx=1,newx,a=aa,b=BIGPRIME,c,q;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
  }
 if (a==1)  {return oldx;}
 if (a==-1) {return -oldx;}
 cout << "invmod0 called with " << a << " -- not invertible!\n";
 abort();
 return 0;
}


#endif // ifndef _XMOD_H
