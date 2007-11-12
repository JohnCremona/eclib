// twoadic.cc: implementation of functions for existence of 2-adic points
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
 
#include <iostream>
#include "mlocsol.h"
#include "twoadic.h"

// The following due to Michael Stoll:

/* try1(poly), with poly a deg 3 polynomial in x, determines if */
/* there is a 2-adic integer a such that poly(a) is a square    */
/* in Q_2. Returns 1 if successful, 0 otherwise.                */
long try1(long poly[4])
{ long c, v, sq, j, k;
  /* see if we can lift to obtain a zero of poly.      */
  /* the condition is  val(poly[0]) > 2*val(poly[1]).  */
#ifdef SHOW_POLY      
 cout<<"("<<poly[3]<<","<<poly[2]<<","<<poly[1]<<","<<poly[0]<<")"<<endl;
#endif
  c = poly[0];
  v = val(c);
  // Hensel test for a root:
  if(c == 0 || v > 2*val(poly[1])) 
    {
      //      cout<<" has a 2-adic root via Hensel"<<endl;
      return 1;
    }
  // Newton polygon test for a root:
  if((val(poly[0]) + val(poly[2]) > 2*val(poly[1])) &&
     (2*val(poly[0]) + val(poly[3]) > 3*val(poly[1])))
    {
      //      cout<<" has a 2-adic root via Newton Polygon"<<endl;
      return 1;
    }
  /* constant term a square -> OK */
  if((v & 0x1) == 0)
  { c >>= v; c &= 0x7;
    if(c == 1) return 1;
    sq = (((c & 0x3) == 3) ? 2 : 3) + v;
  }
  else sq = v+1;
  /* if not, see if we can already read off that it's impossible */
  /* to have a point with integral 2-adic x-coordinate.          */
  for(j = 1; j <= 3; j++) if(val(poly[j]) < sq) goto refine;
  return 0;
refine:
  /* now we must refine */
  { long newpoly[4];
    /* construct poly(2x+1) */
    for(j = 0; j < 4; j++) newpoly[j] = poly[j];
    for(j = 0; j < 4; j++)
    { long h = newpoly[3];
      for(k = 2; k >= j; k--) h = newpoly[k] += h;
      newpoly[j] <<= j;
    }
    if(try1(newpoly)) return 1;  /* success with 2x */
    /* construct poly(2x) */
    for(j = 0; j < 4; j++) newpoly[j] = poly[j]<<j;
    if(try1(newpoly)) return 1; /* success with 2x */
    return 0;
} }
/*
It is possible to be a little bit more sophisticated about the condition
for the existence of a 2-adic root -- it is enough to have
  val(poly[0]) > val(poly[1])  and
  a vertex of the Newton polygon at 1, i.e.,
    val(poly[0]) + val(poly[2]) > 2*val(poly[1]) and
    2*val(poly[0]) + val(poly[3]) > 3*val(poly[1]) .
The condition in the code comes from Hensel's lemma and implies the
conditions given above.
*/

long try1(bigint poly[4])
{ long c, v, sq, j, k;
 bigint mc; 
  /* see if we can lift to obtain a zero of poly.      */
  /* the condition is  val(poly[0]) > 2*val(poly[1]).  */
#ifdef SHOW_POLY      
 cout<<"("<<poly[3]<<","<<poly[2]<<","<<poly[1]<<","<<poly[0]<<")"<<endl;
#endif
  mc = poly[0];
  v = val(mc);
  // Hensel test for a root:
  if(mc == 0 || v > 2*val(poly[1])) 
    {
      //      cout<<" has a 2-adic root via Hensel"<<endl;
      return 1;
    }
  // Newton polygon test for a root:
  if((val(poly[0]) + val(poly[2]) > 2*val(poly[1])) &&
     (2*val(poly[0]) + val(poly[3]) > 3*val(poly[1])))
    {
      //      cout<<" has a 2-adic root via Newton Polygon"<<endl;
      return 1;
    }
  /* constant term a square -> OK */
  if((v & 0x1) == 0)
    { mc >>= v; c  = posmod(mc,8);  //& 0x7;
    if(c == 1) return 1;
    sq = (((c & 0x3) == 3) ? 2 : 3) + v;
  }
  else sq = v+1;
  /* if not, see if we can already read off that it's impossible */
  /* to have a point with integral 2-adic x-coordinate.          */
  for(j = 1; j <= 3; j++) if(val(poly[j]) < sq) goto refine;
  return 0;
refine:
  /* now we must refine */
  { bigint newpoly[4];
    /* construct poly(2x) */
    for(j = 0; j < 4; j++) newpoly[j] = poly[j]<<j;
    if(try1(newpoly)) return 1; /* success with 2x */
    /* construct poly(2x+1) */
    for(j = 0; j < 4; j++) newpoly[j] = poly[j];
    for(j = 0; j < 4; j++)
    { bigint h = newpoly[3];
      for(k = 2; k >= j; k--) h = newpoly[k] += h;
      newpoly[j] <<= j;
    }
    if(try1(newpoly)) return 1;  /* success with 2x */
    return 0;
} }

long case1(long a, long b) // A=4a, B=4b
{
  //  cout<<"In case1() with a="<<a<<", b="<<b<<endl;
  long c=a+3, d=2*a+b+2;
// The polynomial in 4x+2 is (16, 24, 4*c, d)
  long d8=d&7, c2=c&1;
  long d4=d8&3;
  if((d4==2)||(d4==3)) return 0;
  if(d4==1) return (c2||(d8==1));
  if(c2) return 1;
// Now, d = 0 mod 4. Divide by 4 to get (4, 6, c, d/4)
  d>>=2;
  c>>=1; 
  a=b=1;
  while(1)
    {
      // The polynomial is (4*a, 4*b+2, 2*c, d)
      //      cout<<"(a,b,c,d)=("<<a<<","<<b<<","<<c<<","<<d<<")"<<endl;
      d8=d&7; d4=d8&3; c2=c&1;
      if (c2) return ((d4==0)||(d4==1));
      if (d4&1) return ((d8==1)||  (((4*a+4*b+2*c-1)&7)==d8));
      if (d4 == 0)
	{    // any solution x must be even
	  d >>=2; c >>=1; a <<=1;
	}
      else  // any solution x must be odd
	{
	  d = a+b+c/2+(d+2)/4; c = 3*a+2*b+c/2+1; b = 3*a+b; a = 2*a;
	}
    }
}

//JC's version:
/*
long case2(long a, long b) // A=4a+1, B=4b+2
{
  //  cout<<"In case2() with a="<<a<<", b="<<b<<endl;
  long c=a+1, d=a+b+1, temp;
// The polynomial in 4x+1 is (16, 12, 4*c, d)
  long d8=d&7, c2=c&1;
  long d4=d8&3;
  if((d4==2)||(d4==3)) return 0;
  if(d4==1) return (c2==0||(d8==1));
  if(c2==1) return 1;
  // Now, d = 0 mod 4. Divide by 4 and replace c,d by c/2,d/4
  // to get (4, 3, 2*c, d)
  d>>=2;
  c>>=1; 
  a=1; b=0;  // The polynomial is (4*a, 3*(4*b+1), 2*c, d)
  while(1)
    {
      //      cout<<"(a,b,c,d)=("<<a<<","<<b<<","<<c<<","<<d<<")"<<endl;
      d8=d&7; d4=d8&3; c2=c&1;
      if(d8==1) return 1;
      if (d4&1) // d is odd
	{
	  if((c2==0)&&(d8==5)) return 1;
	  temp=2*c+d+3;
	  if(temp&3) return 0;
	  if(c2==0) return 1;
	  // Now any root must be odd, so change poly	  
	  d = a+3*b+temp/4; c = 3*a+6*b+(c+3)/2; b = a+b; a = 2*a;
	}
      else
	{
	  temp=4*(a+b)+2+d;
	  if(((temp+2*c)&7)==0) return 1;
	  if(((temp-2*c)&7)==0) return 1;
	  if(d4==2) return 0;
	  // Now any root must be even, so change poly	  
	  d >>=2; c >>=1; a <<=1;  
	}
    }
}
*/

//MS's version:

long case2(long a, long b) // A=4a+1, B=4b+2
{
  long c=a+1, d=a+b+1;
  // The polynomial in 4x+1 is (16, 12, 4*c, d)
  long d8=d&7, c2=c&1;
  long d4=d8&3;
  if((d4==2)||(d4==3)) return 0;
  if(d4==1) return (c2==0||(d8==1));
  if(c2==1) return 1;
  // Now, d = 0 mod 4. Divide by 4 and replace c,d by c/2,d/4
  // to get (4, 3, 2*c, d)
  d>>=2;
  c>>=1;
  a=1; b=0;  // The polynomial is (4*a, 3*(4*b+1), 2*c, d)
  while(1)
    {
      d8=d&7; d4=d8&3;
      if(c&1) // c odd
      { 
	switch(d4){
	case 0: return 1; break;
	case 2: return 0; break;
	case 1: return (d8==1); break;
	case 3: 
	  { // replace f(x) by f(2x+1)/4 and loop
	    d = a + 3*b + (c+1)/2 + (d+1)/4 ;
	    c = 3*a + 6*b + (c+3)/2;
	    b += a ;
	    a <<= 1; 
	  } 
	}
      }
      else   // c even
      { 
	switch(d4){
	case 1: return 1; break;
	case 3: return 0; break;
	case 2: return (((d8+4*(a+b)+2*c+2)&7) == 0);
	case 0: 
	  {  // replace f(x) by f(2x)/4 and loop
	    d >>= 2; c >>= 1; a <<= 1; 
	  }
	}
      }
    }
}

// bigint versions of case1() and case2():

long case1(bigint a, bigint b) // A=4a, B=4b
{
  //  cout<<"In case1() with a="<<a<<", b="<<b<<endl;
  bigint c=a+3, d=2*a+b+2;
// The polynomial in 4x+2 is (16, 24, 4*c, d)
  long d8=posmod(d,8), c2=posmod(c,2);
  long d4=d8&3;
  //  cout<<"c="<<c<<", d="<<d<<endl;
  //  cout<<"d4="<<d4<<", d8="<<d8<<endl;
  //  cout<<"c2="<<c2<<endl;
  if((d4==2)||(d4==3)) return 0;
  if(d4==1) return (c2||(d8==1));
  if(c2) return 1;
// Now, d = 0 mod 4. Divide by 4 to get (4, 6, c, d/4)
  d>>=2;
  c>>=1; 
  a=b=1;
  while(1)
    {
      // The polynomial is (4*a, 4*b+2, 2*c, d)
      //      cout<<"(a,b,c,d)=("<<a<<","<<b<<","<<c<<","<<d<<")"<<endl;
      d8=posmod(d,8); d4=d8&3; c2=posmod(c,2);
      if (c2) return ((d4==0)||(d4==1));
      if (d4&1) return ((d8==1)||  ((posmod(4*a+4*b+2*c-1,8))==d8));
      if (d4 == 0)
	{    // any solution x must be even
	  d >>=2; c >>=1; a <<=1;
	}
      else  // any solution x must be odd
	{
	  d = a+b+c/2+(d+2)/4; c = 3*a+2*b+c/2+1; b = 3*a+b; a = 2*a;
	}
    }
}

long case2(bigint a, bigint b) // A=4a+1, B=4b+2
{
  bigint c=a+1, d=a+b+1;
  // The polynomial in 4x+1 is (16, 12, 4*c, d)
  long d8=posmod(d,8), c2=posmod(c,2);
  long d4=d8&3;
  if((d4==2)||(d4==3)) return 0;
  if(d4==1) return (c2==0||(d8==1));
  if(c2==1) return 1;
  // Now, d = 0 mod 4. Divide by 4 and replace c,d by c/2,d/4
  // to get (4, 3, 2*c, d)
  d>>=2;
  c>>=1;
  a=1; b=0;  // The polynomial is (4*a, 3*(4*b+1), 2*c, d)
  while(1)
    {
      d8=posmod(d,8); d4=d8&3; c2=posmod(c,2);
      if(c2) // c odd
      { 
	switch(d4){
	case 0: return 1; break;
	case 2: return 0; break;
	case 1: return (d8==1); break;
	case 3: 
	  { // replace f(x) by f(2x+1)/4 and loop
	    d = a + 3*b + (c+1)/2 + (d+1)/4 ;
	    c = 3*a + 6*b + (c+3)/2;
	    b += a ;
	    a <<= 1; 
	  } 
	}
      }
      else   // c even
      { 
	switch(d4){
	case 1: return 1; break;
	case 3: return 0; break;
	case 2: return ((posmod(-4*(a+b)-2*c-2,8)) == d8);
	case 0: 
	  {  // replace f(x) by f(2x)/4 and loop
	    d >>= 2; c >>= 1; a <<= 1; 
	  }
	}
      }
    }
}
