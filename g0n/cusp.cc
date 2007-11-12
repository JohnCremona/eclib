// FILE CUSP.CC
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

#include "moddata.h"
#include "symb.h"
#include "cusp.h"

int cusplist::cuspeq(const rational& c1, const rational& c2) const
{
//cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2<<endl;
   long p1 = num(c1), p2 = num(c2), q1 = den(c1), q2 = den(c2);
   long s1,r1,s2,r2;
   bezout(p1,q1,s1,r1);  s1*=q2;
   bezout(p2,q2,s2,r2);  s2*=q1;
   long q3 = N->gcd(q1*q2);
   int ans = ((s1-s2)%q3==0) || (N->plusflag && ((s1+s2)%q3==0));
//cout<<"Returning "<<ans<<endl;   
   return ans;
}
 
long cusplist::index(const rational& c)
{  // adds c to list if not there already, and return index
   long ans=-1;
   for (long i=0; (i<number) && (ans<0); i++) if (cuspeq(c,list[i]))  ans=i;
   if (ans==-1) {list[number]=c; ans=number; number++;}
   return ans;
}
