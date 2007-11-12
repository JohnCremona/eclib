// divpol.cc: implementations of functions for division polynomials
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
 

// NB the external interface uses a simple vector<bigint> of
// coefficients rather than a polynomial type to simplify the
// NTL/LiDIA interface

#include "points.h"
#include "polys.h"
#include "divpol.h"

vector<bigint> makepdivpol(Curvedata* E, int p)
{
  if(p==2)
    {
      bigint b2,b4,b6,b8;
      E->getbi(b2,b4,b6,b8);
      vector<bigint> ans;
      ans.reserve(4);
      ans.push_back(b6);
      ans.push_back(2*b4);
      ans.push_back(b2);
      ans.push_back(BIGINT(4));
      return ans;
    }

  // default for odd p: use recursive method

  bigint a1,a2,a3,a4,a6;
  E->getai(a1,a2,a3,a4,a6);
  return div_pol_odd(a1,a2,a3,a4,a6,p);
}

// div_pol_odd(a1,a2,a3,a4,a6,n) returns the coefficients of the
// polynomial in x whose zeros are the (x-coordinates of the) non-zero
// points P on E=[a1,a2,a3,a4,a6] satisfying nP=0 (odd n)

// The poly itself is found recursively
 
ZPoly div_pol_odd_rec(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
		    const bigint& a6, int n);

vector<bigint> div_pol_odd(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
			   const bigint& a6, int n)
{
  ZPoly pol = div_pol_odd_rec(a1,a2,a3,a4,a6,n);
  int i, d = Degree(pol);
  vector<bigint> ans; ans.reserve(d+1);
  for(i=0; i<=d; i++) ans.push_back(PolyCoeff(pol,i));
  return ans;
}

ZPoly div_pol_odd_rec(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
		    const bigint& a6, int n)
{
  ZPoly X; ZPolySetX(X);
  ZPoly f1 = X*(X*(X+a2)+a4)+a6;
  ZPoly f2 = a1*X+a3;
  ZPoly psi24=(BIGINT(4)*f1+f2*f2); psi24*=psi24;
  ZPoly ans;
  switch(n) {
  case 0: 
    SetDegree(ans,0); SetCoeff(ans,0,0); return ans;
  case 1: case 2: 
    SetDegree(ans,0); SetCoeff(ans,0,1); return ans;
  case 3:
    SetDegree(ans,4);
    SetCoeff(ans,4,3);
    SetCoeff(ans,3,a1*a1+4*a2);
    SetCoeff(ans,2,3*a1*a3+6*a4);
    SetCoeff(ans,1,3*a3*a3+12*a6);
    SetCoeff(ans,0,a1*a1*a6-a1*a3*a4+a2*a3*a3+4*a2*a6-a4*a4);
    return ans;
  case 4:
    SetDegree(ans,6);
    SetCoeff(ans,6,2);
    SetCoeff(ans,5,a1*a1+4*a2);
    SetCoeff(ans,4,5*a1*a3+10*a4);
    SetCoeff(ans,3,10*a3*a3+40*a6);
    SetCoeff(ans,2,10*a1*a1*a6-10*a1*a3*a4+10*a2*a3*a3+40*a2*a6-10*a4*a4);
    SetCoeff(ans,1,a1*a1*a1*a1*a6-a1*a1*a1*a3*a4+a1*a1*a2*a3*a3+8*a1*a1*a2*a6-a1*a1*a4*a4-4*a1*a2*a3*a4-a1*a3*a3*a3-4*a1*a3*a6+4*a2*a2*a3*a3+16*a2*a2*a6-4*a2*a4*a4-2*a3*a3*a4-8*a4*a6);
    SetCoeff(ans,0,a1*a1*a1*a3*a6-a1*a1*a3*a3*a4+2*a1*a1*a4*a6+a1*a2*a3*a3*a3+4*a1*a2*a3*a6-3*a1*a3*a4*a4+2*a2*a3*a3*a4+8*a2*a4*a6-a3*a3*a3*a3-8*a3*a3*a6-2*a4*a4*a4-16*a6*a6);
    return ans;
  default: // general case, use recursion
    // If n is odd, n=2m+1:
    if(n%2==1)
      {
	int m=(n-1)/2;
	ZPoly t1=div_pol_odd_rec(a1,a2,a3,a4,a6,m);
	t1=div_pol_odd_rec(a1,a2,a3,a4,a6,m+2)*t1*t1*t1;
	ZPoly t2=div_pol_odd_rec(a1,a2,a3,a4,a6,m+1);
	t2=div_pol_odd_rec(a1,a2,a3,a4,a6,m-1)*t2*t2*t2;
	if(m%2==1) return t1-psi24*t2;
	return psi24*t1-t2;
      }
    else // n is even, n=2m:
      {
	int m=n/2;
	ZPoly t1=div_pol_odd_rec(a1,a2,a3,a4,a6,m-1);
	t1=div_pol_odd_rec(a1,a2,a3,a4,a6,m+2)*t1*t1;
	ZPoly t2=div_pol_odd_rec(a1,a2,a3,a4,a6,m+1);
	t2=div_pol_odd_rec(a1,a2,a3,a4,a6,m-2)*t2*t2;
	return div_pol_odd_rec(a1,a2,a3,a4,a6,m)*(t1-t2);
      }
  }
}

