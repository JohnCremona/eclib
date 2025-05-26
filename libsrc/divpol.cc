// divpol.cc: implementations of functions for division polynomials
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
 
// NB the external interface uses a simple vector<bigint> of
// coefficients rather than a polynomial type to simplify the
// NTL interface

#include <eclib/polys.h>
#include <eclib/divpol.h>

// The 2-divison polynomial (cubic in x)

ZPoly div_pol_2(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
                const bigint& a6)
{
  ZPoly ans;
  SetDegree(ans,3);
  SetCoeff(ans,3,4);
  SetCoeff(ans,2,a1*a1+4*a2);
  SetCoeff(ans,1,2*a1*a3+4*a4);
  SetCoeff(ans,0,a3*a3+4*a6);
  return ans;
}

// div_pol_odd(a1,a2,a3,a4,a6,n) returns the coefficients of the
// polynomial in x whose zeros are the (x-coordinates of the) non-zero
// points P on E=[a1,a2,a3,a4,a6] satisfying nP=0 (odd n)

// The poly itself is found recursively

ZPoly div_pol_odd(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
                  const bigint& a6, int n)
{
  static const bigint four(4);
  ZPoly X; ZPolySetX(X);
  ZPoly f1 = X*(X*(X+a2)+a4)+a6;
  ZPoly f2 = a1*X+a3;
  ZPoly psi24=(four*f1+f2*f2); psi24*=psi24;
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
	ZPoly t1=div_pol_odd(a1,a2,a3,a4,a6,m);
	t1=div_pol_odd(a1,a2,a3,a4,a6,m+2)*t1*t1*t1;
	ZPoly t2=div_pol_odd(a1,a2,a3,a4,a6,m+1);
	t2=div_pol_odd(a1,a2,a3,a4,a6,m-1)*t2*t2*t2;
	if(m%2==1) return t1-psi24*t2;
	return psi24*t1-t2;
      }
    else // n is even, n=2m:
      {
	int m=n/2;
	ZPoly t1=div_pol_odd(a1,a2,a3,a4,a6,m-1);
	t1=div_pol_odd(a1,a2,a3,a4,a6,m+2)*t1*t1;
	ZPoly t2=div_pol_odd(a1,a2,a3,a4,a6,m+1);
	t2=div_pol_odd(a1,a2,a3,a4,a6,m-2)*t2*t2;
	return div_pol_odd(a1,a2,a3,a4,a6,m)*(t1-t2);
      }
  }
}

ZPoly div_pol(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
              const bigint& a6,int n)
{
  return (n==2? div_pol_2(a1,a2,a3,a4,a6) :  div_pol_odd(a1,a2,a3,a4,a6,n));
}

ZPoly division_polynomial(Curvedata* E, int p)
{
  bigint a1,a2,a3,a4,a6;
  E->getai(a1,a2,a3,a4,a6);
  return (p==2? div_pol_2(a1,a2,a3,a4,a6) :  div_pol_odd(a1,a2,a3,a4,a6,p));
}

// Numerator of the multiplication-by-n map on the x-coordinate

ZPoly mul_by_n_num(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
                   const bigint& a6, int n)
{
  ZPoly X; ZPolySetX(X);
  ZPoly P_2 = div_pol_2(a1,a2,a3,a4,a6);
  ZPoly P_n = div_pol_odd(a1,a2,a3,a4,a6,n);
  ZPoly P_nplus1 = div_pol_odd(a1,a2,a3,a4,a6,n+1);
  ZPoly P_nminus1 = div_pol_odd(a1,a2,a3,a4,a6,n-1);
  ZPoly A = X * P_n * P_n;
  ZPoly B = P_nminus1 * P_nplus1;
  return (n%2==0? A * P_2 - B : A - P_2 * B);
}

// Denominator of the multiplication-by-n map on the x-coordinate

ZPoly mul_by_n_den(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
                   const bigint& a6, int n)
{
  ZPoly P_n = div_pol_odd(a1,a2,a3,a4,a6,n);
  return (n%2==1? P_n * P_n : P_n * P_n * div_pol_2(a1,a2,a3,a4,a6));
}

// Polynomial whose roots are x(Q) for Q satisfying n*Q=P, where x(P)=xP/zP

ZPoly division_points_X_pol(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
                            const bigint& a6,
                            int n,
                            const bigint& xP, const bigint& zP)
{
  ZPoly numpoly = mul_by_n_num(a1, a2, a3, a4, a6, n);
  ZPoly denpoly = mul_by_n_den(a1, a2, a3, a4, a6, n);
  //cout << "numpoly = " << numpoly << endl;
  //cout << "denpoly = " << numpoly << endl;
  return zP * numpoly - xP * denpoly;
}
