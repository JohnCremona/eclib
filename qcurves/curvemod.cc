// curvemod.cc: implementation of class curvemodq for elliptic curve mod q 
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
 

// NB Under LiDIA this is just an interface to LiDIA's
// elliptic_curve<gf_element>;  under NTL it is self-contained

#include "curve.h"
#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)

curvemodq reduce_curve(const Curve& E, const bigint& q)
{
  galois_field gf(q);
  //  cout<<"In reduce_curve(), q="<<q<<", field = "<<gf<<  endl;
  gf_element A1(gf),A2(gf),A3(gf),A4(gf),A6(gf);
  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);
  A1.assign(a1);  A2.assign(a2);  A3.assign(a3);  
  A4.assign(a4);  A6.assign(a6);
  curvemodq c(A1,A2,A3,A4,A6);
  //  cout<<"Returning "<<c<<endl;
  return c;
}

#else // NTL

void curvemodq::set_group_order_via_legendre()  
{
  // Do NOT make these static as the modulus might change!
  gf_element two=to_ZZ_p(2);
  gf_element four=two+two;
  if(!is_zero(order)) return; // order already set!
  order=BIGINT(1);  // point at infinity
  gf_element b2 = a1*a1 + four*a2; 
  gf_element b4 = two*a4 + a1*a3;
  gf_element b6 = a3*a3 + four*a6; 
  bigint ix; NewGF(*Fq,x); NewGF(*Fq,d);
  for(ix=BIGINT(0); ix<q; ix++)
    {
      GFSetZ(x,ix);
      d = ((four*x+b2)*x+(two*b4))*x+b6;
      order+=(1+legendre(rep(d),q));
    }
}

void curvemodq::set_group_order()  
{
  if(((this->q)<100)||((this->q)==181)||((this->q)==331)||((this->q)==547))
    {
      set_group_order_via_legendre();
      return;
    }
  pointmodq P1, P2;
  bigint n1, n2, n;
  my_isomorphism_type(*this,n1,n2,P1,P2);
  order=n1*n2;
}

#endif // end of LiDIA/NTL split

// Division poly functions:


FqPoly makepdivpol(const curvemodq& C, int p)
{
  if(p==2)
    {
      gf_element a1,a2,a3,a4,a6;
      C.get_ai(a1,a2,a3,a4,a6);
      NewFqPoly(get_field(C),f);
      SetDegree(f,3);
      SetCoeff(f,0,a3*a3 + 4*a6);
      SetCoeff(f,1,2*(2*a4 + a1*a3));
      SetCoeff(f,2,a1*a1 + 4*a2);
      SetCoeff(f,3,ItoGF(get_field(C),4));
      return f;
    }

  // default for odd p: use recursive method

  return div_pol_odd(C,p);
}

// div_pol_odd(C,n) returns the polynomial in x whose zeros are the
// (x-coordinates of the) non-zero points P on C satisfying nP=0 (for
// odd n)

// The poly itself is found recursively
 
FqPoly div_pol_odd_rec(const curvemodq& C, int n);

FqPoly div_pol_odd(const curvemodq& C, int n)
{
  return div_pol_odd_rec(C,n);
}

FqPoly div_pol_odd_rec(const curvemodq& C, int n)
{
  const galois_field Fq=get_field(C);
  NewFqPoly(Fq,X);  FqPolyAssignX(X);
  gf_element a1,a2,a3,a4,a6;
  C.get_ai(a1,a2,a3,a4,a6);
  FqPoly f1 = X*(X*(X+a2)+a4)+a6;
  FqPoly f2 = a1*X+a3;
  FqPoly psi24=(ItoGF(Fq,4)*f1+f2*f2); psi24*=psi24;
  NewFqPoly(Fq,ans);
  switch(n) {
  case 0: 
    FqPolyAssign0(ans); return ans;
  case 1: case 2: 
    FqPolyAssign1(ans); return ans;
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
	FqPoly t1=div_pol_odd_rec(C,m);
	t1=div_pol_odd_rec(C,m+2)*t1*t1*t1;
	FqPoly t2=div_pol_odd_rec(C,m+1);
	t2=div_pol_odd_rec(C,m-1)*t2*t2*t2;
	if(m%2==1) return t1-psi24*t2;
	return psi24*t1-t2;
      }
    else // n is even, n=2m:
      {
	int m=n/2;
	FqPoly t1=div_pol_odd_rec(C,m-1);
	t1=div_pol_odd_rec(C,m+2)*t1*t1;
	FqPoly t2=div_pol_odd_rec(C,m+1);
	t2=div_pol_odd_rec(C,m-2)*t2*t2;
	return div_pol_odd_rec(C,m)*(t1-t2);
      }
  }
}

