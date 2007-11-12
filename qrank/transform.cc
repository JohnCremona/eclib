// transform.cc: implementation of quartic transformation functions
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
 
//
// Notation: g(x,z) is replaced by g(m11*x+m12*z,m21*x+m22*z)/m00^2
//

#include "marith.h"
#include "unimod.h"
#include "transform.h"

bigint g_content(const bigint& ga, const bigint& gb, const bigint& gc, 
		 const bigint& gd, const bigint& ge)
     // returns largest SQUARE which divides all
{
  // first find the content:
  bigint ans=abs(ga);
  if(ans==1) return ans;
  ans=gcd(ans,gb);
  if(ans==1) return ans;
  ans=gcd(ans,gc);
  if(ans==1) return ans;
  ans=gcd(ans,gd);
  if(ans==1) return ans;
  ans=gcd(ans,ge);
  if(ans==1) return ans;
  // if content non-trivial, get its divisors whose square divides
  // (as we alreasy have this function to hand...)
  vector<bigint> cdivs = sqdivs(ans);
  // and return the last in the list, which is the biggest:
  return cdivs[(cdivs.size())-1];
}

void apply_transform(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
		     const unimod& m)
{
  bigint m11=m(1,1), m12=m(1,2), m21=m(2,1), m22=m(2,2);
  bigint m112=sqr(m11); bigint m113=m112*m11; bigint m114=m113*m11;
  bigint m212=sqr(m21); bigint m213=m212*m21; bigint m214=m213*m21;
  bigint m222=sqr(m22); bigint m223=m222*m22; bigint m224=m223*m22;
  bigint m122=sqr(m12); bigint m123=m122*m12; bigint m124=m123*m12;

  bigint newa = m214*e + m11*m213*d + m112*m212*c + m113*m21*b + m114*a;
  bigint newe = m224*e + m12*m223*d + m122*m222*c + m123*m22*b + m124*a;
  bigint newb = 4*m213*m22*e + (3*m11*m212*m22+m12*m213)*d
    + 2*(m112*m21*m22+m11*m12*m212) * c
    + (3*m112*m12*m21+m113*m22)*b + 4*m113*m12*a;
  bigint newd = 4*m21*m223*e + (3*m12*m21*m222+m11*m223)*d
    + 2*(m122*m21*m22+m11*m12*m222)*c
    + (m123*m21+ 3*m11*m122*m22)*b + 4*m11*m123*a;
  bigint newc = 6*m212*m222*e + 3*(m12*m212*m22+m11*m21*m222) * d
    + (m122*m212+ 4*m11*m12*m21*m22+m112*m222) * c
    + 3*(m11*m122*m21+m112*m12*m22) * b + 6*m112*m122*a;

  a=newa; b=newb; c=newc; d=newd; e=newe;
}

void apply_transform(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
		     const scaled_unimod& m)
{
  apply_transform(a,b,c,d,e,(unimod)m);
  bigint u2=sqr(m.scale_factor());
  if(u2>1)
    {
      divide_exact(a,u2,a);
      divide_exact(b,u2,b);
      divide_exact(c,u2,c);
      divide_exact(d,u2,d);
      divide_exact(e,u2,e);
    }
}

void xshift(const bigint& alpha,
	    bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	    unimod& m)
{
  e = e+alpha*(d+alpha*(  c+alpha*(  b+  alpha*a)));
  d =          d+alpha*(2*c+alpha*(3*b+4*alpha*a));
  c =                     c+alpha*(3*b+6*alpha*a);
  b =                                b+4*alpha*a;
//a =                                          a;  
  m.x_shift(alpha);
}

void zshift(const bigint& gamma,
	    bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	    unimod& m)
{
  a = a+gamma*(b+gamma*(  c+gamma*(  d+  gamma*e)));
  b =          b+gamma*(2*c+gamma*(3*d+4*gamma*e));
  c =                     c+gamma*(3*d+6*gamma*e);
  d =                                d+4*gamma*e;
//e =                                          e;
  m.y_shift(gamma);
}

void m_invert(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      unimod& m)
{
  swap(a,e); swap(b,d); ::negate(b); ::negate(d);
  m.invert();
}

void m_invert(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      scaled_unimod& m)
{
  swap(a,e); swap(b,d); ::negate(b); ::negate(d);
  m.invert();
}

int check_transform(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& d, const bigint& e,
		    const unimod& m,
		    const bigint& xa, const bigint& xb, const bigint& xc, 
		    const bigint& xd, const bigint& xe)
{
  bigint aa(a), bb(b), cc(c), dd(d), ee(e);
  apply_transform(aa,bb,cc,dd,ee,m);
  if(aa!=xa) return 0;
  if(bb!=xb) return 0;
  if(cc!=xc) return 0;
  if(dd!=xd) return 0;
  if(ee!=xe) return 0;
  return 1;
}

int check_transform(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& d, const bigint& e,
		    const scaled_unimod& m,
		    const bigint& xa, const bigint& xb, const bigint& xc, 
		    const bigint& xd, const bigint& xe)
{
  bigint aa(a), bb(b), cc(c), dd(d), ee(e);
  apply_transform(aa,bb,cc,dd,ee,m);
  if(aa!=xa) return 0;
  if(bb!=xb) return 0;
  if(cc!=xc) return 0;
  if(dd!=xd) return 0;
  if(ee!=xe) return 0;
  return 1;
}

