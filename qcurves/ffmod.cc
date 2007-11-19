// ffmod.cc: implementation of class ffmodq and Weil pairing functions
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
 

// ffmodq is the function field of an elliptic curve mod a prime q
// (or more precisely the affine coordinate ring Fq[x,y])

#include "curve.h"
#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"

galois_field ffmodq::Fq;
curvemodq ffmodq::E;
FqPoly ffmodq::f1;
FqPoly ffmodq::f2;

//  special constructor to initialize the curve and field only:
ffmodq::ffmodq(const curvemodq& EE)
{
  E=EE; Fq=get_field(EE);
  //  cout<<"In ffmodq constructor"<<endl;
  init_f1f2();
}

void ffmodq::init_f1f2(void)
{
  //  cout<<"In ffmodq::init_f1f1()"<<endl;

  NewGF(Fq,a1);  NewGF(Fq,a2);  NewGF(Fq,a3);
  NewGF(Fq,a4);  NewGF(Fq,a6);
  E.get_ai(a1,a2,a3,a4,a6);
  // set f1, f2:
  NewFqPoly(Fq,X); 
  FqPolyAssignX(X);
  f1 = X*(X*(X+a2)+a4)+a6;
  f2 = a1*X+a3;     
}

int ffmodq::operator==(const ffmodq& b) const
{
  return (h1==b.h1) && (h2==b.h2);
}

ffmodq ffmodq::operator+(const ffmodq& b) const
{
  return ffmodq(h1+b.h1,h2+b.h2);
}

ffmodq ffmodq::operator-(const ffmodq& b) const
{
  return ffmodq(h1-b.h1,h2-b.h2);
}

ffmodq ffmodq::operator*(const ffmodq& b) const
{
  return ffmodq(h1*b.h1 + f1*h2*b.h2 , h1*b.h2 + h2*b.h1 - f2*h2*b.h2);
}

ffmodq ffmodq::operator*(const FqPoly& h) const
{
  return ffmodq(h*h1 , h*h2 );
}

ffmodq ffmodq::operator/(const FqPoly& h) const
{
  FqPoly g1 = h1; g1/=h;
  FqPoly g2 = h2; g2/=h;
  return ffmodq(g1, g2);
}

ffmodq ffmodq::operator/(const ffmodq& b) const
{
  if(Degree(b.h2)==-1) return (*this)/b.h1;
  cout<<"ffmodq error:  division by general elements not implemented!"<<endl;
  abort();
  return ffmodq();
}

// this is because of a one-time bug in gf_polynomial::operator():
gf_element evaluate(const FqPoly& f, const gf_element& value)
{
  int d = Degree(f);
  if (d == 0) return PolyCoeff(f,0);

  NewGF(GetField(f),result);
  GFSetZ(result,0);
	
  if (d < 0) return result;

  result = PolyCoeff(f,d);
  for (int i = d-1; i >= 0; i--) result = result*value + PolyCoeff(f,i);
  return result;
}

gf_element ffmodq::evaluate(const pointmodq& P) const
{
  if(P.is_zero()) 
    {cout<<"ffmodq error: attempt to evaluate at "<<P<<endl; abort();}
  //  cout<<"In ffmodq::operator() with this = "<<(*this)<<", P="<<P<<endl;
  gf_element x=P.get_x(), y=P.get_y();
  //  cout<<"x="<<x<<", y="<<y<<endl;
  //  cout<<"h1(x)="<<::evaluate(h1,x)<<endl;
  //  cout<<"h2(x)="<<::evaluate(h2,x)<<endl;
  //  cout<<"returning "<<::evaluate(h1,x)+y*::evaluate(h2,x)<<endl;
  return  ::evaluate(h1,x)+y*::evaluate(h2,x);
}

// vertical(P) has divisor (P)+(-P)-2(0)
ffmodq vertical(const pointmodq& P)
{
  if(P.is_zero()) 
    {
      ffmodq g(BIGINT(1)); return g;
    }
  NewFqPoly(base_field(P),h1);
  FqPolyAssignX(h1);
  return ffmodq(h1-(P.get_x()));
}

// tangent(P) has divisor 2(P)+(-2P)-3(0)
ffmodq tangent(const pointmodq& P)
{
  if(P.is_zero()) 
    {
      ffmodq g(BIGINT(1)); return g;
    }
  gf_element x=P.get_x(), y=P.get_y();
  gf_element a1,a2,a3,a4,a6;
  P.get_curve().get_ai(a1,a2,a3,a4,a6);
  gf_element dyf=y+y+a1*x+a3;
  NewFqPoly(base_field(P),h1);
  FqPolyAssignX(h1);
  // test for 2P=0:
  if(dyf==0) return ffmodq(h1-x);

  gf_element dxf=a1*y-(3*x*x+2*a2*x+a4);
  gf_element slope=-dxf/dyf;
  h1=-y-slope*(h1-x);
  NewFqPoly(base_field(P),h2);
  FqPolyAssign1(h2);
  return ffmodq(h1,h2);
}

//  chord between points:
// chord(P,Q) has divisor (P)+(Q)+(-P-Q)-3(O)
ffmodq chord(const pointmodq& P, const pointmodq& Q)
{
  if(P.is_zero()) return vertical(Q);
  if(Q.is_zero()) return vertical(P);

  const gf_element& xP = P.get_x();
  const gf_element& yP = P.get_y();
  const gf_element& xQ = Q.get_x();
  const gf_element& yQ = Q.get_y();
  gf_element ydiff=yP-yQ;
  gf_element xdiff=xP-xQ;
  if(xdiff==0)
    {
      if(ydiff==0)
	{
	  return tangent(P);
	}
      else
	{
	  return vertical(P);
	}    
    }
  gf_element slope=ydiff/xdiff;
  NewFqPoly(base_field(P),h1); FqPolyAssignX(h1);
  NewFqPoly(base_field(P),h2); FqPolyAssign1(h2);
  h1=-yP-slope*(h1-xP);
  return ffmodq(h1,h2);
}

ffmodq weil_pol(const pointmodq& T, int m)
{
  ffmodq h(BIGINT(1));
  switch(m) {
  case 2: return vertical(T);
  case 3: return tangent(T);
  default:
    {
      int k;
      pointmodq kT = T+T;
      h = tangent(T);
      for(k=2; k<m-1; k++, kT=kT+T)
	{
	  h=h*chord(T,kT);
	  h=h/vertical(kT);
	}
    }
  }
#if(0)
  // Check:
  FqPoly h1 = h.h1, h2=h.h2;
  FqPoly f1 = h.f1, f2=h.f2;
  FqPoly t = f1*h2*h2+f2*h1*h2-h1*h1;
  // that should equal const*(x-x(T))^m
  FqPoly u;
  power(u,vertical(T).h1,m);
  u = PolyCoeff(t,m)*u;
  if(t==u) 
    ;//	cout<<"weil_pol("<<T<<","<<m<<") = "<<h<<" checks OK"<<endl;
  else
    {
      cout<<"Error: weil_pol("<<T<<","<<m<<") = "<<h<<" fails to check"<<endl;
      cout<<"t = "<<t<<endl;
      cout<<"u = "<<u<<endl;
      abort();
    }
#endif
  return h;
}

void ffmodq::output(ostream& os) const
{
  os<<"("<<h1<<")+Y*("<<h2<<")";
}

  // Evaluation of Weil poly at another point S without actually
  // computing the polynomial -- only guaranteed to work if S is not a
  // multiple of T, so we check that m*S is not 0.  To avoid this,
  // evaluate at S+S' and at S' and divide, where m*(S+S') and m*S'
  // are nonzero!

// "Unsafe" version only called internally when S,T,and m*S are nonzero:

gf_element evaluate_weil_pol_0(const pointmodq& T, int m, const pointmodq& S)
{
  gf_element a = T.get_x(); // just to set the field
  GFSetZ(a,1);

  if(m==2) // easy case
    {
      return S.get_x()-T.get_x();
    }

  // We compute m*(T,1) = (m*T,am) = (0,answer) according to Frey-Ruck 
  // (kT,a) holds the current approximation, starting at (0,1)
  // (T2,a2) holds the repeated doubling of (T,1)

  pointmodq  kT = pointmodq(T.get_curve());
  pointmodq  T2 = T;
  gf_element a2 = a; // =1

  ffmodq h;
  int k=m;
  while(k)
    {
      //      cout<<"k = "<<k<<endl;
      if(k&1) // k odd, so add (T2,a2) to (kT,a)
	{
	  //	  cout<<"k = "<<k<<": odd case"<<endl;
	  k--;
	  if(kT.is_zero()) // at the beginning
	    {
	      kT=T2;
	      a =a2;
	    }
	  else
	    {
	      h = chord(kT,T2);
	      a *= a2;
	      a *=  h.evaluate(S);
	      kT+=T2;
	      if(!kT.is_zero())
		{
		  h = vertical(kT);
		  a /=  h.evaluate(S);
		}
	    }
	  //	  cout<<"(kT,a) = ("<<kT<<","<<a<<")"<<endl;
	}
      // now double (T2,a2) unless finished
      if(k)
	{
	  k/=2;
	  a2 *= a2;
	  h = tangent(T2);
	  a2 *=  h.evaluate(S);
	  T2=2*T2;
	  if(!T2.is_zero())
	    {
	      h = vertical(T2);
	      a2 /=  h.evaluate(S);
	    }
	  //	  cout<<"After doubling, (T2,a2) = ("<<T2<<","<<a2<<")"<<endl;
	}
    }
  //  cout<<"At end, (kT,a) = ("<<kT<<","<<a<<")"<<endl;
  return a;
}

gf_element evaluate_weil_pol(const pointmodq& T, int m, const pointmodq& S)
{
  gf_element a = T.get_x(); // just to set the field
  GFSetZ(a,1);
  if(T.is_zero()||S.is_zero()) return a;

  if(!(m*S).is_zero()) 
    return evaluate_weil_pol_0(T,m,S);

  pointmodq R=T.get_curve().random_point();
  while((m*R).is_zero())
    R=T.get_curve().random_point();
  return evaluate_weil_pol_0(T,m,R+S)/evaluate_weil_pol_0(T,m,R);
}

gf_element weil_pairing(const pointmodq& S, const pointmodq& T, int m)
{
  gf_element a = T.get_x(); // just to set the field
  GFSetZ(a,0); // for return on error condition 
  // cout<<"Evaluating Weil Pairing of order "<<m<<" on "<<S<<" and "<<T<<endl;
  if(!(m*T).is_zero()) {cout<<"error in Weil pairing of "<<S<<" and "<<T<<" and order "<<m<<": m*T is not 0"<<endl; abort(); return a;}
  if(!(m*S).is_zero()) {cout<<"error in Weil pairing of "<<S<<" and "<<T<<" and order "<<m<<": m*S is not 0"<<endl; abort(); return a;}

  GFSetZ(a,1); // for return if trivial
  if(T.is_zero()||S.is_zero()) return a;
  if(T==S) return a;
  if(m==2) return -a;

  // now m>=3

  ffmodq fT = weil_pol(T, m);
  //  cout<<"T = "<<T<<", fT = "<<fT<<endl;
  ffmodq fS = weil_pol(S, m);
  //  cout<<"S = "<<S<<", fS = "<<fS<<endl;

  pointmodq R1 = T.get_curve().random_point();
  pointmodq R2 = T.get_curve().random_point();
  pointmodq R3=R1-R2;
  pointmodq ST=S-T;
  while(R3.is_zero()||(R3==-T)||(R3==S)||(R3==ST))
    {
      R2=T.get_curve().random_point();
      R3=R1-R2;
    }
  //  cout<<"numerator   = "<<(fT.evaluate(S-R3)*fS.evaluate(R3))<<endl;
  //  cout<<"denominator = "<<(fS.evaluate(T+R3)*fT.evaluate(-R3))<<endl;
  return (fT.evaluate(S-R3)*fS.evaluate(R3)) / 
         (fS.evaluate(T+R3)*fT.evaluate(-R3));
  //  return evaluate_weil_pol(T,m,S)/evaluate_weil_pol(S,m,T);
}
