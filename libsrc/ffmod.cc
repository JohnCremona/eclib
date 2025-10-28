// ffmod.cc: implementation of class ffmodq and Weil pairing functions
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
 

// ffmodq is the function field of an elliptic curve mod a prime q
// (or more precisely the affine coordinate ring Fq[x,y])

#include "eclib/ffmod.h"

galois_field ffmodq::Fq;
curvemodq ffmodq::E;
ZZ_pX ffmodq::f1;
ZZ_pX ffmodq::f2;

void ffmodq::init(const curvemodq& EE)
{
  E = EE;
  Fq = get_field(EE);
  ZZ_p a1, a2, a3, a4, a6;
  E.get_ai(a1,a2,a3,a4,a6);
  // set f1, f2:
  ZZ_pX X;
  SetX(X);
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

ffmodq ffmodq::operator*(const ZZ_pX& h) const
{
  return ffmodq(h*h1 , h*h2 );
}

ffmodq ffmodq::operator/(const ZZ_pX& h) const
{
  ZZ_pX g1 = h1; g1/=h;
  ZZ_pX g2 = h2; g2/=h;
  return ffmodq(g1, g2);
}

ffmodq ffmodq::operator/(const ffmodq& b) const
{
  if(deg(b.h2)==-1) return (*this)/b.h1;
  cerr<<"ffmodq error:  division by general elements not implemented!"<<endl;
  return ffmodq();
}

// this is because of a one-time bug in gf_polynomial::operator():
ZZ_p evaluate(const ZZ_pX& f, const ZZ_p& value)
{
  int d = deg(f);
  if (d == 0) return coeff(f,0);

  ZZ_p result = to_ZZ_p(0);

  if (d < 0) return result;

  result = coeff(f,d);
  for (int i = d-1; i >= 0; i--) result = result*value + coeff(f,i);
  return result;
}

ZZ_p ffmodq::evaluate(const pointmodq& P) const
{
  //  cout<<"In ffmodq::operator() with this = "<<(*this)<<", P="<<P<<endl;
  ZZ_p x=P.get_x(), y=P.get_y();
  if(P.is_zero()) 
    {cerr<<"ffmodq error: attempt to evaluate at "<<P<<endl; return x;}
  //  cout<<"x="<<x<<", y="<<y<<endl;
  //  cout<<"h1(x)="<<::evaluate(h1,x)<<endl;
  //  cout<<"h2(x)="<<::evaluate(h2,x)<<endl;
  //  cout<<"returning "<<::evaluate(h1,x)+y*::evaluate(h2,x)<<endl;
  return  ::evaluate(h1,x)+y*::evaluate(h2,x);
}

// vertical(P) has divisor (P)+(-P)-2(0)
ffmodq vertical(const pointmodq& P)
{
  static const ZZ one(1);
  if(P.is_zero())
    {
      ffmodq g(one);
      return g;
    }
  ZZ_pX h1;
  SetX(h1);
  return ffmodq(h1-(P.get_x()));
}

// tangent(P) has divisor 2(P)+(-2P)-3(0)
ffmodq tangent(const pointmodq& P)
{
  static const ZZ one(1);
  if(P.is_zero())
    {
      ffmodq g(one);
      return g;
    }
  ZZ_p x=P.get_x(), y=P.get_y();
  ZZ_p a1,a2,a3,a4,a6;
  P.get_curve().get_ai(a1,a2,a3,a4,a6);
  ZZ_p dyf=y+y+a1*x+a3;
  ZZ_pX h1;
  SetX(h1);
  // test for 2P=0:
  if(dyf==0) return ffmodq(h1-x);

  ZZ_p dxf=a1*y-(3*x*x+2*a2*x+a4);
  ZZ_p slope=-dxf/dyf;
  h1=-y-slope*(h1-x);
  ZZ_pX h2(1);
  return ffmodq(h1,h2);
}

//  chord between points:
// chord(P,Q) has divisor (P)+(Q)+(-P-Q)-3(O)
ffmodq chord(const pointmodq& P, const pointmodq& Q)
{
  if(P.is_zero()) return vertical(Q);
  if(Q.is_zero()) return vertical(P);

  const ZZ_p& xP = P.get_x();
  const ZZ_p& yP = P.get_y();
  const ZZ_p& xQ = Q.get_x();
  const ZZ_p& yQ = Q.get_y();
  ZZ_p ydiff=yP-yQ;
  ZZ_p xdiff=xP-xQ;
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
  ZZ_p slope=ydiff/xdiff;
  ZZ_pX h1;
  SetX(h1);
  ZZ_pX h2(1);
  h1=-yP-slope*(h1-xP);
  return ffmodq(h1,h2);
}

ffmodq weil_pol(const pointmodq& T, int m)
{
  static const ZZ one(1);
  ffmodq h(one);
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
  ZZ_pX h1 = h.h1, h2=h.h2;
  ZZ_pX f1 = h.f1, f2=h.f2;
  ZZ_pX t = f1*h2*h2+f2*h1*h2-h1*h1;
  // that should equal const*(x-x(T))^m
  ZZ_pX u;
  power(u,vertical(T).h1,m);
  u = coeff(t,m)*u;
  if(t==u) 
    ;//	cout<<"weil_pol("<<T<<","<<m<<") = "<<h<<" checks OK"<<endl;
  else
    {
      cerr<<"Error: weil_pol("<<T<<","<<m<<") = "<<h<<" fails to check"<<endl;
      cerr<<"t = "<<t<<endl;
      cerr<<"u = "<<u<<endl;
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

ZZ_p evaluate_weil_pol_0(const pointmodq& T, int m, const pointmodq& S)
{
  if(m==2) // easy case
    {
      return S.get_x()-T.get_x();
    }

  // We compute m*(T,1) = (m*T,am) = (0,answer) according to Frey-Ruck 
  // (kT,a) holds the current approximation, starting at (0,1)
  // (T2,a2) holds the repeated doubling of (T,1)

  pointmodq  kT = pointmodq(T.get_curve());
  pointmodq  T2 = T;
  ZZ_p a = to_ZZ_p(1);
  ZZ_p a2 = a;

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

ZZ_p evaluate_weil_pol(const pointmodq& T, int m, const pointmodq& S)
{
  ZZ_p a = to_ZZ_p(1);
  if(T.is_zero()||S.is_zero()) return a;

  if(!(m*S).is_zero())
    return evaluate_weil_pol_0(T,m,S);

  pointmodq R=T.get_curve().random_point();
  while((m*R).is_zero())
    R=T.get_curve().random_point();
  return evaluate_weil_pol_0(T,m,R+S)/evaluate_weil_pol_0(T,m,R);
}

ZZ_p weil_pairing(const pointmodq& S, const pointmodq& T, int m)
{
  ZZ_p a = to_ZZ_p(0); // for return on error condition
  // cout<<"Evaluating Weil Pairing of order "<<m<<" on "<<S<<" and "<<T<<endl;
  if(!(m*T).is_zero())
    {
      cerr<<"error in Weil pairing of "<<S<<" and "<<T<<" and order "<<m<<": m*T is not 0"<<endl;
      return a;
    }
  if(!(m*S).is_zero())
    {
      cerr<<"error in Weil pairing of "<<S<<" and "<<T<<" and order "<<m<<": m*S is not 0"<<endl;
      return a;
    }

  a = to_ZZ_p(1); // for return if trivial
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
