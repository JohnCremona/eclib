// elog.cc: implementations of elliptic logarithm functions
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
 
#include "compproc.h"
#include "elog.h"

bigfloat ssqrt(const bigfloat& x)
{
  if(x<0) 
    {
      cout<<"Attempts to take real square root of "<<x<<endl; 
      return to_bigfloat(0);
    }
  return sqrt(x);
}

void boundedratapprox(bigfloat x, bigint& a, bigint& b, const bigint& maxden);


// Given an elliptic curve and its (precomputed) periods, and a point
// P=(x,y), returns the unique complex number z such that 

// (1) \wp(z)=x+b2/12, \wp'(z)=2y+a1*x+a3, 

// (2) either z is real and 0\le z\lt w1, or Delta>0, z-w2/2 is real
// and 0\le z-w2/2\le w1.  

// Here, [w1,w2] is the standard period lattice basis

// c.f. Cohen page 399

//#define DEBUG_ELOG

bigcomplex ellpointtoz(const Curvedata& E, const Cperiods& per, const bigfloat& x, const bigfloat& y)
{
  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);
  bigfloat ra1=I2bigfloat(a1);
  bigfloat ra2=I2bigfloat(a2);
  bigfloat ra3=I2bigfloat(a3);
  bigfloat xP(x), yP(y);
  int posdisc = (sign(getdiscr(E))>0);
  
  bigcomplex e1,e2,e3;
  getei(E,e1,e2,e3);
  if (posdisc) reorder1(e1,e2,e3);  else reorder2(e1,e2,e3);
  bigfloat re1=real(e1);

  bigcomplex w1,w2;
  per.getwRI(w1,w2);
#ifdef DEBUG_ELOG
  cout<<"w1 = "<<w1<<endl;
  cout<<"w2 = "<<w2<<endl;
#endif

  if(posdisc) // all roots real, e1>e2>e3
    {
#ifdef DEBUG_ELOG
      cout<<"positive discriminant"<<endl; 
#endif
      bigfloat re2=real(e2);
      bigfloat re3=real(e3);
#ifdef DEBUG_ELOG
      cout<<"Real roots, should be in descending order:\n"
      	  <<re1<<"\n"<<re2<<"\n"<<re3<<endl;
#endif
      bigfloat a1, a = sqrt(re1-re3);
      bigfloat b1, b = sqrt(re1-re2);
      bigfloat c;

      int egg = (xP<re1); // if P is on the "egg", replace it by P+T3
			  // where T3=(e3,y3) is a 2-torsion point on
			  // the egg coming from w2/2 on the lattice
      if(egg) 
	{
	  bigfloat y3 = -(ra1*re3+ra3)/2;
	  bigfloat lambda=(yP-y3)/(xP-re3);
	  bigfloat xP3 = lambda*(lambda+ra1)-ra2-xP-re3;
	  yP = lambda*(xP-xP3)-yP-ra1*xP3-ra3;
	  xP = xP3;
#ifdef DEBUG_ELOG
	  cout<<"Point on egg, replacing by ";
	  cout<<" ("<<xP<<","<<yP<<")"<<endl;
#endif
	}
      c = sqrt(xP-re3);    
      
      while(!is_approx_zero((a-b)/a))
	{
	  a1=(a+b)/2;
	  b1=sqrt(a*b);
	  c=(c+sqrt(c*c+b*b-a*a))/2;
	  a=a1; b=b1;
	}
#ifdef DEBUG_ELOG
      cout<<"After AGM loop, |a-b|=("<<abs(a-b)<<endl;
#endif

      bigcomplex z(asin(a/c)/a);
#ifdef DEBUG_ELOG
      cout<<"Basic z = "<<z<<endl;
#endif
      if((2*yP+ra1*xP+ra3)>0) 
	{
	  z = w1-z;
#ifdef DEBUG_ELOG
	  cout<<"(adjusted) point in upper half, replacing z by "<<z<<endl;
#endif
	}
      if( egg )
	{
	  z = z + w2/to_bigfloat(2);
#ifdef DEBUG_ELOG
	  cout<<"adding half imaginary period since point was on egg, now z = "<<z<<endl;
#endif
	}
      return z;
    }
  else // negative disc
    {
#ifdef DEBUG_ELOG
      cout<<"negative discriminant"<<endl;
      cout<<"Real root = " <<re1<<endl;
#endif
      // Here we use formulae equivalent to those in COhen, but better
      // behaved when roots are close together!
      bigcomplex zz = sqrt(e1-e2);
      bigfloat beta = abs(e1-e2);
      bigfloat a1, b1, a = 2*abs(zz), b = 2*real(zz);
      bigfloat c = (xP-re1+beta)/sqrt(xP-re1);
#ifdef DEBUG_ELOG
      cout<<"a,b,c = "<<a<<", "<<b<<", "<<c<<endl;
#endif
      while(!is_approx_zero((a-b)/a))
	{
	  a1=(a+b)/2;
	  b1=sqrt(a*b);
	  c=(c+sqrt(c*c+b*b-a*a))/2;
	  a=a1; b=b1;
#ifdef DEBUG_ELOG
	  cout<<"a,b,(a-b)/a = "<<a<<", "<<b<<", "<<(a-b)/a<<endl;
#endif
	}
      bigfloat z = asin(a/c);
#ifdef DEBUG_ELOG
      cout<<"Basic z = "<<z<<endl;
#endif
      bigfloat w = (2*yP+ra1*xP+ra3);
      if(w*((xP-re1)*(xP-re1)-beta*beta) >= 0)
	{
	  z=Pi()-z;
#ifdef DEBUG_ELOG
	  cout<<"After first adjustment, z = "<<z<<endl;
#endif
	}
      z/=a;
      if(w>0) 
	{
#ifdef DEBUG_ELOG
	  cout<<"After second adjustment, z = "<<z<<endl;
#endif
	  z+=(Pi()/a);
	}
      return bigcomplex(z);
    }
}


//#define DEBUG_EZP

// Cperiods is a class containing a basis for the period lattice L;
// it knows how to compute points from z mod L; so this function
// effectively does the same as PARI's ellztopoint()
//
// First function:  returns x,y as complex numbers

vector<bigcomplex> ellztopoint(Curvedata& E,  Cperiods& per, const bigcomplex& z)
{
  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);
  bigfloat ra1=I2bigfloat(a1);
  bigfloat ra2=I2bigfloat(a2);
  bigfloat ra3=I2bigfloat(a3);
  bigcomplex cx,cy;
  Cperiods per2 = per;  // since XY_coords changes the normalization
  per2.XY_coords(cx,cy,z);
  cx = cx-(ra1*ra1+4*ra2)/to_bigfloat(12);
  cy = (cy - ra1*cx - ra3)/to_bigfloat(2);
#ifdef DEBUG_EZP
  cout<<"In ellztopoint() with E = "<<(Curve)E<<endl;
  cout<<"periods = "<<per2<<endl;
  cout<<"z       = "<<z<<endl;
  cout<<"point = ("<<cx<<","<<cy<<")"<<endl;
#endif
  vector<bigcomplex> ans;
  ans.push_back(cx);
  ans.push_back(cy);
  return ans;
}

// Second function, expects to return a rational point.
// User supplies a denominator for the point; if it doesn't work, the
// Point returned is 0 on the curve

Point ellztopoint(Curvedata& E, Cperiods& per, const bigcomplex& z, const bigint& den)
{
  if(is_zero(z)) {return Point(E);}
  vector<bigcomplex> CP = ellztopoint(E,per,z);
  bigcomplex cx=CP[0],cy=CP[1];
  bigint nx,ny,dx,dy;
  boundedratapprox(real(cx),nx,dx,den);
  boundedratapprox(real(cy),ny,dy,den);
#ifdef DEBUG_EZP
  cout<<"Rounded x = "<<nx<<"/"<<dx<<endl;
  cout<<"Rounded y = "<<ny<<"/"<<dy<<endl;
#endif
  Point P(E, nx*dy, ny*dx, dx*dy);
  if(P.isvalid()) 
    {
#ifdef DEBUG_EZP
      cout<<"ellztopoint returning valid point "<<P<<endl; 
#endif
      return P;
    }
  return Point(E);
}

// Returns a (possibly empty) vector of solutions to m*Q=P

vector<Point> division_points(Curvedata& E,  const Point& P, int m)
{
  Cperiods cp(E);
  return division_points(E,cp,P,m);
}

//#define DEBUG_DIVPT

vector<Point> division_points(Curvedata& E,  Cperiods& per, const Point& P, int m)
{
#ifdef DEBUG_DIVPT
  cout<<"division_points("<<(Curve)E<<","<<P<<","<<m<<")"<<endl;  
#endif
  vector<Point> ans;
  if(m==0) 
    {
      cout<<"division_points() called with m=0!"<<endl; 
      return ans;
    }
  if(m<0) m=-m;

  bigcomplex w1, w2;
  per.getwRI(w1,w2);
  Cperiods per2 = per;  // since XY_coords changes the normalization

  int posdisc = (sign(getdiscr(E))>0);
  int k, egg;
  Point Q(E);

  bigcomplex z(to_bigfloat(0)), w;
  bigint den;
  int zero_flag = P.iszero();
  if(zero_flag)  
    {
      den=BIGINT(1);
      if(even(m)) ans=two_torsion(E); // computed algebraically
      else        ans.push_back(P);   // (more robust)
    }
  else 
    {
      z = elliptic_logarithm(E,per2,P);
      den=getZ(P);
    }
#ifdef DEBUG_DIVPT
  cout<<"posdisc=  "<<posdisc<<endl;
  cout<<"zero_flag="<<zero_flag<<endl;
  cout<<"den=      "<<den<<endl;
#endif

  if(posdisc)
    {
      egg = !is_real(z);
      bigcomplex half_w2 = w2/to_bigfloat(2);
#ifdef DEBUG_DIVPT
  cout<<"egg_flag= "<<egg<<endl;
#endif
      if(egg)  // P is on the "egg"
	{
	  if(even(m)) return ans;  // no solutions!
	  for(k=0; k<m; k++)       // now m is odd, Q on egg too
	    {
	      w = real(z+to_bigfloat(k)*w1)/m + half_w2;
	      Q = ellztopoint(E,per2,w,den);
	      if(!Q.iszero()
		 &&(m*Q==P)
		 &&(find(ans.begin(),ans.end(),Q)==ans.end()))
		ans.push_back(Q);	      
	    }
	}
      else //  P is on the connected component
	{
	  for(k=0; k<m; k++)
	    {
	      if(zero_flag&&((k==0)||(2*k==m))) 
		continue; // already have 2-torsion
	      w = real(z+to_bigfloat(k)*w1)/m;
	      if((k>0)||(!zero_flag))
		{
		  Q = ellztopoint(E,per2,w,den);
		  if(!Q.iszero()
		     &&(m*Q==P)
		     &&(find(ans.begin(),ans.end(),Q)==ans.end()))
		    ans.push_back(Q);	      
		}
	      if(even(m))
		{
		  Q = ellztopoint(E,per2,w+half_w2,den);
		  if(!Q.iszero()
		     &&(m*Q==P)
		     &&(find(ans.begin(),ans.end(),Q)==ans.end()))
		    ans.push_back(Q);	      
		}
	    }
	}
    }
  else // negative discriminant (so z is real)
    {
      for(k=0; k<m; k++)
	{
	  if(zero_flag&&((k==0)||(2*k==m))) 
	    continue; // already have 2-torsion
	  w = real(z+to_bigfloat(k)*w1)/m;
	  Q = ellztopoint(E,per2,w,den);
	  if(!Q.iszero()
	     &&(m*Q==P)
	     &&(find(ans.begin(),ans.end(),Q)==ans.end()))
	    ans.push_back(Q);	      
	}      
    }
  return ans;
}

// Returns a vector of solutions to m*Q=0 (including Q=0)

// First version will compute the Cperiods itself, so best to use the
// second one if more than one call is to be made for the same curve

vector<Point> torsion_points(Curvedata& E,int m)
{
  Cperiods cp(E);
  return torsion_points(E,cp,m);
}

vector<Point> torsion_points(Curvedata& E,  Cperiods& per, int m)
{
  Point P(E);
  return division_points(E,per,P,m);
}

void boundedratapprox(bigfloat x, bigint& a, bigint& b, const bigint& maxden)
{
  //  cout<<"bounded ratapprox of "<<x<<" (maxden = "<<maxden<<")"<<endl;
  bigint c, x0, x1, x2, y0, y1, y2;
  bigfloat rc, xx, diff, eps = to_bigfloat(1.0e-6);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
  diff = 1; c=x2=y2=0;
  while ((abs(y2)<maxden)&&!is_approx_zero(diff)) // ( diff > eps )
    { c = Iround( xx ); rc=I2bigfloat(c);
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      diff = abs( x - I2bigfloat(x2)/I2bigfloat(y2) );
      //      cout<<"x2 = "<<x2<<",\ty2 = "<<y2<<",\tdiff = "<<diff<<endl;
      if ( abs(xx - rc) < eps ) diff = 0;
      else xx = 1/(xx - rc);
    }
  a = x2; b = y2;
  if ( b < 0 )
    {::negate(a); ::negate(b); }
}

