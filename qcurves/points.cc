// points.cc:  implementations for Point class for points on elliptic curves
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
 
// originally adapted from Elliptic.cc by Oisin McGuiness

#include "points.h"  // which includes curve.h
#include "cperiods.h"
#ifdef NTL_INTS
#include <NTL/RR.h>   // for the realify_point() function
#endif

//#define DEBUG_TORSION

//
// Point member functions
//

void Point::reduce(void)
{
  if(sign(Z)==0){  // point at infinity
    X = 0 ;
    Y = 1 ;
    return ;
  }
  if(Z==1){  // integral point, no work needed
    return ;
  }
  bigint d = gcd(X, Y) ; d = gcd(d, Z) ; // now d = gcd(X, Y, Z)
  if(sign(d)==0){
    cerr << "##Point::reduce:  Bad point: gcd == 0 \n" ;
    X = 0 ; Y = 1 ; Z = 0;
    abort() ;
    return ;
  }
  if(d!=1) {
    X /=  d ;
    Y /=  d ;
    Z /=  d ;
  }
  if(sign(Z)<0){
    ::negate(X) ;
    ::negate(Y) ;
    ::negate(Z) ;
  }
}

// Point input: 3 formats allowed are 
// [x:y:z], [x,y], [x/z,y/z] with any type of brackets

istream& operator>>(istream & is, Point& P)
{
  char c; 
  //  is.flags( is.flags() | ios::dec );  //force decimal input (bug fix)
  is>>c;
  //  cout<<"First char read = "<<c<<"\n";
  bigint x,y,z,z2,z3;
  is >> x >> c;
  //  cout<<"Next char read = "<<c<<"\n";
  switch(c) {
  case ',': is >> y >> c;  P.X=x; P.Y=y; P.Z=1; break;
  case '/': is >> z2 >> c >> y >> c >> z3 >> c; z=z3/z2; 
    //    cout<<"x="<<x<<", z2="<<z2<<", y="<<y<<", z3="<<z3<<", z="<<z<<endl;
    P.X=x*z; P.Y=y; P.Z=z3; break;
  case ':': P.X=x; is >> P.Y >> c >> P.Z >> c; break;
  case ']': default: P.X=P.Y=P.Z=0; // null point
  }

  P.ord = 0; P.height = -1.0; if(P.isvalid()) P.reduce(); return is; 
}

// test of equality of points
int eq(const Point&P, const Point&Q) // assumes same curve
{
  bigint temp = P.X*Q.Z-P.Z*Q.X;
  int ans = (sign(temp)==0);
  if(!ans) return 0;
  temp = P.Y*Q.Z-P.Z*Q.Y;
  ans = (sign(temp)==0); 
  if(!ans) return 0;
  temp = P.Y*Q.X-P.X*Q.Y;
  ans = (sign(temp)==0); 
  return ans;
}

// Shifting coords: 
// N.B. (1) This changes the curve too
//      (2) This uses coords (x/z^2,y/z^3)

void raw_shift_point(const bigint& x, const bigint& y, const bigint& z, 
                const bigint& u, const bigint& r, 
                const bigint& s, const bigint& t, 
                bigint& newx, bigint& newy, bigint& newz, int back)
{
  //cout<<"x="<<x<<", y="<<y<<", z="<<z<<endl;
  //cout<<"\t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
  const bigint& z2 = z*z;
  const bigint& z3 = z*z2;
  if(back)
    {
      const bigint& u2 = u*u;
      const bigint& u3 = u*u2;
      newx = u2*x + r*z2;
      newy = u3*y + u2*s*x*z + t*z3;
      newz = z;
    }
  else
    {
      newx = x - r*z2;
      newy = y - s*x*z + (r*s-t)*z3;
      newz = u*z;
    }
  //cout<<"newx="<<newx<<", newy="<<newy<<", newz="<<newz<<endl;

  const bigint& newz2 = newz*newz;
  const bigint& v2 = gcd(newx, newz2);
  const bigint& v3 = gcd(newy, newz*newz2);
  if(sign(v2)!=0) {newx/=v2;   newy/=v3;   newz/=(v3/v2);}
  //cout<<"newx="<<newx<<", newy="<<newy<<", newz="<<newz<<endl;
}

// Shifts P via T = (u, r, s, t) on to newc.   
// N.B. Assumes that (1) T(P.E) = newc
//                   (2) T(P) is on newc
// "back" means use reverse transform

Point shift(const Point& P, Curvedata* newc, 
	    const bigint& u, const bigint& r, const bigint& s, 
	    const bigint& t, int back)
{
  if(P.iszero()) return Point(newc);
  if(!P.isvalid())
    cout << "Attempting to shift the point " << P 
	 << "which is not a valid point on its curve " << P.getcurve() << "!\n";
  const bigint& oldz = gcd(P.X, P.Z);
  const bigint& oldx = P.X/oldz;
  const bigint& oldy = P.Y;
  bigint newx, newy, newz;
  raw_shift_point(oldx, oldy, oldz, u, r, s, t, newx, newy, newz, back);
  
  Point ans(newc, newx*newz, newy, pow(newz,3));
  if(!ans.isvalid())
    {
      cout << "Result of shifting the point " << P;
      cout << " is " << ans << " which is not a valid point on its curve " 
	   << ans.getcurve() << "!\n";
    }
  return ans;
}

// addition and subtraction
void Point::operator+=(const Point& Q)  // P += Q; 
{
  Point sum = (*this) + Q ;
  E = sum.E; X = sum.X; Y = sum.Y; Z = sum.Z; reduce();
  ord = 0; height = -1;
}

void Point::operator-=(const Point& Q)  // P -= Q ; 
{
  Point sum = -Q;  sum+=(*this);
  E = sum.E; X = sum.X; Y = sum.Y; Z = sum.Z; reduce();
  ord = 0; height = -1;
}
        
Point Point::operator+(const Point& Q) const // P + Q
{
  Point ans(E);
  // make sure that both points are on the same curve
  if(E != Q.E){
    cerr << "## Can't add points on different curves!" << "\n" ;
    return ans;
  }
  // NOTE: we don't do any reductions here, but rely on assignment to do it
  // check special cases first
  if(Z==0) return Q ;
  if(Q.iszero()) return *this ;
  if(eq(*this , Q))      {return Q.twice();}
  Point minusQ=-Q;
  if(eq(*this , minusQ)) {return ans;}       // zero
        
  // we now have genuine work to do
  // let's set up some local variables to avoid repeated references
  // coefficients
  bigint A1,A2,A3,A4,A6;  E->getai(A1,A2,A3,A4,A6);
  // coordinates of P
  const bigint& X1 = X ;
  const bigint& Y1 = Y ;
  const bigint& Z1 = Z ;
  // coordinates of Q
  const bigint& X2 = Q.X ;
  const bigint& Y2 = Q.Y ;
  const bigint& Z2 = Q.Z ;
  const bigint& Z12 = Z1 * Z2 ;
        
  const bigint& L = - Y2 * Z1 + Y1 * Z2 ;    /* lambda */
  const bigint& M = - X2 * Z1 + X1 * Z2 ;    /* mu    */
  const bigint& N = - Y1 * X2 + Y2 * X1 ;    /* nu   */
        
  const bigint& Mz =  M * M * Z12  ; 
  
  const bigint& t = L * L * Z12 + M * ( A1 * L * Z12 - M * (A2 * Z12
                                        + X1 * Z2 + X2 * Z1 ) )  ;
  const bigint& newX =  M * t ;
  const bigint& newY = - (  t * (L + A1 * M) +   Mz * (N + A3 * M )) ;    
  const bigint& newZ = M * Mz ;
  ans.init(E, newX, newY, newZ);
  return ans;
}

Point Point::operator-(const Point& Q) const // P - Q
{
  Point ans = (-Q) ;
  ans += (*this);
  return ans;
}

Point Point::twice(void) const // doubles P
{
  // do trivial cases
  Point ans(E);
  if( Z==0) return ans;
  bigint A1,A2,A3,A4,A6;  E->getai(A1,A2,A3,A4,A6);
  Point minusthis = -(*this);
  if(eq(*this,minusthis)) return ans;   // order 2

  const bigint& Zsq = Z * Z ;
  
  const bigint& L = 3 * X * X  +  2 * A2 * X * Z  +  
                                     A4 * Zsq  -  A1 * Y * Z ;
  const bigint& M = 2 * Y  +  A1 * X  +  A3 * Z ;
        
  const bigint& Mz = M * Z ;
  const bigint& N =  - X * X * X -A3 * Y * Zsq  + A4 * X *Zsq  +  
    2 * A6 * Z *Zsq ;
        
  const bigint& t  = L * L  +  Mz * ( A1 * L  -  M * ( A2 * Z  + 2 * X ) );  

  const bigint& newX = t * Mz ;
  const bigint& newY = - ( L * t +  Mz * ( A1 * t+ M * ( N  + A3 * Mz * Z) ) );
  const bigint& newZ = Mz * Mz * Mz;

  ans.init(E, newX, newY, newZ) ;
  return ans;
}

Point Point::operator-(void) const // -Q
{
  Point negation(*this);
  negation.Y = -Y - (E->a3)*Z - (E->a1)*X;
  return negation;
}

// calculates nP
Point operator*(int n, const Point& P)
{
  Point ans(*(P.E));
  if(P.iszero() || n == 0) return ans;
  int negative = (n < 0) ;
  if(negative) n = - n ;
  if(n == 1) {
    ans = P;
    if(negative) ans =  -P;
    return ans ;
  }
  // now n >= 2
  if(n == 2){
    ans = P.twice() ;
    if(negative) ans =  -ans;
    return ans ;
  }
  // now n >= 3
  if(n&1) ans = P ;    // else ans is ZERO
  Point  Q = P ; 
  while(n > 1){
    Q = Q.twice() ; // 2^k P
    n /= 2 ;
    if(n&1) ans = ans + Q ;
  }
  if(negative) ans = -ans ;
  return ans ;
}

int order(Point& p)
{
  // ASSUME that point is valid; check before calling if unknown
  if (p.ord) {return p.ord;}
  bigint eight, z=p.Z; eight=8;
  if (is_zero(z))  {p.ord = 1; return 1; }
  if (z>eight)     {p.ord =-1; return -1;}
  Point q = p;  long ord=1;
  while ( (sign(z)!=0) && (z<=eight) ) // (worst denom of a torsion point)
    {ord++;   q+=p;  z = q.Z;}
  if (sign(z)!=0) ord = -1;
  p.ord = ord;
  return ord;
}

int order(Point& p, vector<Point>& multiples)
{
#ifdef DEBUG_TORSION
  cout<<"In order() with p = "<<p<<endl;
#endif
  // ASSUME that point is valid; check before calling if unknown
  multiples.clear(); multiples.reserve(13);
  multiples.push_back(Point(*(p.E)));
  if (p.iszero()) {p.ord=1; return 1; }
  multiples.push_back(p);
  Point q = p;
  bigint eight; eight=8;
  while ( (!q.iszero()) && (q.Z<=eight) 
	                 && (multiples.size()<13) ) // 12 is max poss order
  { 
    q+=p;
    if (!q.iszero()) multiples.push_back(q);
  }
  if (q.iszero()) 
    p.ord = multiples.size();
  else
    p.ord = -1;
#ifdef DEBUG_TORSION
  cout<<"Returning ord="<<p.ord<<", multiples = "<<multiples<<endl;
#endif
  return p.ord;
}

int Point::isvalid() const // P on its curve ?
{
  if(E == 0){
    cerr << "## Bad point: null curve pointer!\n" ;
    return 0 ;
  }
  if((sign(X)==0)&&(sign(Y)==0)&&(sign(Z)==0)) return 0;   // Null point, useful for terminating input
  if(sign(Z)==0) return 1 ;                    // Point at infinity is on a curve
  else{
    // should calculate 
    //                       y^2 +a1 x y + a3 y 
    //  and
    //                       x^3 + a2 x^2 + a4 x + a6
    // where 
    //          x = X/Z, y = Y/Z
    // and verify equality.
    //
    // In homogeneous coordinates:
    //
    // Lhs = Y^2 Z + a1 XYZ + a3 YZ^2 = (YZ) *(Y + a1 X + a3 Z)
    //
    //
    // Rhs = X^3 +a2 X^2 Z + a4 X Z^2 + a6 Z^3
    //
    bigint A1,A2,A3,A4,A6;    E->getai(A1,A2,A3,A4,A6);
    const bigint& Lhs = Y*Z*(Y + A1*X + A3*Z) ;
    const bigint& Rhs = A6*pow(Z,3) + X*(A4*Z*Z + X*(A2*Z + X)) ;
    return Lhs == Rhs ;
  }
}

//#define DEBUG_REALIFY

// the real x and y coords of the point
// 
// There are different version for the various arithmetic options so
// that when we are not using multiprecision floating point, we do not
// overflow when the homogeneous coordinates are large
//
// NTL:  temporarily use RR type
// LiDIA without m.p.: temporarily use bigrational type
// LiDIA with m.p.: just do it

void realify_point(const Point& P, bigfloat&x, bigfloat& y)
#ifdef NTL_INTS
{
  RR zp=to_RR(getZ(P));  
#ifdef NTL_ALL
  x=to_RR(getX(P))/zp;
  y=to_RR(getY(P))/zp;
#else
  x=to_double(to_RR(getX(P))/zp);
  y=to_double(to_RR(getY(P))/zp);
#endif
#ifdef DEBUG_REALIFY
  cout<<"realifying P="<<P<<" (NTL version)"<<endl;
  cout<<"Real point = ("<<x<<","<<y<<")"<<endl;
#endif
#ifndef NTL_ALL
  if((abs(x)==INFINITY)||(abs(y)==INFINITY))
    {
      cout<<"After converting P="<<P<<" to doubles, ";
      cout<<"Real point = ("<<x<<","<<y<<")"<<endl;
      cout<<"insufficient precision to continue, aborting"<<endl;
      abort();
    }
#endif
}
#else // not NTL...
#ifdef LiDIA_INTS
#ifndef LiDIA_ALL
{
  bigrational xp(getX(P),getZ(P)), yp(getY(P),getZ(P));
  x=dbl(xp);
  y=dbl(yp);
#ifdef DEBUG_REALIFY
  cout<<"P="<<P<<endl;
  cout<<"Rational point = ("<<xp<<","<<yp<<")"<<endl;
  cout<<"Real point = ("<<x<<","<<y<<")"<<endl;
#endif
  if((abs(x)==INFINITY)||(abs(y)==INFINITY))
    {
      cout<<"After converting P="<<P<<" to doubles, ";
      cout<<"Real point = ("<<x<<","<<y<<")"<<endl;
      cout<<"insufficient precision to continue, aborting"<<endl;
      abort();
    }
}
#else // LiDIA multiprecision
{
  bigfloat z = I2bigfloat(getZ(P));
  x = I2bigfloat(getX(P))/z;
  y = I2bigfloat(getY(P))/z;
#ifdef DEBUG_REALIFY
  cout<<"P="<<P<<endl;
  cout<<"Real point = ("<<x<<","<<y<<")"<<endl;
#endif
}
#endif
#endif
#endif              

// find all the torsion points on a curve (Curvedata member function)

vector<Point> old_torsion_points(Curvedata& E);  // code is below

long Curvedata::get_ntorsion()
{
  if(ntorsion==0)
    {
#ifdef DEBUG_TORSION
      cout<<"Calling torsion_points() on "<<(Curve)(*this)<<endl;
#endif
      vector<Point> ans = torsion_points(*this);
      ntorsion = ans.size();
#ifdef DEBUG_TORSION
      cout<<"torsion_points() returns "<<ans<<" of size "<<ntorsion<<endl;
#endif
    }
#ifdef DEBUG_TORSION
      cout<<"Curvedata::get_ntorsion() returns "<<ntorsion<<endl;
#endif
  return ntorsion;
}

// Cperiods is a class containing a basis for the period lattice L;
// it knows how to compute points from z mod L; so this function
// effectively does the same as PARI's pointell() (called ellzto point
// in the next PARI release)
//
Point make_tor_pt(Curvedata& E, Cperiods& per, 
		  const bigfloat& ra1, const bigfloat& ra2, const bigfloat& ra3, 
		  const bigcomplex& z)
{
  bigcomplex cx,cy;
#ifdef DEBUG_TORSION
  cout<<"In make_tor_pt() with z = " << z << endl;
#endif
  per.XY_coords(cx,cy,z);
  cx = cx-(ra1*ra1+4*ra2)/to_bigfloat(12);
  cy = (cy - ra1*cx - ra3)/to_bigfloat(2);
#ifdef DEBUG_TORSION
  cout<<"(x,y) = ("<<(cx)<<","<<(cy)<<")\n";
#endif
  bigint x=Iround(real(cx)), y=Iround(real(cy));
  Point P(E, x, y);
  return P;
}


vector<Point> two_torsion(Curvedata& E)
{
#ifdef DEBUG_TORSION
  cout<<"\nIn two_torsion() with curve "<<(Curve)E<<"\n";
#endif
  bigint a1, a2, a3, a4, a6, b2, b4, b6, b8;
  E.getai(a1,a2,a3,a4,a6);
  E.getbi(b2,b4,b6,b8);
  int scaled=0;
  if (odd(a1) || odd(a3)) 
    { 
      b4*=BIGINT(8);
      b6*=BIGINT(16);
      scaled=1;
    }
  else 
    {
      b2=a2; b4=a4; b6=a6;
    }
  vector<bigint> xlist = Introotscubic(b2,b4,b6);
  int n, n2tors = xlist.size();

  // If there are 3 points of order 2, we order them for consistency:

  bigint ei;
  if(n2tors==3)
    {
      if(abs(xlist[0])>abs(xlist[1])) 
	{ei=xlist[0];xlist[0]=xlist[1];xlist[1]=ei;}
      if(abs(xlist[0])>abs(xlist[2])) 
	{ei=xlist[0];xlist[0]=xlist[2];xlist[2]=ei;}
      if(abs(xlist[1])>abs(xlist[2])) 
	{ei=xlist[1];xlist[1]=xlist[2];xlist[2]=ei;}
    }
  //  sort(xlist.begin(),xlist.end()); would be unambiguous!
  vector<Point> two_tors;
  two_tors.push_back(Point(E)) ;     // zero point
  for(n=0; n<n2tors; n++)
    {
      ei = xlist[n];
      if(scaled) 
	two_tors.push_back(Point(E,2*ei,-a1*ei-4*a3,BIGINT(8)));
      else 
	two_tors.push_back(Point(E,ei,BIGINT(0),BIGINT(1)));
    }
#ifdef DEBUG_TORSION
  cout<<"\ntwo_torsion() returning "<<two_tors<<"\n";
#endif
  return two_tors;
}

// Returns vector of x values such that x/3 is the rational x-coord of
// a point in E[3], possibly with y not rational.  If y is rational
// then these x will in fact be multiples of 3 since rational
// 3-torsion is integral
vector<bigint> three_torsion_x(Curvedata& E)
{
#ifdef DEBUG_TORSION
  cout<<"\nIn three_torsion_x() with curve "<<(Curve)E<<"\n";
#endif
  bigint b2, b4, b6, b8;
  E.getbi(b2,b4,b6,b8);
  vector<bigint> xlist = Introotsquartic(b2,9*b4,27*b6,27*b8);
  // NB None of the implementation of Introosquartic() in marith.cc
  // fix the order of the roots, which depends on either (1) the order
  // of the factors in NTL's Z[X] factorization routine, or (2) the
  // order of the roots returned by LiDIA's complex roots() function.
  // HENCE the order of the 3-torsion point x-coordinates (when there
  // are two) is not well-defined without the sorting done here.
#ifdef DEBUG_TORSION
  cout<<"\nthree_torsion_x() finds unsorted xlist =  "<<xlist<<"\n";
#endif
  if(xlist.size()==2)
    {
      sort(xlist.begin(),xlist.end());
#ifdef DEBUG_TORSION
      cout<<"\nthree_torsion_x() returns sorted xlist =  "<<xlist<<"\n";
#endif
    }
  return xlist;
}

vector<Point> three_torsion(Curvedata& E)
{
#ifdef DEBUG_TORSION
  cout<<"\nIn three_torsion() with curve "<<(Curve)E<<"\n";
#endif
  bigint a1, a2, a3, a4, a6, b2, b4, b6, b8, xi, d, rd;
  E.getai(a1,a2,a3,a4,a6);
  E.getbi(b2,b4,b6,b8);
  vector<bigint> xlist = three_torsion_x(E);
  vector<Point> three_tors;
  three_tors.push_back(Point(E)) ;     // zero point
  for(unsigned int n=0; n<xlist.size(); n++)
    {
      xi = xlist[n];
      if(xi%3==0) // 3-torsion must be integral
	{
	  xi/=3;
	  d = ((4*xi+b2)*xi+(2*b4))*xi+b6;
	  if(isqrt(d,rd)) 			    
	    {
	      Point P(E,2*xi,rd-(a1*xi+a3),BIGINT(2));
	      three_tors.push_back(P);
	      three_tors.push_back(-P);
	    }
	}
    }
#ifdef DEBUG_TORSION
  cout<<"\nthree_torsion() returning "<<three_tors<<"\n";
#endif
  return three_tors;
}

// New torsion point routine using Mazur's theorem to limit possibilities 
// and computing possible real torsion from the period lattice.  Suggestion 
// of Darrin Doud.
//
vector<Point> torsion_points(Curvedata& E)  // After Darrin Doud, adapted by JC
{
  if ( E.isnull() ) return vector<Point>(0);
//
// table[i][] contains a list of possible maximal orders for a point,
// given that the 2-torsion subgroup has order i
//
  static long table[5][5] = {{},{5,7,9,3},{12,6,8,4,10},{},{8,6,4}};
  static long nt[5] = {0,4,5,0,3};
#ifdef DEBUG_TORSION
  cout<<"\nIn torsion_points() with curve "<<(Curve)E<<"\n";
#endif
  bigint a1,a2,a3,a4,a6,sa2,sa4,sa6, x, y ;
  E.getai(a1,a2,a3,a4,a6);
  bigfloat ra1=I2bigfloat(a1), ra2=I2bigfloat(a2), ra3=I2bigfloat(a3);
  long i,j,ntp=1,nt2;
  vector<Point> points;
  vector<Point> cycle;
  Cperiods per(E);  bigcomplex w1,w2,z,z2; 
  per.getwRI(w1,w2); z2=w2/to_bigfloat(2);
  int ncc = getconncomp(E);
  int found;
#ifdef DEBUG_TORSION
  cout<<"Periods: "<<per<<"\nReal Period = "<<w1<<endl;
  cout<<ncc<<" real component(s)"<<endl;
#endif
  points.push_back(Point(E)) ;       // zero point
  Point p, q;

  // We find the two-torsion algebraically
  
  vector<Point> two_tors = two_torsion(E);
  nt2=two_tors.size();
#ifdef DEBUG_TORSION
  cout<<"Size of 2-torsion subgroup = " << nt2 << endl;
  cout << two_tors << endl;
#endif

  for(i=0, found=0; (i<nt[nt2])&&(!found); i++)
    {
      long ni=table[nt2][i]; 
#ifdef DEBUG_TORSION
      cout<<"Looking for a point of order " << ni << "\n";
#endif
      if(ni==3)
	{
	  p=Point(E);
	  vector<Point> p3 = three_torsion(E);
	  if(p3.size()>1) {p=p3[1];}
	}
      else
	{
	  z=w1/to_bigfloat(ni);
	  p = make_tor_pt(E,per,ra1,ra2,ra3,z);
	}
#ifdef DEBUG_TORSION
      cout<<"p = " << p <<"?\n";
#endif
      found=(p.isvalid())&&(order(p,cycle)==ni);
      if(!found&&(ncc==2)&&even(ni))
	{
	  p = make_tor_pt(E,per,ra1,ra2,ra3,z+z2);
#ifdef DEBUG_TORSION
	  cout<<"p = " << p <<"?\n";
#endif
	  found=(p.isvalid())&&(order(p,cycle)==ni);
	  if(!found&&((ni%4)==2))
	    {
	      p = make_tor_pt(E,per,ra1,ra2,ra3,z+z+z2);
#ifdef DEBUG_TORSION
	      cout<<"p = " << p <<"?\n";
#endif
	      found=(p.isvalid())&&(order(p,cycle)==ni);
	    }
	}
      if(found)
	{
	  ntp=ni;
	  points=cycle;
#ifdef DEBUG_TORSION
	  cout<<"Found a point " << p << " of order " << ni << endl;
	  cout<<"generating subgroup "<<cycle<<endl;
#endif
	}
#ifdef DEBUG_TORSION
      else
	{
	  cout<<"none\n";
	}
#endif
    }

#ifdef DEBUG_TORSION
  cout<<"Number of points in cyclic part = " << ntp << endl;
  cout<<points<<endl;
#endif
  
  if(ntp==1) {return two_tors;}                   // C1, C2 or C2xC2
  if(nt2==4)
  // non-cyclic case, C2xC4, C2xC6 or C2xC8
  // Find a point of order 2 not in the second factor and add it in
    {
      for(i=1; i<4; i++) // 0 is first point in two_tors !
	{
	  p=two_tors[i];
	  if(find(points.begin(),points.end(),p)==points.end())
	    {
#ifdef DEBUG_TORSION
	      cout<<"Using 2-torsion point "<<p<<" as coset rep\n";
#endif
	      for(j=0; j<ntp; j++)
		{
		  points.push_back(points[j]+p);
		}
	      ntp*=2;
	      break;
	    }
#ifdef DEBUG_TORSION
	  cout<<"2-torsion point "<<p<<" is in subgroup\n";
#endif
	}
    }
  return points;
}

vector<Point> old_torsion_points(Curvedata& E)
{
  if ( E.isnull() ) return vector<Point>(0);
  bigint a1,a2,a3,a4,a6,sa2,sa4,sa6, d, x, y ;
  E.getai(a1,a2,a3,a4,a6);
  long i,nroots,ntp;  int scaled_flag;
  vector<Point> points;
  ntp = 1;
  points.push_back(Point(E)) ;     // zero point
  if ( (sign(a1)==0) && (sign(a3)==0) )
    {sa2=a2; sa4=a4; sa6=a6; scaled_flag=0; }
  else
    {sa2 = a1*a1 + 4*a2;
     sa4 = 8*a1*a3 + 16*a4;
     sa6 = 16*a3*a3 + 64*a6;
     scaled_flag=1; 
    }
  d = sa2*sa2*(sa4*sa4 - 4*sa2*sa6) + 18*sa2*sa4*sa6
                    - 4*sa4*sa4*sa4 - 27*sa6*sa6;
  Point p;

// First test y=0 for points of order 2:

  vector<bigint> xlist = Introotscubic(sa2, sa4, sa6);
  nroots = xlist.size();
  for (i=0; i<nroots; i++)
    {
      x = xlist[i];
      if (scaled_flag) p.init(E, 2*x, - a1*x - 4*a3, BIGINT(8));
      else p.init(E, x, BIGINT(0));
      points.push_back(p);
    }

// Now test y such that y^2 divides d:

  vector<bigint> possible_y( sqdivs(d) );
  vector<bigint>::iterator yvar=possible_y.begin();
  while( yvar!=possible_y.end())
  {
    y = *yvar++;
    xlist = Introotscubic(sa2, sa4, sa6-y*y);
    nroots = xlist.size();
    for (i=0; i<nroots; i++)
    {
      x = xlist[i];
      if (scaled_flag) p.init(E, 2*x, y - a1*x - 4*a3, BIGINT(8));
      else p.init(E, x, y);
      if (order(p) > 0)
        {
          points.push_back(p);
          points.push_back(-p); // N.B. order>2 here!
        }
    }
  }
  return points;
}

// end of file: points.cc
