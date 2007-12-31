// p2points.cc:  implementations for P2Point class for points in P^2(Q)
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
 
#include "p2points.h"
#ifdef NTL_INTS
#include <NTL/RR.h>   // for the realify_point() function
#endif

//
// P2Point member functions
//

void P2Point::reduce(void)
{
  if(Z==1) return;   // integral point, no work needed
  bigint d = gcd(gcd(X, Y), Z);
  if(sign(d)==0){
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

// P2Point input: 3 formats allowed are 
// [x:y:z], [x,y], [x/z,y/z] with any type of brackets

istream& operator>>(istream & is, P2Point& P)
{
  char c; 
  is>>c;  // swallow first bracket
  bigint x,y,dx,dy;
  is >> x >> c; // swallow comma or colon
  switch(c) {
  case ',':
    P.X=x; is >> P.Y >> c; P.Z=BIGINT(1); break;
  case '/': is >> dx >> c >> y >> c >> dy >> c; 
    P.X=x*dy; P.Y=y*dx; P.Z=dx*dy; break;
  case ':': P.X=x; is >> P.Y >> c >> P.Z >> c; break;
  default: P.X=P.Y=P.Z=BIGINT(0); // null point
  }
  P.reduce(); 
  return is; 
}

// test of equality of points
int eq(const P2Point&P, const P2Point&Q)
{
  if(sign(P.X*Q.Z-P.Z*Q.X)) return 0;
  if(sign(P.Y*Q.Z-P.Z*Q.Y)) return 0;
  if(sign(P.Y*Q.X-P.X*Q.Y)) return 0;
  return 1;
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

void P2Point::getrealcoordinates(bigfloat&x, bigfloat& y) const
#ifdef NTL_INTS
{
  RR zp=to_RR(Z);  
#ifdef NTL_ALL
  x=to_RR(X)/zp;
  y=to_RR(Y)/zp;
#else
  x=to_double(to_RR(X)/zp);
  y=to_double(to_RR(Y)/zp);
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
  bigrational xp, yp;  getaffinecoordinates(xp,yp);
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

// Coordinate transforms useful for elliptic curve points 
P2Point scale(const P2Point& P, const bigint& u, int back)
{
  if(u==BIGINT(1)) return P;
  bigint u2=u*u;
  bigint u3=u*u2;
  if(back) 
    return P2Point(u2*P.X,u3*P.Y,P.Z);
  else
    return P2Point(u*P.X,P.Y,u3*P.Z);
} 

P2Point scale(const P2Point& P, long u, int back)
{
  if(u==1) return P;
  return scale(P,BIGINT(u),back);
}

P2Point shift(const P2Point& P,
	      const bigint& r, const bigint& s, const bigint& t, 
	      int back)
{
  if(back)
    return P2Point(P.X+r*P.Z, P.Y+s*P.X+t*P.Z, P.Z);
  else
    return P2Point(P.X-r*P.Z, P.Y-s*P.X+(r*s-t)*P.Z, P.Z);
}

P2Point transform(const P2Point& P,
		  const bigint& u, 
		  const bigint& r, const bigint& s, const bigint& t, 
		  int back)
{
  if(back)
    return shift(scale(P,u,1),r,s,t,1);
  else
    return scale(shift(P,r,s,t,0),u,0);
}

// end of file: p2points.cc
