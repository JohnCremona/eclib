// p2points.h:  declarations of P2Point class for points in P^2(Q)
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
 
// allow for multiple includes
#ifndef _P2POINT_
#define _P2POINT_

#include "bigrat.h"

//
// class for  points in P^2(Q), used as base for points on elliptic curves etc.
//

class P2Point{ 
  friend class Point;
  bigint X ; // homogeneous coordinates
  bigint Y ;
  bigint Z ;
  void reduce(); //divide out coordinate gcd
public:
  // constructors 
  P2Point(void)       // not a real point
    { X=0; Y=0; Z=0;} 
  P2Point(const bigint& x, const bigint& y, const bigint& z)
    : X(x), Y(y), Z(z)
    { reduce(); }
  P2Point(const bigint& x, const bigint& y, long z)
    : X(x), Y(y), Z(BIGINT(z))
    { reduce(); }
  P2Point(const bigint& x, long y, long z)
    : X(x), Y(BIGINT(y)), Z(BIGINT(z))
    { reduce(); }
  P2Point(long x, long y, long z)
    : X(BIGINT(x)), Y(BIGINT(y)), Z(BIGINT(z))
    { reduce(); }
  P2Point(const bigint& x, const bigint& y)
    : X(x), Y(y), Z(BIGINT(1))
    { ; } // no need to reduce 
  /* The following creates ambiguities owing to the bigint->bigrational coercion
  P2Point(const bigrational& x, const bigrational& y)
    : X(num(x)*den(y)), Y(num(y)*den(x)), Z(den(x)*den(y))
    { reduce(); }
  */
  P2Point(const P2Point& Q)
    :X(Q.X), Y(Q.Y), Z(Q.Z)
    { ; }
  ~P2Point(void) {;}
                
  // input and output
  friend inline ostream& operator<<(ostream & os, const P2Point& P)
    {return os << "[" << P.X << ":" << P.Y << ":" << P.Z << "]" ;}

  friend inline void output_pari(ostream&os, const P2Point& P)
    {
      bigint xp=P.X, yp=P.Y, zp=P.Z;
      if(is_zero(zp)) {os<<"[0]"; return;}
      if(is_one(zp)) {os<<"["<<xp<<","<<yp<<"]"; return;}
      bigint z=gcd(xp,zp);
      os<<"["<<(xp/z)<<"/"<<(zp/z)<<","<<yp<<"/"<<zp<<"]";
    }

// P2Point input: 3 formats allowed are 
// [x:y:z], [x,y], [x/z,y/z] with any stype of brackets
  friend istream& operator>>(istream & is, P2Point& P);

  // test of equality of points        
  int operator==(const P2Point& Q) const
  {
    return eq(*this,Q);
  }
  int operator!=(const P2Point& Q) const { return !(*this == Q); }
  friend int eq(const P2Point&P, const P2Point&Q);

  // assignment (p.init(.) is quicker than p=P2Point(.) for existing p)
  void init(const bigint& x, const bigint& y, const bigint& z)
  {X=x; Y=y; Z=z; reduce(); }
  void init(const bigint& x, const bigint& y)
  {X=x; Y=y; Z = BIGINT(1); }
  void operator=(const P2Point& Q) // P1 = P2
    { X=Q.X ; Y=Q.Y; Z=Q.Z;  }

  // Coordinate transforms useful for elliptic curve points 
  friend P2Point scale(const P2Point& P, const bigint& u, int back=0); 
  friend P2Point scale(const P2Point& P, long u=1, int back=0); 
  friend P2Point shift(const P2Point& P,
		       const bigint& r, const bigint& s, const bigint& t, 
		       int back=0); 
  friend P2Point transform(const P2Point& P,
			   const bigint& u, 
			   const bigint& r, const bigint& s, const bigint& t, 
			   int back=0); 

  void getcoordinates(bigint& x, bigint& y, bigint& z) const
    {x=X; y=Y; z=Z; }
  void getaffinecoordinates(bigrational& x, bigrational& y) const
  {x=bigrational(X,Z); y=bigrational(Y,Z); }
  void getrealcoordinates(bigfloat&x, bigfloat& y) const;
  friend inline bigint getX(const P2Point& p) {return p.X; }
  friend inline bigint getY(const P2Point& p) {return p.Y; }
  friend inline bigint getZ(const P2Point& p) {return p.Z; }
  int isintegral() const { return Z==BIGINT(1); }
  int isinfinite() const { return Z==BIGINT(0); }

}; // end of p2point class

// the real x and y coords of the point
inline void realify_point(const P2Point& P, bigfloat&x, bigfloat& y)
{
  P.getrealcoordinates(x,y);
}

// end of file: p2points.h

#endif
