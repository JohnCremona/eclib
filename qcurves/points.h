// points.h:  declarations of Point class for points on elliptic curves
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
 
// originally adapted from Elliptic.h by Oisin McGuiness

// allow for multiple includes
#ifndef _PELLIPTIC_
#define _PELLIPTIC_

#include "curve.h"

//
// class for  points on elliptic curves
//

class Point{ 
  bigint X ; // homogeneous coordinates
  bigint Y ;
  bigint Z ;
  Curvedata *E; // pointer to the curve that the point is on
                   // warning: no check is made that client doesn't alter curve
  int ord; // 0 if not calculated yet, -1 if infinite
  bigfloat height; // -1.0 if not calculated yet, 0.0 for torsion point
  void reduce(); //divide out coordinate gcd
public:
  // constructors 
  Point(void)       // not a real point
    : E(0), ord(0), height(to_bigfloat(-1.0))
    { X=0; Y=0; Z=0;} 
  Point(Curvedata &EE)      // set to point at infinity
    : E(&EE), ord(1), height(to_bigfloat(0.0))
    { X=0; Y=1; Z=0;}      
  Point(Curvedata *EE)      // set to point at infinity
    : E(EE), ord(1), height(to_bigfloat(0.0))
    { X=0; Y=1; Z=0;}      
  Point(Curvedata &EE, const bigint& x, const bigint& y, const bigint& z)
    : X(x), Y(y), Z(z), E(&EE), ord(0), height(to_bigfloat(-1.0))
    { reduce(); }
  Point(Curvedata *EE, const bigint& x, const bigint& y, const bigint& z)
    : X(x), Y(y), Z(z), E(EE), ord(0), height(to_bigfloat(-1.0))
    { reduce(); }
  Point(Curvedata &EE, const bigint& x, const bigint& y)
    : X(x), Y(y), E(&EE), ord(0), height(to_bigfloat(-1.0))
    { Z=1;} // no need to reduce already integral
  Point(Curvedata *EE, const bigint& x, const bigint& y)
    : X(x), Y(y), E(EE), ord(0), height(to_bigfloat(-1.0))
    { Z=1;} // no need to reduce already integral
  Point(const Point& Q)
    :X(Q.X), Y(Q.Y), Z(Q.Z), E(Q.E), ord(Q.ord), height(Q.height)
    { ; }
  ~Point(void) {;}
                
  // input and output
  friend inline ostream& operator<<(ostream & os, const Point& P)
    {return os << "[" << P.X << ":" << P.Y << ":" << P.Z << "]" ;}

  friend inline void output_pari(ostream&os, const Point& P)
    {
      bigint xp=P.X, yp=P.Y, zp=P.Z;
      if(is_zero(zp)) {os<<"[0]"; return;}
      if(is_one(zp)) {os<<"["<<xp<<","<<yp<<"]"; return;}
      bigint z=gcd(xp,zp);
      os<<"["<<(xp/z)<<"/"<<(zp/z)<<","<<yp<<"/"<<zp<<"]";
    }

// Point input: 3 formats allowed are 
// [x:y:z], [x,y], [x/z,y/z] with any stype of brackets
  friend istream& operator>>(istream & is, Point& P);

  // test of equality of points        
  int operator==(const Point& Q) const
  {
    if(E != Q.E) return 0 ;      // different curves!
    return eq(*this,Q);
  }
  int operator!=(const Point& Q) const { return !(*this == Q); }
  friend int eq(const Point&P, const Point&Q);

  // assignment (p.init(.) is quicker than p=Point(.) for existing p)
  void init(Curvedata &EE,
            const bigint& x, const bigint& y, const bigint& z)
    {E=&EE; X=x; Y=y; Z=z; reduce(); ord=0; height=-1.0; }
  void init(Curvedata *EE,
            const bigint& x, const bigint& y, const bigint& z)
    {E=EE; X=x; Y=y; Z=z; reduce(); ord=0; height=-1.0; }
  void init(Curvedata &EE,
            const bigint& x, const bigint& y)
    {E=&EE; X=x; Y=y; Z = 1; ord=0; height=-1.0; }
  void init(Curvedata *EE,
            const bigint& x, const bigint& y)
    {E=EE; X=x; Y=y; Z = 1; ord=0; height=-1.0; }
  void operator=(const Point& Q) // P1 = P2
    { E=Q.E; X=Q.X ; Y=Q.Y; Z=Q.Z; ord=Q.ord; height=Q.height; }

  friend Point shift(const Point& p,  Curvedata* newc, 
		     const bigint& u, const bigint& r, 
		     const bigint& s, const bigint& t, int back=0); 

  void operator+=(const Point&) ; // P1 += P2 ; order and height unknown
  void operator-=(const Point&) ; // P1 -= P2 ; ditto

  // addition of points, etc
  Point operator+(const Point &) const ; // P1 + P2
  Point operator-(const Point &) const ; // P1 - P2
  Point operator-(void) const ; // -P
  Point twice(void) const ; // doubles P
  friend Point operator*(int, const Point&) ; // n*P
                
  // access functions
  Curve getcurve() const {return *E;}
  void getcoordinates(bigint& x, bigint& y, bigint& z) const
    {x=X; y=Y; z=Z; }
  friend inline bigint getX(const Point& p) {return p.X; }
  friend inline bigint getY(const Point& p) {return p.Y; }
  friend inline bigint getZ(const Point& p) {return p.Z; }
  friend int order(Point& p);       // calculate and set if not set
  friend int order(Point& p,  vector<Point>&multiples);
                        // also create and return list of multiples
  friend bigfloat height(Point& P);    //calculate and set if not set
  friend bigfloat realheight(const Point& P);
  friend bigfloat pheight(const Point& P, const bigint& p);

  // useful logical tests
  int iszero() const { return Z==0; }
  int isintegral() const { return Z==1; }
  int isvalid() const ; // P on its curve ?

}; // end of point class

// the real x and y coords of the point
void realify_point(const Point& P, bigfloat&x, bigfloat& y);

// the real component of the canonical height
bigfloat realheight(const bigfloat& x, const Curvedata* E);


// the height pairing of two points
bigfloat height_pairing(Point& P, Point& Q);  

// regulator of a list of n points
bigfloat regulator(vector<Point>& points);  // not a const array; heights get set.

void raw_shift_point(const bigint& x, const bigint& y, const bigint& z, 
                const bigint& u, const bigint& r, 
                const bigint& s, const bigint& t, 
                bigint& newx, bigint& newy, bigint& newz, int back=0);


// torsion functions
// N.B. Don't make the params const here

vector<Point> two_torsion(Curvedata& E);
vector<bigint> three_torsion_x(Curvedata& E);
vector<Point> three_torsion(Curvedata& E);
vector<Point> torsion_points(Curvedata& E);

inline long ntorsion(Curvedata& E)
{
  return E.get_ntorsion();
}

// end of file: points.h

#endif
