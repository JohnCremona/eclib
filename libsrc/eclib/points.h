// points.h:  declarations of Point class for points on elliptic curves
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
 
// originally adapted from Elliptic.h by Oisin McGuiness

// allow for multiple includes
#ifndef _ECLIB_POINTS_
#define _ECLIB_POINTS_

#include "polys.h"
#include "p2points.h"
#include "divpol.h"

//
// class for  points on elliptic curves
//

class Point;
Point transform(const Point& p,  Curvedata* newc, 
			 const bigint& u, 
			 const bigint& r, const bigint& s, const bigint& t, 
			 int back=0); 

class Point : public P2Point { 
  Curvedata *E;    // pointer to the curve that the point is on
  int ord;         // order: 0 if not calculated yet, -1 if infinite
  bigfloat height; // -1.0 if not calculated yet, 0.0 for torsion point
public:
  // constructors
  Point(void)
    : P2Point(), E(0), ord(0), height(to_bigfloat(-1.0))
    { ; }
  explicit Point(Curvedata &EE)      // set to point at infinity
    : P2Point(0,1,0), E(&EE), ord(1), height(to_bigfloat(0.0))
    { ; }
  explicit Point(Curvedata *EE)      // set to point at infinity
    : P2Point(0,1,0), E(EE), ord(1), height(to_bigfloat(0.0))
    { ; }
  Point(Curvedata &EE, const bigint& x, const bigint& y, const bigint& z)
    : P2Point(x,y,z), E(&EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(Curvedata &EE, const P2Point& P)
    : P2Point(P), E(&EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(Curvedata *EE, const bigint& x, const bigint& y, const bigint& z)
    : P2Point(x,y,z), E(EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(Curvedata *EE, const P2Point& p)
    : P2Point(p), E(EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(Curvedata &EE, const bigint& x, const bigint& y)
    : P2Point(x,y), E(&EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(Curvedata *EE, const bigint& x, const bigint& y)
    : P2Point(x,y), E(EE), ord(0), height(to_bigfloat(-1.0))
    { ; }
  Point(const Point& Q)
    : P2Point(Q), E(Q.E), ord(Q.ord), height(Q.height)
    { ; }
  ~Point(void) {;}

  // input and output are inherited from P2Point class but the input
  // function must initialize the ord and height fields too
  friend istream& operator>>(istream & is, Point& P)
  {
    is>>(P2Point&)P;
    P.ord=0;
    P.height=to_bigfloat(-1.0);
    // NB P's Curve should have been set when it was constructed.
    return is;
  }

  // test of equality of points        
  int operator==(const Point& Q) const
  {
    if(E != Q.E) return 0 ;      // different curves!
    return eq(*this,Q);
  }
  int operator!=(const Point& Q) const { return !(*this == Q); }

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

  friend Point transform(const Point& p,  Curvedata* newc, 
			 const bigint& u, 
			 const bigint& r, const bigint& s, const bigint& t, 
			 int back); 

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
  friend int order(Point& p);       // calculate and set if not set
  friend int order(Point& p,  vector<Point>&multiples);
                        // also create and return list of multiples
  friend bigfloat height(Point& P);    //calculate and set if not set
  friend bigfloat realheight(const Point& P);
  friend bigfloat pheight(const Point& P, const bigint& p);

  vector<Point> division_points(int m); // list of Q s.t. n*Q=this

  // useful logical tests
  int is_zero() const { return isinfinite(); }
  int isvalid() const ; // P on its curve ?
  int is_torsion() { return order(*this)>0; } // will compute and set the order if needed
  int is_on_real_identity_component() const;
  int is_on_egg() const {return !is_on_real_identity_component();}
// return 1 if P mod p is nonsingular (or for p=0 if it is on the real identity component):
  int has_good_reduction(long p) const;
  int has_good_reduction(const bigint& p) const;
  // return 1 if P has good reduction at all p in list, if not then p0 holds first bad prime
  int has_good_reduction(const vector<bigint>& plist, bigint& p0, int check_real=0) const;

}; // end of point class

// list of 0,1 or 2 points with given x-coordinate:
vector<Point> points_from_x(Curvedata &E, const bigrational& x);

// the real x and y coords of the point
void realify_point(const Point& P, bigfloat&x, bigfloat& y);

// the real component of the canonical height
bigfloat realheight(const bigfloat& x, const Curvedata* E);


// the height pairing of two points
bigfloat height_pairing(Point& P, Point& Q);

// regulator of a list of n points
bigfloat regulator(vector<Point>& points);  // not a const array; heights get set.

// torsion functions returning
// (exact==0) list of points in E[m] for m=2, m=3 and general m
// (exact==1) list of points of exact order m for m=2, m=3 and general m

// N.B. Don't make the params const here

vector<Point> two_torsion(Curvedata& E, int exact=0);
vector<bigint> three_torsion_x(Curvedata& E);
vector<Point> three_torsion(Curvedata& E, int exact=0);

// List m-torsion points, i.e. points in E[m], or points of exact order m if exact==1
vector<Point> m_torsion(Curvedata& E, long m, int exact=0);
vector<Point> torsion_points(Curvedata& E);

inline long ntorsion(Curvedata& E)
{
  return E.get_ntorsion();
}

// end of file: points.h

#endif
