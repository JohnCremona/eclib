// pointsmod.h: declaration of classes pointmodq and curvemodqbasis
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

// curvemodqbasis is derived from curvemodq (see file curvemod.h) and
// contains a Z-basis for the group of points

// The baby-step-giant step algorithm in my_bg_algorithm was
// originally adapted from LiDIA's bg_algorithm(); it has some changes.

// The point-counting and group structure algorithm in
// my_isomorphism_type() provide the same functionality as LiDIA's
// isomorphism_type() but has been rewritten from scratch; a main
// difference from the LiDIA version is the use of Weil pairing when
// the group is not cyclic.  This is only intended for use when q is
// small-medium sized (NOT cryptographic!).  The current
// implementation is only for prime fields, but the same strategy
// would work over arbitrary finite fields.

// allow for multiple includes
#ifndef _ECLIB_POINTSMOD_H
#define _ECLIB_POINTSMOD_H

#include "points.h"
#include "curvemod.h"

class ffmodq;

// Class for points on an elliptic curve mod q
// 

class galois_field;

class pointmodq{ 
  gf_element X ; // inhomogeneous coordinates
  gf_element Y ; //
  int is0flag;  // set iff it's the point at infinity
  ZZ order; // 0 if not set
  curvemodq E;                  //  the curve it's on

public:
  // constructors 
  pointmodq(void) :E() {}
  explicit pointmodq(const curvemodq& EE ) :is0flag(1), order(1), E(EE) {} //  the point at oo
  pointmodq(const gf_element&x, const gf_element&y, const curvemodq& EE) 
    :X(x), Y(y), is0flag(0), order(0), E(EE)
    {
      if(!on_curve())
	cout<<"Error!  ("<<x<<","<<y<<") is not on "<<(EE)<<endl;
    }
  pointmodq(const ZZ&x, const ZZ&y, const curvemodq& EE) 
    :X(to_ZZ_p(x)), Y(to_ZZ_p(y)), is0flag(0), order(0), E(EE)
    {
      if(!on_curve())
	cout<<"Error!  ("<<x<<","<<y<<") is not on "<<(EE)<<endl;
    }
  pointmodq(const gf_element&x, const curvemodq& EE);  // a point with X=x or oo if none
  pointmodq(const pointmodq& P ) :X(P.X), Y(P.Y), is0flag(P.is0flag), order(P.order), E(P.E) {;}

  // assignment
  void operator=(const pointmodq& P) {is0flag=P.is0flag; E=P.E; X=P.X; Y=P.Y; order=P.order;}

  // access
  int is_zero() const { return is0flag;}
  gf_element get_x() const 
    {
      if (is0flag){return to_ZZ_p(0);} 
      return X;
    }
  gf_element get_y() const 
    {
      if (is0flag){return to_ZZ_p(1);} 
      return Y;
    }
  curvemodq get_curve() const {return E;}
  // output
  void output(ostream& os) const;

  // test of equality of points        
  int operator==(const pointmodq& Q) const
  {
    if(E!=(Q.E)) return 0; // different curve!
    int fl=Q.is0flag;
    if(is0flag) return fl;
    if(fl) return 0;
    return (X==Q.X)  && (Y==Q.Y);
  }
  int operator!=(const pointmodq& Q) const { return !(*this == Q); }

  // test of validity:
  int on_curve() const
    {
      if(is0flag) return 1;
      return (Y*(Y+(E.a1)*X+(E.a3))-(X*(X*(X+(E.a2))+(E.a4))+(E.a6)))==to_ZZ_p(0);
    }
  // make a point with given x & return true, or return false if none
  int set_x_coordinate(const gf_element& x);

  // order: get_order() computes if not yet set
  void set_order(const ZZ& n) {order=n;} // use with caution!
  ZZ get_order();
  ZZ get_order(const ZZ& lower, const ZZ& upper); //if bounds known
  ZZ get_order(const ZZ& mult); // use if multiple of order known

  // addition of points, etc
  pointmodq operator+(const pointmodq & Q) const ; // add Q to this
  pointmodq operator-(const pointmodq & Q) const ; // sub Q from this
  pointmodq operator-(void) const ; // -P
  pointmodq negate(void) const ; // negates P
  pointmodq twice(void) const ; // doubles P

  void operator+=(const pointmodq & P)
    {
      *this =  (*this)+P;
    }
  void operator-=(const pointmodq & P)
    {
      *this =  (*this)-P;
    }
  friend pointmodq operator*(long, const pointmodq&) ; // n*P
  friend pointmodq operator*(const ZZ&, const pointmodq&) ; // n*P
  friend ZZ order_point(pointmodq& P); // not const as may set the order
  friend galois_field base_field(const pointmodq& P);

  friend class ffmodq;
};

inline ostream& operator<<(ostream& os, const pointmodq& P)
{
  P.output(os);
  return os;
}

inline galois_field base_field(const pointmodq& P)
{
  return galois_field((P.get_curve()).get_modulus());
}

inline ZZ order_point(pointmodq& P) // not const as may set the order
{ return P.get_order();}

pointmodq reduce_point(const Point& P,  const curvemodq& Emodq);

class curvemodqbasis : public curvemodq { 
  ZZ n1,n2,n;           // n=n1*n2 = #E(Fq)
  pointmodq P1,P2;          // basis for E(F_q)
  void set_basis();         // computes basis  
  int lazy_flag;            // if 1, only computes a "lazy basis"
                            // with P2=0 and P1 of "large" order
 public:

  curvemodqbasis(void) :curvemodq(), n(0), n1(0), n2(0) {}
  explicit curvemodqbasis(const curvemodq& C, int lazy=0)
    :curvemodq(C), lazy_flag(lazy)
  {
    set_basis();
  }
  curvemodqbasis(const Curve& E, const ZZ& q, int lazy=0) 
    :curvemodq(reduce_curve(E,q)) 
  {
    lazy_flag=lazy;
    set_basis();
  }

  //  ~curvemodqbasis(void) {;}

  ZZ get_order() {return n;}
  ZZ get_exponent() {return n1;}
  pointmodq get_gen(int i);

  vector<pointmodq> get_pbasis(int p);
  vector<pointmodq> get_pbasis_from_roots(int p,  const vector<gf_element>& xi);
  vector<pointmodq> get_pbasis_via_divpol(int p);
  vector<pointmodq> get_pbasis_via_divpol(int p, const ZZX& pdivpol);

  friend class TLSS;
};

ZZ my_bg_algorithm(const pointmodq& PP,
		    const pointmodq& QQ,
                    const ZZ& lower,
		    const ZZ& upper,
		    bool info=false);

void set_hasse_bounds(const ZZ& q, ZZ& l, ZZ& u);
ZZ my_order_point(const pointmodq& PP);
ZZ my_order_point(const pointmodq& PP, 
		   const ZZ& lower, const ZZ& upper);
ZZ my_order_point(const pointmodq& PP, const ZZ& mult);

// returns minimal m>0 s.t. m*Q is in <P> with m*Q=a*P.  Special case:
// if <Q> and <P> are disjoint, then m=order(Q) and a=0.
ZZ linear_relation( pointmodq& P, pointmodq& Q, ZZ& a);

// Replace P (of order ordP) with a point whose order is lcm(ordP,order(Q))
void merge_points_1(pointmodq& PP, ZZ& ordP, pointmodq& Q);

// Given independent generators P1,P2 with orders n1, n2 and n2|n1,
// and a new point Q: 
//
// (1) If ord(Q)|ord(P1) -- the normal case -- replace P2 with a point
// whose order mod <P1> is lcm of ord(P2) and ord(Q) mod <P1> 
//
// (2) Else replace P1 as with merge_points_1 and reset P2

void merge_points_2(pointmodq& P1, ZZ& n1, pointmodq& P2, ZZ& n2, 
		    const ZZ& n2target, pointmodq& Q);

inline bool less(const gf_element& a, const gf_element& b)
{
  return LiftGF(a)<LiftGF(b);
}

// find a point of "large" order
void one_generator(curvemodq& Cq, ZZ& n1, pointmodq& P1);

// find full Z-basis
void my_isomorphism_type(curvemodq& C, 
			 ZZ& n1, ZZ& n2, pointmodq& P1, pointmodq& P2);
void my_isomorphism_type_new(curvemodq& Cq, 
    		     ZZ& n1, ZZ& n2, pointmodq& P1, pointmodq& P2);

void set_order_point(pointmodq& P, const ZZ& n);


#endif // #define _POINTSMOD_
