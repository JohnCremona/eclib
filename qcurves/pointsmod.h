// pointsmod.h: declaration of classes pointmodq and curvemodqbasis
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
 // and functions for point counting and elliptic curve discrete log


// Under LiDIA pointmodq is a wrapper for point<gf_element>; under NTL
// it is self-contained.  We provided a common interface.

// curvemodqbasis is derived from curvemodq (see file curvemod.h) and
// contains a Z-basis for the group of points

// The baby-step-giant step algorithm in my_bg_algorithm is adapted
// from LiDIA's bg_algorithm() with few changes.

// The point-counting and group structure algorithm in
// my_isomorphism_type() provide the same functionality as LiDIA's
// isomorphism_type() but has been rewritten from scratch by JEC; a
// main difference from the LiDIA version is the use of Weil pairing
// when the group is not cyclic.  This is only intended for use when q
// is small-medium sized (NOT cryptographic!) -- as is also true of
// LiDIA's isomorphism_type().  The current implementation is only for
// prime fields, but the same strategy would work over arbitrary
// finite fields.

// allow for multiple includes
#ifndef _POINTSMOD_
#define _POINTSMOD_

class ffmodq;

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)

#define pointmodq point<gf_element>

inline galois_field base_field(const pointmodq& P)
{
  return P.get_x().get_field();
}

#else // NTL

// Class for points on an elliptic curve mod q
// 

class galois_field;

class pointmodq{ 
  gf_element X ; // inhomogeneous coordinates
  gf_element Y ; //
  int is0flag;  // set iff it's the point at infinity
  bigint order; // 0 if not set
  curvemodq E;                  //  the curve it's on

public:
  // constructors 
  pointmodq(void) :E() {;}
  pointmodq(const curvemodq& EE ) :is0flag(1), order(BIGINT(1)), E(EE) {;} //  the point at oo
  pointmodq(const gf_element&x, const gf_element&y, const curvemodq& EE) 
    :X(x), Y(y), is0flag(0), order(BIGINT(0)), E(EE)
    {
      if(!on_curve())
	cout<<"Error!  ("<<x<<","<<y<<") is not on "<<(EE)<<endl;
    }
  pointmodq(const bigint&x, const bigint&y, const curvemodq& EE) 
    :is0flag(0), order(BIGINT(0)), E(EE)
    {
      X=to_ZZ_p(x);
      Y=to_ZZ_p(y);
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
      if (is0flag){cerr<<"error in get_x();"<<endl; return to_ZZ_p(0);} 
      return X;
    }
  gf_element get_y() const 
    {
      if (is0flag){cerr<<"error in get_x();"<<endl; return to_ZZ_p(0);} 
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
  void set_order(const bigint& n) {order=n;} // use with caution!
  bigint get_order();
  bigint get_order(const bigint& lower, const bigint& upper); //if bounds known
  bigint get_order(const bigint& mult); // use if multiple of order known

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
  friend pointmodq operator*(const bigint&, const pointmodq&) ; // n*P
  friend bigint order_point(pointmodq& P); // not const as may set the order
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

inline bigint order_point(pointmodq& P) // not const as may set the order
{ return P.get_order();}

#endif // end of LiDIA/NTL split

pointmodq reduce_point(const Point& P,  const curvemodq& Emodq);

class curvemodqbasis : public curvemodq { 
  bigint n1,n2,n;           // n=n1*n2 = #E(Fq)
  pointmodq P1,P2;          // basis for E(F_q)
  void set_basis();         // computes basis  
  int lazy_flag;            // if 1, only computes a "lazy basis"
                            // with P2=0 and P1 of "large" order
 public:

  curvemodqbasis(void) :curvemodq(){n=n1=n2=0;}
  curvemodqbasis(const curvemodq& C, int lazy=0) 
    :curvemodq(C) 
  {
    lazy_flag=lazy;
    set_basis();
  }
  curvemodqbasis(const Curve& E, const bigint& q, int lazy=0) 
    :curvemodq(reduce_curve(E,q)) 
  {
    lazy_flag=lazy;
    set_basis();
  }

  //  ~curvemodqbasis(void) {;}

  bigint get_order() {return n;}
  pointmodq get_gen(int i);

  vector<pointmodq> get_pbasis(int p);
  vector<pointmodq> get_pbasis_from_roots(int p,  const vector<gf_element>& xi);
  vector<pointmodq> get_pbasis_via_divpol(int p);
  vector<pointmodq> get_pbasis_via_divpol(int p, const vector<bigint>& pdivpol);

  friend class TLSS;
};

bigint my_bg_algorithm(const pointmodq& PP,
		    const pointmodq& QQ,
                    const bigint& lower,
		    const bigint& upper,
		    bool info=false);

void set_hasse_bounds(const bigint& q, bigint& l, bigint& u);
bigint my_order_point(const pointmodq& PP);
bigint my_order_point(const pointmodq& PP, 
		   const bigint& lower, const bigint& upper);
bigint my_order_point(const pointmodq& PP, const bigint& mult);

// returns minimal m>0 s.t. m*Q is in <P> with m*Q=a*P.  Special case:
// if <Q> and <P> are disjoint, then m=order(Q) and a=0.
bigint linear_relation( pointmodq& P, pointmodq& Q, bigint& a);

// Replace P (of order ordP) with a point whose order is lcm(ordP,order(Q))
void merge_points_1(pointmodq& PP, bigint& ordP, pointmodq& Q);

// Given independent generators P1,P2 with orders n1, n2 and n2|n1,
// and a new point Q: 
//
// (1) If ord(Q)|ord(P1) -- the normal case -- replace P2 with a point
// whose order mod <P1> is lcm of ord(P2) and ord(Q) mod <P1> 
//
// (2) Else replace P1 as with merge_points_1 and reset P2

void merge_points_2(pointmodq& P1, bigint& n1, pointmodq& P2, bigint& n2, 
		    const bigint& n2target, pointmodq& Q);

inline bool less(const gf_element& a, const gf_element& b)
{
  return LiftGF(a)<LiftGF(b);
}

// find a point of "large" order
void one_generator(curvemodq& Cq, bigint& n1, pointmodq& P1);

// find full Z-basis
void my_isomorphism_type(curvemodq& C, 
			 bigint& n1, bigint& n2, pointmodq& P1, pointmodq& P2);
void my_isomorphism_type_new(curvemodq& Cq, 
    		     bigint& n1, bigint& n2, pointmodq& P1, pointmodq& P2);

void set_order_point(pointmodq& P, const bigint& n);


#endif
