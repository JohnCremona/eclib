// curvemod.h: declaration of class curvemodq for elliptic curve mod (prime) q 
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
 

// NB Under LiDIA this is just an interface to LiDIA's
// elliptic_curve<gf_element>;  under NTL it is self-contained


// allow for multiple includes
#ifndef _CURVEMOD_
#define _CURVEMOD_


// Class for an elliptic curve mod q
// 
// The q is fixed -- must be set before use!
//

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)

#include "LiDIA/elliptic_curve.h"
#include "LiDIA/point.h"
class elliptic_curve<gf_element>;
class point<gf_element>;
#define pointmodq point<gf_element>
#define curvemodq elliptic_curve<gf_element>

// Two class methods not in LiDIA but provided below for NTL version:
#define set_group_order group_order
#define set_group_order_via_legendre group_order

curvemodq reduce_curve(const Curve& E, const bigint& q);

inline bigint get_modulus(const curvemodq& E)
{
  return E.discriminant().get_field().characteristic();
}

inline galois_field get_field(const curvemodq& E)
{
  return E.discriminant().get_field();
}

#else // NTL

class pointmodq;

class curvemodq  { 

protected:
  galois_field* Fq;           // pointer to ground field
  bigint q;                   //  the modulus
  gf_element a1,a2,a3,a4,a6;  //  the coefficients mod q
  bigint order;               // number of points (0 if not set)

public:
  // constructors 
  curvemodq(void);
  curvemodq(const Curve& E, const bigint& qq);
  ~curvemodq();
  curvemodq(const curvemodq& C); // copy constructor
  void operator=(const curvemodq& C); // assignment
  
  // access
  void get_ai(gf_element& aa1, gf_element& aa2, gf_element& aa3, gf_element& aa4, gf_element& aa6) const
    {
      aa1=a1; aa2=a2; aa3=a3; aa4=a4; aa6=a6;
    }

  // output
  void output(ostream& os) const {os<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"] mod "<<q;}
  
  // equality test
  int operator==(const curvemodq& C) const
    { return (q==C.q)&&(a1==C.a1)&&(a2==C.a2)&&(a3==C.a3)&&(a4==C.a4)&&(a6==C.a6);}
  int operator!=(const curvemodq& C) const
    { return (q!=C.q)||(a1!=C.a1)||(a2!=C.a2)||(a3!=C.a3)||(a4!=C.a4)||(a6!=C.a6);}

  // random point
  pointmodq random_point(); // defined in pointsmod.cc

  bigint get_modulus() const { return q; }

  void set_group_order_via_legendre(); 
  void set_group_order();
  bigint group_order() {if(is_zero(order)) set_group_order(); return order;}
  
  friend class galois_field;
  friend class ffmodq;
  friend class pointmodq;
  friend class TLSS;
  friend galois_field get_field(const curvemodq& C);
  friend bigint get_modulus(const curvemodq& E);
  };

inline ostream& operator<<(ostream& os, const curvemodq& C)
{
  C.output(os);
  return os;
}

  // just a wrapper round the constructor, for LiDIA-compatibility:

inline curvemodq reduce_curve(const Curve& E, const bigint& q)
{
  return curvemodq(E,q);
}

inline galois_field get_field(const curvemodq& C) { return galois_field(C.q); }
inline bigint get_modulus(const curvemodq& C) { return C.q; }

#endif // LiDIA/NTL split

FqPoly div_pol_odd(const curvemodq& C, int n); 
FqPoly makepdivpol(const curvemodq& C, int p);


#endif
