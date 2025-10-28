// curvemod.h: declaration of class curvemodq for elliptic curve mod (prime) q 
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
 

// allow for multiple includes
#ifndef _ECLIB_CURVEMOD_
#define _ECLIB_CURVEMOD_

#include "polys.h"
#include "curve.h"

// Class for an elliptic curve mod q
// 
// The q is fixed -- must be set before use!
//


class pointmodq;

class curvemodq  { 

protected:
  galois_field* Fq;           // pointer to ground field
  ZZ q;                   //  the modulus
  ZZ_p a1,a2,a3,a4,a6;  //  the coefficients mod q
  ZZ order;               // number of points (0 if not set)

public:
  // constructors 
  curvemodq(void);
  curvemodq(const Curve& E, const ZZ& qq);
  ~curvemodq();
  curvemodq(const curvemodq& C); // copy constructor
  void operator=(const curvemodq& C); // assignment
  
  // access
  void get_ai(ZZ_p& aa1, ZZ_p& aa2, ZZ_p& aa3, ZZ_p& aa4, ZZ_p& aa6) const
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

  ZZ get_modulus() const { return q; }

  void set_group_order_via_legendre(); 
  void set_group_order();
  ZZ group_order() {if(is_zero(order)) set_group_order(); return order;}
  
  friend class galois_field;
  friend class ffmodq;
  friend class pointmodq;
  friend class TLSS;
  friend galois_field get_field(const curvemodq& C);
  friend ZZ get_modulus(const curvemodq& E);
  };

inline ostream& operator<<(ostream& os, const curvemodq& C)
{
  C.output(os);
  return os;
}

  // just a wrapper round the constructor, (originally for for LiDIA-compatibility):

inline curvemodq reduce_curve(const Curve& E, const ZZ& q)
{
  return curvemodq(E,q);
}

inline galois_field get_field(const curvemodq& C) { return galois_field(C.q); }
inline ZZ get_modulus(const curvemodq& C) { return C.q; }

ZZ_pX div_pol_odd(const curvemodq& C, int n); 
ZZ_pX makepdivpol(const curvemodq& C, int p);


#endif // #define _CURVEMOD_
