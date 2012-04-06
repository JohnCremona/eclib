// ffmod.h: declaration of class ffmodq and Weil pairing functions
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
 

// ffmodq is the function field of an elliptic curve mod a prime q
// (or more precisely the affine coordinate ring Fq[x,y])

// allow for multiple includes
#ifndef _FFMOD_
#define _FFMOD_

#include "curve.h"
#include "curvemod.h"
#include "pointsmod.h"

class ffmodq{ 

  public:
  static galois_field Fq;        //  the constant field
  static curvemodq E;            //  the curve mod q
  static FqPoly f1, f2;          //  f2=a1*x+a3, f2=x^3+a2*x^2+a4*x+a6
  FqPoly h1, h2;                 //  for h1+y*h2

public:

  //  constructors
  
  //  special one to initialize the curve and field only:
  ffmodq(const curvemodq& EE);

  //  normal ones:
  ffmodq(void) 
    {
      init_h1h2();
      FqPolyAssign0(h1);
      FqPolyAssign0(h2);
    }

  ffmodq(const gf_element& c) 
    {
      init_h1h2();
      FqPolyAssignGF(h1,c);
      FqPolyAssign0(h2);
    }

  ffmodq(const bigint& c) 
    {
      init_h1h2();
      FqPolyAssignZ(h1,c);
      FqPolyAssign0(h2);
    }

  ffmodq(const FqPoly& hh1) 
    {
      init_h1h2();
      h1=hh1; 
      FqPolyAssign0(h2);
    }

  ffmodq(const FqPoly& hh1, const FqPoly& hh2) {h1=hh1; h2=hh2;}

  //  initialization
  void init_f1f2(void);

  void init_h1h2(void) 
    {    
      FqPolySetField(h1,Fq);
      FqPolySetField(h2,Fq);
    }

  // assignment
  void operator=(const ffmodq& a) {h1=a.h1; h2=a.h2;}

  // equality test
  int operator==(const ffmodq& b) const;
  int operator!=(const ffmodq& b) const {return !((*this)==b);}

  // output
  void output(ostream& os) const;

  // addition, subtraction, multiplication
  ffmodq operator+(const ffmodq& b) const;
  ffmodq operator-(const ffmodq& b) const;
  ffmodq operator*(const ffmodq& b) const;
  ffmodq operator*(const FqPoly& h) const;

  //  division
  ffmodq operator/(const FqPoly& h) const;
  ffmodq operator/(const ffmodq& b) const;

  //  evaluation at a point:
  gf_element evaluate(const pointmodq& P) const;
  gf_element operator()(const pointmodq& P) const {return this->evaluate(P);}

  //  vertical line through a point:
  friend ffmodq vertical(const pointmodq& P);

  //  tangent at a point:
  friend ffmodq tangent(const pointmodq& P);

  //  chord between points:
  friend ffmodq chord(const pointmodq& P, const pointmodq& Q);


};

// weil_pol(T,m): T is a point of finite order m; returns a function
// f_T whose divisor is m(T)-m(0).

// The second version evaluates that at another point S without
// actually computing the polynomial

ffmodq weil_pol(const pointmodq& T, int m);

gf_element evaluate_weil_pol(const pointmodq& T, int m, const pointmodq& S);

gf_element weil_pairing(const pointmodq& S, const pointmodq& T, int m);

inline ostream& operator<<(ostream& os, const ffmodq& f)
{
  f.output(os);
  return os;
}


#endif
