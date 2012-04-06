// quadratic.h: declaration of class for handling integer quadratics
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
 

// Stored as bigint* arrays q of size 3 representing q[0]*x^2+q[1]*x+q[2]

#ifndef _QUADRATIC_H
#define _QUADRATIC_H      1 //flags that this file has been included

#include "unimod.h"

class quadratic {
  friend class unimod;
private: 
  bigint * coeffs;  // will always have length 3
  // init just allocates memory
  void init();
public:
  void set(long a, long b, long c) 
    {coeffs[0]=a; coeffs[1]=b; coeffs[2]=c;}
  void set(const  bigint& a, const bigint& b, const bigint& c) 
    {coeffs[0]=a; coeffs[1]=b; coeffs[2]=c;}
  void set(const  bigint* q)
    {coeffs[0]=q[0]; coeffs[1]=q[1]; coeffs[2]=q[2];}
  void set(const  quadratic& q)
    {coeffs[0]=q.coeffs[0]; coeffs[1]=q.coeffs[1]; coeffs[2]=q.coeffs[2];}
  quadratic() {init(); set(0,0,0);}
  ~quadratic();
  quadratic(const  bigint& a, const bigint& b, const bigint& c) {init(); set(a,b,c);}
  quadratic(long a, long b, long c) {init(); set(a,b,c);}
  quadratic(const  bigint* q) {init(); set(q);}
  quadratic(const  quadratic& q) {init(); set(q);}
  void operator=(const quadratic& q) {set(q);}
  bigint eval(const bigint& x, const bigint& z) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x*z + coeffs[2]*sqr(z);}
  bigint operator()(const bigint& x, const bigint& z) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x*z + coeffs[2]*sqr(z);}
  bigint eval(const bigint& x) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x + coeffs[2];}
  bigint operator()(const bigint& x) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x + coeffs[2];}
  bigint coeff(int i); 
  bigint operator[](int i) const;
  void set_coeff(int i, const bigint& a)
    {if((i>=0)&&(i<=2)) coeffs[i]=a;}
  friend bigint resultant(const quadratic& q1, const quadratic& q2);
  bigint disc() const  {return sqr(coeffs[1])-4*coeffs[0]*coeffs[2];}
  void output(ostream& os=cout) const
    {
      os<<"["<<coeffs[0]<<","<<coeffs[1]<<","<<coeffs[2]<<"]";
    }
  friend ostream& operator<<(ostream& os, const quadratic& q);
  void transform(const unimod& m);
  // In the next 4 functions, m already holds a unimod and is updated:
  void x_shift(const bigint& a, unimod& m);
  void y_shift(const bigint& a, unimod& m);
  void invert(unimod& m);
  void reduce(unimod& m);
};

bigint resultant(const quadratic& q1, const quadratic& q2);

inline ostream& operator<<(ostream& os, const quadratic& q)
{
  return os<<"["<<q.coeffs[0]<<","<<q.coeffs[1]<<","<<q.coeffs[2]<<"]";
}

#endif
