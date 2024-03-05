// quadratic.h: declaration of class for handling integer quadratics
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
 

// Stored as bigint* arrays q of size 3 representing q[0]*x^2+q[1]*x+q[2]

#ifndef _ECLIB_QUADRATIC_H
#define _ECLIB_QUADRATIC_H      1 //flags that this file has been included

#include "unimod.h"

class quadratic {
  friend class unimod;
private:
  vector<bigint> coeffs;
public:
  void set(const  bigint& a, const bigint& b, const bigint& c) {coeffs = {a, b, c};}
  quadratic() { bigint zero(0); coeffs={zero,zero,zero};}
  quadratic(const  bigint& a, const bigint& b, const bigint& c) {coeffs = {a, b, c};}
  quadratic(long a, long b, long c)  {coeffs = {BIGINT(a), BIGINT(b), BIGINT(c)};}
  quadratic(const  bigint* q) {coeffs = {q[0], q[1], q[2]};}
  quadratic(const  quadratic& q) {coeffs = q.coeffs;}
  void operator=(const quadratic& q)  {coeffs = q.coeffs;}
  bigint eval(const bigint& x, const bigint& z) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x*z + coeffs[2]*sqr(z);}
  bigint operator()(const bigint& x, const bigint& z) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x*z + coeffs[2]*sqr(z);}
  bigint eval(const bigint& x) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x + coeffs[2];}
  bigint operator()(const bigint& x) const
    {return coeffs[0]*sqr(x) + coeffs[1]*x + coeffs[2];}
  bigint coeff(int i) const {return coeffs[i];}
  bigint operator[](int i) const  {return coeffs[i];}
  void set_coeff(int i, const bigint& a) {if((i>=0)&&(i<=2)) coeffs[i]=a;}
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
