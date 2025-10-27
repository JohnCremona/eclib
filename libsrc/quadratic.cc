// quadratic.cc: implementation of class for handling integer quadratics
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
 
#include "eclib/quadratic.h"

void quadratic::transform(const unimod& m)
{
  ZZ newq0 = eval(m.m11, m.m21);
  ZZ newq2 = eval(m.m12, m.m22);
  coeffs[1] = 2 * (coeffs[0]*m.m11*m.m12 + coeffs[2]*m.m21*m.m22) 
               + coeffs[1] * (m.m11*m.m22+m.m12*m.m21);
  coeffs[0] = newq0;
  coeffs[2] = newq2;
}

void quadratic::x_shift(const ZZ& a, unimod& m)
{
  const ZZ& aq0=a*coeffs[0];
  coeffs[2]+= (aq0+coeffs[1])*a;
  coeffs[1]+= 2*aq0;
  m.x_shift(a);
}

void quadratic::y_shift(const ZZ& a, unimod& m)
{
  const ZZ& aq2=a*coeffs[2];
  coeffs[0]+= (aq2+coeffs[1])*a;
  coeffs[1]+= 2*aq2;
  m.y_shift(a);
}

void quadratic::invert(unimod& m)
{
  swap(coeffs[0],coeffs[2]);  coeffs[1]=-coeffs[1];
  m.invert();
}

void quadratic::reduce(unimod& m)
{
  //  cout<<"Reducing quadratic  "<<(*this)<<endl;
  if(coeffs[0]<0) 
    {
      coeffs[0]=-coeffs[0]; 
      coeffs[2]=-coeffs[2]; 
      coeffs[1]=-coeffs[1];
    }
  ZZ a = roundover(-coeffs[1],2*coeffs[0]);
  x_shift(a,m);
  int reduced = (coeffs[0]<=coeffs[2]);
  while(!reduced)
    {
      invert(m);
      a = roundover(-coeffs[1],2*coeffs[0]);
      x_shift(a,m);
      reduced = (coeffs[0]<=coeffs[2]);
    }
  //  cout<<"Reduced quadratic = "<<(*this)<<endl;
  //  cout<<"  via matrix "<<m<<endl;
  //  cout<<"  of determinant "<<m.det()<<endl;
}

ZZ resultant(const quadratic& q1, const quadratic& q2)
{
  const ZZ&
    a1=q1.coeffs[0], b1=q1.coeffs[1], c1=q1.coeffs[2],
    a2=q2.coeffs[0], b2=q2.coeffs[1], c2=q2.coeffs[2];
  return
    sqr(c2*a1) -a1*c2*b2*b1 + (-2*c2*a2 + sqr(b2))*c1*a1 + c2*a2*sqr(b1) - b2*a2*c1*b1 + sqr(a2*c1);
}
