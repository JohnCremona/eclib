// quadratic.cc: implementation of class for handling integer quadratics
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
 
#include "marith.h"
#include "unimod.h"
#include "quadratic.h"

quadratic::~quadratic() {delete [] coeffs;}

void quadratic::init()
{
  coeffs = new bigint[3];
}

bigint quadratic::coeff(int i) 
    {if((i>=0)&&(i<=2)) return coeffs[i]; else {bigint ans; return ans;}}

bigint quadratic::operator[](int i) const
    {if((i>=0)&&(i<=2)) return coeffs[i]; else {bigint ans; return ans;}}

void quadratic::transform(const unimod& m)
{
  bigint newq0 = eval(m.m11, m.m21);
  bigint newq2 = eval(m.m12, m.m22);
  coeffs[1] = 2 * (coeffs[0]*m.m11*m.m12 + coeffs[2]*m.m21*m.m22) 
               + coeffs[1] * (m.m11*m.m22+m.m12*m.m21);
  coeffs[0] = newq0;
  coeffs[2] = newq2;
}

void quadratic::x_shift(const bigint& a, unimod& m)
{
  const bigint& aq0=a*coeffs[0];
  coeffs[2]+= (aq0+coeffs[1])*a;
  coeffs[1]+= 2*aq0;
  m.x_shift(a);
}

void quadratic::y_shift(const bigint& a, unimod& m)
{
  const bigint& aq2=a*coeffs[2];
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
  bigint a = roundover(-coeffs[1],2*coeffs[0]);
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

bigint resultant(const quadratic& q1, const quadratic& q2)
{return sqr(q2.coeffs[0]*q1.coeffs[2]) + sqr(q1.coeffs[0]*q2.coeffs[2]) +
   q2.coeffs[2]*q2.coeffs[0]*sqr(q1.coeffs[1]) - q2.coeffs[1]*q1.coeffs[1]*(q2.coeffs[0]*q1.coeffs[2] + q1.coeffs[0]*q2.coeffs[2]) +
   (sqr(q2.coeffs[1])-2*q2.coeffs[0]*q2.coeffs[2])*q1.coeffs[0]*q1.coeffs[2];
}
