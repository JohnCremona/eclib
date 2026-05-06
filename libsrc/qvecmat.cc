// File QVECMAT.CC: classes for working with rational vectors and matrices
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

// NB These classes do not provide all the functionality one might
// think of, just what is needed for the classes Field, FieldElement,
// FieldIso, Newspace and Newform.  Also, the base integer type is
// currently fixed to be ZZ (NTL integers) but could easily be
// templated to work also with int, long and INT (wrapping FLINT
// integers).

#include "eclib/qvecmat.h"
#include "eclib/polred.h"

void Qvec::cancel() // divides through by gcd(content(numerator, denom))
{
  if (denom<0)
    {
      denom = -denom;
      numerator = -numerator;
    }
  ZZ g = gcd(content(numerator), denom);
  if (IsOne(g))
    return;
  denom /= g;
  numerator /= g;
}

// String for pretty printing, used in default <<, or (if raw) raw
// output, suitable for re-input:
string Qvec::str(int raw) const
{
  ostringstream s;
  if (raw)
    s << numerator << " " << denom;
  else
    {
      s << numerator;
      if (denom>1)
        s << "/" << denom;
    }
  return s.str();
}

Qvec Qvec::operator+(const Qvec& b) const
{
  Qvec a = *this;
  if (!b.is_zero())
    a += b;
  return a;
}

void Qvec::operator+=(const Qvec& b)
{
  if (b.is_zero())
    return;
  numerator = b.denom*numerator + denom*b.numerator;
  denom *= b.denom;
  cancel();
}

Qvec Qvec::operator-(const Qvec& b) const
{
  Qvec a = *this;
  if (!b.is_zero())
    a -= b;
  return a;
}

void Qvec::operator-=(const Qvec& b)
{
  if (b.is_zero())
    return;
  numerator = b.denom*numerator - denom*b.numerator;
  denom *= b.denom;
  cancel();
}

void Qmat::cancel() // divides through by gcd(content(numerator, denom))
{
  if (denom<0)
    {
      denom = -denom;
      numerator = -numerator;
    }
  ZZ g = gcd(numerator.content(), denom);
  if (IsOne(g))
    return;
  denom /= g;
  numerator /= g;
}

// String for pretty printing, used in default <<, or (if raw) raw
// output, suitable for re-input:
string Qmat::str(int raw) const
{
  ostringstream s;
  if (raw)
    s << numerator << " " << denom;
  else
    {
      s << numerator;
      if (denom>1)
        s << "/" << denom;
    }
  return s.str();
}

Qmat Qmat::inverse() const
{
  mat_m inverse_numerator;
  ZZ d = ::inverse(numerator, inverse_numerator);
  return Qmat(denom * inverse_numerator, d);
}

void Qmat::setcol(int i, const Qvec& v)
{
  ZZ new_d = lcm(denom, v.denom);
  numerator *= (new_d/denom);
  denom = new_d;
  numerator.setcol(i, (new_d/v.denom) * v.numerator);
  cancel();
}

void Qmat::setrow(int i, const Qvec& v)
{
  ZZ new_d = lcm(denom, v.denom);
  numerator *= (new_d/denom);
  denom = new_d;
  numerator.setrow(i, (new_d/v.denom) * v.numerator);
  cancel();
}

// Compute matrix B, with column j equal to A^(j-1)v for j from 1 to
// d, where v is the first unit vector, so that A*B=B*C where C is the
// companion matrix of A's minpoly. We assume that A has integral char
// poly and return onyl the 'numerator' of C whose denominator is 1.
mat_m Qmat::companion_transform(Qmat& B, Qmat& Binv) const
{
  int d = nrows();
  Qvec v(d);
  v.set_unit_vector(1);
  B = Qmat(d,d);
  B.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = (*this)*v;
      B.setcol(i,v);
    }
  Binv = B.inverse();
  mat_m C = (Binv*(*this)*B).get_numerator();
  assert (CompanionMatrix(this->charpoly()) == mat_to_mat_ZZ(C));
  return C;
}

Qmat Qmat::operator+(const Qmat& b) const
{
  Qmat a = *this;
  if (!b.is_zero())
    a += b;
  return a;
}

void Qmat::operator+=(const Qmat& b)
{
  if (b.is_zero())
    return;
  numerator = b.denom*numerator + denom*b.numerator;
  denom *= b.denom;
  cancel();
}

Qmat Qmat::operator-(const Qmat& b) const
{
  Qmat a = *this;
  if (!b.is_zero())
    a -= b;
  return a;
}

void Qmat::operator-=(const Qmat& b)
{
  if (b.is_zero())
    return;
  numerator = b.denom*numerator - denom*b.numerator;
  denom *= b.denom;
  cancel();
}
