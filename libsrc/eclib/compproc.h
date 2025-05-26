// compproc.h: declarations of functions using complex numbers
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
 
#ifndef _ECLIB_COMPPROC_H_
#define _ECLIB_COMPPROC_H_

#include <eclib/interface.h>

int is_small(bigfloat x);
int is_small(const bigcomplex& z);
int is_real(const bigcomplex& z);

void orderreal(bigfloat& e1, bigfloat& e2, bigfloat& e3);  // puts in decreasing order
bigcomplex root(const bigcomplex& z, int n);
bigcomplex cagm(const bigcomplex& a, const bigcomplex& b);
bigcomplex normalize(bigcomplex& w1, bigcomplex& w2);
void getc4c6(const bigcomplex& w1, const bigcomplex& w2, bigcomplex& c4, bigcomplex &c6);
bigcomplex discriminant(const bigcomplex& b, const bigcomplex& c, const bigcomplex& d);

vector<bigcomplex> solvecubic(const bigcomplex& c1, const bigcomplex& c2, const bigcomplex& c3);
inline vector<bigcomplex> solvecubic(const bigfloat& c1, const bigfloat& c2, const bigfloat& c3)
{return solvecubic(bigcomplex(c1), bigcomplex(c2), bigcomplex(c3));}

vector<bigcomplex> solvequartic(const bigcomplex& a, const bigcomplex& b, const bigcomplex& c, const bigcomplex& d);
inline vector<bigcomplex> solvequartic(const bigfloat& a, const bigfloat& b, const bigfloat& c, const bigfloat& d)
{return solvequartic(bigcomplex(a), bigcomplex(b), bigcomplex(c), bigcomplex(d));}

vector<bigcomplex> solverealquartic(const bigfloat& a, const bigfloat& b, const bigfloat& c, const bigfloat& d, const bigfloat& e);

void quadsolve(const bigfloat& p, const bigfloat& q, bigcomplex& root1,bigcomplex& root2);

#endif
