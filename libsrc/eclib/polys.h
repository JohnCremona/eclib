// polys.h : defines interface to NTL polynomials
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
#ifndef _ECLIB_POLYS_H
#define _ECLIB_POLYS_H

// NB --warning!  In NTL there is a universal modulus set by a single
// call to ZZ_p::init(q).

#include "gf.h"
#include "bigrat.h"

// #define SetDegree(f,d)
#define Degree(f) deg((f))

#define FqPoly ZZ_pX
#define NewFqPoly(field,name) FqPoly name
#define FqPolySetField(f,field) // do nothing
#define FqPolyAssignGF(f,c)  f=((c))
#define FqPolyAssignZ(f,c)  f=(to_ZZ_p(c))
#define FqPolyAssign0(f)  f=to_ZZ_p(0)
#define FqPolyAssign1(f)  f=to_ZZ_p(1)
#define FqPolyEval(f,c)  eval((f),to_ZZ_p(c))
#define FqPolyAssignX(f) SetX(f)
#define GetField(f) (galois_field(ZZ_p::modulus()))

vector<ZZ> rootsmod(const vector<ZZ>& coeffs, ZZ p);
vector<gf_element> roots(const FqPoly& f);
vector<bigrational> roots(const vector<ZZ>& coeffs);
vector<bigrational> roots(const ZZX& f);
vector<ZZ> introots(const ZZX& f);

vector<ZZ> Introotscubic(const ZZ& a, const ZZ& b, const ZZ& c);
vector<ZZ> Introotsquartic(const ZZ& a, const ZZ& b, const ZZ& c,
                            const ZZ& d);



// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
int nrootscubic(const ZZ& bb, const ZZ& cc, const ZZ& dd, const ZZ& p);

FqPoly reduce(const ZZX& f, const galois_field& Fq);

inline FqPoly reduce(const ZZX& f, const ZZ& q)
{return reduce(f,galois_field(q));}



#endif // #define _POLYS_
