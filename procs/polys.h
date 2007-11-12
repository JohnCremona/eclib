// polys.h : contains includes and macros for uniform interface to
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
 //           polynomials (including polynomials mod q) in NTL and LiDIA

// allow for multiple includes
#ifndef _POLYS_
#define _POLYS_

// NB --warning!  In NTL there is a universal modulus set by a single
// call to ZZ_p::init(q), while in LiDIA each Fp_polynomial has a
// pointer to its own modulus, which must be set on creation.


#include "gf.h"  // scalars

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)
#include <LiDIA/factorization.h>
#include "LiDIA/gf_element.h"
#include "LiDIA/gf_polynomial.h"

#define ZPoly polynomial<bigint>
#define PolyCoeff(f,i) (f)[(i)]
#define ZPolySetX(f) f.assign_x()
#define SetCoeff(f,i,c) (f)[(i)]=(c)
#define SetDegree(f,d) (f).set_degree((d))
#define Degree(f) (f).degree()

#define FqPoly polynomial<gf_element>
#define NewFqPoly(field,name) FqPoly name(field)
#define FqPolySetField(f,field) f.set_field(field)
#define FqPolyAssignGF(f,c)  f.assign((c))
#define FqPolyAssignZ(f,c)  f.assign(ZtoGF(f.get_field(),c))
#define FqPolyAssign0(f)  f.assign_zero()
#define FqPolyAssign1(f)  f.assign_one()
#define FqPolyEval(f,c)  (f)((c))
#define FqPolyAssignX(f) (f.assign_x(f.get_field()))
#define GetField(f) ((f).get_field())

#endif

#ifdef NTL_INTS
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#define ZPoly ZZX
#define PolyCoeff(f,i) coeff((f),(i))
#define ZPolySetX(f) SetX(f);
#define SetDegree(f,d)
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

#endif

vector<bigint> rootsmod(const vector<bigint>& coeffs, bigint p);
vector<gf_element> roots(const FqPoly& f);

// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)
int nrootscubic(const bigint& bb, const bigint& cc, const bigint& dd, const bigint& p);

FqPoly reduce(const ZPoly& f, const galois_field& Fq);

inline FqPoly reduce(const ZPoly& f, const bigint& q)
{return reduce(f,galois_field(q));}



#endif
