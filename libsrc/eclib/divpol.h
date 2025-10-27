// divpol.h: declaration of functions for division polynomials
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
#ifndef _ECLIB_DIVPOL_
#define _ECLIB_DIVPOL_

#include "curve.h"

// These functions return the vector of coefficients (starting with
// the constant term) rather than polynomial types

// div_pol_odd(a1,a2,a3,a4,a6,n) returns the coefficients of the
// polynomial in x whose zeros are the (x-coordinates of the) non-zero
// points P on E=[a1,a2,a3,a4,a6] satisfying nP=0 (odd n)

// The poly itself is found recursively

// Despite the name, for even n this returns a correct n-division
// polynomial without the 2-torsion factor, i.e. the polynomial whoe
// roots are the x-coordinates of the points P satisfying nP=0, 2P!=0.

ZPoly div_pol_odd(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
                  const ZZ& a6,int n);

ZPoly div_pol_2(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
                const ZZ& a6);

ZPoly div_pol(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
              const ZZ& a6,int n);

ZPoly division_polynomial(Curvedata* EE, int p);

// Numerator and denominator of the multiplication-by-n map on the x-coordinate

ZPoly mul_by_n_num(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
                   const ZZ& a6, int n);

ZPoly mul_by_n_den(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
                   const ZZ& a6, int n);

// Polynomial whose roots are x(Q) for Q satisfying n*Q=P, where x(P)=xP/zP

ZPoly division_points_X_pol(const ZZ& a1,const ZZ& a2,const ZZ& a3,const ZZ& a4,
                            const ZZ& a6,
                            int n,
                            const ZZ& xP, const ZZ& zP);

#endif // #define _DIVPOL_
