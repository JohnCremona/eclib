// divpol.h: declaration of functions for division polynomials
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
 
// allow for multiple includes
#ifndef _DIVPOL_
#define _DIVPOL_

// These functions return the vector of coefficients (starting with
// the constant term rather than polynomial types (which depend on
// whether NTL or LiDIA is being used)

// div_pol_odd(a1,a2,a3,a4,a6,n) returns the coefficients of the
// polynomial in x whose zeros are the (x-coordinates of the) non-zero
// points P on E=[a1,a2,a3,a4,a6] satisfying nP=0 (odd n)

// The poly itself is found recursively

vector<bigint> div_pol_odd(const bigint& a1,const bigint& a2,const bigint& a3,const bigint& a4,
			   const bigint& a6,int n); 


vector<bigint> makepdivpol(Curvedata* EE, int p);

#endif
