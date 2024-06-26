// transform.h: declaration of quartic transformation functions
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
 //
// Notation: g(x,z) is replaced by g(m11*x+m12*z,m21*x+m22*z)/m00^2
//

#ifndef _ECLIB_TRANSFORM_H
#define _ECLIB_TRANSFORM_H      1
                           //flags that this file has been included

#include <eclib/unimod.h>

void apply_transform(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
		     const unimod& m);

void apply_transform(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
		     const scaled_unimod& m);

int check_transform(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& d, const bigint& e,
		    const unimod& m,
		    const bigint& xa, const bigint& xb, const bigint& xc, 
		    const bigint& xd, const bigint& xe);

int check_transform(const bigint& a, const bigint& b, const bigint& c, 
		    const bigint& d, const bigint& e,
		    const scaled_unimod& m,
		    const bigint& xa, const bigint& xb, const bigint& xc, 
		    const bigint& xd, const bigint& xe);

void xshift(const bigint& alpha,
	    const bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	    unimod& m);

void zshift(const bigint& gamma,
	    bigint& a, bigint& b, bigint& c, bigint& d, const bigint& e,
	    unimod& m);

void m_invert(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      unimod& m);

void m_invert(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      scaled_unimod& m);

#endif
