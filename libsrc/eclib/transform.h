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

#include "unimod.h"

void apply_transform(ZZ& a, ZZ& b, ZZ& c, ZZ& d, ZZ& e,
		     const unimod& m);

void apply_transform(ZZ& a, ZZ& b, ZZ& c, ZZ& d, ZZ& e,
		     const scaled_unimod& m);

int check_transform(const ZZ& a, const ZZ& b, const ZZ& c, 
		    const ZZ& d, const ZZ& e,
		    const unimod& m,
		    const ZZ& xa, const ZZ& xb, const ZZ& xc, 
		    const ZZ& xd, const ZZ& xe);

int check_transform(const ZZ& a, const ZZ& b, const ZZ& c, 
		    const ZZ& d, const ZZ& e,
		    const scaled_unimod& m,
		    const ZZ& xa, const ZZ& xb, const ZZ& xc, 
		    const ZZ& xd, const ZZ& xe);

void xshift(const ZZ& alpha,
	    const ZZ& a, ZZ& b, ZZ& c, ZZ& d, ZZ& e,
	    unimod& m);

void zshift(const ZZ& gamma,
	    ZZ& a, ZZ& b, ZZ& c, ZZ& d, const ZZ& e,
	    unimod& m);

void m_invert(ZZ& a, ZZ& b, ZZ& c, ZZ& d, ZZ& e,
	      unimod& m);

void m_invert(ZZ& a, ZZ& b, ZZ& c, ZZ& d, ZZ& e,
	      scaled_unimod& m);

#endif
