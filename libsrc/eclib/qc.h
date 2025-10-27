// qc.h: declaration of function for mapping quartic point to curve
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
 
#ifndef _ECLIB_QC_H
#define _ECLIB_QC_H      1
                           //flags that this file has been included

// Given a quartic g with a point (x0:y0:z0) on it,
// constructs a point P on the corresponding minimal elliptic curve.
// This is supposed to be on the given curve E: error printed if not.
// The point is returned indirectly.

#include "mquartic.h"

void qc(quartic& g,
        const ZZ& x0,  const ZZ& y0,  const ZZ& z0,
        Curvedata * E,  
	Curvedata* IJ_curve, 
	const ZZ& tr_u, const ZZ& tr_r, 
	const ZZ& tr_s, const ZZ& tr_t,  
	Point& P, int verbose=0);

#endif
