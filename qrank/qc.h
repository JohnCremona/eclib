// qc.h: declaration of function for mapping quartic point to curve
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
 
// Given a quartic g with a point (x0:y0:z0) on it,
// constructs a point P on the corresponding minimal elliptic curve.
// This is supposed to be on the given curve E: error printed if not.
// The point is returned indirectly.

void qc(quartic& g,
        const bigint& x0,  const bigint& y0,  const bigint& z0,
        Curvedata * E,  
	Curvedata* IJ_curve, 
	const bigint& tr_u, const bigint& tr_r, 
	const bigint& tr_s, const bigint& tr_t,  
	Point& P, int verbose=0);
