// minim.h: declaration of quartic minimization functions
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
 
#if     !defined(_ECLIB_MINIM_H)
#define _ECLIB_MINIM_H      1       //flags that this file has been included

#include "unimod.h"

ZZ g_content(const ZZ& ga, const ZZ& gb, const ZZ& gc, 
		 const ZZ& gd, const ZZ& ge);
// returns largest f s.t. f^2 divides all coeffs

ZZ root_p(const ZZ& a, const ZZ& b, const ZZ& c, 
		 const ZZ& d, const ZZ& e, const ZZ& p);
// assuming p|I, p|J, returns the unique alpha mod p 
// modulo which quartic has a root of multiplicity at least 3
// returns -1 if multiple root is at infinity (if a=b=0 mod p)
// (program does not actaully use this dubious feature)

int minim_p(ZZ& ga, ZZ& gb, ZZ& gc, 
	      ZZ& gd, ZZ& ge, const ZZ& p,
	      scaled_unimod& m);
// assuming p^4|I, p^6|J, (or stronger conditions when p=2 or p=3)
// returns an equivalent quartic with invariants divided by p^4, p^6;
// m holds the transformation matrix, must be initialized (say with identity)
// returns success, can be 0 only for p=2

int is_nonmin(int smallp, long vpi, long vpj, long vpd, int assume_locsol);
// Given vpi = val(p,I) and vpj=val(p,J) returns 1 if non-minimal
// smallp = p if p=2,3 else =1.
// p=3: needs also vpd=val(p,disc)
// p=2: may or may not be minimizable, but worth a try
// (The commented out condition is sufficient but NOT necessary)

void minim_all(ZZ& ga, ZZ& gb, ZZ& gc, ZZ& gd, ZZ& ge, 
	       ZZ& I, ZZ& J, const vector<ZZ>& plist, 
	       scaled_unimod& m,
	       int assume_locsol, int verb=0);

#endif
