// minim.h: declaration of quartic minimization functions
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
 
bigint g_content(const bigint& ga, const bigint& gb, const bigint& gc, 
		 const bigint& gd, const bigint& ge);
// returns largest f s.t. f^2 divides all coeffs

bigint root_p(const bigint& a, const bigint& b, const bigint& c, 
		 const bigint& d, const bigint& e, const bigint& p);
// assuming p|I, p|J, returns the unique alpha mod p 
// modulo which quartic has a root of multiplicity at least 3
// returns -1 if multiple root is at infinity (if a=b=0 mod p)
// (program does not actaully use this dubious feature)

int minim_p(bigint& ga, bigint& gb, bigint& gc, 
	      bigint& gd, bigint& ge, const bigint& p,
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

void minim_all(bigint& ga, bigint& gb, bigint& gc, bigint& gd, bigint& ge, 
	       bigint& I, bigint& J, const vector<bigint>& plist, 
	       scaled_unimod& m,
	       int assume_locsol, int verb=0);

