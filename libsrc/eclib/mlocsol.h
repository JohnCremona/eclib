// mlocsol.h: declaration of functions for local solubility of quartics 
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
 
#if     !defined(_ECLIB_MLOCSOL_H)
#define _ECLIB_MLOCSOL_H      1       //flags that this file has been included

#include "mquartic.h"

// Checks for solublility in Qp:

int qpsoluble(const quartic& g, const bigint& p); 
int qpsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d, 
	      const bigint& e, const bigint& p);
int qpsoluble(const bigint& a, const bigint& c, const bigint& e, const bigint& p);
// latter assumes b=d=0

int Rsoluble(const quartic& g);
int Rsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,  
	      const bigint& e);


// Checks for local solubility in Qp for all p in plist; 
//if not, badp will hold the first p for which NOT soluble in Qp:

int locallysoluble(const quartic& g, const vector<bigint>& plist, bigint& badp);
int locallysoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d, 
	      const bigint& e, const vector<bigint>& plist, bigint& badp);
int locallysoluble(const bigint& a, const bigint& c, const bigint& e, 
		   const vector<bigint>& plist, bigint& badp);
// latter assumes b=d=0


/* Samir Siksek's Local Solubility Test for odd p */

int local_sol(const bigint& p, vector<bigint> c, int verbose=0);

 // Checks for solublility in Qp 
int new_qpsoluble(const quartic& g, const bigint& p, int verbose=0);
int new_qpsoluble(const bigint& a, const bigint& b, const bigint& c, 
		  const bigint& d, const bigint& e, 
		  const bigint& p, int verbose=0);
int new_qpsoluble(const bigint& a, const bigint& c, const bigint& e, 
		  const bigint& p, int verbose=0);

int new_zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,
	      const bigint& e, const bigint& p, int verbose=0);

#endif
