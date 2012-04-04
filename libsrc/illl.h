// illl.h: declarations of functions for integer LLL
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
 
// Integral LLL algorithm (Cohen's book page 94)
// 
// b is an array of n+1 vectors indexed from 0 to n.  
// b[1]...b[n] are the lattice basis, while b[0] holds the coefficients 
// of the (diagonal) Gram matrix, so the inner product of b[i] and b[j] 
// is sum(k,b[0][k]b[i][k]*b[j][k]). 
//

#ifndef _ILLL_H
#define _ILLL_H      1     //flags that this file has been included

#include "mvector.h"

bigint sdot(const vec_m* b, int i, int j);

void lll_reduce(const int n, vec_m* b);

//
// Uses Pohst-Zassenhaus Algorithm (page 190) to find all vectors of
// length < c, where the quadratic form is again given by b[0].
//
// NB The following DOES NOT WORK: I blindly implemented P-Z without
// doing the necessary preliminary completing of the square.
// So DO NOT USE
//
//void list_short_vecs(const int n, vec_m* b, const bigint& c);

#endif
