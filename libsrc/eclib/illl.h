// illl.h: declarations of functions for integer LLL
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
 
// Integral LLL algorithm (Cohen's book page 94)
// 
// b is an array of n+1 vectors indexed from 0 to n.  
// b[1]...b[n] are the lattice basis, while b[0] holds the coefficients 
// of the (diagonal) Gram matrix, so the inner product of b[i] and b[j] 
// is sum(k,b[0][k]b[i][k]*b[j][k]). 
//

#ifndef _ECLIB_ILLL_H
#define _ECLIB_ILLL_H      1     //flags that this file has been included

#include "vector.h"

// dot product of b[i] and b[j] with weights from b[0]:
bigint sdot(const vector<vec_m>& b, int i, int j);

// LLL-reduce the vectors b[1],...,b[n] in place.

// b is an array of n+1 vectors indexed from 0 to n.
// b[1]...b[n] are the lattice basis, while b[0] holds the coefficients
// of the (diagonal) Gram matrix, so the inner product of b[i] and b[j]
// is sdot(b, i, j) = sum(k,b[0][k]*b[i][k]*b[j][k]).
//

void lll_reduce(const int n, vector<vec_m>& b);

#endif
