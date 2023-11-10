// reduce.h:  declaration of quartic reduction functions
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
 
#ifndef _ECLIB_REDUCE_H
#define _ECLIB_REDUCE_H      1
                           //flags that this file has been included

#include <eclib/unimod.h>

// NB In both functions, the unimod m must hold the cumulative
// transformation on entry. It will be updated to show the new
// cumulative transformation

// Full reduction via covariant quadratic:

void reduce(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	    unimod& m);

// Simple shift to minimise b:

void reduce_b(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e,
	      unimod& m);

// Compute the quadratic covariant of a real quartic:
bigfloat* quadratic_covariant(bigint& a, bigint& b, bigint& c, bigint& d, bigint& e);

// Given a pos. def. quadratic x^2+b*x+c, returns a unimod which
// reduces it (whose inverse takes its root into the fundamental
// region).
unimod reduce_quad(const bigfloat& b, const bigfloat& c);

#endif
