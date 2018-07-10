// getcurve.h: declaration of function getcurve() for curve input
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#if     !defined(_ECLIB_GETCURVE_H)
#define _ECLIB_GETCURVE_H      1       //flags that this file has been included

#include <eclib/curve.h>

// Read in a curve as [a1,a2,a3,a4,a6] with ai integers:
int getcurve(Curvedata& CD, int verb);

// Read in a curve as [a1,a2,a3,a4,a6] with ai rational:
int getcurve(vector<bigrational>& ai, int verb);

#endif
