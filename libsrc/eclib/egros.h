// egros.h - declarations of functions for curves with good reduction outside S
//////////////////////////////////////////////////////////////////////////
//
// Copyright 2025 John Cremona
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
 

#include <set>
#include "curve.h"

// Test whether a curve with good reduction outside S and this j-invariant could exist
// (using criteria from Cremona-Lingham)
int is_j_possible(const bigrational& j, const vector<bigint>& S);

// Return integers representing QQ(S,n)
vector<bigint> twist_factors(const vector<bigint>& S, int n); // only intended for n=2,4,6

// Return list of curves with good reduction outside S and j=1728
vector<CurveRed> egros_from_j_1728(const vector<bigint>& S);
// Return list of curves with good reduction outside S and j=0
vector<CurveRed> egros_from_j_0(const vector<bigint>& S);
// Return list of curves with good reduction outside S and any fixed j
vector<CurveRed> egros_from_j(const bigrational& j, const vector<bigint>& S);

// Test whether N is a possible conductor for j=0
int is_N_possible_j_0(const bigint& N, const vector<bigint>& support);

// Test whether N is a possible conductor for j=1728
int is_N_possible_j_1728(const bigint& N, const vector<bigint>& support);
